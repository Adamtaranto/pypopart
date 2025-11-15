"""
Median-Joining Network (MJN) algorithm for haplotype network construction.

Implements the method from Bandelt, Forster & Rohl (1999):
"Median-joining networks for inferring intraspecific phylogenies"
Molecular Biology and Evolution 16: 37-48
"""

import itertools
from typing import List, Optional, Set, Tuple

from ..core.alignment import Alignment
from ..core.distance import DistanceMatrix, hamming_distance
from ..core.graph import HaplotypeNetwork
from ..core.haplotype import (
    Haplotype,
)
from ..core.haplotype import identify_haplotypes_from_alignment as identify_haplotypes
from ..core.sequence import Sequence
from .msn import MinimumSpanningNetwork


class MedianJoiningNetwork(MinimumSpanningNetwork):
    """
    Construct haplotype network using Median-Joining algorithm.

    The MJN method extends MSN by inferring median vectors (ancestral or
    unsampled haplotypes) that simplify the network structure. It combines
    minimum spanning network principles with median vector inference to
    create a more parsimonious representation of haplotype relationships.

    This is particularly useful for datasets with missing intermediate
    haplotypes and complex reticulation patterns.
    """

    def __init__(
        self,
        distance_method: str = 'hamming',
        epsilon: float = 0.0,
        max_median_vectors: Optional[int] = None,
        simplify: bool = True,
        **kwargs,
    ):
        """
        Initialize MJN algorithm.

        Parameters
        ----------
        distance_method :
            Method for calculating distances.
        epsilon :
            Weight parameter for controlling network complexity.
        max_median_vectors :
            Maximum number of median vectors to add.
        simplify :
            Whether to simplify network after median vector addition.
        **kwargs :
            Additional parameters.
        """
        super().__init__(distance_method, epsilon=epsilon, **kwargs)
        self.max_median_vectors = max_median_vectors
        self.simplify = simplify
        self._median_counter = 0

    def construct_network(
        self, alignment: Alignment, distance_matrix: Optional[DistanceMatrix] = None
    ) -> HaplotypeNetwork:
        """
            Construct MJN from sequence alignment with iterative refinement.

        Parameters
        ----------
            alignment :
                Multiple sequence alignment.
            distance_matrix :
                Optional pre-computed distance matrix.

        Returns
        -------
            Haplotype network with inferred median vectors.
        """
        # Identify unique haplotypes
        haplotypes = identify_haplotypes(alignment)

        if len(haplotypes) <= 2:
            # Too few haplotypes for median vector inference
            return super().construct_network(alignment, distance_matrix)

        # Calculate or use provided distance matrix
        if distance_matrix is None:
            distance_matrix = self.calculate_distances(alignment)

        self._distance_matrix = distance_matrix

        # Iteratively refine network with median vectors (matches C++ behavior)
        network = self._iterative_median_joining(haplotypes, alignment, distance_matrix)

        return network

    def _iterative_median_joining(
        self, haplotypes: List, alignment: Alignment, distance_matrix: DistanceMatrix
    ) -> HaplotypeNetwork:
        """
            Refine network iteratively by adding median vectors (C++ algorithm).

            Repeatedly:
            1. Build MSN from current haplotypes
            2. Find quasi-median vectors for triplets
            3. Add medians that reduce network cost
            4. Remove obsolete vertices (degree < 2)
            5. Repeat until no new medians are added

        Parameters
        ----------
            haplotypes :
                Initial list of observed haplotypes.
            alignment :
                Sequence alignment.
            distance_matrix :
                Distance matrix.

        Returns
        -------
            Refined network with median vectors.
        """
        # Track all sequences seen (observed + inferred)
        all_sequences = {h.sequence.data for h in haplotypes}
        current_haplotypes = list(haplotypes)

        iteration = 0
        max_iterations = 50  # Prevent infinite loops
        changed = True
        old_msn_length = -1

        while changed and iteration < max_iterations:
            iteration += 1
            changed = False

            # Build MSN from current haplotypes
            msn = self._build_msn_for_iteration(current_haplotypes, distance_matrix)

            # Calculate MSN total length
            msn_length = sum(
                msn.get_edge_distance(u, v) or 0 for u, v in msn.graph.edges()
            )

            # Check for convergence
            if old_msn_length > 0 and msn_length >= old_msn_length:
                # Network is not improving, stop
                break

            old_msn_length = msn_length

            # Remove obsolete median vectors (degree < 2)
            current_haplotypes = self._remove_obsolete_medians(msn, current_haplotypes)

            # Find minimum cost for median vectors
            min_cost = float('inf')
            candidate_medians = []

            # Check all triplets in current MSN
            for triplet in self._find_all_triplets_in_msn(msn):
                h1_id, h2_id, h3_id = triplet

                h1 = msn.get_haplotype(h1_id)
                h2 = msn.get_haplotype(h2_id)
                h3 = msn.get_haplotype(h3_id)

                # Compute quasi-median sequences
                quasi_medians = self._compute_quasi_medians(
                    h1.sequence, h2.sequence, h3.sequence
                )

                # Check each quasi-median
                for median_seq in quasi_medians:
                    if median_seq not in all_sequences:
                        # Calculate cost
                        cost = self._compute_median_cost(
                            h1.sequence, h2.sequence, h3.sequence, median_seq
                        )

                        if cost < min_cost:
                            min_cost = cost

                        candidate_medians.append((median_seq, cost, triplet))

            # Add medians within epsilon of minimum cost
            for median_seq, cost, _triplet in candidate_medians:
                if cost <= min_cost + self.epsilon:
                    # Create median haplotype
                    median_hap = Haplotype(
                        sequence=Sequence(
                            id=f'Median_{self._median_counter}', data=median_seq
                        ),
                        sample_ids=[],
                    )
                    self._median_counter += 1

                    # Add to working set
                    current_haplotypes.append(median_hap)
                    all_sequences.add(median_seq)
                    changed = True

                    if (
                        self.max_median_vectors
                        and self._median_counter >= self.max_median_vectors
                    ):
                        break

            if not changed:
                break

        # Final MSN construction
        final_network = self._build_msn_for_iteration(
            current_haplotypes, distance_matrix
        )

        # Final cleanup
        if self.simplify:
            final_network = self._simplify_network(final_network)

        return final_network

    def _build_msn_for_iteration(
        self, haplotypes: List, distance_matrix: DistanceMatrix
    ) -> HaplotypeNetwork:
        """
            Build MSN from current haplotypes for one iteration.

        Parameters
        ----------
            haplotypes :
                Current list of haplotypes.
            distance_matrix :
                Distance matrix (may need recalculation).

        Returns
        -------
            MSN network.
        """
        # Recalculate distances including any new median vectors
        n = len(haplotypes)
        labels = [h.id for h in haplotypes]

        import numpy as np

        matrix = np.zeros((n, n))

        for i in range(n):
            for j in range(i + 1, n):
                dist = hamming_distance(
                    haplotypes[i].sequence, haplotypes[j].sequence, ignore_gaps=True
                )
                matrix[i, j] = matrix[j, i] = dist

        new_dist_matrix = DistanceMatrix(labels, matrix)

        # Build MSN using parent class method
        # Create temporary alignment
        sequences = [h.sequence for h in haplotypes]
        from ..core.alignment import Alignment

        temp_alignment = Alignment(sequences)

        # Use parent MSN construction
        network = super().construct_network(temp_alignment, new_dist_matrix)

        return network

    def _find_all_triplets_in_msn(
        self, network: HaplotypeNetwork
    ) -> List[Tuple[str, str, str]]:
        """
            Find all triplets (connected triples) in MSN.

            A triplet consists of three nodes where at least two edges exist
            between them (not necessarily a triangle).

        Parameters
        ----------
            network :
                Current MSN.

        Returns
        -------
            List of triplets.
        """
        triplets = []
        hap_ids = list({h.id for h in network.haplotypes})

        for i, h1 in enumerate(hap_ids):
            neighbors_h1 = set(network.get_neighbors(h1))

            for j in range(i + 1, len(hap_ids)):
                h2 = hap_ids[j]

                if h2 in neighbors_h1:
                    # h1 and h2 are connected
                    neighbors_h2 = set(network.get_neighbors(h2))

                    # Find common neighbors
                    for h3 in neighbors_h1.intersection(neighbors_h2):
                        if h3 != h1 and h3 != h2:
                            triplets.append((h1, h2, h3))

        return triplets

    def _compute_quasi_medians(
        self, seq1: Sequence, seq2: Sequence, seq3: Sequence
    ) -> Set[str]:
        """
            Compute quasi-median sequences (Steiner tree approach).

            For positions where all three sequences differ, generates all
            possible combinations (creating a set of quasi-medians).

            This matches the C++ computeQuasiMedianSeqs implementation.

        Parameters
        ----------
                seq1, seq2, seq3: Three sequences

        Returns
        -------
            Set of quasi-median sequence strings.
        """
        if len(seq1) != len(seq2) or len(seq1) != len(seq3):
            return set()

        # Build initial quasi-median with '*' at ambiguous positions
        qm_seq = []
        has_star = False

        for i in range(len(seq1)):
            c1, c2, c3 = seq1.data[i], seq2.data[i], seq3.data[i]

            if c1 == c2 or c1 == c3:
                qm_seq.append(c1)
            elif c2 == c3:
                qm_seq.append(c2)
            else:
                # All three differ
                qm_seq.append('*')
                has_star = True

        if not has_star:
            # Simple median exists
            return {''.join(qm_seq)}

        # Resolve '*' positions by generating all combinations
        medians = set()
        stack = [''.join(qm_seq)]

        while stack:
            current = stack.pop()
            first_star = current.find('*')

            if first_star == -1:
                # No more stars, add to result
                medians.add(current)
            else:
                # Replace star with each of the three bases
                for base in [
                    seq1.data[first_star],
                    seq2.data[first_star],
                    seq3.data[first_star],
                ]:
                    new_seq = current[:first_star] + base + current[first_star + 1 :]
                    stack.append(new_seq)

        return medians

    def _compute_median_cost(
        self, seq1: Sequence, seq2: Sequence, seq3: Sequence, median_seq: str
    ) -> int:
        """
            Compute cost of a median vector.

            Cost is sum of distances from median to the three sequences.

        Parameters
        ----------
                seq1, seq2, seq3: Three sequences forming the triplet
            median_seq :
                Candidate median sequence string.

        Returns
        -------
            Total cost (sum of Hamming distances).
        """
        median = Sequence(id='temp', data=median_seq)

        dist1 = hamming_distance(seq1, median, ignore_gaps=True)
        dist2 = hamming_distance(seq2, median, ignore_gaps=True)
        dist3 = hamming_distance(seq3, median, ignore_gaps=True)

        return dist1 + dist2 + dist3

    def _remove_obsolete_medians(
        self, network: HaplotypeNetwork, haplotypes: List
    ) -> List:
        """
            Remove median vectors with degree < 2 (obsolete).

            Matches C++ removeObsoleteVerts behavior.

        Parameters
        ----------
            network :
                Current network.
            haplotypes :
                List of all haplotypes.

        Returns
        -------
            Updated list with obsolete medians removed.
        """
        changed = True

        while changed:
            changed = False
            to_remove = []

            for haplotype in haplotypes:
                if haplotype.frequency == 0:  # Only check inferred medians
                    degree = network.get_degree(haplotype.id)
                    if degree < 2:
                        to_remove.append(haplotype)
                        changed = True

            # Remove obsolete medians
            for hap in to_remove:
                haplotypes = [h for h in haplotypes if h.id != hap.id]

        return haplotypes

    def _add_median_vectors(
        self, network: HaplotypeNetwork, sequence_length: int
    ) -> HaplotypeNetwork:
        """
            Infer and add median vectors to the network.

        Parameters
        ----------
            network :
                Initial haplotype network.
            sequence_length :
                Length of sequences.

        Returns
        -------
            Network with median vectors added.
        """
        self._median_counter = 0
        median_vectors_added = 0

        # Get all triangles (3-cliques) in the network
        triangles = self._find_triplets(network)

        for hap1_id, hap2_id, hap3_id in triangles:
            if (
                self.max_median_vectors
                and median_vectors_added >= self.max_median_vectors
            ):
                break

            # Get the three haplotypes
            hap1 = network.get_haplotype(hap1_id)
            hap2 = network.get_haplotype(hap2_id)
            hap3 = network.get_haplotype(hap3_id)

            # Calculate median vector
            median_seq = self._calculate_median(
                hap1.sequence, hap2.sequence, hap3.sequence
            )

            if median_seq is None:
                continue

            # Check if median already exists in network
            if self._sequence_in_network(median_seq, network):
                continue

            # Calculate distances from median to the three haplotypes
            median_hap = Haplotype(
                sequence=median_seq, sample_ids=[]
            )
            self._median_counter += 1

            dist1 = hamming_distance(
                median_hap.sequence, hap1.sequence, ignore_gaps=True
            )
            dist2 = hamming_distance(
                median_hap.sequence, hap2.sequence, ignore_gaps=True
            )
            dist3 = hamming_distance(
                median_hap.sequence, hap3.sequence, ignore_gaps=True
            )

            # Check if median vector simplifies the network
            # (reduces total edge weight)
            original_edges = [
                network.get_edge_distance(hap1_id, hap2_id) or 0,
                network.get_edge_distance(hap1_id, hap3_id) or 0,
                network.get_edge_distance(hap2_id, hap3_id) or 0,
            ]
            original_weight = sum(d for d in original_edges if d > 0)
            new_weight = dist1 + dist2 + dist3

            if new_weight < original_weight:
                # Add median vector to network
                network.add_haplotype(median_hap)

                # Connect median to the three haplotypes
                if dist1 > 0:
                    network.add_edge(median_hap.id, hap1_id, distance=dist1)
                if dist2 > 0:
                    network.add_edge(median_hap.id, hap2_id, distance=dist2)
                if dist3 > 0:
                    network.add_edge(median_hap.id, hap3_id, distance=dist3)

                median_vectors_added += 1

        return network

    def _find_triplets(self, network: HaplotypeNetwork) -> List[Tuple[str, str, str]]:
        """
            Find all triplets (triangles) in the network.

        Parameters
        ----------
            network :
                Haplotype network.

        Returns
        -------
            List of triplets (id1, id2, id3).
        """
        triplets = []
        hap_ids = list({h.id for h in network.haplotypes})

        # Check all combinations of 3 haplotypes
        for hap1, hap2, hap3 in itertools.combinations(hap_ids, 3):
            # Check if they form a triangle (all three edges exist)
            has_12 = network.has_edge(hap1, hap2)
            has_13 = network.has_edge(hap1, hap3)
            has_23 = network.has_edge(hap2, hap3)

            if has_12 and has_13 and has_23:
                triplets.append((hap1, hap2, hap3))

        return triplets

    def _calculate_median(
        self, seq1: Sequence, seq2: Sequence, seq3: Sequence
    ) -> Optional[Sequence]:
        """
            Calculate median sequence of three sequences.

            For each position, the median is the most common nucleotide
            among the three sequences. If all three are different, no
            clear median exists for that position.

        Parameters
        ----------
                seq1, seq2, seq3: Three Sequence objects

        Returns
        -------
            Median Sequence, or None if no clear median exists.
        """
        if len(seq1) != len(seq2) or len(seq1) != len(seq3):
            return None

        median_data = []

        for i in range(len(seq1)):
            c1, c2, c3 = seq1.data[i], seq2.data[i], seq3.data[i]

            # Find most common nucleotide
            if c1 == c2:
                median_data.append(c1)
            elif c1 == c3:
                median_data.append(c1)
            elif c2 == c3:
                median_data.append(c2)
            else:
                # All three different - no clear median
                # Use majority rule or first sequence's base
                # For simplicity, we'll say no median exists for this triplet
                return None

        median_seq = Sequence(id='median_temp', data=''.join(median_data))

        return median_seq

    def _sequence_in_network(
        self, sequence: Sequence, network: HaplotypeNetwork
    ) -> bool:
        """
            Check if a sequence already exists in the network.

        Parameters
        ----------
            sequence :
                Sequence to check.
            network :
                Haplotype network.

        Returns
        -------
            True if sequence exists in network.
        """
        seq_data = sequence.data

        for haplotype in network.haplotypes:
            if haplotype.sequence.data == seq_data:
                return True

        return False

    def _simplify_network(self, network: HaplotypeNetwork) -> HaplotypeNetwork:
        """
            Simplify network by removing unnecessary median vectors.

            A median vector is unnecessary if:
            - It has degree 2 (only connects two nodes)
            - It can be replaced by a direct edge without increasing total weight

        Parameters
        ----------
            network :
                Network with median vectors.

        Returns
        -------
            Simplified network.
        """
        # Find median vectors with degree 2
        to_remove = []

        for hap_id in list({h.id for h in network.haplotypes}):
            if not hap_id.startswith('Median_'):
                continue

            degree = network.get_degree(hap_id)

            if degree == 2:
                # Get the two neighbors
                neighbors = network.get_neighbors(hap_id)
                if len(neighbors) == 2:
                    n1, n2 = neighbors

                    # Get distances
                    d1 = network.get_edge_distance(hap_id, n1) or 0
                    d2 = network.get_edge_distance(hap_id, n2) or 0

                    # Check if direct connection exists
                    direct_dist = network.get_edge_distance(n1, n2)

                    if direct_dist is None:
                        # No direct connection - add it with combined distance
                        network.add_edge(n1, n2, distance=d1 + d2)

                    # Remove the median vector
                    to_remove.append(hap_id)

        # Remove marked median vectors
        for hap_id in to_remove:
            network.remove_haplotype(hap_id)

        return network

    def get_parameters(self) -> dict:
        """Get algorithm parameters."""
        params = super().get_parameters()
        params['max_median_vectors'] = self.max_median_vectors
        params['simplify'] = self.simplify
        return params

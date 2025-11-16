"""
Parsimony Network algorithm for haplotype network construction.

This module implements the Parsimony Network algorithm, which constructs
haplotype networks by sampling edges from multiple parsimony trees. The
algorithm creates a consensus network that includes edges that appear
frequently across random parsimony tree topologies.

The approach provides a more robust network than a single tree by capturing
uncertainty in phylogenetic reconstruction and potential reticulation events.

References
----------
.. [1] Excoffier, L. & Smouse, P. E. (1994). Using allele frequencies and
       geographic subdivision to reconstruct gene trees within a species:
       molecular variance parsimony. Genetics, 136(1), 343-359.
.. [2] Templeton, A. R., Boerwinkle, E., & Sing, C. F. (1987). A cladistic
       analysis of phenotypic associations with haplotypes inferred from
       restriction endonuclease mapping. I. Basic theory and an analysis of
       alcohol dehydrogenase activity in Drosophila. Genetics, 117(2), 343-351.
"""

import random
from collections import defaultdict
from typing import Dict, List, Optional, Set, Tuple

import numpy as np

from ..core.alignment import Alignment
from ..core.distance import DistanceMatrix
from ..core.graph import HaplotypeNetwork
from ..core.haplotype import Haplotype, identify_haplotypes_from_alignment
from ..core.sequence import Sequence
from .base import NetworkAlgorithm


class ParsimonyNetwork(NetworkAlgorithm):
    """
    Construct a haplotype network using the Parsimony Network algorithm.

    The Parsimony Network algorithm constructs a consensus network by sampling
    edges from multiple random parsimony tree topologies. This approach captures
    phylogenetic uncertainty and can represent reticulation events.

    The algorithm works by:
    1. Generating multiple random parsimony trees
    2. Sampling edges from these trees with frequency-based weighting
    3. Building a network that includes frequently occurring edges
    4. Inferring ancestral sequences at internal nodes

    Parameters
    ----------
    distance_method : str, default='hamming'
        Method for calculating pairwise distances between sequences.
        Options: 'hamming', 'jukes_cantor', 'kimura_2p', 'tamura_nei'.
    n_trees : int, default=100
        Number of random parsimony trees to sample.
    min_edge_frequency : float, default=0.05
        Minimum frequency (0-1) for an edge to be included in the network.
        Lower values create more reticulate networks.
    random_seed : int, optional
        Random seed for reproducibility.
    **kwargs : dict
        Additional parameters passed to base NetworkAlgorithm.

    Attributes
    ----------
    n_trees : int
        Number of trees to sample
    min_edge_frequency : float
        Minimum edge frequency threshold
    _random_seed : int
        Random seed for reproducibility
    _median_counter : int
        Counter for generating unique median vertex IDs

    Examples
    --------
    >>> from pypopart.algorithms import ParsimonyNetwork
    >>> from pypopart.io import load_alignment
    >>>
    >>> # Load alignment
    >>> alignment = load_alignment('sequences.fasta')
    >>>
    >>> # Construct Parsimony Network
    >>> pn = ParsimonyNetwork(n_trees=100, min_edge_frequency=0.1)
    >>> network = pn.build_network(alignment)

    Notes
    -----
    The Parsimony Network algorithm is computationally moderate, suitable for
    small to medium datasets (< 200 sequences). The computational cost scales
    with the number of trees sampled.

    The algorithm can create reticulate networks when multiple edges have similar
    frequencies, representing alternative evolutionary pathways.

    See Also
    --------
    MinimumSpanningTree : Simpler tree-based method
    TCS : Statistical parsimony network
    TightSpanWalker : Exact tight span computation
    """

    def __init__(
        self,
        distance_method: str = 'hamming',
        n_trees: int = 100,
        min_edge_frequency: float = 0.05,
        random_seed: Optional[int] = None,
        **kwargs,
    ):
        """
        Initialize Parsimony Network algorithm.

        Parameters
        ----------
        distance_method : str, default='hamming'
            Method for calculating distances.
        n_trees : int, default=100
            Number of random parsimony trees to sample.
        min_edge_frequency : float, default=0.05
            Minimum frequency for edge inclusion.
        random_seed : int, optional
            Random seed for reproducibility.
        **kwargs : dict
            Additional algorithm parameters.
        """
        super().__init__(distance_method=distance_method, **kwargs)
        self.n_trees = n_trees
        self.min_edge_frequency = min_edge_frequency
        self._random_seed = random_seed
        self._median_counter = 0

        if random_seed is not None:
            random.seed(random_seed)
            np.random.seed(random_seed)

    def construct_network(
        self, alignment: Alignment, distance_matrix: Optional[DistanceMatrix] = None
    ) -> HaplotypeNetwork:
        """
        Construct haplotype network using Parsimony Network algorithm.

        Parameters
        ----------
        alignment : Alignment
            Multiple sequence alignment.
        distance_matrix : DistanceMatrix, optional
            Pre-computed distance matrix. If None, will be calculated.

        Returns
        -------
        HaplotypeNetwork
            Constructed haplotype network.
        """
        # Handle empty or single sequence alignments
        if len(alignment) == 0:
            return HaplotypeNetwork()

        if len(alignment) == 1:
            network = HaplotypeNetwork()
            seq = alignment[0]
            haplotype = Haplotype(sequence=seq, sample_ids=[seq.id])
            network.add_haplotype(haplotype)
            return network

        # Calculate distances if not provided
        if distance_matrix is None:
            distance_matrix = self.calculate_distances(alignment)

        # Identify unique haplotypes
        haplotypes = identify_haplotypes_from_alignment(alignment)

        # Initialize network
        network = HaplotypeNetwork()

        # Add all observed haplotypes as nodes
        for haplotype in haplotypes:
            network.add_haplotype(haplotype)

        # Get distance matrix as numpy array
        dist_array = distance_matrix.matrix

        # Sample edges from random parsimony trees
        edge_counts = self._sample_edges_from_trees(
            haplotypes, dist_array, alignment.length
        )

        # Add edges that meet frequency threshold
        self._add_consensus_edges(network, edge_counts, haplotypes, alignment)

        return network

    def _sample_edges_from_trees(
        self, haplotypes: List[Haplotype], dist_array: np.ndarray, seq_length: int
    ) -> Dict[Tuple[str, str], int]:
        """
        Sample edges from multiple random parsimony trees.

        For each tree iteration, we build an approximate parsimony tree using
        neighbor joining with random tie-breaking, then extract all edges.

        Parameters
        ----------
        haplotypes : List[Haplotype]
            List of haplotypes.
        dist_array : np.ndarray
            Distance matrix.
        seq_length : int
            Length of sequences.

        Returns
        -------
        Dict[Tuple[str, str], int]
            Dictionary mapping edge (as tuple of IDs) to count across trees.
        """
        edge_counts: Dict[Tuple[str, str], int] = defaultdict(int)
        n = len(haplotypes)

        for tree_idx in range(self.n_trees):
            # Build a random parsimony tree using modified neighbor joining
            tree_edges = self._build_random_parsimony_tree(
                haplotypes, dist_array.copy()
            )

            # Count edges
            for edge in tree_edges:
                # Sort to make undirected
                sorted_edge = tuple(sorted(edge))
                edge_counts[sorted_edge] += 1

        return edge_counts

    def _build_random_parsimony_tree(
        self, haplotypes: List[Haplotype], dist_matrix: np.ndarray
    ) -> List[Tuple[str, str]]:
        """
        Build a random parsimony tree using neighbor joining with random tie-breaking.

        This is a simplified parsimony tree construction that uses distance-based
        neighbor joining with randomized selection when multiple pairs have equal scores.

        Parameters
        ----------
        haplotypes : List[Haplotype]
            List of haplotypes.
        dist_matrix : np.ndarray
            Distance matrix.

        Returns
        -------
        List[Tuple[str, str]]
            List of edges (pairs of haplotype IDs).
        """
        n = len(haplotypes)
        if n < 2:
            return []

        # Track active clusters (initially each haplotype is a cluster)
        active_clusters = {i: [haplotypes[i].id] for i in range(n)}
        edges = []

        # Keep joining until only one cluster remains
        while len(active_clusters) > 1:
            # Find pair of clusters to join
            best_pairs = []
            best_score = float('inf')

            cluster_ids = list(active_clusters.keys())
            for i in range(len(cluster_ids)):
                for j in range(i + 1, len(cluster_ids)):
                    idx_i = cluster_ids[i]
                    idx_j = cluster_ids[j]

                    # Get distance between clusters (use minimum for simplicity)
                    min_dist = float('inf')
                    for id_i in active_clusters[idx_i]:
                        i_pos = next(k for k, h in enumerate(haplotypes) if h.id == id_i)
                        for id_j in active_clusters[idx_j]:
                            j_pos = next(k for k, h in enumerate(haplotypes) if h.id == id_j)
                            if i_pos < n and j_pos < n:
                                min_dist = min(min_dist, dist_matrix[i_pos, j_pos])

                    # Track best pairs (with ties)
                    if min_dist < best_score:
                        best_score = min_dist
                        best_pairs = [(idx_i, idx_j)]
                    elif min_dist == best_score:
                        best_pairs.append((idx_i, idx_j))

            # Randomly select among tied pairs
            if best_pairs:
                idx_i, idx_j = random.choice(best_pairs)

                # Add edges between clusters
                # For simplicity, connect all pairs between clusters
                for id_i in active_clusters[idx_i]:
                    for id_j in active_clusters[idx_j]:
                        edges.append((id_i, id_j))
                        break  # Only add one edge per cluster pair
                    break

                # Merge clusters
                active_clusters[idx_i].extend(active_clusters[idx_j])
                del active_clusters[idx_j]

        return edges

    def _add_consensus_edges(
        self,
        network: HaplotypeNetwork,
        edge_counts: Dict[Tuple[str, str], int],
        haplotypes: List[Haplotype],
        alignment: Alignment,
    ) -> None:
        """
        Add consensus edges that meet frequency threshold.

        Parameters
        ----------
        network : HaplotypeNetwork
            Network to add edges to.
        edge_counts : Dict[Tuple[str, str], int]
            Edge counts from tree sampling.
        haplotypes : List[Haplotype]
            List of haplotypes.
        alignment : Alignment
            Sequence alignment for distance calculation.

        Returns
        -------
        None
            Network is modified in place.
        """
        # Calculate frequency threshold
        min_count = int(self.n_trees * self.min_edge_frequency)

        # Create haplotype ID to sequence mapping
        hap_seqs = {h.id: h.sequence for h in haplotypes}

        # Add edges that meet threshold
        for edge, count in edge_counts.items():
            if count >= min_count:
                id1, id2 = edge

                # Calculate distance between sequences
                seq1 = hap_seqs[id1]
                seq2 = hap_seqs[id2]
                distance = self._calculate_pairwise_distance(seq1, seq2)

                # Add edge if not already present
                if not network.has_edge(id1, id2):
                    network.add_edge(id1, id2, distance=distance)

                    # Optionally add median vertices for edges > 1 mutation
                    if distance > 1:
                        self._add_median_vertices_along_edge(
                            network, id1, id2, seq1, seq2, int(distance)
                        )

    def _calculate_pairwise_distance(self, seq1: Sequence, seq2: Sequence) -> float:
        """
        Calculate distance between two sequences.

        Parameters
        ----------
        seq1 : Sequence
            First sequence.
        seq2 : Sequence
            Second sequence.

        Returns
        -------
        float
            Distance between sequences.
        """
        if len(seq1.data) != len(seq2.data):
            raise ValueError("Sequences must have equal length")

        # Simple Hamming distance
        distance = sum(c1 != c2 for c1, c2 in zip(seq1.data, seq2.data))
        return float(distance)

    def _add_median_vertices_along_edge(
        self,
        network: HaplotypeNetwork,
        id1: str,
        id2: str,
        seq1: Sequence,
        seq2: Sequence,
        num_mutations: int,
    ) -> None:
        """
        Add median vertices along an edge with multiple mutations.

        Creates intermediate nodes for edges representing more than one mutation,
        which is required for proper network topology.

        Parameters
        ----------
        network : HaplotypeNetwork
            Network to modify.
        id1 : str
            First node ID.
        id2 : str
            Second node ID.
        seq1 : Sequence
            First sequence.
        seq2 : Sequence
            Second sequence.
        num_mutations : int
            Number of mutations between sequences.

        Returns
        -------
        None
            Network is modified in place.
        """
        if num_mutations <= 1:
            return

        # Remove the direct edge
        if network.has_edge(id1, id2):
            network.remove_edge(id1, id2)

        # Find positions that differ
        diff_positions = [
            i for i, (c1, c2) in enumerate(zip(seq1.data, seq2.data)) if c1 != c2
        ]

        # Create intermediate sequences
        prev_id = id1
        for i in range(1, num_mutations):
            # Create median sequence by changing i positions
            median_data = list(seq1.data)
            for j in range(i):
                median_data[diff_positions[j]] = seq2.data[diff_positions[j]]

            median_seq_str = ''.join(median_data)

            # Check if this sequence already exists in network
            existing_id = None
            for hap_id in network.haplotype_ids:
                hap = network.get_haplotype(hap_id)
                if hap.data == median_seq_str:
                    existing_id = hap_id
                    break

            if existing_id:
                median_id = existing_id
            else:
                # Create new median vertex
                median_id = f"Median_{self._median_counter}"
                self._median_counter += 1

                median_seq = Sequence(id=median_id, data=median_seq_str)
                median_haplotype = Haplotype(sequence=median_seq, sample_ids=[])
                network.add_haplotype(median_haplotype, median_vector=True)

            # Add edge from previous to median
            if not network.has_edge(prev_id, median_id):
                network.add_edge(prev_id, median_id, distance=1.0)

            prev_id = median_id

        # Add final edge to id2
        if not network.has_edge(prev_id, id2):
            network.add_edge(prev_id, id2, distance=1.0)

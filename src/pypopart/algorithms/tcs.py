"""
TCS (Statistical Parsimony) algorithm for haplotype network construction.

Implements the method from Clement, Posada & Crandall (2000):
"TCS: a computer program to estimate gene genealogies"
Molecular Ecology 9: 1657-1659
"""

import math
from typing import List, Optional, Tuple

from ..core.alignment import Alignment
from ..core.distance import DistanceMatrix, hamming_distance
from ..core.graph import HaplotypeNetwork
from ..core.haplotype import identify_haplotypes_from_alignment as identify_haplotypes
from .base import NetworkAlgorithm


class TCS(NetworkAlgorithm):
    """
    Construct haplotype network using Statistical Parsimony (TCS algorithm).

    The TCS method uses a 95% parsimony criterion to determine which
    haplotypes can be connected with statistical confidence. It estimates
    the maximum number of mutational differences (connection limit) that
    can be explained by parsimony at the 95% confidence level.

    This method is particularly useful for intraspecific data and closely
    related sequences where parsimony is a reasonable assumption.
    """

    def __init__(
        self,
        distance_method: str = 'hamming',
        confidence: float = 0.95,
        connection_limit: Optional[int] = None,
        **kwargs,
    ):
        """
        Initialize TCS algorithm.

        Args:
            distance_method: Method for calculating distances (should be hamming)
            confidence: Confidence level for parsimony criterion (default 0.95)
            connection_limit: Maximum connection distance (auto-calculated if None)
            **kwargs: Additional parameters
        """
        super().__init__(distance_method, **kwargs)
        self.confidence = confidence
        self.connection_limit = connection_limit

        if distance_method not in ('hamming',):
            # TCS is designed for discrete mutation counts
            import warnings

            warnings.warn(
                f"TCS is designed for hamming distance. Using '{distance_method}' "
                'may not produce theoretically correct results.', stacklevel=2
            )

    def construct_network(
        self, alignment: Alignment, distance_matrix: Optional[DistanceMatrix] = None
    ) -> HaplotypeNetwork:
        """
        Construct TCS network from sequence alignment.

        Args:
            alignment: Multiple sequence alignment
            distance_matrix: Optional pre-computed distance matrix

        Returns:
            Haplotype network constructed using statistical parsimony
        """
        # Identify unique haplotypes
        haplotypes = identify_haplotypes(alignment)

        if len(haplotypes) == 0:
            return HaplotypeNetwork()

        if len(haplotypes) == 1:
            network = HaplotypeNetwork()
            network.add_haplotype(haplotypes[0])
            return network

        # Calculate distances between haplotypes
        haplotype_dist_matrix = self._calculate_haplotype_distances(haplotypes)
        self._distance_matrix = haplotype_dist_matrix

        # Calculate connection limit if not provided
        if self.connection_limit is None:
            self.connection_limit = self._calculate_connection_limit(
                alignment.length, len(haplotypes)
            )

        # Build network using parsimony criterion
        edges = self._build_parsimony_network(haplotypes, haplotype_dist_matrix)

        # Construct network
        network = self._build_network_from_edges(haplotypes, edges)

        return network

    def _calculate_connection_limit(
        self, sequence_length: int, num_haplotypes: int
    ) -> int:
        """
        Calculate maximum parsimony connection limit.

        Uses the formula from Templeton et al. (1992) to estimate the
        maximum number of mutational differences that can be explained
        by parsimony at the specified confidence level.

        Args:
            sequence_length: Length of aligned sequences
            num_haplotypes: Number of unique haplotypes

        Returns:
            Maximum connection distance for parsimony criterion
        """
        # Calculate probability of finding n or more differences by chance
        # Using cumulative Poisson distribution

        # Expected number of differences under infinite sites model
        # Adjusted for sample size
        n = num_haplotypes

        # Calculate connection limit using Poisson distribution
        # P(X >= k) >= confidence
        # We find the largest k where P(X <= k-1) < (1 - confidence)

        connection_limit = 1
        cumulative_prob = 0.0

        for k in range(1, sequence_length + 1):
            # Calculate probability using Poisson approximation
            # For each possible number of mutations

            # Probability that k or fewer mutations occurred by chance
            prob_k = self._poisson_probability(k, sequence_length, n)
            cumulative_prob += prob_k

            if cumulative_prob >= (1 - self.confidence):
                connection_limit = k
                break

        return max(1, connection_limit)

    def _poisson_probability(self, k: int, seq_length: int, sample_size: int) -> float:
        """
        Calculate probability using Poisson distribution.

        Args:
            k: Number of mutations
            seq_length: Sequence length
            sample_size: Number of sequences

        Returns:
            Probability
        """
        # Simple Poisson probability: P(X = k) = (λ^k * e^(-λ)) / k!
        # where λ is the expected number of mutations

        # Expected mutations based on coalescent theory
        # Simplified: λ ≈ 2 * sample_size for neutral markers
        lambda_param = 2.0 * math.log(sample_size) if sample_size > 1 else 1.0

        # Poisson probability
        try:
            prob = (lambda_param**k) * math.exp(-lambda_param) / math.factorial(k)
        except (OverflowError, ValueError):
            # For large k, use approximation
            prob = 0.0

        return prob

    def _build_parsimony_network(
        self, haplotypes: List, distance_matrix: DistanceMatrix
    ) -> List[Tuple[str, str, float]]:
        """
        Build network using statistical parsimony criterion.

        Connects haplotypes in order of increasing distance, up to
        the connection limit. More frequent haplotypes are connected first.

        Args:
            haplotypes: List of Haplotype objects
            distance_matrix: Distance matrix

        Returns:
            List of edges (id1, id2, distance)
        """
        # Sort haplotypes by frequency (most frequent first)
        sorted_haps = sorted(haplotypes, key=lambda h: h.frequency, reverse=True)
        hap_ids = [h.id for h in sorted_haps]

        # Track which haplotypes are in the network
        in_network = {hap_ids[0]}  # Start with most frequent
        not_in_network = set(hap_ids[1:])

        edges = []

        # For each distance level up to connection limit
        for dist_level in range(1, self.connection_limit + 1):
            if not not_in_network:
                break

            # Find all connections at this distance level
            connections = []

            for hap_in in in_network:
                for hap_out in not_in_network:
                    dist = distance_matrix.get_distance(hap_in, hap_out)

                    if abs(dist - dist_level) < 0.5:  # Allow for rounding
                        connections.append((hap_in, hap_out, dist))

            # Add connections, preferring those to more frequent haplotypes
            # Sort by frequency of the haplotype being added
            hap_freq = {h.id: h.frequency for h in haplotypes}
            connections.sort(key=lambda c: hap_freq[c[1]], reverse=True)

            for hap_in, hap_out, dist in connections:
                if hap_out in not_in_network:
                    edges.append((hap_in, hap_out, dist))
                    in_network.add(hap_out)
                    not_in_network.remove(hap_out)

        # Handle any remaining unconnected haplotypes
        # These form separate networks (sub-networks)
        if not_in_network:
            # Connect remaining haplotypes among themselves
            remaining = list(not_in_network)
            for i, hap1 in enumerate(remaining):
                for hap2 in remaining[i + 1 :]:
                    dist = distance_matrix.get_distance(hap1, hap2)
                    if dist <= self.connection_limit:
                        edges.append((hap1, hap2, dist))

        return edges

    def _build_network_from_edges(
        self, haplotypes: List, edges: List[Tuple[str, str, float]]
    ) -> HaplotypeNetwork:
        """
        Build HaplotypeNetwork from haplotypes and edges.

        Args:
            haplotypes: List of Haplotype objects
            edges: List of edges (id1, id2, distance)

        Returns:
            Constructed haplotype network
        """
        network = HaplotypeNetwork()

        # Add all haplotypes
        for haplotype in haplotypes:
            network.add_haplotype(haplotype)

        # Add edges
        for id1, id2, dist in edges:
            network.add_edge(id1, id2, distance=int(round(dist)))

        return network

    def _calculate_haplotype_distances(self, haplotypes: List) -> DistanceMatrix:
        """
        Calculate pairwise distances between haplotypes.

        Args:
            haplotypes: List of Haplotype objects

        Returns:
            DistanceMatrix with distances between haplotypes
        """
        import numpy as np

        n = len(haplotypes)
        labels = [h.id for h in haplotypes]
        matrix = np.zeros((n, n))

        for i in range(n):
            for j in range(i + 1, n):
                dist = hamming_distance(
                    haplotypes[i].sequence,
                    haplotypes[j].sequence,
                    ignore_gaps=self.params.get('ignore_gaps', True),
                )
                matrix[i, j] = matrix[j, i] = dist

        return DistanceMatrix(labels, matrix)

    def get_parameters(self) -> dict:
        """Get algorithm parameters."""
        params = super().get_parameters()
        params['confidence'] = self.confidence
        params['connection_limit'] = self.connection_limit
        return params

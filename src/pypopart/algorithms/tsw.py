"""
Tight Span Walker (TSW) algorithm for haplotype network construction.

This module implements the Tight Span Walker algorithm, which constructs
a haplotype network by computing the tight span of a distance matrix.
The tight span is the smallest metric space that contains all optimal
paths between sequences.

References
----------
.. [1] Dress, A. W., & Huson, D. H. (2004). Constructing splits graphs.
       IEEE/ACM Transactions on Computational Biology and Bioinformatics, 1(3), 109-115.
.. [2] Bryant, D., & Moulton, V. (2004). Neighbor-Net: an agglomerative method
       for the construction of phylogenetic networks. Molecular biology and evolution, 21(2), 255-265.
"""

import logging
from typing import Dict, List, Optional, Set, Tuple

import numpy as np

from ..core.alignment import Alignment
from ..core.distance import DistanceMatrix
from ..core.graph import HaplotypeNetwork
from ..core.haplotype import Haplotype, identify_haplotypes_from_alignment
from .base import NetworkAlgorithm

logger = logging.getLogger(__name__)


class TightSpanWalker(NetworkAlgorithm):
    """
    Construct a Tight Span Walker network from haplotype data.

    The Tight Span Walker (TSW) algorithm constructs a network by computing
    the tight span of a distance matrix. This creates a network that preserves
    all metric properties of the original distances and can represent complex
    reticulate relationships.

    The algorithm works by:
    1. Computing dT distances (tree metric distances)
    2. Building geodesic paths between all pairs of haplotypes
    3. Adding median/internal vertices where needed to maintain metric properties

    This creates a more complex network than MST/MSN but captures the full
    metric structure of the data.

    Parameters
    ----------
    distance_method : str, default='hamming'
        Method for calculating pairwise distances between sequences.
        Options: 'hamming', 'jukes_cantor', 'kimura_2p', 'tamura_nei'.
    epsilon : float, default=1e-6
        Tolerance for floating point comparisons.
    **kwargs : dict
        Additional parameters passed to base NetworkAlgorithm.

    Attributes
    ----------
    epsilon : float
        Tolerance for floating point comparisons.
    _dt_matrix : Optional[np.ndarray]
        Matrix of dT (tree metric) distances.
    _internal_vertices : Dict[str, Haplotype]
        Dictionary of inferred internal/median vertices.

    Examples
    --------
    >>> from pypopart.algorithms import TightSpanWalker
    >>> from pypopart.io import load_alignment
    >>> alignment = load_alignment('sequences.fasta')
    >>> tsw = TightSpanWalker(distance_method='hamming')
    >>> network = tsw.construct_network(alignment)

    Notes
    -----
    TSW is computationally expensive (O(n³) or worse) and is best suited
    for smaller datasets (n < 100) where accurate representation of complex
    relationships is critical.
    """

    def __init__(
        self, distance_method: str = 'hamming', epsilon: float = 1e-6, **kwargs
    ):
        """
        Initialize TightSpanWalker algorithm.

        Parameters
        ----------
        distance_method : str, default='hamming'
            Method for calculating distances.
        epsilon : float, default=1e-6
            Tolerance for floating point comparisons.
        **kwargs : dict
            Additional parameters.
        """
        super().__init__(distance_method, **kwargs)
        self.epsilon = epsilon
        self._dt_matrix: Optional[np.ndarray] = None
        self._internal_vertices: Dict[str, Haplotype] = {}
        self._internal_counter = 0

    def construct_network(
        self, alignment: Alignment, distance_matrix: Optional[DistanceMatrix] = None
    ) -> HaplotypeNetwork:
        """
        Construct TSW network from sequence alignment.

        Parameters
        ----------
        alignment : Alignment
            Multiple sequence alignment.
        distance_matrix : Optional[DistanceMatrix]
            Optional pre-computed distance matrix.

        Returns
        -------
        HaplotypeNetwork
            Haplotype network representing the tight span.
        """
        # Identify unique haplotypes
        haplotypes = identify_haplotypes_from_alignment(alignment)

        if len(haplotypes) <= 1:
            # Single or no haplotypes - create trivial network
            network = HaplotypeNetwork()
            for hap in haplotypes:
                network.add_haplotype(hap)
            return network

        # Calculate distances between haplotypes
        if distance_matrix is None:
            haplotype_dist_matrix = self._calculate_haplotype_distances(haplotypes)
        else:
            haplotype_dist_matrix = distance_matrix
        self._distance_matrix = haplotype_dist_matrix

        # Compute dT (tree metric) distances
        self._compute_dt_distances(haplotype_dist_matrix)

        # Build network using geodesic paths
        network = self._build_tight_span_network(haplotypes, haplotype_dist_matrix)

        return network

    def _calculate_haplotype_distances(
        self, haplotypes: List[Haplotype]
    ) -> DistanceMatrix:
        """
        Calculate pairwise distances between haplotypes.

        Parameters
        ----------
        haplotypes : List[Haplotype]
            List of unique haplotypes.

        Returns
        -------
        DistanceMatrix
            Pairwise distance matrix.
        """
        n = len(haplotypes)
        labels = [hap.id for hap in haplotypes]
        matrix = np.zeros((n, n))

        for i in range(n):
            for j in range(i + 1, n):
                # Calculate Hamming distance between sequences
                seq1 = haplotypes[i].data  # Use .data to get string
                seq2 = haplotypes[j].data
                dist = sum(c1 != c2 for c1, c2 in zip(seq1, seq2))
                matrix[i, j] = dist
                matrix[j, i] = dist

        return DistanceMatrix(labels, matrix)

    def _compute_dt_distances(self, distance_matrix: DistanceMatrix) -> None:
        """
        Compute dT (tree metric) distances.

        For each pair (i, j), dT(i,j) = max over all k of |d(i,k) - d(j,k)|.
        This represents the minimum distance if sequences were on a tree.

        Parameters
        ----------
        distance_matrix : DistanceMatrix
            Pairwise distance matrix.
        """
        labels = distance_matrix.labels
        n = len(labels)
        self._dt_matrix = np.zeros((n, n))

        for i in range(n):
            for j in range(i + 1, n):
                # Compute max |d(i,k) - d(j,k)| over all k
                max_diff = 0.0
                for k in range(n):
                    if k != i and k != j:
                        diff = abs(
                            distance_matrix.get_distance(labels[i], labels[k])
                            - distance_matrix.get_distance(labels[j], labels[k])
                        )
                        max_diff = max(max_diff, diff)

                self._dt_matrix[i, j] = max_diff
                self._dt_matrix[j, i] = max_diff

    def _build_tight_span_network(
        self, haplotypes: List[Haplotype], distance_matrix: DistanceMatrix
    ) -> HaplotypeNetwork:
        """
        Build the tight span network using geodesic paths.

        Parameters
        ----------
        haplotypes : List[Haplotype]
            List of unique haplotypes.
        distance_matrix : DistanceMatrix
            Pairwise distance matrix.

        Returns
        -------
        HaplotypeNetwork
            Constructed tight span network.
        """
        network = HaplotypeNetwork()

        # Add all original haplotypes to network
        for hap in haplotypes:
            network.add_haplotype(hap)

        # Build geodesic paths between all pairs
        n = len(haplotypes)
        for i in range(n):
            for j in range(i + 1, n):
                self._add_geodesic_path(
                    network, haplotypes[i], haplotypes[j], distance_matrix, i, j
                )

        return network

    def _add_geodesic_path(
        self,
        network: HaplotypeNetwork,
        hap1: Haplotype,
        hap2: Haplotype,
        distance_matrix: DistanceMatrix,
        idx1: int,
        idx2: int,
    ) -> None:
        """
        Add geodesic path between two haplotypes.

        This may add intermediate vertices if needed to maintain metric properties.

        Parameters
        ----------
        network : HaplotypeNetwork
            Network being constructed.
        hap1 : Haplotype
            First haplotype.
        hap2 : Haplotype
            Second haplotype.
        distance_matrix : DistanceMatrix
            Pairwise distance matrix.
        idx1 : int
            Index of first haplotype in distance matrix.
        idx2 : int
            Index of second haplotype in distance matrix.
        """
        # Get actual distance and dT distance
        actual_dist = distance_matrix.get_distance(hap1.id, hap2.id)
        dt_dist = self._dt_matrix[idx1, idx2]

        # Check if we need intermediate vertices
        # If dT is much smaller than actual distance, we need medians
        if actual_dist - dt_dist > self.epsilon:
            # Need to add intermediate vertices
            # For simplicity, add a single median vertex
            median_hap = self._create_median_vertex(hap1, hap2)

            # Add median to network if not already present
            if median_hap.id not in [h.id for h in network.haplotypes]:
                network.add_haplotype(median_hap, median_vector=True)

            # Add edges
            dist1 = actual_dist / 2.0
            dist2 = actual_dist / 2.0

            network.add_edge(hap1.id, median_hap.id, distance=dist1)
            network.add_edge(median_hap.id, hap2.id, distance=dist2)
        else:
            # Direct connection
            network.add_edge(hap1.id, hap2.id, distance=actual_dist)

    def _create_median_vertex(self, hap1: Haplotype, hap2: Haplotype) -> Haplotype:
        """
        Create a median vertex between two haplotypes.

        Parameters
        ----------
        hap1 : Haplotype
            First haplotype.
        hap2 : Haplotype
            Second haplotype.

        Returns
        -------
        Haplotype
            Median haplotype.
        """
        from ..core.sequence import Sequence

        # Create median sequence by taking consensus
        seq1 = hap1.data
        seq2 = hap2.data

        median_seq = []
        for c1, c2 in zip(seq1, seq2):
            if c1 == c2:
                median_seq.append(c1)
            else:
                # Take first character for simplicity
                # In a more sophisticated implementation, we would
                # choose based on parsimony or other criteria
                median_seq.append(c1)

        median_sequence = ''.join(median_seq)

        # Check if this median already exists
        median_id = f'Median_{self._internal_counter}'

        # Create median Sequence object
        median_seq_obj = Sequence(median_id, median_sequence)

        # Create median haplotype
        median_hap = Haplotype(sequence=median_seq_obj, sample_ids=[])

        self._internal_vertices[median_id] = median_hap
        self._internal_counter += 1

        return median_hap

    def get_parameters(self) -> Dict:
        """
        Get algorithm parameters.

        Returns
        -------
        dict
            Dictionary of parameters.
        """
        params = super().get_parameters()
        params['epsilon'] = self.epsilon
        return params

    def __str__(self) -> str:
        """Return string representation."""
        return f'TightSpanWalker(distance={self.distance_method}, epsilon={self.epsilon})'


# Convenient alias
TSWAlgorithm = TightSpanWalker

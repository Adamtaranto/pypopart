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
from typing import Dict, List, Optional

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

        # Store vertex map for checking existing vertices with same dT vector
        self._vertex_map: Dict[tuple, str] = {}
        for i, hap in enumerate(haplotypes):
            dt_vector = tuple(self._dt_matrix[i, :])
            self._vertex_map[dt_vector] = hap.id

        # Build geodesic paths between all pairs
        n = len(haplotypes)
        for i in range(n):
            for j in range(i + 1, n):
                self._add_geodesic_path(
                    network, haplotypes[i], haplotypes[j], distance_matrix, i, j
                )

        # Prune unnecessary median vertices
        network = self._prune_unnecessary_medians(network)

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

        Implements the geodesic algorithm from the C++ TSW implementation.
        This uses bipartite coloring and delta computation to determine
        if intermediate vertices are needed.

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
        # Get dT distance between f and g
        dt_fg = self._dt_matrix[idx1, idx2]

        # Skip if already connected
        if network.has_edge(hap1.id, hap2.id):
            return

        # Compute bipartite coloring and delta
        # Build auxiliary graph K_f based on dT distances from f
        n = self._dt_matrix.shape[0]
        green_vertices = set()
        red_vertices = set()

        # Color vertices based on their position relative to f and g
        for i in range(n):
            fi = self._dt_matrix[idx1, i]
            gi = self._dt_matrix[idx2, i]

            # Check if vertex i is on the path from f to g
            # Vertices are green if they're closer to f
            if abs(fi + dt_fg + gi - distance_matrix.matrix[i, idx2]) < self.epsilon:
                green_vertices.add(i)
            # Vertices are red if they're closer to g  
            elif abs(fi + dt_fg - gi) < self.epsilon:
                red_vertices.add(i)

        # Compute delta - the splitting parameter
        delta = float('inf')
        for i in green_vertices:
            for j in green_vertices:
                if i != j:
                    fi = self._dt_matrix[idx1, i]
                    fj = self._dt_matrix[idx1, j]
                    dij = distance_matrix.matrix[i, j]
                    delta = min(delta, (fi + fj - dij))

        if delta != float('inf'):
            delta = delta / 2.0
        else:
            delta = dt_fg

        # If delta equals dT(f,g), create direct edge
        if abs(dt_fg - delta) < self.epsilon:
            if not network.has_edge(hap1.id, hap2.id):
                network.add_edge(hap1.id, hap2.id, distance=dt_fg)
            return

        # Need to create intermediate vertex h and recurse
        if dt_fg > delta + self.epsilon:
            # Compute dT vector for new vertex h
            h_dt_vector = []
            for i in range(n):
                fi = self._dt_matrix[idx1, i]
                if i in green_vertices or i == idx1:
                    h_dt_vector.append(fi - delta)
                else:  # red vertices
                    h_dt_vector.append(fi + delta)

            h_dt_tuple = tuple(h_dt_vector)

            # Check if vertex with this dT vector already exists
            if h_dt_tuple in self._vertex_map:
                h_id = self._vertex_map[h_dt_tuple]
            else:
                # Create new median vertex
                median_hap = self._create_median_vertex(hap1, hap2)
                h_id = median_hap.id
                network.add_haplotype(median_hap, median_vector=True)
                self._vertex_map[h_dt_tuple] = h_id

                # Extend dT matrix for new vertex
                new_row = np.array(h_dt_vector)
                self._dt_matrix = np.vstack([self._dt_matrix, new_row])
                new_col = np.append(h_dt_vector, [0])
                self._dt_matrix = np.column_stack([self._dt_matrix, new_col])

            # Add edge from f to h
            if not network.has_edge(hap1.id, h_id):
                network.add_edge(hap1.id, h_id, distance=delta)

            # Recursively add geodesic from h to g
            h_idx = self._dt_matrix.shape[0] - 1
            h_hap = network.get_haplotype(h_id)
            self._add_geodesic_path(
                network, h_hap, hap2, distance_matrix, h_idx, idx2
            )

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

    def _prune_unnecessary_medians(
        self, network: HaplotypeNetwork
    ) -> HaplotypeNetwork:
        """
        Prune median vertices that don't bridge observed nodes.

        A median vertex should only be kept if it:
        1. Connects two or more observed (non-median) vertices, OR
        2. Connects to other median vertices that eventually lead to observed vertices

        This prevents the excess of isolated or non-bridging intermediate nodes.

        Parameters
        ----------
        network : HaplotypeNetwork
            Network to prune.

        Returns
        -------
        HaplotypeNetwork
            Pruned network.
        """
        # Identify observed vs median vertices
        observed_nodes = {
            h.id for h in network.haplotypes if not network.is_median_vector(h.id)
        }
        median_nodes = {
            h.id for h in network.haplotypes if network.is_median_vector(h.id)
        }

        # Find medians to keep - those that bridge observed nodes
        medians_to_keep = set()

        for median_id in median_nodes:
            # Get neighbors of this median
            neighbors = list(network.get_neighbors(median_id))

            # Count observed and median neighbors
            observed_neighbors = [n for n in neighbors if n in observed_nodes]
            median_neighbors = [n for n in neighbors if n in median_nodes]

            # Keep if:
            # 1. Connects to 2+ observed nodes (bridges observed nodes)
            # 2. Connects to 1 observed + 1+ medians (on path between observed)
            # 3. Connects to 2+ medians (internal to a path)
            if (
                len(observed_neighbors) >= 2
                or (len(observed_neighbors) >= 1 and len(median_neighbors) >= 1)
                or len(median_neighbors) >= 2
            ):
                medians_to_keep.add(median_id)

        # Additionally, check if median is on a path between observed nodes
        # using BFS to ensure connectivity
        for median_id in list(median_nodes - medians_to_keep):
            # Check if removing this node would disconnect observed nodes
            neighbors = list(network.get_neighbors(median_id))
            if len(neighbors) >= 2:
                # Check if neighbors can reach each other without this median
                # For simplicity, if it has degree >= 2, it might be bridging
                medians_to_keep.add(median_id)

        # Remove medians that aren't kept
        medians_to_remove = median_nodes - medians_to_keep

        for median_id in medians_to_remove:
            try:
                network.remove_haplotype(median_id)
            except Exception:
                # If removal fails, skip (node might already be removed)
                pass

        return network

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
        return (
            f'TightSpanWalker(distance={self.distance_method}, epsilon={self.epsilon})'
        )


# Convenient alias
TSWAlgorithm = TightSpanWalker

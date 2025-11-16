"""
Tight Span Walker algorithm for haplotype network construction.

This module implements the Tight Span Walker (TSW) algorithm, which constructs
haplotype networks by computing the tight span of a distance matrix. The tight
span is a geometric structure that captures the optimal metric space representation
of distance data.

The algorithm creates reticulate networks that can represent evolutionary scenarios
with recombination/hybridization by ensuring all distances in the network respect
the original distance matrix through addition of hypothetical ancestral/intermediate
sequences where needed.

References
----------
.. [1] Dress, A. W., Huber, K. T., Koolen, J., Moulton, V., & Spillner, A. (2012).
       Basic Phylogenetic Combinatorics. Cambridge University Press.
.. [2] Bandelt, H. J., Forster, P., & Röhl, A. (1999). Median-joining networks for
       inferring intraspecific phylogenies. Molecular Biology and Evolution, 16(1), 37-48.
"""

import sys
from collections import deque
from typing import Dict, List, Optional, Set, Tuple

import numpy as np

from ..core.alignment import Alignment
from ..core.distance import DistanceMatrix
from ..core.graph import HaplotypeNetwork
from ..core.haplotype import Haplotype, identify_haplotypes_from_alignment
from ..core.sequence import Sequence
from .base import NetworkAlgorithm


class TightSpanWalker(NetworkAlgorithm):
    """
    Construct a haplotype network using the Tight Span Walker algorithm.

    The Tight Span Walker algorithm computes the tight span of a distance matrix,
    which is the smallest convex set containing all optimal paths between sequences.
    This creates a reticulate network that exactly represents the metric properties
    of the distance data.

    The algorithm works by:
    1. Computing dT distances (tree metric) for all sequence pairs
    2. Building geodesic paths between all pairs of sequences
    3. Adding median vertices where needed to maintain metric properties
    4. Deduplicating vertices based on their dT vectors

    Parameters
    ----------
    distance_method : str, default='hamming'
        Method for calculating pairwise distances between sequences.
        Options: 'hamming', 'jukes_cantor', 'kimura_2p', 'tamura_nei'.
    **kwargs : dict
        Additional parameters passed to base NetworkAlgorithm.

    Attributes
    ----------
    _dT_matrix : np.ndarray
        Matrix of dT (tree metric) distances between vertices
    _vertex_map : Dict[tuple, str]
        Maps dT vectors to vertex IDs for deduplication
    _median_counter : int
        Counter for generating unique median vertex IDs

    Examples
    --------
    >>> from pypopart.algorithms import TightSpanWalker
    >>> from pypopart.io import load_alignment
    >>>
    >>> # Load alignment
    >>> alignment = load_alignment('sequences.fasta')
    >>>
    >>> # Construct Tight Span Walker network
    >>> tsw = TightSpanWalker()
    >>> network = tsw.build_network(alignment)

    Notes
    -----
    The Tight Span Walker algorithm is computationally expensive (O(n³) or worse)
    and is best suited for small to medium datasets (< 100 sequences). For larger
    datasets, consider using faster approximation methods like MSN or TCS.

    Unlike MST which produces a tree, TSW can create reticulate networks with
    cycles, making it suitable for representing complex evolutionary relationships
    involving recombination or hybridization.

    See Also
    --------
    MinimumSpanningNetwork : Faster approximation for network construction
    TCS : Statistical parsimony network construction
    MedianJoiningNetwork : Alternative reticulation method
    """

    def __init__(self, distance_method: str = 'hamming', **kwargs):
        """
        Initialize Tight Span Walker algorithm.

        Parameters
        ----------
        distance_method : str, default='hamming'
            Method for calculating distances.
        **kwargs : dict
            Additional algorithm parameters.
        """
        super().__init__(distance_method=distance_method, **kwargs)
        self._dT_matrix: Optional[np.ndarray] = None
        self._vertex_map: Dict[tuple, str] = {}
        self._median_counter: int = 0
        # Floating point comparison tolerance
        self._epsilon = np.finfo(float).eps

    def construct_network(
        self, alignment: Alignment, distance_matrix: Optional[DistanceMatrix] = None
    ) -> HaplotypeNetwork:
        """
        Construct haplotype network using Tight Span Walker algorithm.

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
        n_samples = len(haplotypes)

        # Compute dT matrix (tree metric distances)
        self._dT_matrix = self._compute_dT(dist_array, n_samples)

        # Initialize vertex map with observed haplotypes
        self._vertex_map = {}
        for i, haplotype in enumerate(haplotypes):
            dT_vector = tuple(self._dT_matrix[i, :n_samples])
            self._vertex_map[dT_vector] = haplotype.id

        # Build geodesic paths between all pairs of sequences
        self._median_counter = 0
        haplotype_ids = [h.id for h in haplotypes]

        for i in range(len(haplotype_ids)):
            for j in range(i):
                self._geodesic(
                    network,
                    haplotype_ids[i],
                    haplotype_ids[j],
                    dist_array,
                    n_samples,
                )

        return network

    def _compute_dT(self, dist_matrix: np.ndarray, n_samples: int) -> np.ndarray:
        """
        Compute dT (tree metric) distances.

        For each pair of sequences (i, j), dT(i,j) = max over all k of |d(i,k) - d(j,k)|.
        This represents the minimum distance if sequences were on a tree.

        Parameters
        ----------
        dist_matrix : np.ndarray
            Pairwise distance matrix.
        n_samples : int
            Number of sample sequences.

        Returns
        -------
        np.ndarray
            Matrix of dT distances (will grow as median vertices are added).
        """
        # Initialize dT matrix (will expand as we add median vertices)
        dT = np.zeros((n_samples, n_samples), dtype=float)

        for i in range(n_samples):
            for j in range(i):
                # dT(i,j) = max_k |d(i,k) - d(j,k)|
                dT_ij = np.max(np.abs(dist_matrix[i, :n_samples] - dist_matrix[j, :n_samples]))
                dT[i, j] = dT[j, i] = dT_ij

        return dT

    def _about_equal(self, a: float, b: float) -> bool:
        """
        Compare two floats for approximate equality.

        Uses relative tolerance based on the magnitude of the values.

        Parameters
        ----------
        a : float
            First value.
        b : float
            Second value.

        Returns
        -------
        bool
            True if values are approximately equal.
        """
        # Use numpy's isclose with reasonable tolerances
        # rtol=1e-9, atol=1e-12 (more permissive than machine epsilon)
        return bool(np.isclose(a, b, rtol=1e-9, atol=1e-12))

    def _get_dT(self, i: int, j: int) -> float:
        """
        Get dT distance between two vertices.

        Parameters
        ----------
        i : int
            Index of first vertex.
        j : int
            Index of second vertex.

        Returns
        -------
        float
            dT distance.
        """
        if i >= self._dT_matrix.shape[0] or j >= self._dT_matrix.shape[1]:
            raise ValueError(f"Invalid index for dT distance: ({i}, {j})")
        return self._dT_matrix[i, j]

    def _geodesic(
        self,
        network: HaplotypeNetwork,
        f_id: str,
        g_id: str,
        dist_matrix: np.ndarray,
        n_samples: int,
    ) -> None:
        """
        Compute geodesic path between two vertices f and g.

        This is the core recursive function that builds the tight span by adding
        median vertices where needed to maintain metric properties.

        Parameters
        ----------
        network : HaplotypeNetwork
            Network being constructed.
        f_id : str
            ID of first vertex.
        g_id : str
            ID of second vertex.
        dist_matrix : np.ndarray
            Distance matrix (for original samples).
        n_samples : int
            Number of original sample sequences.

        Returns
        -------
        None
            Network is modified in place.
        """
        # Get indices of f and g in the vertex list
        all_haplotype_ids = network.nodes
        f_idx = all_haplotype_ids.index(f_id)
        g_idx = all_haplotype_ids.index(g_id)

        # Get dT distance between f and g
        dT_fg = self._get_dT(f_idx, g_idx)

        # Check if already connected with correct distance
        if network.has_edge(f_id, g_id):
            existing_dist = network.get_edge_distance(f_id, g_id)
            if self._about_equal(existing_dist, dT_fg):
                return

        # Compute bipartite coloring to determine which vertices are closer to f vs g
        green_set, red_set = self._compute_bipartite_coloring(
            f_idx, g_idx, all_haplotype_ids, n_samples, dist_matrix
        )

        # Compute delta (splitting parameter)
        delta = self._compute_delta(green_set, f_idx, dist_matrix, n_samples)

        # Case 1: f and g can be connected directly
        if self._about_equal(dT_fg, delta):
            if not network.has_edge(f_id, g_id):
                network.add_edge(f_id, g_id, distance=delta)
            return

        # Case 2: Need to create intermediate vertex h
        if dT_fg > delta:
            # Compute dT vector for new median vertex
            new_dT_vector = self._compute_new_vertex_dT(
                f_idx, delta, green_set, red_set, all_haplotype_ids
            )

            # Check if this vertex already exists
            dT_tuple = tuple(new_dT_vector[: len(all_haplotype_ids)])
            if dT_tuple in self._vertex_map:
                h_id = self._vertex_map[dT_tuple]
            else:
                # Create new median vertex
                h_id = f"Median_{self._median_counter}"
                self._median_counter += 1

                # Create empty sequence for median (actual sequence reconstruction would go here)
                median_seq = Sequence(id=h_id, data="")
                median_haplotype = Haplotype(sequence=median_seq, sample_ids=[])

                network.add_haplotype(median_haplotype, median_vector=True)

                # Update dT matrix to include new vertex
                self._expand_dT_matrix(new_dT_vector, n_samples)

                # Store in vertex map
                self._vertex_map[dT_tuple] = h_id

            # Connect f to h
            if not network.has_edge(f_id, h_id):
                network.add_edge(f_id, h_id, distance=delta)

            # Recursively connect h to g
            self._geodesic(network, h_id, g_id, dist_matrix, n_samples)

        # Case 3: Would create negative edge length - shouldn't happen
        elif dT_fg < delta:
            raise ValueError(
                f"Negative edge length detected between {f_id} and {g_id}. "
                f"dT={dT_fg}, delta={delta}"
            )

    def _compute_bipartite_coloring(
        self,
        f_idx: int,
        g_idx: int,
        all_ids: List[str],
        n_samples: int,
        dist_matrix: np.ndarray,
    ) -> Tuple[Set[int], Set[int]]:
        """
        Compute bipartite coloring of vertices based on proximity to f and g.

        Green vertices are closer to f, red vertices are closer to g.

        Parameters
        ----------
        f_idx : int
            Index of vertex f.
        g_idx : int
            Index of vertex g.
        all_ids : List[str]
            List of all vertex IDs.
        n_samples : int
            Number of original samples.
        dist_matrix : np.ndarray
            Distance matrix.

        Returns
        -------
        Tuple[Set[int], Set[int]]
            Sets of green and red vertex indices.
        """
        green = set()
        red = set()

        dT_fg = self._get_dT(f_idx, g_idx)

        for i in range(len(all_ids)):
            if i == f_idx or i == g_idx:
                continue

            dT_fi = self._get_dT(f_idx, i)
            dT_gi = self._get_dT(g_idx, i)

            # Check which side of the split this vertex is on
            # Vertex i is on the f side (green) if dT(f,i) + dT(f,g) + dT(g,i) equals path through
            if dT_fi < dT_gi:
                # Check if i is on optimal path from f to g
                if self._about_equal(dT_fi + dT_fg + dT_gi, dT_fg + 2 * dT_gi):
                    green.add(i)
            elif dT_gi < dT_fi:
                if self._about_equal(dT_gi + dT_fg + dT_fi, dT_fg + 2 * dT_fi):
                    red.add(i)

        # f is always green, g is always red (conceptually)
        green.add(f_idx)
        red.add(g_idx)

        return green, red

    def _compute_delta(
        self, green_set: Set[int], f_idx: int, dist_matrix: np.ndarray, n_samples: int
    ) -> float:
        """
        Compute delta (splitting parameter) for geodesic computation.

        Delta is half the minimum over all pairs of green vertices of
        (dT(f,i) + dT(f,j) - d(i,j))

        Parameters
        ----------
        green_set : Set[int]
            Set of green vertex indices.
        f_idx : int
            Index of vertex f.
        dist_matrix : np.ndarray
            Distance matrix.
        n_samples : int
            Number of original samples.

        Returns
        -------
        float
            Delta value.
        """
        min_delta = float('inf')

        green_list = list(green_set)
        for i in green_list:
            dT_fi = self._get_dT(f_idx, i)
            for j in green_list:
                if i == j:
                    continue
                dT_fj = self._get_dT(f_idx, j)

                # Get distance d(i,j)
                if i < n_samples and j < n_samples:
                    d_ij = dist_matrix[i, j]
                else:
                    # If one is a median, use dT distance
                    d_ij = self._get_dT(i, j)

                delta_ij = dT_fi + dT_fj - d_ij
                min_delta = min(min_delta, delta_ij)

        return min_delta / 2.0

    def _compute_new_vertex_dT(
        self,
        f_idx: int,
        delta: float,
        green_set: Set[int],
        red_set: Set[int],
        all_ids: List[str],
    ) -> np.ndarray:
        """
        Compute dT vector for a new median vertex.

        For vertices on the green side (closer to f), dT decreases by delta.
        For vertices on the red side (closer to g), dT increases by delta.

        Parameters
        ----------
        f_idx : int
            Index of vertex f.
        delta : float
            Splitting parameter.
        green_set : Set[int]
            Green vertex indices.
        red_set : Set[int]
            Red vertex indices.
        all_ids : List[str]
            All vertex IDs.

        Returns
        -------
        np.ndarray
            dT vector for new vertex.
        """
        new_dT = np.zeros(len(all_ids) + 1)

        for i in range(len(all_ids)):
            dT_fi = self._get_dT(f_idx, i)

            if i in green_set:
                new_dT[i] = dT_fi - delta
            elif i in red_set:
                new_dT[i] = dT_fi + delta
            else:
                # Shouldn't happen - all vertices should be colored
                new_dT[i] = dT_fi

        # dT to itself is 0
        new_dT[len(all_ids)] = 0.0

        return new_dT

    def _expand_dT_matrix(self, new_dT_vector: np.ndarray, n_samples: int) -> None:
        """
        Expand dT matrix to include a new median vertex.

        Parameters
        ----------
        new_dT_vector : np.ndarray
            dT distances from all vertices to new vertex.
        n_samples : int
            Number of original samples.

        Returns
        -------
        None
            Updates self._dT_matrix in place.
        """
        n_current = self._dT_matrix.shape[0]
        n_new = n_current + 1

        # Create expanded matrix
        new_matrix = np.zeros((n_new, n_new))

        # Copy existing values
        new_matrix[:n_current, :n_current] = self._dT_matrix

        # Add new row and column
        new_matrix[n_current, :] = new_dT_vector[:n_new]
        new_matrix[:, n_current] = new_dT_vector[:n_new]

        self._dT_matrix = new_matrix

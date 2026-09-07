"""
Tight Span Walker (TSW) algorithm, ported from PopART's TightSpanWalker.

Constructs a metric-preserving haplotype network from the tight span of
the distance matrix (Dress et al. 2012). Like PopART, it operates on
condensed site patterns with site-weighted distances, computes the
tight-span metric dT(i,j) = max_k |d(i,k) - d(j,k)|, and walks a
geodesic between every pair of sample vertices, materialising internal
(median) vertices along the way. Internal vertices carry no sequence,
only a dT vector, and are memoised so shared vertices are reused.
"""

import math
from typing import Dict, List, Optional, Tuple

import networkx as nx

from ..core.alignment import Alignment
from ..core.distance import DistanceMatrix
from ..core.graph import HaplotypeNetwork
from ..core.haplotype import Haplotype
from ..core.haplotype import identify_haplotypes_from_alignment
from ..core.sequence import Sequence
from ..core.site_patterns import condense_site_patterns, is_ambiguous
from .base import NetworkAlgorithm


class TightSpanWalker(NetworkAlgorithm):
    """
    Construct a haplotype network using the Tight Span Walker algorithm.

    Parameters
    ----------
    distance_method : str, default='hamming'
        Kept for interface compatibility; TSW always uses PopART's
        site-weighted Hamming distances over condensed site patterns.
    tolerance : float, default=1e-6
        Relative tolerance for floating point comparisons (PopART's
        aboutEqual uses FLT_EPSILON).
    **kwargs : dict
        Additional parameters passed to the base class.

    Notes
    -----
    Internal vertices are inferred median points of the tight span; they
    have no sequence and are flagged with median markers. Best suited to
    small and medium datasets (the walk is O(n^2) geodesics, each
    scanning O(n^2) sample pairs).
    """

    def __init__(
        self, distance_method: str = 'hamming', tolerance: float = 1e-6, **kwargs
    ):
        """
        Initialize TightSpanWalker algorithm.

        Parameters
        ----------
        distance_method : str, default='hamming'
            Method for calculating distances.
        tolerance : float, default=1e-6
            Tolerance for floating point comparisons.
        **kwargs : dict
            Additional parameters.
        """
        super().__init__(distance_method, **kwargs)
        self.tolerance = tolerance

        # Algorithm state (reset per construct_network call)
        self._n_samples = 0
        self._d: List[List[float]] = []
        self._dt: Dict[Tuple[int, int], float] = {}
        self._vertex_map: Dict[Tuple[float, ...], int] = {}
        self._graph = nx.Graph()
        self._n_vertices = 0

    def construct_network(
        self, alignment: Alignment, distance_matrix: Optional[DistanceMatrix] = None
    ) -> HaplotypeNetwork:
        """
        Construct TSW network from sequence alignment.

        Parameters
        ----------
        alignment : Alignment
            Multiple sequence alignment.
        distance_matrix : DistanceMatrix, optional
            Ignored; TSW computes site-weighted distances over condensed
            site patterns, as PopART does.

        Returns
        -------
        HaplotypeNetwork
            Metric-preserving network with inferred internal vertices.
        """
        haplotypes = identify_haplotypes_from_alignment(alignment)

        if len(haplotypes) == 0:
            return HaplotypeNetwork()
        if len(haplotypes) == 1:
            network = HaplotypeNetwork()
            network.add_haplotype(haplotypes[0])
            return network

        condensed, weights = condense_site_patterns([h.data for h in haplotypes])
        n = len(haplotypes)

        self._n_samples = n
        self._d = [
            [
                self._weighted_distance(condensed[i], condensed[j], weights)
                for j in range(n)
            ]
            for i in range(n)
        ]
        self._compute_dt()
        self._vertex_map = {}
        self._graph = nx.Graph()
        self._graph.add_nodes_from(range(n))
        self._n_vertices = n

        for i in range(n):
            self._vertex_map[tuple(self._dt_row(i))] = i

        for i in range(n):
            for j in range(i):
                self._geodesic(i, j)

        # Build the HaplotypeNetwork: samples keep their haplotypes,
        # internal tight-span vertices become sequence-less medians
        network = HaplotypeNetwork()
        for haplotype in haplotypes:
            network.add_haplotype(haplotype)
        index_to_label = {i: haplotypes[i].id for i in range(n)}
        for idx in range(n, self._n_vertices):
            label = f'Median_{idx - n}'
            index_to_label[idx] = label
            median_seq = Sequence(
                id=label, data='', description='Inferred tight-span vertex'
            )
            network.add_haplotype(
                Haplotype(sequence=median_seq, sample_ids=[]), median_vector=True
            )
            network.graph.nodes[label]['is_median'] = True
        for u, v, attrs in self._graph.edges(data=True):
            network.add_edge(
                index_to_label[u], index_to_label[v], distance=attrs['distance']
            )

        return network

    # ------------------------------------------------------------------
    # Distance and dT bookkeeping
    # ------------------------------------------------------------------

    @staticmethod
    def _weighted_distance(s1: str, s2: str, weights: List[int]) -> float:
        """
        Site-weighted Hamming distance over condensed sequences.

        Parameters
        ----------
        s1 : str
            First condensed sequence.
        s2 : str
            Second condensed sequence.
        weights : list of int
            Site weights per condensed column.

        Returns
        -------
        float
            Weighted distance (ambiguous positions skipped).
        """
        total = 0
        for c1, c2, w in zip(s1, s2, weights):
            if is_ambiguous(c1) or is_ambiguous(c2):
                continue
            if c1 != c2:
                total += w
        return float(total)

    def _compute_dt(self) -> None:
        """
        Compute the tight-span metric over sample vertices.

        dT(i, j) = max over all samples k (including i and j) of
        abs(d(i,k) - d(j,k)), so dT(i, j) >= d(i, j) always holds.
        """
        n = self._n_samples
        self._dt = {}
        for i in range(n):
            for j in range(i):
                value = max(abs(self._d[i][k] - self._d[j][k]) for k in range(n))
                self._dt[(i, j)] = self._dt[(j, i)] = value

    def _dt_get(self, i: int, j: int) -> float:
        """
        Look up dT between two vertices.

        Parameters
        ----------
        i : int
            First vertex index.
        j : int
            Second vertex index.

        Returns
        -------
        float
            The tight-span distance.
        """
        if i == j:
            return 0.0
        return self._dt[(i, j)]

    def _dt_row(self, i: int) -> List[float]:
        """
        Return the dT vector from a vertex to every sample vertex.

        Parameters
        ----------
        i : int
            Vertex index.

        Returns
        -------
        list of float
            dT(i, k) for each sample k.
        """
        return [self._dt_get(i, k) for k in range(self._n_samples)]

    def _about_equal(self, a: float, b: float) -> bool:
        """
        Compare floats with relative tolerance (PopART's aboutEqual).

        Parameters
        ----------
        a : float
            First value.
        b : float
            Second value.

        Returns
        -------
        bool
            True when the values are equal within tolerance.
        """
        return math.isclose(a, b, rel_tol=self.tolerance, abs_tol=self.tolerance)

    # ------------------------------------------------------------------
    # Geodesic walk (port of TightSpanWalker::geodesic)
    # ------------------------------------------------------------------

    def _geodesic(self, f: int, g: int) -> None:
        """
        Walk the tight-span geodesic between two vertices.

        Iterative version of the C++ tail recursion: each step either
        connects f directly to g (when dT(f, g) equals delta) or
        materialises/reuses the internal vertex h one delta-step from f
        and continues from h.

        Parameters
        ----------
        f : int
            Start vertex index.
        g : int
            End vertex index.

        Raises
        ------
        RuntimeError
            On the same inconsistencies the C++ code reports: an
            uncoloured auxiliary vertex, a negative delta, or an
            apparent negative edge length.
        """
        n = self._n_samples

        while True:
            # Build the auxiliary graph K over sample vertices: an edge
            # (i, j) whenever dT(f,i) + dT(f,j) == d(i,j)
            adjacency: Dict[int, List[int]] = {i: [] for i in range(n)}
            edge_order: List[Tuple[int, int]] = []
            for i in range(n):
                f_i = self._dt_get(f, i)
                for j in range(i):
                    if self._about_equal(f_i + self._dt_get(f, j), self._d[i][j]):
                        adjacency[i].append(j)
                        adjacency[j].append(i)
                        edge_order.append((i, j))

            BLACK, GREEN, RED = 0, 1, 2
            colour = [BLACK] * n

            # Seed colours from edges lying on the f-g geodesic
            dt_fg = self._dt_get(f, g)
            for i, j in edge_order:
                # C++ names: edge (from=v, to=u) = our (i, j) -> v=i, u=j
                v, u = i, j
                fv, gv = self._dt_get(v, f), self._dt_get(v, g)
                fu, gu = self._dt_get(u, f), self._dt_get(u, g)
                weight = self._d[i][j]
                if fv < gv:
                    if (
                        self._about_equal(fv + dt_fg + gu, weight)
                        and colour[u] == BLACK
                    ):
                        colour[u] = GREEN
                        for neighbour in adjacency[u]:
                            colour[neighbour] = RED
                elif fu < gu:
                    if (
                        self._about_equal(fu + dt_fg + gv, weight)
                        and colour[v] == BLACK
                    ):
                        colour[v] = GREEN
                        for neighbour in adjacency[v]:
                            colour[neighbour] = RED

            # Propagate 2-colouring by BFS from sample vertex 0
            marked = [False] * n
            queue = [0]
            while queue:
                v = queue.pop(0)
                marked[v] = True
                if colour[v] == RED:
                    for u in adjacency[v]:
                        if not marked[u]:
                            queue.append(u)
                else:  # black or green
                    colour[v] = GREEN
                    for u in adjacency[v]:
                        colour[u] = RED
                        if not marked[u]:
                            queue.append(u)

            # delta = half the minimum over green pairs (i == j included)
            delta = float('inf')
            green = [i for i in range(n) if colour[i] == GREEN]
            for a_pos, i in enumerate(green):
                f_i = self._dt_get(f, i)
                for j in green[: a_pos + 1]:
                    delta = min(delta, f_i + self._dt_get(f, j) - self._d[i][j])
            delta /= 2

            if delta < 0:
                raise RuntimeError('Something is wrong, delta should be positive.')

            if self._about_equal(dt_fg, delta):
                if not self._graph.has_edge(f, g):
                    self._graph.add_edge(f, g, distance=delta)
                return

            if dt_fg > delta:
                new_dt_vector: List[float] = []
                for i in range(n):
                    f_i = self._dt_get(f, i)
                    if colour[i] == GREEN:
                        new_dt_vector.append(f_i - delta)
                    elif colour[i] == RED:
                        new_dt_vector.append(f_i + delta)
                    else:
                        raise RuntimeError('Uncoloured vertex!')

                key = tuple(new_dt_vector)
                existing = self._vertex_map.get(key)
                if existing is None:
                    h = self._n_vertices
                    self._n_vertices += 1
                    self._graph.add_node(h)
                    self._vertex_map[key] = h

                    for i in range(n):
                        self._dt[(h, i)] = self._dt[(i, h)] = new_dt_vector[i]
                    # dT between h and previously created internal
                    # vertices: max over sample pairs (both orders) of
                    # d(j,k) - dT(h,j) - dT(i,k)
                    for i in range(n, h):
                        dt_ih = -float('inf')
                        for j in range(n):
                            for k in range(n):
                                dt_ih = max(
                                    dt_ih,
                                    self._d[j][k]
                                    - self._dt_get(h, j)
                                    - self._dt_get(i, k),
                                    self._d[j][k]
                                    - self._dt_get(h, k)
                                    - self._dt_get(i, j),
                                )
                        self._dt[(h, i)] = self._dt[(i, h)] = dt_ih

                    self._graph.add_edge(f, h, distance=delta)
                else:
                    h = existing
                    if not self._graph.has_edge(f, h):
                        self._graph.add_edge(f, h, distance=delta)

                f = h
                continue

            raise RuntimeError('Apparent negative edge length between vertices g and h')

    def get_parameters(self) -> dict:
        """
        Get algorithm parameters.

        Returns
        -------
        dict
            Parameters including tolerance.
        """
        params = super().get_parameters()
        params['tolerance'] = self.tolerance
        return params

    def __repr__(self) -> str:
        """Detailed representation."""
        return (
            f'TightSpanWalker(distance={self.distance_method}, '
            f'tolerance={self.tolerance})'
        )

    def __str__(self) -> str:
        """Return string representation."""
        return self.__repr__()

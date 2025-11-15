"""
Minimum Spanning Network (MSN) algorithm for haplotype network construction.
"""

from typing import Dict, List, Optional, Tuple

from ..core.alignment import Alignment
from ..core.distance import DistanceMatrix
from ..core.graph import HaplotypeNetwork
from ..core.haplotype import identify_haplotypes_from_alignment as identify_haplotypes
from .mst import MinimumSpanningTree


class MinimumSpanningNetwork(MinimumSpanningTree):
    """
    Construct a Minimum Spanning Network from haplotype data.

    MSN extends MST by adding alternative connections at the same distance
    level, creating a network that shows all equally parsimonious relationships
    between haplotypes while removing redundant edges.

    This creates a more realistic representation of genetic relationships than
    a simple tree, as it can represent reticulation events and uncertainty in
    phylogenetic relationships.
    """

    def __init__(
        self,
        distance_method: str = 'hamming',
        epsilon: float = 0.0,
        max_connections: Optional[int] = None,
        **kwargs,
    ):
        """
        Initialize MSN algorithm.

        Parameters
        ----------
        distance_method :
            Method for calculating distances.
        epsilon :
            Tolerance for considering distances equal (default 0.0).
        max_connections :
            Maximum number of alternative connections per node.
        **kwargs :
            Additional parameters.
        """
        super().__init__(distance_method, algorithm='prim', **kwargs)
        self.epsilon = epsilon
        self.max_connections = max_connections

    def construct_network(
        self, alignment: Alignment, distance_matrix: Optional[DistanceMatrix] = None
    ) -> HaplotypeNetwork:
        """
            Construct MSN from sequence alignment.

            Parameters
            ----------
            alignment :
                Multiple sequence alignment.
            distance_matrix :
                Optional pre-computed distance matrix.

        Returns
        -------
            Haplotype network representing the MSN.
        """
        # Identify unique haplotypes
        haplotypes = identify_haplotypes(alignment)

        if len(haplotypes) <= 1:
            return super().construct_network(alignment, distance_matrix)

        # Calculate distances between haplotypes
        haplotype_dist_matrix = self._calculate_haplotype_distances(haplotypes)
        self._distance_matrix = haplotype_dist_matrix

        # Build initial MST
        mst_edges = self._prim_mst(haplotypes, haplotype_dist_matrix)

        # Add alternative connections at same distance
        msn_edges = self._add_alternative_connections(
            haplotypes, mst_edges, haplotype_dist_matrix
        )

        # Remove redundant edges
        final_edges = self._remove_redundant_edges(haplotypes, msn_edges)

        # Construct network
        network = self._build_network(haplotypes, final_edges)

        return network

    def _add_alternative_connections(
        self,
        haplotypes: List,
        mst_edges: List[Tuple[str, str, float]],
        distance_matrix: DistanceMatrix,
    ) -> List[Tuple[str, str, float]]:
        """
            Add alternative connections at the same distance level.

            For each distance level in the MST, add all edges at that distance
            (or within epsilon) that don't create redundancy.

            Parameters
            ----------
            haplotypes :
                List of Haplotype objects.
            mst_edges :
                MST edges from Prim's algorithm.
            distance_matrix :
                Distance matrix.

        Returns
        -------
            Extended list of edges including alternatives.
        """
        hap_ids = [h.id for h in haplotypes]

        # Track which edges are already in the network
        existing_edges = set()
        for id1, id2, dist in mst_edges:
            existing_edges.add((min(id1, id2), max(id1, id2)))

        # Get unique distances from MST
        mst_distances = sorted({dist for _, _, dist in mst_edges})

        all_edges = list(mst_edges)

        # For each distance level, add alternative edges
        for target_dist in mst_distances:
            # Find all possible edges at this distance (within epsilon)
            candidate_edges = []

            for i, id1 in enumerate(hap_ids):
                for id2 in hap_ids[i + 1 :]:
                    edge_key = (min(id1, id2), max(id1, id2))
                    if edge_key in existing_edges:
                        continue

                    dist = distance_matrix.get_distance(id1, id2)

                    # Check if distance matches target (within epsilon)
                    if abs(dist - target_dist) <= self.epsilon:
                        candidate_edges.append((id1, id2, dist))

            # Add candidate edges that create useful connections
            for id1, id2, dist in candidate_edges:
                # Check if adding this edge would be useful
                # (connects nodes that aren't already directly connected)
                edge_key = (min(id1, id2), max(id1, id2))

                # Add the edge
                all_edges.append((id1, id2, dist))
                existing_edges.add(edge_key)

                # Respect max_connections limit if specified
                if self.max_connections is not None:
                    conn_count1 = sum(1 for e in all_edges if id1 in (e[0], e[1]))
                    conn_count2 = sum(1 for e in all_edges if id2 in (e[0], e[1]))

                    if (
                        conn_count1 > self.max_connections
                        or conn_count2 > self.max_connections
                    ):
                        # Remove this edge if it violates max_connections
                        all_edges.pop()
                        existing_edges.remove(edge_key)

        return all_edges

    def _remove_redundant_edges(
        self, haplotypes: List, edges: List[Tuple[str, str, float]]
    ) -> List[Tuple[str, str, float]]:
        """
            Remove redundant edges from the network.

            An edge is redundant if removing it doesn't disconnect the network
            and there exists an alternative path of the same or shorter total length.

            Parameters
            ----------
            haplotypes :
                List of Haplotype objects.
            edges :
                List of edges.

        Returns
        -------
            List of non-redundant edges.
        """
        if len(edges) <= len(haplotypes) - 1:
            # Already minimal - can't remove any edges without disconnecting
            return edges

        # Build adjacency list
        adjacency: Dict[str, List[Tuple[str, float]]] = {}
        for id1, id2, dist in edges:
            if id1 not in adjacency:
                adjacency[id1] = []
            if id2 not in adjacency:
                adjacency[id2] = []
            adjacency[id1].append((id2, dist))
            adjacency[id2].append((id1, dist))

        # Try to remove each edge and check if network remains connected
        non_redundant = []

        for edge in edges:
            id1, id2, dist = edge

            # Temporarily remove edge
            adjacency[id1] = [(n, d) for n, d in adjacency[id1] if n != id2]
            adjacency[id2] = [(n, d) for n, d in adjacency[id2] if n != id1]

            # Check if still connected using BFS
            if self._is_connected(adjacency, id1, id2):
                # Check if alternative path exists with same or shorter length
                alt_path_length = self._shortest_path_length(adjacency, id1, id2)
                if alt_path_length is not None and alt_path_length <= dist:
                    # Edge is redundant - don't add it back
                    continue

            # Edge is not redundant - add it back
            adjacency[id1].append((id2, dist))
            adjacency[id2].append((id1, dist))
            non_redundant.append(edge)

        return non_redundant

    def _is_connected(
        self, adjacency: Dict[str, List[Tuple[str, float]]], start: str, end: str
    ) -> bool:
        """
            Check if two nodes are connected using BFS.

            Parameters
            ----------
            adjacency :
                Adjacency list representation.
            start :
                Start node ID.
            end :
                End node ID.

        Returns
        -------
            True if connected, False otherwise.
        """
        if start == end:
            return True

        if start not in adjacency or end not in adjacency:
            return False

        visited = {start}
        queue = [start]

        while queue:
            current = queue.pop(0)

            if current == end:
                return True

            for neighbor, _ in adjacency.get(current, []):
                if neighbor not in visited:
                    visited.add(neighbor)
                    queue.append(neighbor)

        return False

    def _shortest_path_length(
        self, adjacency: Dict[str, List[Tuple[str, float]]], start: str, end: str
    ) -> Optional[float]:
        """
            Find shortest path length between two nodes using Dijkstra's algorithm.

            Parameters
            ----------
            adjacency :
                Adjacency list representation.
            start :
                Start node ID.
            end :
                End node ID.

        Returns
        -------
            Shortest path length, or None if no path exists.
        """
        import heapq

        if start not in adjacency or end not in adjacency:
            return None

        # Priority queue: (distance, node)
        pq = [(0, start)]
        distances = {start: 0}
        visited = set()

        while pq:
            current_dist, current = heapq.heappop(pq)

            if current in visited:
                continue

            visited.add(current)

            if current == end:
                return current_dist

            for neighbor, edge_dist in adjacency.get(current, []):
                if neighbor not in visited:
                    new_dist = current_dist + edge_dist
                    if neighbor not in distances or new_dist < distances[neighbor]:
                        distances[neighbor] = new_dist
                        heapq.heappush(pq, (new_dist, neighbor))

        return distances.get(end)

    def get_parameters(self) -> dict:
        """Get algorithm parameters."""
        params = super().get_parameters()
        params['epsilon'] = self.epsilon
        params['max_connections'] = self.max_connections
        return params

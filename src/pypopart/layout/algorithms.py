"""
Layout algorithms for network visualization in PyPopART.

Provides various layout algorithms for positioning nodes in haplotype networks,
including force-directed, hierarchical, and custom layouts.
"""

import json
from typing import Dict, List, Optional, Tuple

import networkx as nx
import numpy as np

from ..core.graph import HaplotypeNetwork


class LayoutAlgorithm:
    """
    Base class for layout algorithms.

    Provides interface for computing node positions in network visualizations.
    """

    def __init__(self, network: HaplotypeNetwork):
        """
        Initialize layout algorithm with a network.

        Args:
            network: HaplotypeNetwork object
        """
        self.network = network
        self.graph = network._graph

    def compute(self, **kwargs) -> Dict[str, Tuple[float, float]]:
        """
        Compute node positions.

        Args:
            **kwargs: Algorithm-specific parameters

        Returns:
            Dictionary mapping node IDs to (x, y) positions
        """
        raise NotImplementedError('Subclasses must implement compute()')

    def save_layout(
        self, layout: Dict[str, Tuple[float, float]], filename: str
    ) -> None:
        """
        Save layout to a JSON file.

        Args:
            layout: Node positions dictionary
            filename: Output filename
        """
        # Convert tuples to lists for JSON serialization
        layout_serializable = {node: list(pos) for node, pos in layout.items()}

        with open(filename, 'w') as f:
            json.dump(layout_serializable, f, indent=2)

    @staticmethod
    def load_layout(filename: str) -> Dict[str, Tuple[float, float]]:
        """
        Load layout from a JSON file.

        Args:
            filename: Input filename

        Returns:
            Dictionary mapping node IDs to (x, y) positions
        """
        with open(filename, 'r') as f:
            layout_data = json.load(f)

        # Convert lists back to tuples
        return {node: tuple(pos) for node, pos in layout_data.items()}


class ForceDirectedLayout(LayoutAlgorithm):
    """
    Force-directed layout using spring algorithm.

    Simulates physical spring forces between connected nodes to create
    aesthetically pleasing layouts.
    """

    def compute(
        self,
        k: Optional[float] = None,
        iterations: int = 50,
        seed: Optional[int] = None,
        **kwargs,
    ) -> Dict[str, Tuple[float, float]]:
        """
        Compute force-directed layout.

        Args:
            k: Optimal distance between nodes (None for auto)
            iterations: Number of iterations for optimization
            seed: Random seed for reproducibility
            **kwargs: Additional parameters passed to spring_layout

        Returns:
            Node positions dictionary
        """
        layout = nx.spring_layout(
            self.graph, k=k, iterations=iterations, seed=seed, **kwargs
        )
        # Convert numpy arrays to tuples
        return {node: tuple(pos) for node, pos in layout.items()}


class CircularLayout(LayoutAlgorithm):
    """
    Circular layout arranging nodes in a circle.

    Places nodes evenly spaced around a circle, useful for showing
    connectivity patterns.
    """

    def compute(
        self, scale: float = 1.0, center: Optional[Tuple[float, float]] = None, **kwargs
    ) -> Dict[str, Tuple[float, float]]:
        """
        Compute circular layout.

        Args:
            scale: Scale factor for the layout
            center: Center position (x, y)
            **kwargs: Additional parameters passed to circular_layout

        Returns:
            Node positions dictionary
        """
        layout = nx.circular_layout(self.graph, scale=scale, center=center, **kwargs)
        # Convert numpy arrays to tuples
        return {node: tuple(pos) for node, pos in layout.items()}


class RadialLayout(LayoutAlgorithm):
    """
    Radial layout with center node and concentric rings.

    Places a central node at the origin and arranges other nodes
    in concentric circles based on distance from center.
    """

    def compute(
        self, center_node: Optional[str] = None, scale: float = 1.0, **kwargs
    ) -> Dict[str, Tuple[float, float]]:
        """
        Compute radial layout.

        Args:
            center_node: Node to place at center (most connected if None)
            scale: Scale factor for the layout
            **kwargs: Additional parameters

        Returns:
            Node positions dictionary
        """
        if not self.graph.nodes():
            return {}

        # Find center node if not specified
        if center_node is None:
            # Use node with highest degree
            center_node = max(self.graph.nodes(), key=lambda n: self.graph.degree(n))

        # Calculate distances from center using shortest path
        try:
            distances = nx.single_source_shortest_path_length(self.graph, center_node)
            # Add disconnected nodes at max distance + 1
            max_dist = max(distances.values()) if distances else 0
            for node in self.graph.nodes():
                if node not in distances:
                    distances[node] = max_dist + 1
        except nx.NetworkXError:
            # If graph is disconnected, use all nodes at max distance
            distances = dict.fromkeys(self.graph.nodes(), 1)
            distances[center_node] = 0

        # Group nodes by distance
        max_distance = max(distances.values()) if distances else 0
        rings: Dict[int, List[str]] = {i: [] for i in range(max_distance + 1)}

        for node, dist in distances.items():
            rings[dist].append(node)

        # Position nodes
        positions = {}
        positions[center_node] = (0.0, 0.0)

        for ring_idx, nodes in rings.items():
            if ring_idx == 0:
                continue

            radius = scale * ring_idx
            n_nodes = len(nodes)

            for i, node in enumerate(nodes):
                angle = 2 * np.pi * i / n_nodes
                x = radius * np.cos(angle)
                y = radius * np.sin(angle)
                positions[node] = (x, y)

        return positions


class HierarchicalLayout(LayoutAlgorithm):
    """
    Hierarchical layout arranging nodes in levels.

    Creates a tree-like structure with nodes arranged in horizontal
    levels based on distance from a root node.
    """

    def compute(
        self,
        root_node: Optional[str] = None,
        vertical: bool = True,
        width: float = 2.0,
        height: float = 2.0,
        **kwargs,
    ) -> Dict[str, Tuple[float, float]]:
        """
        Compute hierarchical layout.

        Args:
            root_node: Root node for hierarchy (most connected if None)
            vertical: If True, levels are horizontal; if False, levels are vertical
            width: Total width of the layout
            height: Total height of the layout
            **kwargs: Additional parameters

        Returns:
            Node positions dictionary
        """
        if not self.graph.nodes():
            return {}

        # Find root node if not specified
        if root_node is None:
            root_node = max(self.graph.nodes(), key=lambda n: self.graph.degree(n))

        # Calculate distances from root using BFS
        try:
            distances = nx.single_source_shortest_path_length(self.graph, root_node)
        except nx.NetworkXError:
            # If graph is disconnected, use default distances
            distances = dict.fromkeys(self.graph.nodes(), 1)
            distances[root_node] = 0

        # Group nodes by level
        max_level = max(distances.values()) if distances else 0
        levels: Dict[int, List[str]] = {i: [] for i in range(max_level + 1)}

        for node, level in distances.items():
            levels[level].append(node)

        # Position nodes
        positions = {}

        for level_idx, nodes in levels.items():
            n_nodes = len(nodes)

            if vertical:
                # Horizontal levels (top to bottom)
                y = height * (1 - level_idx / max(max_level, 1))

                for i, node in enumerate(nodes):
                    if n_nodes == 1:
                        x = width / 2
                    else:
                        x = width * i / (n_nodes - 1)
                    positions[node] = (x, y)
            else:
                # Vertical levels (left to right)
                x = width * level_idx / max(max_level, 1)

                for i, node in enumerate(nodes):
                    if n_nodes == 1:
                        y = height / 2
                    else:
                        y = height * (1 - i / (n_nodes - 1))
                    positions[node] = (x, y)

        return positions


class KamadaKawaiLayout(LayoutAlgorithm):
    """
    Kamada-Kawai layout algorithm.

    Uses energy minimization to position nodes based on graph-theoretic
    distances.
    """

    def compute(
        self, scale: float = 1.0, center: Optional[Tuple[float, float]] = None, **kwargs
    ) -> Dict[str, Tuple[float, float]]:
        """
        Compute Kamada-Kawai layout.

        Args:
            scale: Scale factor for the layout
            center: Center position (x, y)
            **kwargs: Additional parameters passed to kamada_kawai_layout

        Returns:
            Node positions dictionary
        """
        layout = nx.kamada_kawai_layout(
            self.graph, scale=scale, center=center, **kwargs
        )
        # Convert numpy arrays to tuples
        return {node: tuple(pos) for node, pos in layout.items()}


class ManualLayout(LayoutAlgorithm):
    """
    Manual layout with user-specified positions.

    Allows manual positioning of nodes or adjustment of existing layouts.
    """

    def __init__(
        self,
        network: HaplotypeNetwork,
        initial_positions: Optional[Dict[str, Tuple[float, float]]] = None,
    ):
        """
        Initialize manual layout.

        Args:
            network: HaplotypeNetwork object
            initial_positions: Starting positions for nodes
        """
        super().__init__(network)
        self.positions = initial_positions or {}

    def compute(self, **kwargs) -> Dict[str, Tuple[float, float]]:
        """
        Return current manual positions.

        Args:
            **kwargs: Ignored

        Returns:
            Node positions dictionary
        """
        # Fill in missing nodes with default layout
        if len(self.positions) < len(self.graph.nodes()):
            default_layout = nx.spring_layout(self.graph)
            for node in self.graph.nodes():
                if node not in self.positions:
                    self.positions[node] = default_layout[node]

        return self.positions

    def set_position(self, node: str, position: Tuple[float, float]) -> None:
        """
        Set position for a specific node.

        Args:
            node: Node ID
            position: (x, y) coordinates
        """
        if node not in self.graph.nodes():
            raise ValueError(f"Node '{node}' not in network")

        self.positions[node] = position

    def move_node(self, node: str, dx: float, dy: float) -> None:
        """
        Move a node by a relative offset.

        Args:
            node: Node ID
            dx: X offset
            dy: Y offset
        """
        if node not in self.positions:
            raise ValueError(f"Node '{node}' has no position set")

        x, y = self.positions[node]
        self.positions[node] = (x + dx, y + dy)


class LayoutManager:
    """
    Manager for network layout computation and persistence.

    Provides high-level interface for computing, saving, and loading
    network layouts.
    """

    def __init__(self, network: HaplotypeNetwork):
        """
        Initialize layout manager.

        Args:
            network: HaplotypeNetwork object
        """
        self.network = network
        self._algorithms = {
            'force_directed': ForceDirectedLayout,
            'spring': ForceDirectedLayout,  # Alias
            'circular': CircularLayout,
            'radial': RadialLayout,
            'hierarchical': HierarchicalLayout,
            'kamada_kawai': KamadaKawaiLayout,
            'manual': ManualLayout,
        }

    def compute_layout(
        self, algorithm: str = 'force_directed', **kwargs
    ) -> Dict[str, Tuple[float, float]]:
        """
        Compute network layout using specified algorithm.

        Args:
            algorithm: Layout algorithm name
            **kwargs: Algorithm-specific parameters

        Returns:
            Node positions dictionary

        Raises:
            ValueError: If algorithm not recognized
        """
        if algorithm not in self._algorithms:
            available = ', '.join(self._algorithms.keys())
            raise ValueError(f"Unknown algorithm '{algorithm}'. Available: {available}")

        layout_class = self._algorithms[algorithm]
        layout_algo = layout_class(self.network)

        return layout_algo.compute(**kwargs)

    def save_layout(
        self, layout: Dict[str, Tuple[float, float]], filename: str
    ) -> None:
        """
        Save layout to file.

        Args:
            layout: Node positions dictionary
            filename: Output filename (JSON format)
        """
        algo = LayoutAlgorithm(self.network)
        algo.save_layout(layout, filename)

    def load_layout(self, filename: str) -> Dict[str, Tuple[float, float]]:
        """
        Load layout from file.

        Args:
            filename: Input filename (JSON format)

        Returns:
            Node positions dictionary
        """
        return LayoutAlgorithm.load_layout(filename)

    def get_available_algorithms(self) -> List[str]:
        """
        Get list of available layout algorithms.

        Returns:
            List of algorithm names
        """
        return list(self._algorithms.keys())


class GeographicLayout(LayoutAlgorithm):
    """
    Geographic layout based on latitude/longitude coordinates.

    Places nodes at their geographic positions, optionally handling
    multiple locations for a single haplotype.
    """

    def __init__(self, network: HaplotypeNetwork):
        """
        Initialize geographic layout algorithm.

        Args:
            network: HaplotypeNetwork object
        """
        super().__init__(network)
        self._projection = None

    def compute(
        self,
        coordinates: Optional[Dict[str, Tuple[float, float]]] = None,
        projection: str = 'mercator',
        handle_multiple_locations: str = 'centroid',
        jitter_amount: float = 0.0,
        **kwargs,
    ) -> Dict[str, Tuple[float, float]]:
        """
        Compute geographic layout from coordinates.

        Args:
            coordinates: Dictionary mapping node IDs to (latitude, longitude) tuples
            projection: Map projection to use ('mercator', 'platecarree', 'orthographic')
            handle_multiple_locations: How to handle nodes in multiple locations:
                - 'centroid': Use centroid of all locations
                - 'first': Use first location
                - 'jitter': Place at first location with small random offset
            jitter_amount: Amount of random jitter to add (in degrees) when multiple
                          nodes share the same location
            **kwargs: Additional projection parameters

        Returns:
            Dictionary mapping node IDs to (x, y) positions in projected coordinates

        Raises:
            ValueError: If coordinates are not provided or invalid
        """
        if coordinates is None:
            # Try to extract coordinates from node metadata
            coordinates = self._extract_coordinates_from_metadata()

        if not coordinates:
            raise ValueError(
                'No geographic coordinates provided. '
                'Coordinates must be provided as argument or in node metadata '
                'with "latitude" and "longitude" attributes.'
            )

        # Validate coordinates
        for node_id, (lat, lon) in coordinates.items():
            if not -90 <= lat <= 90:
                raise ValueError(f'Invalid latitude for node {node_id}: {lat}')
            if not -180 <= lon <= 180:
                raise ValueError(f'Invalid longitude for node {node_id}: {lon}')

        # Apply projection
        layout = {}
        for node_id in self.graph.nodes():
            if node_id not in coordinates:
                # Node has no coordinates, will need to place it later
                continue

            lat, lon = coordinates[node_id]
            x, y = self._project_coordinates(lat, lon, projection)

            # Add jitter if requested
            if jitter_amount > 0:
                x += np.random.uniform(-jitter_amount, jitter_amount)
                y += np.random.uniform(-jitter_amount, jitter_amount)

            layout[node_id] = (x, y)

        # Handle nodes without coordinates using default layout
        if len(layout) < len(self.graph.nodes()):
            # Use spring layout for nodes without coordinates
            missing_nodes = [n for n in self.graph.nodes() if n not in layout]

            if len(missing_nodes) == len(self.graph.nodes()):
                # No nodes have coordinates, use full spring layout
                return nx.spring_layout(self.graph)
            else:
                # Create subgraph of missing nodes and layout separately
                subgraph = self.graph.subgraph(missing_nodes)
                sub_layout = nx.spring_layout(subgraph)

                # Scale and offset to place in unused space
                if layout:
                    # Find bounds of geographic layout
                    xs, ys = zip(*layout.values())
                    x_min, x_max = min(xs), max(xs)
                    y_min, y_max = min(ys), max(ys)
                    x_range = x_max - x_min

                    # Place missing nodes to the right of the main layout
                    for node, (x, y) in sub_layout.items():
                        layout[node] = (
                            x_max + x_range * 0.2 + x * x_range * 0.3,
                            y_min + y * (y_max - y_min),
                        )

        return layout

    def _extract_coordinates_from_metadata(self) -> Dict[str, Tuple[float, float]]:
        """
        Extract coordinates from node metadata.

        Returns:
            Dictionary mapping node IDs to (latitude, longitude) tuples
        """
        from ..io.metadata import extract_coordinates

        coordinates = {}
        for node_id in self.graph.nodes():
            node_data = self.graph.nodes[node_id]

            # Try to extract coordinates from node attributes
            if 'latitude' in node_data and 'longitude' in node_data:
                try:
                    coords = extract_coordinates(
                        node_data,
                        lat_column='latitude',
                        lon_column='longitude',
                        validate=True,
                    )
                    if coords:
                        coordinates[node_id] = coords
                except (ValueError, KeyError):
                    # Skip nodes with invalid coordinates
                    pass

        return coordinates

    def _project_coordinates(
        self, lat: float, lon: float, projection: str
    ) -> Tuple[float, float]:
        """
        Project geographic coordinates to map coordinates.

        Args:
            lat: Latitude
            lon: Longitude
            projection: Projection name

        Returns:
            Tuple of (x, y) projected coordinates
        """
        if projection.lower() == 'mercator':
            # Web Mercator projection
            x = lon
            # Handle poles (avoid division by zero)
            if abs(lat) >= 85:
                lat = 85 if lat > 0 else -85
            y = np.degrees(np.log(np.tan(np.radians(45 + lat / 2))))
        elif projection.lower() == 'platecarree':
            # Equirectangular projection (no transformation)
            x = lon
            y = lat
        elif projection.lower() == 'orthographic':
            # Orthographic projection (3D sphere view)
            # For simplicity, use simple perspective projection
            lat_rad = np.radians(lat)
            lon_rad = np.radians(lon)
            x = np.cos(lat_rad) * np.sin(lon_rad)
            y = np.sin(lat_rad)
        else:
            # Default to plate carrée
            x = lon
            y = lat

        return (x, y)

"""
Geographic visualization for haplotype networks.

Provides map-based visualizations using cartopy for static maps
and folium for interactive maps.
"""

from typing import Any, Dict, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np

from ..core.graph import HaplotypeNetwork
from ..layout.algorithms import GeographicLayout


class GeoVisualizer:
    """
    Geographic network visualizer using matplotlib and cartopy.

    Creates static map-based visualizations with haplotype networks
    overlaid on geographic maps.
    """

    def __init__(self, network: HaplotypeNetwork):
        """
        Initialize geographic visualizer.

        Args:
            network: HaplotypeNetwork object to visualize
        """
        self.network = network
        self.figure = None
        self.ax = None

    def plot(
        self,
        coordinates: Optional[Dict[str, Tuple[float, float]]] = None,
        projection: str = 'mercator',
        extent: Optional[Tuple[float, float, float, float]] = None,
        node_size_scale: float = 300.0,
        node_color_map: Optional[Dict[str, str]] = None,
        population_colors: Optional[Dict[str, str]] = None,
        edge_width_scale: float = 1.0,
        show_labels: bool = True,
        show_land: bool = True,
        show_coastlines: bool = True,
        show_borders: bool = False,
        land_color: str = 'lightgray',
        ocean_color: str = 'lightblue',
        figsize: Tuple[float, float] = (16, 12),
        title: Optional[str] = None,
        output_file: Optional[str] = None,
        **kwargs,
    ) -> Tuple[plt.Figure, plt.Axes]:
        """
        Create a geographic network visualization.

        Args:
            coordinates: Dictionary mapping node IDs to (latitude, longitude) tuples
            projection: Map projection ('mercator', 'platecarree', 'orthographic')
            extent: Map extent as (lon_min, lon_max, lat_min, lat_max)
            node_size_scale: Scaling factor for node sizes
            node_color_map: Custom color mapping {node_id: color}
            population_colors: Color mapping for populations {pop_name: color}
            edge_width_scale: Scaling factor for edge widths
            show_labels: Whether to show node labels
            show_land: Whether to show land features
            show_coastlines: Whether to show coastlines
            show_borders: Whether to show country borders
            land_color: Color for land areas
            ocean_color: Color for ocean areas
            figsize: Figure size (width, height) in inches
            title: Plot title
            output_file: Output file path (if provided, saves figure)
            **kwargs: Additional arguments

        Returns:
            Figure and axes objects

        Raises:
            ImportError: If cartopy is not installed
        """
        try:
            import cartopy.crs as ccrs
            import cartopy.feature as cfeature
        except ImportError as e:
            raise ImportError(
                'cartopy is required for geographic visualization. '
                'Install it with: pip install cartopy'
            ) from e

        # Set up projection
        if projection.lower() == 'mercator':
            proj = ccrs.Mercator()
        elif projection.lower() == 'platecarree':
            proj = ccrs.PlateCarree()
        elif projection.lower() == 'orthographic':
            proj = ccrs.Orthographic()
        else:
            proj = ccrs.PlateCarree()

        # Create figure and axes
        self.figure = plt.figure(figsize=figsize)
        self.ax = plt.axes(projection=proj)

        # Set extent if provided
        if extent:
            self.ax.set_extent(extent, crs=ccrs.PlateCarree())

        # Add map features
        if show_land:
            self.ax.add_feature(
                cfeature.LAND, facecolor=land_color, edgecolor='none', zorder=0
            )
        if show_coastlines:
            self.ax.add_feature(
                cfeature.COASTLINE, linewidth=0.5, edgecolor='black', zorder=1
            )
        if show_borders:
            self.ax.add_feature(
                cfeature.BORDERS, linewidth=0.3, edgecolor='gray', zorder=1
            )

        # Extract coordinates if not provided
        if coordinates is None:
            geo_layout = GeographicLayout(self.network)
            coordinates = geo_layout._extract_coordinates_from_metadata()

        # Prepare node attributes
        node_sizes = self._compute_node_sizes(node_size_scale)
        node_colors = self._compute_node_colors(node_color_map, population_colors)

        # Prepare edge attributes
        edge_widths = self._compute_edge_widths(edge_width_scale)

        # Draw edges
        graph = self.network._graph
        for u, v in graph.edges():
            if u in coordinates and v in coordinates:
                lat1, lon1 = coordinates[u]
                lat2, lon2 = coordinates[v]

                width = edge_widths.get((u, v), 1.0)

                self.ax.plot(
                    [lon1, lon2],
                    [lat1, lat2],
                    color='gray',
                    linewidth=width,
                    alpha=0.7,
                    transform=ccrs.PlateCarree(),
                    zorder=2,
                )

        # Draw nodes
        for node in graph.nodes():
            if node in coordinates:
                lat, lon = coordinates[node]
                size = node_sizes.get(node, 100)
                color = node_colors.get(node, 'lightblue')

                # Check if median vector
                is_median = self.network.is_median_vector(node)
                marker = 's' if is_median else 'o'

                self.ax.plot(
                    lon,
                    lat,
                    marker=marker,
                    markersize=np.sqrt(size) / 5,
                    color=color,
                    markeredgecolor='black',
                    markeredgewidth=1,
                    transform=ccrs.PlateCarree(),
                    zorder=3,
                )

                # Add label if requested
                if show_labels:
                    self.ax.text(
                        lon,
                        lat,
                        node,
                        fontsize=8,
                        ha='center',
                        va='bottom',
                        transform=ccrs.PlateCarree(),
                        zorder=4,
                    )

        # Add title
        if title:
            self.ax.set_title(title, fontsize=14, fontweight='bold')

        # Add gridlines
        gl = self.ax.gridlines(draw_labels=True, alpha=0.3)
        gl.top_labels = False
        gl.right_labels = False

        plt.tight_layout()

        # Save if output file provided
        if output_file:
            self.figure.savefig(output_file, dpi=300, bbox_inches='tight')

        return self.figure, self.ax

    def _compute_node_sizes(self, scale: float) -> Dict[str, float]:
        """
        Compute node sizes based on frequency/count.

        Args:
            scale: Scaling factor

        Returns:
            Dictionary mapping node IDs to sizes
        """
        sizes = {}
        graph = self.network._graph

        for node in graph.nodes():
            # Get frequency from node attributes
            freq = graph.nodes[node].get('frequency', 1)
            sizes[node] = scale * np.sqrt(freq)

        return sizes

    def _compute_node_colors(
        self,
        color_map: Optional[Dict[str, str]] = None,
        population_colors: Optional[Dict[str, str]] = None,
    ) -> Dict[str, str]:
        """
        Compute node colors.

        Args:
            color_map: Custom color mapping
            population_colors: Population-based colors

        Returns:
            Dictionary mapping node IDs to colors
        """
        if color_map:
            return color_map

        colors = {}
        graph = self.network._graph

        for node in graph.nodes():
            if self.network.is_median_vector(node):
                colors[node] = 'lightgray'
            elif population_colors:
                pop = graph.nodes[node].get('population', 'unknown')
                colors[node] = population_colors.get(pop, 'lightblue')
            else:
                colors[node] = 'lightblue'

        return colors

    def _compute_edge_widths(self, scale: float) -> Dict[Tuple[str, str], float]:
        """
        Compute edge widths based on mutation distance.

        Args:
            scale: Scaling factor

        Returns:
            Dictionary mapping edge tuples to widths
        """
        widths = {}
        graph = self.network._graph

        for u, v in graph.edges():
            # Inverse relationship: fewer mutations = thicker line
            mutations = graph[u][v].get('mutations', 1)
            width = scale * (1.0 / max(mutations, 0.5))
            widths[(u, v)] = width

        return widths


class InteractiveGeoVisualizer:
    """
    Interactive geographic network visualizer using folium.

    Creates interactive HTML maps with haplotype networks overlaid
    on web-based maps (OpenStreetMap, etc.).
    """

    def __init__(self, network: HaplotypeNetwork):
        """
        Initialize interactive geographic visualizer.

        Args:
            network: HaplotypeNetwork object to visualize
        """
        self.network = network

    def plot(
        self,
        coordinates: Optional[Dict[str, Tuple[float, float]]] = None,
        base_map: str = 'OpenStreetMap',
        zoom_start: int = 4,
        node_size_scale: float = 10.0,
        node_color_map: Optional[Dict[str, str]] = None,
        population_colors: Optional[Dict[str, str]] = None,
        show_edges: bool = True,
        show_labels: bool = True,
        output_file: Optional[str] = None,
        **kwargs,
    ) -> Any:
        """
        Create an interactive geographic network visualization.

        Args:
            coordinates: Dictionary mapping node IDs to (latitude, longitude) tuples
            base_map: Base map style ('OpenStreetMap', 'Stamen Terrain', 'CartoDB positron')
            zoom_start: Initial zoom level (1-18)
            node_size_scale: Scaling factor for node sizes
            node_color_map: Custom color mapping {node_id: color}
            population_colors: Color mapping for populations {pop_name: color}
            show_edges: Whether to show edges
            show_labels: Whether to show labels in popups
            output_file: Output HTML file path (if provided, saves map)
            **kwargs: Additional arguments

        Returns:
            Folium Map object

        Raises:
            ImportError: If folium is not installed
        """
        try:
            import folium
        except ImportError as e:
            raise ImportError(
                'folium is required for interactive geographic visualization. '
                'Install it with: pip install folium'
            ) from e

        # Extract coordinates if not provided
        if coordinates is None:
            geo_layout = GeographicLayout(self.network)
            coordinates = geo_layout._extract_coordinates_from_metadata()

        if not coordinates:
            raise ValueError(
                'No geographic coordinates available. '
                'Add latitude/longitude metadata to nodes.'
            )

        # Calculate center of map
        lats, lons = zip(*coordinates.values())
        center_lat = sum(lats) / len(lats)
        center_lon = sum(lons) / len(lons)

        # Create map
        m = folium.Map(
            location=[center_lat, center_lon],
            zoom_start=zoom_start,
            tiles=base_map,
        )

        # Prepare node attributes
        node_sizes = self._compute_node_sizes(node_size_scale)
        node_colors = self._compute_node_colors(node_color_map, population_colors)

        # Draw edges
        if show_edges:
            graph = self.network._graph
            for u, v in graph.edges():
                if u in coordinates and v in coordinates:
                    lat1, lon1 = coordinates[u]
                    lat2, lon2 = coordinates[v]

                    mutations = graph[u][v].get('mutations', 1)

                    folium.PolyLine(
                        locations=[[lat1, lon1], [lat2, lon2]],
                        color='gray',
                        weight=2,
                        opacity=0.5,
                        popup=f'{u} - {v}: {mutations} mutations',
                    ).add_to(m)

        # Draw nodes
        graph = self.network._graph
        for node in graph.nodes():
            if node in coordinates:
                lat, lon = coordinates[node]
                size = node_sizes.get(node, 10)
                color = node_colors.get(node, 'blue')

                # Get node info for popup
                freq = graph.nodes[node].get('frequency', 1)
                pop = graph.nodes[node].get('population', 'unknown')
                is_median = self.network.is_median_vector(node)

                popup_text = f'<b>{node}</b><br>'
                popup_text += f'Frequency: {freq}<br>'
                if not is_median:
                    popup_text += f'Population: {pop}<br>'
                else:
                    popup_text += 'Type: Median Vector<br>'
                popup_text += f'Location: {lat:.4f}, {lon:.4f}'

                folium.CircleMarker(
                    location=[lat, lon],
                    radius=size,
                    color='black',
                    fillColor=color,
                    fillOpacity=0.7,
                    weight=2,
                    popup=folium.Popup(popup_text, max_width=250),
                ).add_to(m)

        # Save if output file provided
        if output_file:
            m.save(output_file)

        return m

    def _compute_node_sizes(self, scale: float) -> Dict[str, float]:
        """Compute node sizes based on frequency."""
        sizes = {}
        graph = self.network._graph

        for node in graph.nodes():
            freq = graph.nodes[node].get('frequency', 1)
            sizes[node] = scale * np.sqrt(freq)

        return sizes

    def _compute_node_colors(
        self,
        color_map: Optional[Dict[str, str]] = None,
        population_colors: Optional[Dict[str, str]] = None,
    ) -> Dict[str, str]:
        """Compute node colors."""
        if color_map:
            return color_map

        colors = {}
        graph = self.network._graph

        for node in graph.nodes():
            if self.network.is_median_vector(node):
                colors[node] = 'gray'
            elif population_colors:
                pop = graph.nodes[node].get('population', 'unknown')
                colors[node] = population_colors.get(pop, 'blue')
            else:
                colors[node] = 'blue'

        return colors

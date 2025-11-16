"""
Visualization module for PyPopART.

Provides functions for creating static and interactive visualizations
of haplotype networks, including geographic visualizations.
"""

from .cytoscape_plot import (
    InteractiveCytoscapePlotter,
    create_cytoscape_network,
)
from .geo_visualization import GeoVisualizer, InteractiveGeoVisualizer
from .interactive_plot import (
    InteractiveNetworkPlotter,
    create_interactive_figure,
    plot_interactive_network,
)
from .static_plot import StaticNetworkPlotter, create_publication_figure, plot_network

__all__ = [
    # Static plotting
    'StaticNetworkPlotter',
    'plot_network',
    'create_publication_figure',
    # Interactive plotting
    'InteractiveNetworkPlotter',
    'plot_interactive_network',
    'create_interactive_figure',
    # Cytoscape plotting
    'InteractiveCytoscapePlotter',
    'create_cytoscape_network',
    # Geographic plotting
    'GeoVisualizer',
    'InteractiveGeoVisualizer',
]

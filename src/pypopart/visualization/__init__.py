"""
Visualization module for PyPopART.

Provides functions for creating static and interactive visualizations
of haplotype networks, including geographic visualizations.
"""

from .interactive_plot import (
    InteractiveNetworkPlotter,
    create_interactive_figure,
    plot_interactive_network,
)
from .static_plot import StaticNetworkPlotter, create_publication_figure, plot_network
from .geo_visualization import GeoVisualizer, InteractiveGeoVisualizer

__all__ = [
    # Static plotting
    'StaticNetworkPlotter',
    'plot_network',
    'create_publication_figure',
    # Interactive plotting
    'InteractiveNetworkPlotter',
    'plot_interactive_network',
    'create_interactive_figure',
    # Geographic plotting
    'GeoVisualizer',
    'InteractiveGeoVisualizer',
]

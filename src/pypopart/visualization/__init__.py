"""
Visualization module for PyPopART.

Provides functions for creating static and interactive visualizations
of haplotype networks.
"""

from .static_plot import (
    StaticNetworkPlotter,
    plot_network,
    create_publication_figure
)

from .interactive_plot import (
    InteractiveNetworkPlotter,
    plot_interactive_network,
    create_interactive_figure
)

__all__ = [
    # Static plotting
    'StaticNetworkPlotter',
    'plot_network',
    'create_publication_figure',
    # Interactive plotting
    'InteractiveNetworkPlotter',
    'plot_interactive_network',
    'create_interactive_figure',
]

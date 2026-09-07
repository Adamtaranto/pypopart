"""
Static network visualization using matplotlib for PyPopART.

Provides matplotlib-based plotting functions for haplotype networks
with customizable node sizes, colors, edge styles, and layouts.
"""

from typing import Any, Dict, List, Optional, Tuple

from matplotlib.lines import Line2D
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

from ..core.graph import HaplotypeNetwork


class StaticNetworkPlotter:
    """
    Static network plotter using matplotlib.

    Generates publication-quality static plots of haplotype networks
    with customizable styling for nodes, edges, labels, and legends.

    Parameters
    ----------
    network : HaplotypeNetwork
        HaplotypeNetwork object to visualize.
    """

    def __init__(self, network: HaplotypeNetwork):
        """
        Initialize plotter with a haplotype network.

        Parameters
        ----------
        network : HaplotypeNetwork
            HaplotypeNetwork object to visualize.
        """
        self.network = network
        self.figure = None
        self.ax = None

    def plot(
        self,
        layout: Optional[Dict[str, Tuple[float, float]]] = None,
        layout_algorithm: str = 'spring',
        node_size_scale: float = 300.0,
        node_color_map: Optional[Dict[str, str]] = None,
        population_colors: Optional[Dict[str, str]] = None,
        edge_width_scale: float = 1.0,
        show_labels: bool = True,
        show_mutations: bool = True,
        median_vector_color: str = 'lightgray',
        median_vector_marker: str = 's',
        figsize: Tuple[float, float] = (12, 10),
        title: Optional[str] = None,
        **kwargs,
    ) -> Tuple[plt.Figure, plt.Axes]:
        """
        Create a static network plot.

        Parameters
        ----------
        layout : Dict[str, Tuple[float, float]], optional
            Pre-computed node positions {node_id: (x, y)}.
        layout_algorithm : str, default='spring'
            NetworkX layout algorithm ('spring', 'circular', 'kamada_kawai').
        node_size_scale : float, default=300.0
            Scaling factor for node sizes.
        node_color_map : Dict[str, str], optional
            Custom color mapping {node_id: color}.
        population_colors : Dict[str, str], optional
            Color mapping for populations {pop_name: color}.
        edge_width_scale : float, default=1.0
            Scaling factor for edge widths.
        show_labels : bool, default=True
            Whether to show node labels.
        show_mutations : bool, default=True
            Whether to show mutation counts on edges.
        median_vector_color : str, default='lightgray'
            Color for median vector nodes.
        median_vector_marker : str, default='s'
            Marker shape for median vectors ('s'=square, 'o'=circle).
        figsize : Tuple[float, float], default=(12, 10)
            Figure size (width, height) in inches.
        title : str, optional
            Plot title.
        **kwargs : dict
            Additional arguments passed to networkx drawing functions.

        Returns
        -------
        Tuple[plt.Figure, plt.Axes]
            Figure and axes objects.
        """
        # Create figure and axes
        self.figure, self.ax = plt.subplots(figsize=figsize)

        # Get graph and compute layout if not provided
        graph = self.network._graph
        if layout is None:
            layout = self._compute_layout(graph, layout_algorithm)

        # Prepare node attributes
        node_sizes = self._compute_node_sizes(node_size_scale)
        node_colors = self._compute_node_colors(
            node_color_map, population_colors, median_vector_color
        )

        # Prepare edge attributes
        edge_widths = self._compute_edge_widths(edge_width_scale)

        # Separate haplotypes and median vectors
        haplotype_nodes = [
            n for n in graph.nodes() if not self.network.is_median_vector(n)
        ]
        median_nodes = self.network.median_vector_ids

        # Draw haplotype nodes
        if haplotype_nodes:
            hap_sizes = [node_sizes[n] for n in haplotype_nodes]
            hap_colors = [node_colors[n] for n in haplotype_nodes]
            nx.draw_networkx_nodes(
                graph,
                layout,
                nodelist=haplotype_nodes,
                node_size=hap_sizes,
                node_color=hap_colors,
                node_shape='o',
                edgecolors='black',
                linewidths=1.5,
                ax=self.ax,
                **{k: v for k, v in kwargs.items() if k.startswith('node_')},
            )

        # Draw median vector nodes
        if median_nodes:
            med_sizes = [node_sizes[n] for n in median_nodes]
            med_colors = [node_colors[n] for n in median_nodes]
            nx.draw_networkx_nodes(
                graph,
                layout,
                nodelist=median_nodes,
                node_size=med_sizes,
                node_color=med_colors,
                node_shape=median_vector_marker,
                edgecolors='black',
                linewidths=1.5,
                ax=self.ax,
                **{k: v for k, v in kwargs.items() if k.startswith('node_')},
            )

        # Draw edges
        nx.draw_networkx_edges(
            graph,
            layout,
            width=edge_widths,
            edge_color='gray',
            alpha=0.6,
            ax=self.ax,
            **{k: v for k, v in kwargs.items() if k.startswith('edge_')},
        )

        # Draw labels if requested
        if show_labels:
            labels = {n: n for n in graph.nodes()}
            nx.draw_networkx_labels(
                graph,
                layout,
                labels=labels,
                font_size=8,
                font_weight='bold',
                ax=self.ax,
            )

        # Draw mutation counts on edges if requested
        if show_mutations:
            self._draw_edge_labels(graph, layout)

        # Add title
        if title:
            self.ax.set_title(title, fontsize=14, fontweight='bold', pad=20)
        elif self.network.name:
            self.ax.set_title(self.network.name, fontsize=14, fontweight='bold', pad=20)

        # Remove axes
        self.ax.axis('off')

        # Tight layout
        plt.tight_layout()

        return self.figure, self.ax

    def add_legend(
        self,
        population_colors: Optional[Dict[str, str]] = None,
        show_median_vectors: bool = True,
        show_size_scale: bool = True,
        loc: str = 'best',
        **kwargs,
    ) -> None:
        """
        Add a legend to the plot.

        Parameters
        ----------
        population_colors : Dict[str, str], optional
            Population color mapping {pop_name: color}.
        show_median_vectors : bool, default=True
            Whether to include median vectors in legend.
        show_size_scale : bool, default=True
            Whether to show node size scale.
        loc : str, default='best'
            Legend location.
        **kwargs : dict
            Additional arguments passed to plt.legend().
        """
        if self.ax is None:
            raise ValueError('No plot exists. Call plot() first.')

        legend_elements = []

        # Add population colors
        if population_colors:
            for pop_name, color in sorted(population_colors.items()):
                legend_elements.append(mpatches.Patch(color=color, label=pop_name))

        # Add median vectors
        if show_median_vectors and len(self.network.median_vector_ids) > 0:
            legend_elements.append(
                Line2D(
                    [0],
                    [0],
                    marker='s',
                    color='w',
                    markerfacecolor='lightgray',
                    markersize=10,
                    markeredgecolor='black',
                    markeredgewidth=1.5,
                    label='Median Vector',
                    linestyle='None',
                )
            )

        # Add size scale examples if requested
        if show_size_scale:
            # Find range of frequencies
            frequencies = []
            for node in self.network._graph.nodes():
                if not self.network.is_median_vector(node):
                    hap = self.network.get_haplotype(node)
                    if hap:
                        frequencies.append(hap.frequency)

            if frequencies:
                min_freq = min(frequencies)
                max_freq = max(frequencies)

                # Add size legend for min and max
                if min_freq != max_freq:
                    legend_elements.append(
                        Line2D(
                            [0],
                            [0],
                            marker='o',
                            color='w',
                            markerfacecolor='gray',
                            markersize=5,
                            markeredgecolor='black',
                            markeredgewidth=1,
                            label=f'n={min_freq}',
                            linestyle='None',
                        )
                    )
                    legend_elements.append(
                        Line2D(
                            [0],
                            [0],
                            marker='o',
                            color='w',
                            markerfacecolor='gray',
                            markersize=12,
                            markeredgecolor='black',
                            markeredgewidth=1,
                            label=f'n={max_freq}',
                            linestyle='None',
                        )
                    )

        if legend_elements:
            self.ax.legend(
                handles=legend_elements,
                loc=loc,
                frameon=True,
                fancybox=True,
                shadow=True,
                **kwargs,
            )

    def add_scale_bar(
        self,
        num_mutations: int = 1,
        position: Tuple[float, float] = (0.05, 0.05),
        length: float = 0.1,
        **kwargs,
    ) -> None:
        """
        Add a scale bar showing mutation distance.

        Parameters
        ----------
        num_mutations : int, default=1
            Number of mutations represented by scale bar.
        position : Tuple[float, float], default=(0.05, 0.05)
            Position as fraction of axes (x, y).
        length : float, default=0.1
            Length of scale bar as fraction of axes width.
        **kwargs : dict
            Additional arguments for the line and text.
        """
        if self.ax is None:
            raise ValueError('No plot exists. Call plot() first.')

        # Get axes limits
        xlim = self.ax.get_xlim()
        ylim = self.ax.get_ylim()

        # Calculate absolute position and length
        x_start = xlim[0] + (xlim[1] - xlim[0]) * position[0]
        y_pos = ylim[0] + (ylim[1] - ylim[0]) * position[1]
        bar_length = (xlim[1] - xlim[0]) * length

        # Draw scale bar
        self.ax.plot(
            [x_start, x_start + bar_length],
            [y_pos, y_pos],
            'k-',
            linewidth=2,
            solid_capstyle='butt',
        )

        # Add text label
        label = f'{num_mutations} mutation{"s" if num_mutations != 1 else ""}'
        self.ax.text(
            x_start + bar_length / 2,
            y_pos - (ylim[1] - ylim[0]) * 0.02,
            label,
            ha='center',
            va='top',
            fontsize=10,
            fontweight='bold',
        )

    def add_statistics_annotation(
        self,
        stats: Optional[Dict[str, Any]] = None,
        position: Tuple[float, float] = (0.02, 0.98),
        **kwargs,
    ) -> None:
        """
        Add network statistics as text annotation.

        Parameters
        ----------
        stats : Dict[str, Any], optional
            Dictionary of statistics to display.
        position : Tuple[float, float], default=(0.02, 0.98)
            Position as fraction of axes (x, y).
        **kwargs : dict
            Additional arguments for the text box.
        """
        if self.ax is None:
            raise ValueError('No plot exists. Call plot() first.')

        if stats is None:
            # Get basic network stats
            net_stats = self.network.calculate_stats()
            stats = {
                'Haplotypes': net_stats.num_haplotypes,
                'Samples': net_stats.total_samples,
                'Median Vectors': net_stats.num_median_vectors,
                'Edges': net_stats.num_edges,
            }

        # Format statistics text
        text_lines = []
        for key, value in stats.items():
            if isinstance(value, float):
                text_lines.append(f'{key}: {value:.2f}')
            else:
                text_lines.append(f'{key}: {value}')
        text = '\n'.join(text_lines)

        # Add text box
        props = {
            'boxstyle': 'round',
            'facecolor': 'white',
            'alpha': 0.8,
            'edgecolor': 'black',
        }
        props.update(kwargs.get('bbox', {}))

        self.ax.text(
            position[0],
            position[1],
            text,
            transform=self.ax.transAxes,
            fontsize=10,
            verticalalignment='top',
            bbox=props,
        )

    def save(
        self, filename: str, dpi: int = 300, bbox_inches: str = 'tight', **kwargs
    ) -> None:
        """
        Save the plot to a file.

        Parameters
        ----------
        filename : str
            Output filename (extension determines format: .png, .pdf, .svg).
        dpi : int, default=300
            Resolution in dots per inch.
        bbox_inches : str, default='tight'
            Bounding box setting.
        **kwargs : dict
            Additional arguments passed to plt.savefig().
        """
        if self.figure is None:
            raise ValueError('No plot exists. Call plot() first.')

        self.figure.savefig(filename, dpi=dpi, bbox_inches=bbox_inches, **kwargs)

    def _compute_layout(
        self, graph: nx.Graph, algorithm: str
    ) -> Dict[str, Tuple[float, float]]:
        """
        Compute node layout via the shared LayoutManager.

        Parameters
        ----------
        graph : nx.Graph
            NetworkX graph (unused; the manager works on self.network).
        algorithm : str
            Layout algorithm name.

        Returns
        -------
        dict
            Dictionary mapping node IDs to (x, y) positions.
        """
        from ..layout.algorithms import LayoutManager

        return LayoutManager(self.network).compute_layout(algorithm)

    def _compute_node_sizes(self, scale: float) -> Dict[str, float]:
        """
        Compute node sizes based on haplotype frequencies.

        Parameters
        ----------
        scale : float
            Scaling factor for node sizes.

        Returns
        -------
        Dict[str, float]
            Dictionary mapping node IDs to sizes.
        """
        sizes = {}
        for node in self.network._graph.nodes():
            if self.network.is_median_vector(node):
                # Median vectors get a fixed small size
                sizes[node] = scale * 0.3
            else:
                hap = self.network.get_haplotype(node)
                if hap:
                    # Size proportional to square root of frequency for better visual scaling
                    sizes[node] = scale * np.sqrt(hap.frequency)
                else:
                    sizes[node] = scale * 0.5

        return sizes

    def _compute_node_colors(
        self,
        node_color_map: Optional[Dict[str, str]],
        population_colors: Optional[Dict[str, str]],
        median_vector_color: str,
    ) -> Dict[str, str]:
        """
        Compute node colors based on population or custom mapping.

        Parameters
        ----------
        node_color_map : Dict[str, str], optional
            Custom node color mapping.
        population_colors : Dict[str, str], optional
            Population color mapping.
        median_vector_color : str
            Color for median vectors.

        Returns
        -------
        Dict[str, str]
            Dictionary mapping node IDs to colors.
        """
        from .style import node_color

        colors = {}
        for node in self.network._graph.nodes():
            hap = self.network.get_haplotype(node)
            colors[node] = node_color(
                node,
                hap,
                self.network.is_median_vector(node),
                node_color_map=node_color_map,
                population_colors=population_colors,
                median_vector_color=median_vector_color,
            )
        return colors

    def _compute_edge_widths(self, scale: float) -> List[float]:
        """
        Compute edge widths based on mutation distances.

        Parameters
        ----------
        scale : float
            Scaling factor for edge widths.

        Returns
        -------
        List[float]
            List of edge widths.
        """
        widths = []
        graph = self.network._graph

        for u, v in graph.edges():
            # Get edge weight (distance/mutations)
            weight = graph[u][v].get('weight', 1)
            # Inverse relationship: fewer mutations = thicker line
            width = scale * max(0.5, 3.0 / max(weight, 1))
            widths.append(width)

        return widths

    def _draw_edge_labels(
        self, graph: nx.Graph, layout: Dict[str, Tuple[float, float]]
    ) -> None:
        """
        Draw mutation counts on edges.

        Parameters
        ----------
        graph : nx.Graph
            NetworkX graph.
        layout : Dict[str, Tuple[float, float]]
            Node positions.
        """
        edge_labels = {}
        for u, v in graph.edges():
            weight = graph[u][v].get('weight', 1)
            if weight > 0:
                edge_labels[(u, v)] = int(weight)

        if edge_labels:
            nx.draw_networkx_edge_labels(
                graph,
                layout,
                edge_labels=edge_labels,
                font_size=7,
                bbox={
                    'boxstyle': 'round',
                    'facecolor': 'white',
                    'alpha': 0.7,
                    'edgecolor': 'none',
                },
                ax=self.ax,
            )


def plot_network(network: HaplotypeNetwork, **kwargs) -> Tuple[plt.Figure, plt.Axes]:
    """
    Plot a haplotype network.

    Parameters
    ----------
    network : HaplotypeNetwork
        HaplotypeNetwork object to visualize.
    **kwargs : dict
        Arguments passed to StaticNetworkPlotter.plot().

    Returns
    -------
    Tuple[plt.Figure, plt.Axes]
        Figure and axes objects.
        Example:
        >>> from pypopart.core.graph import HaplotypeNetwork
        >>> from pypopart.visualization.static_plot import plot_network
        >>> network = HaplotypeNetwork()
        >>> # ... build network ...
        >>> fig, ax = plot_network(network, layout_algorithm='spring')
        >>> plt.show()
    """
    plotter = StaticNetworkPlotter(network)
    return plotter.plot(**kwargs)


def create_publication_figure(
    network: HaplotypeNetwork,
    population_colors: Optional[Dict[str, str]] = None,
    filename: Optional[str] = None,
    **kwargs,
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Create a publication-ready figure with legend and scale bar.

    Parameters
    ----------
    network : HaplotypeNetwork
        HaplotypeNetwork object to visualize.
    population_colors : Dict[str, str], optional
        Color mapping for populations.
    filename : str, optional
        Optional filename to save figure.
    **kwargs : dict
        Additional arguments passed to plot().

    Returns
    -------
    Tuple[plt.Figure, plt.Axes]
        Figure and axes objects.
        Example:
        >>> fig, ax = create_publication_figure(
        ...     network,
        ...     population_colors={'PopA': 'red', 'PopB': 'blue'},
        ...     filename='network.pdf'
        ... )
    """
    plotter = StaticNetworkPlotter(network)

    # Create plot with defaults optimized for publication
    fig, ax = plotter.plot(
        population_colors=population_colors,
        show_labels=True,
        show_mutations=True,
        figsize=(10, 8),
        **kwargs,
    )

    # Add legend if population colors provided
    if population_colors:
        plotter.add_legend(
            population_colors=population_colors,
            show_median_vectors=True,
            show_size_scale=True,
            loc='upper right',
        )

    # Add scale bar
    plotter.add_scale_bar(num_mutations=1)

    # Add statistics
    plotter.add_statistics_annotation()

    # Save if filename provided
    if filename:
        plotter.save(filename)

    return fig, ax

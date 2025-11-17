"""
Interactive network visualization using Dash Cytoscape for PyPopART.

Provides Cytoscape-based interactive plotting for haplotype networks
with manual node repositioning, pie chart nodes, and legend support.
"""

import base64
import colorsys
import math
from typing import Dict, List, Optional, Tuple

import numpy as np

from ..core.graph import HaplotypeNetwork


class InteractiveCytoscapePlotter:
    """
    Interactive network plotter using Dash Cytoscape.

    Creates interactive visualizations of haplotype networks with
    manual node repositioning, pie chart nodes for population data,
    and customizable legends.
    """

    def __init__(self, network: HaplotypeNetwork):
        """
        Initialize interactive Cytoscape plotter with a haplotype network.

        Parameters
        ----------
        network :
            HaplotypeNetwork object to visualize.
        """
        self.network = network
        self.elements = None
        self.stylesheet = None

    @staticmethod
    def generate_pie_chart_svg(pie_data: List[Dict]) -> str:
        """
        Generate SVG pie chart as Data URI for node background.

        Parameters
        ----------
        pie_data :
            List of dicts with 'percent' and 'color' keys for each segment.

        Returns
        -------
            Data URI string with base64-encoded SVG.
        """
        size = 100
        center = size / 2
        radius = size / 2

        svg_parts = [
            f'<svg width="{size}" height="{size}" xmlns="http://www.w3.org/2000/svg">'
        ]

        # Start angle for pie slices (in radians)
        current_angle = -math.pi / 2  # Start at top (12 o'clock position)

        for segment in pie_data:
            percent = segment['percent']
            color = segment['color']

            # Calculate angle for this segment
            angle_size = (percent / 100) * 2 * math.pi

            # Calculate start and end points
            start_x = center + radius * math.cos(current_angle)
            start_y = center + radius * math.sin(current_angle)

            current_angle += angle_size

            end_x = center + radius * math.cos(current_angle)
            end_y = center + radius * math.sin(current_angle)

            # Use large-arc-flag if angle > 180 degrees
            large_arc = 1 if angle_size > math.pi else 0

            # Create pie slice path
            if percent == 100:
                # Full circle
                svg_parts.append(
                    f'<circle cx="{center}" cy="{center}" r="{radius}" fill="{color}"/>'
                )
            else:
                # Pie slice
                path = (
                    f'M {center},{center} '
                    f'L {start_x},{start_y} '
                    f'A {radius},{radius} 0 {large_arc},1 {end_x},{end_y} '
                    f'Z'
                )
                svg_parts.append(f'<path d="{path}" fill="{color}"/>')

        svg_parts.append('</svg>')
        svg_str = ''.join(svg_parts)

        # Encode as Data URI
        encoded = base64.b64encode(svg_str.encode('utf-8')).decode('utf-8')
        return f'data:image/svg+xml;base64,{encoded}'

    def create_elements(
        self,
        layout: Optional[Dict[str, Tuple[float, float]]] = None,
        node_size_scale: float = 20.0,
        population_colors: Optional[Dict[str, str]] = None,
        show_labels: bool = True,
        show_edge_labels: bool = True,
        median_vector_color: str = '#D3D3D3',
        node_labels: Optional[Dict[str, str]] = None,
    ) -> List[Dict]:
        """
        Create Cytoscape elements from network data.

        Parameters
        ----------
        layout :
            Pre-computed node positions {node_id: (x, y)}.
        node_size_scale :
            Scaling factor for node sizes.
        population_colors :
            Color mapping for populations {pop_name: color}.
        show_labels :
            Whether to show node labels.
        show_edge_labels :
            Whether to show edge labels with mutation counts.
        median_vector_color :
            Color for median vector nodes.
        node_labels :
            Custom labels for nodes {node_id: label}.

        Returns
        -------
            List of Cytoscape element dictionaries.
        """
        elements = []
        graph = self.network._graph

        # Create nodes
        for node in graph.nodes():
            hap = self.network.get_haplotype(node)
            is_median = self.network.is_median_vector(node)

            # Get position
            pos = layout.get(node, (0, 0)) if layout else (0, 0)

            # Calculate size
            if is_median:
                size = node_size_scale * 0.8
            elif hap:
                # Use sqrt of frequency, but ensure minimum size
                size = max(
                    node_size_scale * 0.5, node_size_scale * np.sqrt(hap.frequency)
                )
            else:
                size = node_size_scale * 0.5

            # Determine label
            if show_labels:
                if node_labels and node in node_labels:
                    label = node_labels[node]
                else:
                    label = node
            else:
                label = ''

            # Build node data
            node_data = {
                'id': node,
                'label': label,
                'size': size,
                'is_median': is_median,
            }

            # Add population pie chart data if available
            if not is_median and hap and population_colors:
                pop_counts = hap.get_frequency_by_population()
                if pop_counts and len(pop_counts) > 1:
                    # Multiple populations - prepare for pie chart display
                    total = sum(pop_counts.values())
                    pie_data = []
                    pie_colors = []
                    pie_sizes = []

                    for pop, count in sorted(pop_counts.items()):
                        if count > 0:
                            percent = (count / total) * 100
                            pie_data.append(
                                {
                                    'population': pop,
                                    'value': count,
                                    'percent': percent,
                                    'color': population_colors.get(pop, '#cccccc'),
                                }
                            )
                            pie_colors.append(population_colors.get(pop, '#cccccc'))
                            pie_sizes.append(percent)

                    # Store pie chart data for custom rendering
                    node_data['pie_data'] = pie_data
                    node_data['has_pie'] = True
                    node_data['pie_colors'] = pie_colors
                    node_data['pie_sizes'] = pie_sizes

                    # Generate SVG pie chart as Data URI
                    node_data['pie_svg'] = self.generate_pie_chart_svg(pie_data)

                    # Use transparent background to show pie chart
                    node_data['color'] = 'transparent'
                elif pop_counts:
                    # Single population
                    node_data['has_pie'] = False
                    pop = list(pop_counts.keys())[0]
                    node_data['color'] = population_colors.get(pop, '#87CEEB')
                else:
                    node_data['has_pie'] = False
                    node_data['color'] = '#87CEEB'
            else:
                node_data['has_pie'] = False
                if is_median:
                    node_data['color'] = median_vector_color
                else:
                    node_data['color'] = '#87CEEB'  # lightblue

            # Add hover information
            if is_median:
                node_data['hover'] = f'{node} (Median Vector)'
            elif hap:
                hover_lines = [f'{node}', f'Frequency: {hap.frequency}']
                pop_counts = hap.get_frequency_by_population()
                if pop_counts:
                    hover_lines.append('Populations:')
                    for pop, count in sorted(pop_counts.items()):
                        hover_lines.append(f'  {pop}: {count}')
                node_data['hover'] = '\n'.join(hover_lines)
            else:
                node_data['hover'] = node

            # Create element with position
            element = {
                'data': node_data,
                'position': {
                    'x': float(pos[0] * 100),
                    'y': float(pos[1] * 100),
                },  # Scale for visibility
                'grabbable': True,
            }

            elements.append(element)

        # Create edges
        for u, v in graph.edges():
            # Get mutation count from distance attribute
            distance = graph[u][v].get('distance', 1)
            # Get weight for layout purposes (should be uniform)
            weight = graph[u][v].get('weight', 1.0)

            edge_data = {
                'id': f'{u}-{v}',
                'source': u,
                'target': v,
                'distance': distance,  # Mutation count for labels and proportional layout
                'weight': weight,  # Uniform weight for standard layouts
                'label': str(int(distance))
                if show_edge_labels and distance > 0
                else '',
            }

            elements.append({'data': edge_data})

        self.elements = elements
        return elements

    def create_stylesheet(
        self,
        population_colors: Optional[Dict[str, str]] = None,
        median_vector_color: str = '#D3D3D3',
    ) -> List[Dict]:
        """
        Create Cytoscape stylesheet for network visualization.

        Parameters
        ----------
        population_colors :
            Color mapping for populations.
        median_vector_color :
            Color for median vector nodes.

        Returns
        -------
            List of stylesheet dictionaries.
        """
        stylesheet = [
            # Default node style - circular markers
            {
                'selector': 'node',
                'style': {
                    'content': 'data(label)',
                    'text-valign': 'center',
                    'text-halign': 'center',
                    'background-color': 'data(color)',
                    'width': 'data(size)',
                    'height': 'data(size)',
                    'shape': 'ellipse',
                    'border-width': 2,
                    'border-color': '#000000',
                    'font-size': '10px',
                    'font-weight': 'bold',
                    'text-outline-width': 2,
                    'text-outline-color': '#ffffff',
                },
            },
            # Median vector style - also circular but distinguished by color
            {
                'selector': 'node[is_median = true]',
                'style': {
                    'shape': 'ellipse',
                    'background-color': median_vector_color,
                },
            },
            # Default edge style
            {
                'selector': 'edge',
                'style': {
                    'width': 'mapData(weight, 1, 10, 3, 1)',
                    'line-color': '#969696',
                    'target-arrow-color': '#969696',
                    'curve-style': 'bezier',
                    'opacity': 0.6,
                },
            },
            # Edge label style
            {
                'selector': 'edge[label]',
                'style': {
                    'label': 'data(label)',
                    'font-size': '10px',
                    'text-background-color': '#ffffff',
                    'text-background-opacity': 0.7,
                    'text-background-padding': '3px',
                    'color': '#333333',
                },
            },
            # Highlighted/selected node style
            {
                'selector': 'node:selected',
                'style': {
                    'border-width': 4,
                    'border-color': '#ff0000',
                },
            },
        ]

        self.stylesheet = stylesheet
        return stylesheet

    def create_pie_stylesheet(self, population_colors: Dict[str, str]) -> List[Dict]:
        """
        Create stylesheet with pie chart support for nodes.

        Parameters
        ----------
        population_colors :
            Color mapping for populations.

        Returns
        -------
            List of stylesheet rules for pie chart nodes.
        """
        pie_styles = []

        # Style for pie chart nodes - use SVG background image
        pie_styles.append(
            {
                'selector': 'node[pie_svg]',
                'style': {
                    'background-color': 'transparent',
                    'background-image': 'data(pie_svg)',
                    'background-fit': 'contain',
                    'background-clip': 'node',
                    'border-width': 2,
                    'border-color': '#000000',
                },
            }
        )

        return pie_styles

    def generate_population_colors(self, populations: List[str]) -> Dict[str, str]:
        """
        Generate distinct colors for populations using HSV color space.

        Parameters
        ----------
        populations :
            List of population names.

        Returns
        -------
            Dictionary mapping population names to hex colors.
        """
        n = len(populations)
        colors = {}

        for i, pop in enumerate(sorted(populations)):
            # Generate evenly spaced hues
            hue = i / n
            # Use high saturation and value for vivid colors
            saturation = 0.7
            value = 0.9
            # Convert to RGB
            r, g, b = colorsys.hsv_to_rgb(hue, saturation, value)
            # Convert to hex
            hex_color = '#{:02x}{:02x}{:02x}'.format(
                int(r * 255), int(g * 255), int(b * 255)
            )
            colors[pop] = hex_color

        return colors


def create_cytoscape_network(
    network: HaplotypeNetwork,
    layout: Optional[Dict[str, Tuple[float, float]]] = None,
    population_colors: Optional[Dict[str, str]] = None,
    node_size_scale: float = 20.0,
    show_labels: bool = True,
    show_edge_labels: bool = True,
    median_vector_color: str = '#D3D3D3',
    node_labels: Optional[Dict[str, str]] = None,
) -> Tuple[List[Dict], List[Dict]]:
    """
    Create Cytoscape elements and stylesheet for a haplotype network.

    Parameters
    ----------
    network :
        HaplotypeNetwork object to visualize.
    layout :
        Pre-computed node positions {node_id: (x, y)}.
    population_colors :
        Color mapping for populations.
    node_size_scale :
        Scaling factor for node sizes.
    show_labels :
        Whether to show node labels.
    show_edge_labels :
        Whether to show edge labels with mutation counts.
    median_vector_color :
        Color for median vector nodes.
    node_labels :
        Custom labels for nodes {node_id: label}.

    Returns
    -------
        Tuple of (elements, stylesheet) for Cytoscape component.
    """
    plotter = InteractiveCytoscapePlotter(network)

    # Generate population colors if needed and populations exist
    if population_colors is None:
        # Check if any haplotypes have population data
        populations = set()
        for node in network._graph.nodes():
            if not network.is_median_vector(node):
                hap = network.get_haplotype(node)
                if hap:
                    pop_counts = hap.get_frequency_by_population()
                    if pop_counts:
                        populations.update(pop_counts.keys())

        if populations:
            population_colors = plotter.generate_population_colors(list(populations))

    elements = plotter.create_elements(
        layout=layout,
        node_size_scale=node_size_scale,
        population_colors=population_colors,
        show_labels=show_labels,
        show_edge_labels=show_edge_labels,
        median_vector_color=median_vector_color,
        node_labels=node_labels,
    )

    stylesheet = plotter.create_stylesheet(
        population_colors=population_colors,
        median_vector_color=median_vector_color,
    )

    # Add pie chart styles if we have population colors
    if population_colors:
        stylesheet.extend(plotter.create_pie_stylesheet(population_colors))

    return elements, stylesheet

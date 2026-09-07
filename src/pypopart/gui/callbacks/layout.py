"""Layout application, graph rendering, and node positioning."""

import traceback
from typing import Dict, List, Optional, Tuple

from dash import Input, Output, State, html
from dash.exceptions import PreventUpdate

from pypopart.core.graph import HaplotypeNetwork
from pypopart.layout.algorithms import LayoutManager
from pypopart.visualization.cytoscape_plot import (
    InteractiveCytoscapePlotter,
    create_cytoscape_network,
)


def register(app, logger) -> None:
    """
    Register layout callbacks on the Dash app.

    Parameters
    ----------
    app : dash.Dash
        The Dash application.
    logger : logging.Logger
        Application logger.
    """

    @app.callback(
        [
            Output('geographic-options', 'style'),
            Output('geographic-mode', 'data'),
        ],
        Input('layout-select', 'value'),
    )
    def toggle_geographic_options(layout: str) -> Tuple[Dict, bool]:
        """
        Show/hide geographic options based on layout selection.

        Parameters
        ----------
        layout : str
            Selected layout name.

        Returns
        -------
        Tuple[Dict, bool]
            Style for the geographic options panel, and whether geographic
            mode is active.
        """
        if layout == 'geographic':
            return {'display': 'block'}, True
        else:
            return {'display': 'none'}, False

    @app.callback(
        Output('layout-store', 'data'),
        [
            Input('apply-layout-button', 'n_clicks'),
            Input('network-store', 'data'),
            Input('spacing-slider', 'value'),
        ],
        [
            State('layout-select', 'value'),
            State('metadata-store', 'data'),
            State('map-projection', 'value'),
        ],
        prevent_initial_call=False,
    )
    def apply_layout(
        n_clicks: Optional[int],
        network_data: Optional[Dict],
        spacing_factor: float,
        layout_method: str,
        metadata_data: Optional[Dict],
        projection: str,
    ) -> Optional[Dict]:
        """
        Apply layout algorithm to network.

        Parameters
        ----------
        n_clicks : int, optional
            Button click count from Dash.
        network_data : Dict, optional
            Serialized network from the network store.
        spacing_factor : float
            Layout spacing multiplier.
        layout_method : str
            Selected layout algorithm name.
        metadata_data : Dict, optional
            Serialized metadata from the metadata store.
        projection : str
            Map projection name.

        Returns
        -------
        Optional[Dict]
            Node positions for the layout store, or None on failure.
        """
        if not network_data:
            return None

        try:
            # Reconstruct network
            network = HaplotypeNetwork.from_serialized(network_data)

            # Check if this is a proportional edge layout
            use_proportional_edges = layout_method in [
                'spring_proportional',
                'kamada_kawai_proportional',
            ]

            # Apply proportional edge length layouts
            if use_proportional_edges:
                import networkx as nx

                # Create a copy of the graph
                G = network.graph.copy()

                # Get all mutation distances to calculate relative scaling
                all_distances = [G[u][v].get('distance', 1) for u, v in G.edges()]
                min_dist = min(all_distances) if all_distances else 1
                max_dist = max(all_distances) if all_distances else 1

                logger.info(f'Edge distances range: {min_dist} to {max_dist} mutations')

                # Set edge length proportional to mutation count in NORMALIZED space
                for u, v in G.edges():
                    mutations = G[u][v].get('distance', 1)

                    # Normalize to 0.1-1.0 range based on min/max distances
                    # This keeps proportions but prevents huge layouts
                    if max_dist > min_dist:
                        normalized_length = 0.1 + 0.9 * (
                            (mutations - min_dist) / (max_dist - min_dist)
                        )
                    else:
                        normalized_length = 0.5  # All edges same length

                    G[u][v]['ideal_length'] = normalized_length

                logger.info(
                    f'Applying {layout_method} layout with proportional edge lengths'
                )

                if layout_method == 'spring_proportional':
                    # Spring layout with edge lengths
                    # NetworkX spring_layout uses 'weight' inversely
                    # We need to invert our ideal_length
                    for u, v in G.edges():
                        ideal_len = G[u][v]['ideal_length']
                        # Invert: short edges get high weight (attract more)
                        # Add small epsilon to avoid division by zero
                        G[u][v]['spring_weight'] = 1.0 / (ideal_len + 0.01)

                    positions = nx.spring_layout(
                        G,
                        weight='spring_weight',
                        k=None,  # Let NetworkX calculate optimal k
                        iterations=100,
                        scale=spacing_factor
                        * 1.0,  # Scale is applied to normalized coords
                        seed=42,  # For reproducibility
                    )

                elif layout_method == 'kamada_kawai_proportional':
                    # Kamada-Kawai layout respects edge distances directly
                    # Use ideal_length as the weight (desired distance)
                    positions = nx.kamada_kawai_layout(
                        G, weight='ideal_length', scale=spacing_factor * 1.0
                    )

            # Apply standard layouts (non-proportional)
            else:
                # Geographic layouts were removed with the geo feature;
                # fall back to spring for any stale 'geographic' value.
                layout_manager = LayoutManager(network)
                positions = layout_manager.compute_layout(
                    'spring' if layout_method == 'geographic' else layout_method
                )

                # Apply spacing factor to expand/contract the layout
                if spacing_factor and spacing_factor != 1.0:
                    positions = {
                        node: (pos[0] * spacing_factor, pos[1] * spacing_factor)
                        for node, pos in positions.items()
                    }

            # Convert to serializable format
            layout_data = {node: list(pos) for node, pos in positions.items()}

            return layout_data

        except Exception as e:
            logger.error(f'Error applying layout: {e}')
            logger.error(traceback.format_exc())
            return None

    @app.callback(
        [
            Output('network-graph', 'elements'),
            Output('network-graph', 'stylesheet'),
            Output('network-legend', 'children'),
        ],
        [
            Input('layout-store', 'data'),
            Input('network-store', 'data'),
            Input('geographic-mode', 'data'),
        ],
        [State('metadata-store', 'data'), State('h-number-mapping-store', 'data')],
    )
    def update_network_graph(
        layout_data: Optional[Dict],
        network_data: Optional[Dict],
        geographic_mode: bool,
        metadata_data: Optional[Dict],
        h_number_mapping: Optional[Dict],
    ) -> Tuple[List[Dict], List[Dict], html.Div]:
        """
        Update network visualization with Cytoscape.

        Parameters
        ----------
        layout_data : Dict, optional
            Node positions from the layout store.
        network_data : Dict, optional
            Serialized network from the network store.
        geographic_mode : bool
            Whether the geographic layout mode is active.
        metadata_data : Dict, optional
            Serialized metadata from the metadata store.
        h_number_mapping : Dict, optional
            Custom haplotype label mapping, if uploaded.

        Returns
        -------
        Tuple[List[Dict], List[Dict], html.Div]
            Cytoscape elements, the stylesheet, and the legend content.
        """
        if not network_data or not layout_data:
            # Return empty elements
            return [], [], html.Div('Upload data and compute network to visualize')

        try:
            # Reconstruct network
            network = HaplotypeNetwork.from_serialized(network_data)

            # Convert layout data
            positions = {node: tuple(pos) for node, pos in layout_data.items()}

            # Generate H number labels for nodes
            # Use custom mapping if available, otherwise generate default H numbers
            node_labels = {}
            if h_number_mapping:
                node_labels = h_number_mapping
            else:
                for i, node_id in enumerate(sorted(network.graph.nodes()), start=1):
                    node_labels[node_id] = f'H{i}'

            # Extract population colors from metadata if available
            population_colors = None
            if metadata_data and metadata_data.get('populations'):
                population_colors = metadata_data.get('population_colors', {})

            # If no colors provided, generate them
            if population_colors is None or not population_colors:
                plotter = InteractiveCytoscapePlotter(network)
                populations = set()
                for node in network._graph.nodes():
                    if not network.is_median_vector(node):
                        hap = network.get_haplotype(node)
                        if hap:
                            pop_counts = hap.get_frequency_by_population()
                            if pop_counts:
                                populations.update(pop_counts.keys())
                if populations:
                    population_colors = plotter.generate_population_colors(
                        list(populations)
                    )

            # Extract population mapping from metadata if available
            population_mapping = None
            if metadata_data and metadata_data.get('populations'):
                population_mapping = metadata_data['populations']

            # Create Cytoscape elements and stylesheet
            elements, stylesheet = create_cytoscape_network(
                network,
                layout=positions,
                population_colors=population_colors,
                population_mapping=population_mapping,
                show_labels=True,
                show_edge_labels=True,
                node_labels=node_labels,
            )

            # Add geographic styling if in geographic mode
            if geographic_mode and metadata_data and metadata_data.get('coordinates'):
                # Update stylesheet for geographic mode (add ocean background effect)
                pass

            # Create legend
            legend_items = []
            if population_colors:
                legend_items.append(
                    html.H6('Populations', style={'marginBottom': '5px'})
                )
                for pop, color in sorted(population_colors.items()):
                    legend_items.append(
                        html.Div(
                            [
                                html.Span(
                                    '●',
                                    style={
                                        'color': color,
                                        'fontSize': '20px',
                                        'marginRight': '5px',
                                    },
                                ),
                                html.Span(pop),
                            ],
                            style={'marginBottom': '3px'},
                        )
                    )

                # Add additional legend items
                legend_items.append(html.Hr(style={'margin': '5px 0'}))

                # Mixed population indicator (pie chart)
                legend_items.append(
                    html.Div(
                        [
                            html.Span(
                                '◕',
                                style={
                                    'color': '#000000',
                                    'fontSize': '20px',
                                    'marginRight': '5px',
                                },
                            ),
                            html.Span('Mixed Populations (Pie Chart)'),
                        ],
                        style={'marginBottom': '3px'},
                    )
                )

                # Median vector
                legend_items.append(
                    html.Div(
                        [
                            html.Span(
                                '■',
                                style={
                                    'color': '#D3D3D3',
                                    'fontSize': '20px',
                                    'marginRight': '5px',
                                },
                            ),
                            html.Span('Median Vector'),
                        ]
                    )
                )

            legend = html.Div(legend_items) if legend_items else html.Div()

            return elements, stylesheet, legend

        except Exception as e:
            logger.error(f'Error creating visualization: {e}')
            logger.error(traceback.format_exc())
            error_msg = html.Div(
                [
                    html.Strong('Error creating visualization'),
                    html.Br(),
                    str(e),
                ],
                style={'color': 'red'},
            )
            return [], [], error_msg

    @app.callback(
        [
            Output('layout-store', 'data', allow_duplicate=True),
            Output('manual-edit-flag', 'data', allow_duplicate=True),
        ],
        Input('network-graph', 'elements'),
        State('layout-store', 'data'),
        prevent_initial_call=True,
    )
    def update_node_positions(
        elements: Optional[List[Dict]],
        current_layout: Optional[Dict],
    ) -> Tuple[Optional[Dict], bool]:
        """
        Update node positions when user drags nodes in Cytoscape.

        Parameters
        ----------
        elements : List[Dict], optional
            Current Cytoscape elements.
        current_layout : Dict, optional
            Current layout positions.

        Returns
        -------
        Tuple[Optional[Dict], bool]
            Updated layout positions, and whether the update was applied.
        """
        if not elements or not current_layout:
            raise PreventUpdate

        try:
            updated_layout = current_layout.copy()
            position_changed = False

            # Extract positions from Cytoscape elements
            for element in elements:
                if 'position' in element and 'data' in element:
                    node_id = element['data'].get('id')
                    if node_id:
                        # Cytoscape positions are scaled by 100
                        x = element['position']['x'] / 100
                        y = element['position']['y'] / 100
                        # Check if position actually changed
                        if node_id in current_layout:
                            old_pos = current_layout[node_id]
                            if abs(old_pos[0] - x) > 0.01 or abs(old_pos[1] - y) > 0.01:
                                position_changed = True
                        updated_layout[node_id] = [x, y]

            # Set manual edit flag to True if positions changed
            return updated_layout, position_changed

        except Exception as e:
            logger.error(f'Error updating node positions: {e}')
            raise PreventUpdate from e

    @app.callback(
        Output('network-graph', 'stylesheet', allow_duplicate=True),
        [
            Input('node-size-slider', 'value'),
            Input('edge-width-slider', 'value'),
        ],
        State('network-graph', 'stylesheet'),
        prevent_initial_call=True,
    )
    def update_node_edge_sizes(
        node_size: int,
        edge_width: float,
        current_stylesheet: List[Dict],
    ) -> List[Dict]:
        """
        Update node size and edge width in stylesheet.

        Parameters
        ----------
        node_size : int
            Node size setting from the slider.
        edge_width : float
            Edge width setting from the slider.
        current_stylesheet : List[Dict]
            Current Cytoscape stylesheet.

        Returns
        -------
        List[Dict]
            The stylesheet with the new node and edge sizing.
        """
        if not current_stylesheet:
            raise PreventUpdate

        # Create a copy of the stylesheet
        new_stylesheet = []
        for style in current_stylesheet:
            new_style = style.copy()

            # Update node size - scale proportionally based on data(size)
            if style.get('selector') == 'node':
                if 'style' not in new_style:
                    new_style['style'] = {}
                new_style['style'] = {**new_style['style']}
                # Scale nodes proportionally: multiply data(size) by scale factor
                # Default node size is 40, so scale factor is node_size/40
                scale_factor = node_size / 40.0
                # Use mapData to scale the size attribute proportionally
                # This preserves the relative size differences between nodes
                new_style['style']['width'] = (
                    f'mapData(size, 0, 100, 0, {100 * scale_factor})'
                )
                new_style['style']['height'] = (
                    f'mapData(size, 0, 100, 0, {100 * scale_factor})'
                )

            # Update edge width
            elif style.get('selector') == 'edge':
                if 'style' not in new_style:
                    new_style['style'] = {}
                new_style['style'] = {**new_style['style']}
                new_style['style']['width'] = edge_width

            new_stylesheet.append(new_style)

        return new_stylesheet

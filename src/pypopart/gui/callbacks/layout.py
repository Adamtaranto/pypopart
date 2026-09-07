"""Layout application, graph rendering, and node positioning."""

import traceback
from typing import Dict, List, Optional, Tuple

from dash import Input, Output, State, html
from dash.exceptions import PreventUpdate

from pypopart.core.graph import HaplotypeNetwork
from pypopart.gui.serialization import merge_node_positions
from pypopart.layout.algorithms import LayoutManager, snap_to_grid
from pypopart.visualization.cytoscape_plot import (
    DEFAULT_TICK_THRESHOLD,
    InteractiveCytoscapePlotter,
    create_cytoscape_network,
    create_edge_tick_stylesheet,
)
from pypopart.visualization.style import DEFAULT_MEDIAN_COLOR

#: Node size the base stylesheet is authored against, used to turn the
#: slider value into a proportional scale factor.
BASE_NODE_SIZE = 40.0

#: Cytoscape renders stored layout coordinates multiplied by this, so a grid
#: expressed in Cytoscape pixels divides by it to reach stored units.
CYTOSCAPE_POSITION_SCALE = 100.0

#: How far a node must move, in stored units, to count as dragged. Below
#: this the change is float noise from the position scaling round trip.
POSITION_EPSILON = 0.01


def _grid_step(snap_enabled: Optional[bool], grid_size: Optional[float]) -> float:
    """
    Convert the grid controls into a step in stored layout units.

    Parameters
    ----------
    snap_enabled : bool, optional
        Whether the snap-to-grid switch is on.
    grid_size : float, optional
        Grid spacing in Cytoscape pixels.

    Returns
    -------
    float
        Grid step in stored units, or 0.0 when snapping is off.
    """
    if not snap_enabled or not grid_size:
        return 0.0
    return float(grid_size) / CYTOSCAPE_POSITION_SCALE


def _apply_size_overrides(
    stylesheet: List[Dict], node_size: Optional[float], edge_width: Optional[float]
) -> List[Dict]:
    """
    Re-apply the node size and edge width slider values to a stylesheet.

    The base stylesheet is regenerated whenever the figure is rebuilt, so
    without this the sliders would silently snap back to their defaults.

    Parameters
    ----------
    stylesheet : List[Dict]
        Cytoscape stylesheet to override.
    node_size : float, optional
        Node size from the slider; ``None`` leaves node sizing alone.
    edge_width : float, optional
        Edge width from the slider; ``None`` leaves edge width alone.

    Returns
    -------
    List[Dict]
        A new stylesheet with the slider values applied.
    """
    new_stylesheet = []
    for style in stylesheet:
        new_style = style.copy()

        if style.get('selector') == 'node' and node_size is not None:
            new_style['style'] = {**new_style.get('style', {})}
            # Scale proportionally via mapData so relative node sizes hold.
            scale_factor = node_size / BASE_NODE_SIZE
            mapping = f'mapData(size, 0, 100, 0, {100 * scale_factor})'
            new_style['style']['width'] = mapping
            new_style['style']['height'] = mapping

        elif style.get('selector') == 'edge' and edge_width is not None:
            new_style['style'] = {**new_style.get('style', {})}
            new_style['style']['width'] = edge_width

        new_stylesheet.append(new_style)

    return new_stylesheet


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
        [
            Output('layout-store', 'data'),
            # A freshly computed layout supersedes any manual drags, and the
            # view should re-fit to it.
            Output('node-positions-store', 'data', allow_duplicate=True),
            Output('manual-edit-flag', 'data', allow_duplicate=True),
        ],
        [
            Input('apply-layout-button', 'n_clicks'),
            Input('network-store', 'data'),
            Input('spacing-slider', 'value'),
        ],
        [
            State('layout-select', 'value'),
            State('metadata-store', 'data'),
            State('map-projection', 'value'),
            State('snap-to-grid-toggle', 'value'),
            State('grid-size', 'value'),
        ],
        prevent_initial_call='initial_duplicate',
    )
    def apply_layout(
        n_clicks: Optional[int],
        network_data: Optional[Dict],
        spacing_factor: float,
        layout_method: str,
        metadata_data: Optional[Dict],
        projection: str,
        snap_enabled: Optional[bool],
        grid_size: Optional[float],
    ) -> Tuple[Optional[Dict], None, bool]:
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
        snap_enabled : bool, optional
            Whether computed positions are quantised to the grid.
        grid_size : float, optional
            Grid spacing in Cytoscape pixels.

        Returns
        -------
        Tuple[Optional[Dict], None, bool]
            Node positions for the layout store, a cleared drag store, and
            a cleared manual-edit flag.
        """
        if not network_data:
            return None, None, False

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

            # Snap after the spacing multiply, so the grid is what the user
            # sees rather than what the algorithm happened to produce.
            positions = snap_to_grid(positions, _grid_step(snap_enabled, grid_size))

            # Convert to serializable format
            layout_data = {node: list(pos) for node, pos in positions.items()}

            return layout_data, None, False

        except Exception as e:
            logger.error(f'Error applying layout: {e}')
            logger.error(traceback.format_exc())
            return None, None, False

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
        [
            State('metadata-store', 'data'),
            State('h-number-mapping-store', 'data'),
            State('node-size-slider', 'value'),
            State('edge-width-slider', 'value'),
            State('edge-tick-toggle', 'value'),
            State('edge-tick-threshold', 'value'),
            State('node-positions-store', 'data'),
        ],
    )
    def update_network_graph(
        layout_data: Optional[Dict],
        network_data: Optional[Dict],
        geographic_mode: bool,
        metadata_data: Optional[Dict],
        h_number_mapping: Optional[Dict],
        node_size: Optional[float],
        edge_width: Optional[float],
        show_edge_ticks: Optional[bool],
        edge_tick_threshold: Optional[int],
        dragged_positions: Optional[Dict],
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
        node_size : float, optional
            Current node size slider value, re-applied to the new stylesheet.
        edge_width : float, optional
            Current edge width slider value, re-applied to the new stylesheet.
        show_edge_ticks : bool, optional
            Whether mutation counts render as tick marks.
        edge_tick_threshold : int, optional
            Mutation count above which an edge shows a numeral.
        dragged_positions : Dict, optional
            Manually dragged node positions, which win over the layout.

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

            # Manual drags win over the computed layout.
            positions = merge_node_positions(layout_data, dragged_positions)

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
                show_edge_ticks=bool(show_edge_ticks),
                edge_tick_threshold=edge_tick_threshold or DEFAULT_TICK_THRESHOLD,
            )

            # The base stylesheet is fresh, so the sliders have to be re-applied.
            stylesheet = _apply_size_overrides(stylesheet, node_size, edge_width)

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
                                    'color': DEFAULT_MEDIAN_COLOR,
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
                className='pp-error',
            )
            return [], [], error_msg

    @app.callback(
        [
            # Deliberately NOT layout-store. That store is an Input of
            # update_network_graph, so writing drags there made every drag
            # rebuild all elements, which re-fit the view and pushed the
            # position through a lossy /100 -> *100 float round trip.
            Output('node-positions-store', 'data', allow_duplicate=True),
            Output('manual-edit-flag', 'data', allow_duplicate=True),
        ],
        Input('network-graph', 'elements'),
        [
            State('node-positions-store', 'data'),
            State('layout-store', 'data'),
            State('snap-to-grid-toggle', 'value'),
            State('grid-size', 'value'),
        ],
        prevent_initial_call=True,
    )
    def update_node_positions(
        elements: Optional[List[Dict]],
        current_positions: Optional[Dict],
        current_layout: Optional[Dict],
        snap_enabled: Optional[bool],
        grid_size: Optional[float],
    ) -> Tuple[Optional[Dict], bool]:
        """
        Persist node positions when the user drags nodes in Cytoscape.

        Snapping is applied here as well as in the browser: the clientside
        handler moves the node so the user sees it land on the grid, and
        this guarantees the *stored* position is quantised even if that
        reposition does not make it back into the elements prop.

        Parameters
        ----------
        elements : List[Dict], optional
            Current Cytoscape elements.
        current_positions : Dict, optional
            Previously stored manual positions.
        current_layout : Dict, optional
            Positions from the computed layout, used as the baseline for
            deciding whether anything actually moved.
        snap_enabled : bool, optional
            Whether dragged positions are quantised to the grid.
        grid_size : float, optional
            Grid spacing in Cytoscape pixels.

        Returns
        -------
        Tuple[Optional[Dict], bool]
            The stored positions, and whether a node actually moved.
        """
        if not elements or not current_layout:
            raise PreventUpdate

        try:
            step = _grid_step(snap_enabled, grid_size)
            baseline = merge_node_positions(current_layout, current_positions)

            dragged = {}
            for element in elements:
                if 'position' in element and 'data' in element:
                    node_id = element['data'].get('id')
                    if node_id:
                        # Cytoscape positions are scaled by 100
                        dragged[node_id] = (
                            element['position']['x'] / CYTOSCAPE_POSITION_SCALE,
                            element['position']['y'] / CYTOSCAPE_POSITION_SCALE,
                        )

            # Idempotent, so re-running it on already-snapped positions is a
            # no-op -- this callback fires on every elements change, not only
            # at the end of a drag.
            dragged = snap_to_grid(dragged, step)

            updated = dict(current_positions or {})
            position_changed = False
            for node_id, (x, y) in dragged.items():
                old_pos = baseline.get(node_id)
                if old_pos is not None and (
                    abs(old_pos[0] - x) > POSITION_EPSILON
                    or abs(old_pos[1] - y) > POSITION_EPSILON
                ):
                    position_changed = True
                # Rounded so repeated round trips cannot accumulate noise.
                updated[node_id] = [round(x, 6), round(y, 6)]

            if not position_changed:
                raise PreventUpdate

            return updated, True

        except PreventUpdate:
            raise
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

        return _apply_size_overrides(current_stylesheet, node_size, edge_width)

    @app.callback(
        Output('network-graph', 'stylesheet', allow_duplicate=True),
        [
            Input('edge-tick-toggle', 'value'),
            Input('edge-tick-threshold', 'value'),
        ],
        State('network-graph', 'stylesheet'),
        prevent_initial_call=True,
    )
    def update_edge_tick_style(
        show_ticks: bool,
        threshold: int,
        current_stylesheet: List[Dict],
    ) -> List[Dict]:
        """
        Swap the edge mutation-count rules between tick marks and numerals.

        Tick strings are always present in the element data, so this is a
        stylesheet-only update: no re-layout and no flicker.

        Parameters
        ----------
        show_ticks : bool
            Whether to draw tick marks.
        threshold : int
            Mutation count above which an edge shows a numeral.
        current_stylesheet : List[Dict]
            Current Cytoscape stylesheet.

        Returns
        -------
        List[Dict]
            The stylesheet with regenerated edge mutation-count rules.
        """
        if not current_stylesheet:
            raise PreventUpdate

        # Every generated rule is prefixed 'edge[distance', so the whole
        # group can be dropped and rebuilt in one pass. Order matters:
        # these must stay last to override the catch-all edge[label] rule.
        kept = [
            style
            for style in current_stylesheet
            if not str(style.get('selector', '')).startswith('edge[distance')
        ]
        return kept + create_edge_tick_stylesheet(
            threshold or DEFAULT_TICK_THRESHOLD, bool(show_ticks)
        )

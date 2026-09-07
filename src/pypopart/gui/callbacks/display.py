"""Statistics, alignment, haplotype-summary, metadata and tooltip displays."""

import functools
import json
import traceback
from typing import Dict, List, Optional, Tuple

from dash import Input, Output, State, html
import dash_bootstrap_components as dbc

from pypopart.core.graph import HaplotypeNetwork
from pypopart.stats import (
    calculate_diversity_metrics,
    calculate_network_metrics,
    identify_central_haplotypes,
)


def register(app, logger) -> None:
    """
    Register display callbacks on the Dash app.

    Parameters
    ----------
    app : dash.Dash
        The Dash application.
    logger : logging.Logger
        Application logger.
    """

    @app.callback(
        Output('statistics-display', 'children'), Input('network-store', 'data')
    )
    def update_statistics(network_data: Optional[Dict]) -> html.Div:
        """
        Update statistics display (cached per network payload).

        Parameters
        ----------
        network_data : Dict, optional
            Serialized network from the network store.

        Returns
        -------
        html.Div
            The rendered statistics panel.
        """
        if not network_data:
            return html.Div(
                'Compute a network to see statistics',
                style={'color': 'gray', 'padding': '20px'},
            )
        return _statistics_for(json.dumps(network_data, sort_keys=True))

    @functools.lru_cache(maxsize=4)
    def _statistics_for(network_json: str) -> html.Div:
        """
        Build the statistics view for one serialized network.

        Parameters
        ----------
        network_json : str
            Serialized network payload used as the cache key.

        Returns
        -------
        html.Div
            The rendered statistics panel.
        """
        network_data = json.loads(network_json)

        try:
            # Reconstruct network
            network = HaplotypeNetwork.from_serialized(network_data)

            # Calculate statistics
            network_metrics = calculate_network_metrics(network)
            diversity_metrics = calculate_diversity_metrics(network)
            central_haps = identify_central_haplotypes(network)

            # Create statistics display
            return html.Div(
                [
                    html.H4('Network Statistics'),
                    html.Hr(),
                    html.H5('Basic Metrics'),
                    html.Ul(
                        [
                            html.Li(f'Number of Nodes: {len(network.graph.nodes)}'),
                            html.Li(f'Number of Edges: {len(network.graph.edges)}'),
                            html.Li(f'Network Diameter: {network_metrics.diameter}'),
                            html.Li(
                                f'Average Clustering Coefficient: '
                                f'{network_metrics.clustering_coefficient:.3f}'
                            ),
                            html.Li(
                                f'Reticulation Index: '
                                f'{network_metrics.reticulation_index:.3f}'
                            ),
                        ]
                    ),
                    html.Hr(),
                    html.H5('Diversity Metrics'),
                    html.Ul(
                        [
                            html.Li(
                                f'Haplotype Diversity: '
                                f'{diversity_metrics.haplotype_diversity:.3f}'
                            ),
                            html.Li(
                                f'Shannon Index: {diversity_metrics.shannon_index:.3f}'
                            ),
                        ]
                    ),
                    html.Hr(),
                    html.H5('Central Haplotypes'),
                    _format_central_haplotypes(central_haps),
                ]
            )

        except Exception as e:
            logger.error(f'Error calculating statistics: {e}')
            logger.error(traceback.format_exc())
            return html.Div(
                [
                    html.H5('Error calculating statistics', style={'color': 'red'}),
                    html.P(str(e)),
                    html.Details(
                        [
                            html.Summary('Show traceback'),
                            html.Pre(
                                traceback.format_exc(),
                                style={
                                    'background': '#f5f5f5',
                                    'padding': '10px',
                                    'overflow': 'auto',
                                    'font-size': '12px',
                                },
                            ),
                        ]
                    ),
                ],
                style={'color': 'red', 'padding': '20px'},
            )

    @app.callback(
        Output('alignment-display', 'children'), Input('alignment-store', 'data')
    )
    def update_alignment_display(alignment_data: Optional[Dict]):
        """
        Update alignment viewer with colored nucleotides for polymorphic sites.

        Parameters
        ----------
        alignment_data : Dict, optional
            Serialized alignment from the alignment store.

        Returns
        -------
        html.Div or str
            The rendered alignment view, or a prompt when no data is
            loaded.
        """
        if not alignment_data:
            return 'Upload data to view alignment'

        try:
            # Standard nucleotide color scheme
            nuc_colors = {
                'A': '#64F73F',  # Green
                'a': '#64F73F',
                'C': '#3C88EE',  # Blue
                'c': '#3C88EE',
                'G': '#FFB340',  # Orange/Yellow
                'g': '#FFB340',
                'T': '#EB413E',  # Red
                't': '#EB413E',
                'U': '#EB413E',  # Red (for RNA)
                'u': '#EB413E',
                '-': '#CCCCCC',  # Gray for gaps
                'N': '#999999',  # Dark gray for N
                'n': '#999999',
            }

            sequences = alignment_data['sequences'][:50]  # Limit to first 50
            if not sequences:
                return 'No sequences to display'

            max_id_len = max(len(seq['id']) for seq in sequences)
            seq_length = len(sequences[0]['data'])

            # Identify polymorphic sites (positions with >1 base type)
            polymorphic_sites = set()
            for pos in range(seq_length):
                bases_at_pos = set()
                for seq in sequences:
                    if pos < len(seq['data']):
                        base = seq['data'][pos].upper()
                        bases_at_pos.add(base)
                # Position is polymorphic if it has more than one unique base
                if len(bases_at_pos) > 1:
                    polymorphic_sites.add(pos)

            # Build HTML rows
            rows = []
            for seq in sequences:
                seq_id = seq['id']
                seq_data = seq['data']

                # Create sequence ID span
                id_span = html.Span(
                    f'{seq_id:<{max_id_len}}  ',
                    style={'color': 'black', 'fontWeight': 'bold'},
                )

                # Create colored nucleotide spans
                seq_spans = []
                for pos, base in enumerate(seq_data):
                    if pos in polymorphic_sites:
                        # Color polymorphic positions
                        color = nuc_colors.get(base, '#000000')
                        seq_spans.append(
                            html.Span(
                                base,
                                style={
                                    'backgroundColor': color,
                                    'color': 'white'
                                    if base.upper() not in ['-', 'N']
                                    else 'black',
                                    'padding': '0 1px',
                                },
                            )
                        )
                    else:
                        # Keep invariant positions as plain text
                        seq_spans.append(html.Span(base, style={'color': 'black'}))

                # Combine ID and sequence
                rows.append(
                    html.Div([id_span] + seq_spans, style={'whiteSpace': 'pre'})
                )

            # Add note if sequences were truncated
            if len(alignment_data['sequences']) > 50:
                rows.append(
                    html.Div(
                        f'\n... ({len(alignment_data["sequences"]) - 50} more sequences)',
                        style={
                            'color': 'gray',
                            'fontStyle': 'italic',
                            'marginTop': '10px',
                        },
                    )
                )

            return html.Div(rows)

        except Exception as e:
            logger.error(f'Error displaying alignment: {e}')
            logger.error(traceback.format_exc())
            return f'Error displaying alignment: {str(e)}'

    @app.callback(
        Output('haplotype-summary-display', 'children'),
        [
            Input('network-store', 'data'),
            Input('alignment-store', 'data'),
            Input('metadata-store', 'data'),
            Input('h-number-mapping-store', 'data'),
        ],
    )
    def update_haplotype_summary(
        network_data: Optional[Dict],
        alignment_data: Optional[Dict],
        metadata_data: Optional[Dict],
        h_number_mapping: Optional[Dict],
    ) -> html.Div:
        """
        Update haplotype summary showing H number to sequence name mapping.

        Parameters
        ----------
        network_data : Dict, optional
            Serialized network from the network store.
        alignment_data : Dict, optional
            Serialized alignment from the alignment store.
        metadata_data : Dict, optional
            Serialized metadata from the metadata store.
        h_number_mapping : Dict, optional
            Custom haplotype label mapping, if uploaded.

        Returns
        -------
        html.Div
            The rendered haplotype summary table.
        """
        if not network_data or not alignment_data:
            return html.Div('Compute network to view haplotype summary')

        try:
            # Reconstruct network
            network = HaplotypeNetwork.from_serialized(network_data)

            # Check if we have population data
            has_populations = metadata_data and metadata_data.get('populations')

            # Create mapping of H numbers to sequence names
            haplotype_mapping = []
            for i, node_id in enumerate(sorted(network.graph.nodes()), start=1):
                node_data = network.graph.nodes[node_id]
                sample_ids = node_data.get('sample_ids', [])
                is_median = node_data.get('median_vector', False)
                frequency = node_data.get('frequency', len(sample_ids))

                # Use custom label if available, otherwise default H number
                if h_number_mapping and node_id in h_number_mapping:
                    h_label = h_number_mapping[node_id]
                else:
                    h_label = f'H{i}'

                # Determine if this is an inferred haplotype
                if is_median or len(sample_ids) == 0:
                    haplotype_type = '🔵 Inferred'
                    sample_display = 'None (inferred ancestral haplotype)'
                    populations_display = ''
                else:
                    haplotype_type = '🟢 Observed'
                    sample_display = ', '.join(sample_ids) if sample_ids else 'Unknown'

                    # Collect populations for this haplotype
                    if has_populations:
                        populations = set()
                        for sid in sample_ids:
                            if sid in metadata_data['populations']:
                                populations.add(metadata_data['populations'][sid])
                        populations_display = (
                            ', '.join(sorted(populations)) if populations else ''
                        )
                    else:
                        populations_display = ''

                haplotype_mapping.append(
                    {
                        'h_label': h_label,
                        'node_id': node_id,
                        'type': haplotype_type,
                        'frequency': frequency,
                        'samples': sample_display,
                        'populations': populations_display,
                    }
                )

            # Create table headers (conditionally include populations)
            headers = [
                html.Th('H Number'),
                html.Th('Type'),
                html.Th('Frequency'),
                html.Th('Sample IDs'),
            ]
            if has_populations:
                headers.append(html.Th('Populations'))

            table_header = [html.Thead(html.Tr(headers))]

            # Create table rows
            table_rows = []
            for hap in haplotype_mapping:
                row_cells = [
                    html.Td(hap['h_label'], style={'fontWeight': 'bold'}),
                    html.Td(hap['type']),
                    html.Td(hap['frequency']),
                    html.Td(
                        hap['samples'],
                        style={
                            'maxWidth': '600px',
                            'overflow': 'auto',
                            'whiteSpace': 'normal',
                        },
                    ),
                ]
                if has_populations:
                    row_cells.append(html.Td(hap['populations']))
                table_rows.append(html.Tr(row_cells))

            table_body = [html.Tbody(table_rows)]

            # Count statistics
            n_observed = sum(1 for h in haplotype_mapping if '🟢' in h['type'])
            n_inferred = sum(1 for h in haplotype_mapping if '🔵' in h['type'])

            return html.Div(
                [
                    html.H4('Haplotype Summary'),
                    html.P(
                        [
                            f'Total haplotypes: {len(haplotype_mapping)} ',
                            f'(🟢 {n_observed} observed, 🔵 {n_inferred} inferred)',
                        ]
                    ),
                    html.Hr(),
                    dbc.Table(
                        table_header + table_body,
                        bordered=True,
                        hover=True,
                        responsive=True,
                        striped=True,
                        style={'fontSize': '14px'},
                    ),
                ]
            )

        except Exception as e:
            logger.error(f'Error creating haplotype summary: {e}')
            logger.error(traceback.format_exc())
            return html.Div(
                [
                    html.H5('Error creating haplotype summary', style={'color': 'red'}),
                    html.P(str(e)),
                    html.Details(
                        [
                            html.Summary('Show traceback'),
                            html.Pre(
                                traceback.format_exc(),
                                style={
                                    'background': '#f5f5f5',
                                    'padding': '10px',
                                    'overflow': 'auto',
                                    'font-size': '12px',
                                },
                            ),
                        ]
                    ),
                ],
                style={'color': 'red', 'padding': '20px'},
            )

    @app.callback(
        [
            Output('haplotype-search', 'options'),
            Output('network-graph', 'stylesheet', allow_duplicate=True),
        ],
        [Input('network-store', 'data'), Input('haplotype-search', 'value')],
        [
            State('network-graph', 'stylesheet'),
            State('h-number-mapping-store', 'data'),
        ],
        prevent_initial_call=True,
    )
    def update_search_and_highlight(
        network_data: Optional[Dict],
        selected_h_list: Optional[List[str]],
        current_stylesheet: List[Dict],
        h_number_mapping: Optional[Dict],
    ) -> Tuple[List[Dict], List[Dict]]:
        """
        Update search dropdown options and highlight selected nodes.

        Parameters
        ----------
        network_data : Dict, optional
            Serialized network from the network store.
        selected_h_list : List[str], optional
            Haplotype labels selected in the search box.
        current_stylesheet : List[Dict]
            Current Cytoscape stylesheet.
        h_number_mapping : Dict, optional
            Custom haplotype label mapping, if uploaded.

        Returns
        -------
        Tuple[List[Dict], List[Dict]]
            Dropdown options for the haplotype search, and the stylesheet
            highlighting the selection.
        """
        if not network_data:
            return [], current_stylesheet or []

        try:
            # Reconstruct network
            network = HaplotypeNetwork.from_serialized(network_data)

            # Build H number options - mapping H numbers to node IDs
            h_numbers = []
            h_to_node = {}
            for i, node_id in enumerate(sorted(network.graph.nodes()), start=1):
                # Use custom label if available, otherwise default H number
                if h_number_mapping and node_id in h_number_mapping:
                    h_label = h_number_mapping[node_id]
                else:
                    h_label = f'H{i}'
                h_numbers.append({'label': h_label, 'value': h_label})
                h_to_node[h_label] = node_id

            # Create a clean stylesheet without any highlight styles
            if not current_stylesheet:
                current_stylesheet = []

            # Remove ALL existing search highlight styles (for specific nodes)
            # But preserve the :selected pseudo-selector style for click selection
            base_stylesheet = [
                s
                for s in current_stylesheet
                if not (
                    s.get('selector', '').startswith('node[id = "')
                    and 'border-color' in s.get('style', {})
                    and s.get('style', {}).get('border-color') == '#FF0000'
                    and s.get('style', {}).get('border-width')
                    == '5px'  # Only remove search highlights (5px)
                )
            ]

            # Ensure the :selected pseudo-selector styles are always present for click selection
            has_selected_style = any(
                s.get('selector') == 'node:selected' for s in base_stylesheet
            )
            has_pie_selected_style = any(
                s.get('selector') == 'node[pie_svg]:selected' for s in base_stylesheet
            )

            if not has_selected_style:
                # Re-add the :selected style if it was somehow removed
                base_stylesheet.append(
                    {
                        'selector': 'node:selected',
                        'style': {
                            'border-width': 4,
                            'border-color': '#ff0000',
                            'z-index': 999,
                        },
                    }
                )

            if not has_pie_selected_style:
                # Re-add the pie chart selected style if it was somehow removed
                # This is more specific than node:selected and overrides node[pie_svg] border
                base_stylesheet.append(
                    {
                        'selector': 'node[pie_svg]:selected',
                        'style': {
                            'border-width': 4,
                            'border-color': '#ff0000',
                            'z-index': 999,
                        },
                    }
                )

            # If nodes are selected, add highlight styles
            if selected_h_list:
                # Ensure it's a list (might be single value in some cases)
                if not isinstance(selected_h_list, list):
                    selected_h_list = [selected_h_list]

                # Add highlight style for each selected H number
                for selected_h in selected_h_list:
                    # Map H number to node ID
                    node_id = h_to_node.get(selected_h)
                    if node_id:
                        base_stylesheet.append(
                            {
                                'selector': f'node[id = "{node_id}"]',
                                'style': {
                                    'border-width': '5px',
                                    'border-color': '#FF0000',
                                    'border-style': 'solid',
                                },
                            }
                        )

                return h_numbers, base_stylesheet

            # No selection - return clean stylesheet
            return h_numbers, base_stylesheet

        except Exception as e:
            logger.error(f'Error updating search: {e}')
            return [], current_stylesheet or []

    @app.callback(
        Output('metadata-display', 'children'),
        Output('metadata-warnings', 'children'),
        [Input('alignment-store', 'data'), Input('metadata-store', 'data')],
    )
    def update_metadata_tab(
        alignment_data: Optional[Dict],
        metadata_data: Optional[Dict],
    ) -> Tuple[html.Div, html.Div]:
        """
        Display metadata with all sequence IDs.

        Parameters
        ----------
        alignment_data : Dict, optional
            Serialized alignment from the alignment store.
        metadata_data : Dict, optional
            Serialized metadata from the metadata store.

        Returns
        -------
        Tuple[html.Div, html.Div]
            The metadata table and its summary panel.
        """
        if not alignment_data:
            return html.Div('Upload alignment to view metadata'), html.Div()

        try:
            # Get sequence IDs from alignment
            alignment_ids = {seq['id'] for seq in alignment_data['sequences']}

            # Get metadata IDs if available
            metadata_ids = set()
            metadata_records = {}
            if metadata_data:
                metadata_ids = set(metadata_data.get('sequence_ids', []))
                # Build metadata records
                for sid in metadata_data.get('sequence_ids', []):
                    metadata_records[sid] = {
                        'population': metadata_data.get('populations', {}).get(sid, ''),
                        'latitude': metadata_data.get('coordinates', {})
                        .get(sid, {})
                        .get('lat', ''),
                        'longitude': metadata_data.get('coordinates', {})
                        .get(sid, {})
                        .get('lon', ''),
                    }

            # Union of all IDs
            all_ids = alignment_ids.union(metadata_ids)

            # Check for duplicates in alignment
            alignment_id_list = [seq['id'] for seq in alignment_data['sequences']]
            alignment_duplicates = [
                sid
                for sid in set(alignment_id_list)
                if alignment_id_list.count(sid) > 1
            ]

            # Check for mismatches
            only_in_alignment = alignment_ids - metadata_ids
            only_in_metadata = metadata_ids - alignment_ids

            # Build warnings
            warnings = []
            if alignment_duplicates:
                warnings.append(
                    dbc.Alert(
                        f'⚠️ Duplicate IDs found in alignment: {", ".join(alignment_duplicates)}',
                        color='warning',
                    )
                )
            if only_in_alignment and metadata_data:
                warnings.append(
                    dbc.Alert(
                        f'⚠️ {len(only_in_alignment)} IDs only in alignment (not in metadata)',
                        color='info',
                    )
                )
            if only_in_metadata:
                warnings.append(
                    dbc.Alert(
                        f'⚠️ {len(only_in_metadata)} IDs only in metadata (not in alignment)',
                        color='info',
                    )
                )

            # Get population colors if available
            population_colors = (
                metadata_data.get('population_colors', {}) if metadata_data else {}
            )

            # Build table
            table_header = [
                html.Thead(
                    html.Tr(
                        [
                            html.Th('Sequence ID'),
                            html.Th('In Alignment'),
                            html.Th('In Metadata'),
                            html.Th('Population'),
                            html.Th('Color'),
                            html.Th('Latitude'),
                            html.Th('Longitude'),
                        ]
                    )
                )
            ]

            table_rows = []
            for sid in sorted(all_ids):
                in_alignment = '✓' if sid in alignment_ids else '✗'
                in_metadata = '✓' if sid in metadata_ids else '✗'

                meta = metadata_records.get(sid, {})
                pop = meta.get('population', '')

                # Get color for this population
                color_display = ''
                if pop and population_colors and pop in population_colors:
                    color_hex = population_colors[pop]
                    color_display = html.Div(
                        [
                            html.Span(
                                '●',
                                style={
                                    'color': color_hex,
                                    'fontSize': '16px',
                                    'marginRight': '5px',
                                },
                            ),
                            html.Span(color_hex, style={'fontSize': '12px'}),
                        ]
                    )

                table_rows.append(
                    html.Tr(
                        [
                            html.Td(sid),
                            html.Td(in_alignment, style={'textAlign': 'center'}),
                            html.Td(in_metadata, style={'textAlign': 'center'}),
                            html.Td(pop),
                            html.Td(color_display),
                            html.Td(meta.get('latitude', '')),
                            html.Td(meta.get('longitude', '')),
                        ]
                    )
                )

            table_body = [html.Tbody(table_rows)]

            table = dbc.Table(
                table_header + table_body,
                bordered=True,
                hover=True,
                responsive=True,
                striped=True,
                style={'fontSize': '14px'},
            )

            return table, html.Div(warnings)

        except Exception as e:
            logger.error(f'Error creating metadata display: {e}')
            logger.error(traceback.format_exc())
            return html.Div(f'Error: {str(e)}'), html.Div()

    @app.callback(
        Output('node-tooltip', 'children'),
        [
            Input('network-graph', 'mouseoverNodeData'),
            Input('network-graph', 'mouseoverEdgeData'),
        ],
        [State('network-store', 'data'), State('h-number-mapping-store', 'data')],
    )
    def update_tooltip_content(
        hover_data: Optional[Dict],
        edge_hover_data: Optional[Dict],
        network_data: Optional[Dict],
        h_number_mapping: Optional[Dict],
    ) -> html.Div:
        """
        Update tooltip content based on hovered node.

        Parameters
        ----------
        hover_data : Dict, optional
            Hovered node data from Cytoscape.
        edge_hover_data : Dict, optional
            Hovered edge data from Cytoscape.
        network_data : Dict, optional
            Serialized network from the network store.
        h_number_mapping : Dict, optional
            Custom haplotype label mapping, if uploaded.

        Returns
        -------
        html.Div
            The tooltip content for the hovered node or edge.
        """
        # Hide tooltip if hovering over edge instead of node
        if edge_hover_data and not hover_data:
            return html.Div()

        if not hover_data or not network_data:
            return html.Div()

        try:
            # Reconstruct network
            network = HaplotypeNetwork.from_serialized(network_data)

            # Get node data
            node_id = hover_data.get('id')
            if not node_id:
                return html.Div()

            node_data = network.graph.nodes.get(node_id, {})
            sample_ids = node_data.get('sample_ids', [])
            is_median = node_data.get('median_vector', False)

            # Find H number for this node
            if h_number_mapping and node_id in h_number_mapping:
                h_label = h_number_mapping[node_id]
            else:
                h_label = None
                for i, nid in enumerate(sorted(network.graph.nodes()), start=1):
                    if nid == node_id:
                        h_label = f'H{i}'
                        break

            # Build tooltip content
            if is_median or len(sample_ids) == 0:
                content = html.Div(
                    [
                        html.Strong(h_label or 'Unknown'),
                        html.Br(),
                        html.Span('Inferred median vector'),
                    ]
                )
            else:
                content = html.Div(
                    [
                        html.Strong(h_label or 'Unknown'),
                        html.Br(),
                        html.Span(f'Sequences ({len(sample_ids)}):'),
                        html.Br(),
                        html.Span(
                            ', '.join(sample_ids[:10])
                            + ('...' if len(sample_ids) > 10 else '')
                        ),
                    ]
                )

            return content

        except Exception as e:
            logger.error(f'Error showing tooltip: {e}')
            return html.Div()

    # Use clientside callback for tooltip positioning
    # This gets the actual rendered position from Cytoscape
    app.clientside_callback(
        """
        function() {
            // Set up event listeners once when page loads
            if (window.tooltipSetup) {
                return window.dash_clientside.no_update;
            }
            window.tooltipSetup = true;

            setTimeout(function() {
                try {
                    const cy = document.getElementById('network-graph')._cyreg.cy;
                    const tooltip = document.getElementById('node-tooltip');

                    if (!cy || !tooltip) {
                        console.log('Could not find cytoscape or tooltip element');
                        return;
                    }

                    // Hide tooltip on mouseover edge or background
                    cy.on('mouseover', 'edge', function(evt) {
                        tooltip.style.display = 'none';
                    });

                    cy.on('mouseover', function(evt) {
                        // If target is cy (background), hide tooltip
                        if (evt.target === cy) {
                            tooltip.style.display = 'none';
                        }
                    });

                    // Show tooltip on node hover
                    cy.on('mouseover', 'node', function(evt) {
                        const node = evt.target;
                        const renderedPos = node.renderedPosition();

                        tooltip.style.display = 'block';
                        tooltip.style.left = (renderedPos.x + 15) + 'px';
                        tooltip.style.top = (renderedPos.y - 40) + 'px';
                    });

                    // Hide tooltip when mouse leaves node
                    cy.on('mouseout', 'node', function(evt) {
                        tooltip.style.display = 'none';
                    });

                    console.log('Tooltip event listeners installed');

                } catch (e) {
                    console.log('Error setting up tooltip:', e);
                }
            }, 500);

            return window.dash_clientside.no_update;
        }
        """,
        Output('node-tooltip', 'style', allow_duplicate=True),
        Input('network-graph', 'elements'),
        prevent_initial_call=True,
    )


def _format_central_haplotypes(central: Dict) -> html.Div:
    """
    Format central haplotypes for display.

    Parameters
    ----------
    central : Dict
        Central haplotype metrics.

    Returns
    -------
    html.Div
        The rendered list of central haplotypes.
    """
    try:
        return html.Ul(
            [
                html.Li(
                    f'Degree Centrality: {central["degree_centrality"]} '
                    f'(degree: {central["degree"]})'
                ),
                html.Li(f'Betweenness Centrality: {central["betweenness_centrality"]}'),
                html.Li(f'Closeness Centrality: {central["closeness_centrality"]}'),
            ]
        )
    except Exception:
        return html.Div('Unable to identify central haplotypes')

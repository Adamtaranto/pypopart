"""
Dash-based GUI for PyPopART haplotype network analysis.

This module provides a web-based graphical user interface for PyPopART,
allowing users to upload sequence data, configure network algorithms,
visualize results, and export outputs.
"""

import base64
import logging
import tempfile
import traceback
from typing import Dict, List, Optional, Tuple

import dash
from dash import Input, Output, State, dcc, html
from dash.exceptions import PreventUpdate
import dash_bootstrap_components as dbc
import networkx as nx
import plotly.graph_objects as go

from pypopart.algorithms import (
    TCS,
    MedianJoiningNetwork,
    MinimumSpanningNetwork,
    MinimumSpanningTree,
)
from pypopart.core.alignment import Alignment
from pypopart.core.graph import HaplotypeNetwork
from pypopart.io import FastaReader, NexusReader, PhylipReader
from pypopart.io.network_export import GraphMLExporter, JSONExporter
from pypopart.layout.algorithms import LayoutManager
from pypopart.stats import (
    calculate_diversity_metrics,
    calculate_network_metrics,
    identify_central_haplotypes,
)
from pypopart.visualization.interactive_plot import InteractiveNetworkPlotter


class PyPopARTApp:
    """
    Main PyPopART Dash application class.

    Provides a web-based interface for haplotype network analysis.
    """

    def __init__(self, debug: bool = False, port: int = 8050):
        """
        Initialize PyPopART Dash application.

        Parameters
        ----------
        debug : bool, default=False
            Enable debug mode for development
        port : int, default=8050
            Port number for the web server
        """
        self.debug = debug
        self.port = port

        # Configure logging
        log_level = logging.DEBUG if debug else logging.INFO
        logging.basicConfig(
            level=log_level,
            format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        )
        self.logger = logging.getLogger(__name__)

        self.app = dash.Dash(
            __name__,
            external_stylesheets=[dbc.themes.BOOTSTRAP],
            suppress_callback_exceptions=True,
        )
        self.app.title = 'PyPopART - Haplotype Network Analysis'

        # Storage for application state
        self.alignment = None
        self.network = None
        self.layout_positions = None

        self._setup_layout()
        self._setup_callbacks()

    def _setup_layout(self) -> None:
        """Set up the application layout with all components."""
        self.app.layout = dbc.Container(
            [
                dbc.Row(
                    [
                        dbc.Col(
                            html.H1(
                                'PyPopART: Haplotype Network Analysis',
                                className='text-center mb-4',
                            ),
                            width=12,
                        )
                    ]
                ),
                dbc.Row(
                    [
                        # Left panel - Controls
                        dbc.Col(
                            [
                                self._create_upload_card(),
                                html.Br(),
                                self._create_algorithm_card(),
                                html.Br(),
                                self._create_layout_card(),
                                html.Br(),
                                self._create_export_card(),
                            ],
                            width=3,
                            style={'height': '90vh', 'overflow-y': 'auto'},
                        ),
                        # Right panel - Visualization
                        dbc.Col(
                            [
                                dbc.Tabs(
                                    [
                                        dbc.Tab(
                                            self._create_network_tab(),
                                            label='Network',
                                        ),
                                        dbc.Tab(
                                            self._create_statistics_tab(),
                                            label='Statistics',
                                        ),
                                        dbc.Tab(
                                            self._create_alignment_tab(),
                                            label='Alignment',
                                        ),
                                    ]
                                )
                            ],
                            width=9,
                        ),
                    ]
                ),
                # Hidden stores for data
                dcc.Store(id='alignment-store'),
                dcc.Store(id='network-store'),
                dcc.Store(id='layout-store'),
                dcc.Store(id='computation-status'),
            ],
            fluid=True,
            style={'padding': '20px'},
        )

    def _create_upload_card(self) -> dbc.Card:
        """Create file upload card."""
        return dbc.Card(
            [
                dbc.CardHeader(html.H5('1. Upload Data')),
                dbc.CardBody(
                    [
                        dcc.Upload(
                            id='upload-data',
                            children=dbc.Button(
                                'Select File', color='primary', className='w-100'
                            ),
                            multiple=False,
                        ),
                        html.Div(id='upload-status', className='mt-2'),
                        html.Hr(),
                        html.Small(
                            'Supported formats: FASTA, NEXUS, PHYLIP',
                            className='text-muted',
                        ),
                    ]
                ),
            ]
        )

    def _create_algorithm_card(self) -> dbc.Card:
        """Create algorithm selection and parameter card."""
        return dbc.Card(
            [
                dbc.CardHeader(html.H5('2. Configure Algorithm')),
                dbc.CardBody(
                    [
                        dbc.Label('Algorithm'),
                        dcc.Dropdown(
                            id='algorithm-select',
                            options=[
                                {
                                    'label': 'Minimum Spanning Tree (MST)',
                                    'value': 'mst',
                                },
                                {
                                    'label': 'Minimum Spanning Network (MSN)',
                                    'value': 'msn',
                                },
                                {
                                    'label': 'TCS (Statistical Parsimony)',
                                    'value': 'tcs',
                                },
                                {
                                    'label': 'Median-Joining Network (MJN)',
                                    'value': 'mjn',
                                },
                            ],
                            value='msn',
                        ),
                        html.Br(),
                        html.Div(id='algorithm-parameters'),
                        html.Br(),
                        dbc.Button(
                            'Compute Network',
                            id='compute-button',
                            color='success',
                            className='w-100',
                            disabled=True,
                        ),
                        html.Div(id='computation-feedback', className='mt-2'),
                    ]
                ),
            ]
        )

    def _create_layout_card(self) -> dbc.Card:
        """Create layout configuration card."""
        return dbc.Card(
            [
                dbc.CardHeader(html.H5('3. Layout Options')),
                dbc.CardBody(
                    [
                        dbc.Label('Layout Algorithm'),
                        dcc.Dropdown(
                            id='layout-select',
                            options=[
                                {'label': 'Spring (Force-Directed)', 'value': 'spring'},
                                {'label': 'Circular', 'value': 'circular'},
                                {'label': 'Radial', 'value': 'radial'},
                                {'label': 'Hierarchical', 'value': 'hierarchical'},
                                {'label': 'Kamada-Kawai', 'value': 'kamada_kawai'},
                            ],
                            value='spring',
                        ),
                        html.Br(),
                        dbc.Checklist(
                            id='snap-to-grid',
                            options=[{'label': ' Snap to Grid', 'value': 'snap'}],
                            value=[],
                            inline=True,
                        ),
                        html.Br(),
                        dbc.Button(
                            'Apply Layout',
                            id='apply-layout-button',
                            color='info',
                            className='w-100',
                            disabled=True,
                        ),
                    ]
                ),
            ]
        )

    def _create_export_card(self) -> dbc.Card:
        """Create export options card."""
        return dbc.Card(
            [
                dbc.CardHeader(html.H5('4. Export')),
                dbc.CardBody(
                    [
                        dbc.Label('Format'),
                        dcc.Dropdown(
                            id='export-format',
                            options=[
                                {'label': 'GraphML', 'value': 'graphml'},
                                {'label': 'GML', 'value': 'gml'},
                                {'label': 'JSON', 'value': 'json'},
                                {'label': 'PNG Image', 'value': 'png'},
                                {'label': 'SVG Image', 'value': 'svg'},
                            ],
                            value='graphml',
                        ),
                        html.Br(),
                        dbc.Button(
                            'Download',
                            id='export-button',
                            color='secondary',
                            className='w-100',
                            disabled=True,
                        ),
                        dcc.Download(id='download-data'),
                    ]
                ),
            ]
        )

    def _create_network_tab(self) -> html.Div:
        """Create network visualization tab."""
        return html.Div(
            [
                dcc.Loading(
                    id='loading-network',
                    type='default',
                    children=[
                        dcc.Graph(
                            id='network-graph',
                            style={'height': '85vh'},
                            config={'displayModeBar': True, 'displaylogo': False},
                        )
                    ],
                )
            ]
        )

    def _create_statistics_tab(self) -> html.Div:
        """Create statistics display tab."""
        return html.Div(
            [
                html.Div(
                    id='statistics-display',
                    style={'padding': '20px', 'height': '85vh', 'overflow-y': 'auto'},
                )
            ]
        )

    def _create_alignment_tab(self) -> html.Div:
        """Create alignment viewer tab."""
        return html.Div(
            [
                html.Div(
                    id='alignment-display',
                    style={
                        'padding': '20px',
                        'height': '85vh',
                        'overflow': 'auto',
                        'font-family': 'monospace',
                        'white-space': 'pre',
                    },
                )
            ]
        )

    def _setup_callbacks(self) -> None:
        """Set up all Dash callbacks for interactivity."""

        @self.app.callback(
            [
                Output('upload-status', 'children'),
                Output('alignment-store', 'data'),
                Output('compute-button', 'disabled'),
            ],
            Input('upload-data', 'contents'),
            State('upload-data', 'filename'),
        )
        def handle_file_upload(
            contents: Optional[str], filename: Optional[str]
        ) -> Tuple[html.Div, Optional[Dict], bool]:
            """Handle file upload and parse alignment."""
            if contents is None:
                return html.Div(), None, True

            try:
                content_type, content_string = contents.split(',')
                decoded = base64.b64decode(content_string)

                # Determine file format and parse
                alignment = None
                if filename.endswith('.fasta') or filename.endswith('.fa'):
                    with tempfile.NamedTemporaryFile(
                        mode='wb', suffix='.fasta', delete=False
                    ) as tmp:
                        tmp.write(decoded)
                        tmp.flush()
                        reader = FastaReader(tmp.name)
                        alignment = reader.read_alignment()
                elif filename.endswith('.nex') or filename.endswith('.nexus'):
                    with tempfile.NamedTemporaryFile(
                        mode='wb', suffix='.nexus', delete=False
                    ) as tmp:
                        tmp.write(decoded)
                        tmp.flush()
                        reader = NexusReader(tmp.name)
                        alignment = reader.read_alignment()
                elif filename.endswith('.phy') or filename.endswith('.phylip'):
                    with tempfile.NamedTemporaryFile(
                        mode='wb', suffix='.phylip', delete=False
                    ) as tmp:
                        tmp.write(decoded)
                        tmp.flush()
                        reader = PhylipReader(tmp.name)
                        alignment = reader.read_alignment()
                else:
                    return (
                        dbc.Alert(
                            'Unsupported file format. Use FASTA, NEXUS, or PHYLIP.',
                            color='danger',
                        ),
                        None,
                        True,
                    )

                # Store alignment data
                alignment_data = {
                    'sequences': [
                        {
                            'id': seq.id,
                            'data': seq.data,
                            'metadata': seq.metadata,
                            'description': seq.description or '',
                        }
                        for seq in alignment
                    ],
                    'length': alignment.length,
                    'num_sequences': len(alignment),
                }

                status = dbc.Alert(
                    [
                        html.Strong('Success! '),
                        f'Loaded {len(alignment)} sequences '
                        f'of length {alignment.length}',
                    ],
                    color='success',
                )

                return status, alignment_data, False

            except Exception as e:
                return (
                    dbc.Alert(f'Error parsing file: {str(e)}', color='danger'),
                    None,
                    True,
                )

        @self.app.callback(
            Output('algorithm-parameters', 'children'),
            Input('algorithm-select', 'value'),
        )
        def update_algorithm_parameters(algorithm: str) -> html.Div:
            """Update parameter controls based on selected algorithm."""
            if algorithm == 'mst':
                return html.Div(
                    [
                        dbc.Label('Distance Metric'),
                        dcc.Dropdown(
                            id={'type': 'algorithm-param', 'name': 'distance'},
                            options=[
                                {'label': 'Hamming', 'value': 'hamming'},
                                {'label': 'Jukes-Cantor', 'value': 'jc'},
                                {'label': 'Kimura 2-parameter', 'value': 'k2p'},
                            ],
                            value='hamming',
                        ),
                    ]
                )
            elif algorithm == 'msn':
                return html.Div(
                    [
                        dbc.Label('Distance Metric'),
                        dcc.Dropdown(
                            id={'type': 'algorithm-param', 'name': 'distance'},
                            options=[
                                {'label': 'Hamming', 'value': 'hamming'},
                                {'label': 'Jukes-Cantor', 'value': 'jc'},
                                {'label': 'Kimura 2-parameter', 'value': 'k2p'},
                            ],
                            value='hamming',
                        ),
                    ]
                )
            elif algorithm == 'tcs':
                return html.Div(
                    [
                        dbc.Label('Connection Limit'),
                        dcc.Slider(
                            id={'type': 'algorithm-param', 'name': 'connection_limit'},
                            min=1,
                            max=20,
                            step=1,
                            value=10,
                            marks={i: str(i) for i in range(1, 21, 2)},
                        ),
                        html.Small(
                            'Maximum mutations between connected haplotypes',
                            className='text-muted',
                        ),
                    ]
                )
            elif algorithm == 'mjn':
                return html.Div(
                    [
                        dbc.Label('Epsilon'),
                        dcc.Input(
                            id={'type': 'algorithm-param', 'name': 'epsilon'},
                            type='number',
                            value=0,
                            min=0,
                            step=1,
                            className='form-control',
                        ),
                        html.Br(),
                        html.Small(
                            'Parameter for median vector inference (0 = automatic)',
                            className='text-muted',
                        ),
                    ]
                )
            return html.Div()

        @self.app.callback(
            [
                Output('network-store', 'data'),
                Output('computation-feedback', 'children'),
                Output('apply-layout-button', 'disabled'),
                Output('export-button', 'disabled'),
            ],
            Input('compute-button', 'n_clicks'),
            [
                State('alignment-store', 'data'),
                State('algorithm-select', 'value'),
                State({'type': 'algorithm-param', 'name': dash.ALL}, 'value'),
            ],
            prevent_initial_call=True,
        )
        def compute_network(
            n_clicks: int,
            alignment_data: Dict,
            algorithm: str,
            param_values: List,
        ) -> Tuple[Optional[Dict], html.Div, bool, bool]:
            """Compute haplotype network using selected algorithm."""
            if not alignment_data:
                raise PreventUpdate

            try:
                # Reconstruct alignment from stored data
                from pypopart.core.sequence import Sequence

                sequences = [
                    Sequence(
                        id=seq['id'],
                        data=seq['data'],
                        metadata=seq['metadata'],
                        description=seq['description'],
                    )
                    for seq in alignment_data['sequences']
                ]
                alignment = Alignment(sequences)

                # Extract parameter values - they come in as a list
                # The pattern-matching callback returns values in order
                param_value = param_values[0] if param_values else None

                # Select and configure algorithm
                if algorithm == 'mst':
                    algo = MinimumSpanningTree(distance_metric=param_value or 'hamming')
                elif algorithm == 'msn':
                    algo = MinimumSpanningNetwork(
                        distance_metric=param_value or 'hamming'
                    )
                elif algorithm == 'tcs':
                    algo = TCS(connection_limit=param_value or 10)
                elif algorithm == 'mjn':
                    algo = MedianJoiningNetwork(epsilon=param_value or 0)
                else:
                    raise ValueError(f'Unknown algorithm: {algorithm}')

                # Build network
                network = algo.build_network(alignment)

                # Convert to serializable format
                network_data = {
                    'nodes': [
                        {
                            'id': node,
                            'sequence': network.graph.nodes[node].get('sequence', ''),
                            'frequency': network.graph.nodes[node].get('frequency', 1),
                            'is_median': network.graph.nodes[node].get(
                                'is_median', False
                            ),
                            'samples': network.graph.nodes[node].get('samples', []),
                        }
                        for node in network.graph.nodes()
                    ],
                    'edges': [
                        {
                            'source': u,
                            'target': v,
                            'weight': network.graph[u][v].get('weight', 1),
                        }
                        for u, v in network.graph.edges()
                    ],
                }

                feedback = dbc.Alert(
                    [
                        html.Strong('Network computed! '),
                        f'{len(network.graph.nodes)} nodes, '
                        f'{len(network.graph.edges)} edges',
                    ],
                    color='success',
                )

                return network_data, feedback, False, False

            except Exception as e:
                logging.error(f'Error computing network: {e}')
                logging.error(traceback.format_exc())
                return (
                    None,
                    dbc.Alert(f'Error computing network: {str(e)}', color='danger'),
                    True,
                    True,
                )

        @self.app.callback(
            Output('layout-store', 'data'),
            [Input('apply-layout-button', 'n_clicks'), Input('network-store', 'data')],
            [State('layout-select', 'value'), State('snap-to-grid', 'value')],
            prevent_initial_call=False,
        )
        def apply_layout(
            n_clicks: Optional[int],
            network_data: Optional[Dict],
            layout_method: str,
            snap_to_grid: List[str],
        ) -> Optional[Dict]:
            """Apply layout algorithm to network."""
            if not network_data:
                return None

            try:
                # Reconstruct network graph
                G = nx.Graph()
                for node in network_data['nodes']:
                    G.add_node(
                        node['id'],
                        sequence=node['sequence'],
                        frequency=node['frequency'],
                        is_median=node['is_median'],
                    )
                for edge in network_data['edges']:
                    G.add_edge(edge['source'], edge['target'], weight=edge['weight'])

                # Apply layout
                layout_manager = LayoutManager()
                if layout_method == 'spring':
                    positions = layout_manager.spring_layout(G)
                elif layout_method == 'circular':
                    positions = layout_manager.circular_layout(G)
                elif layout_method == 'radial':
                    positions = layout_manager.radial_layout(G)
                elif layout_method == 'hierarchical':
                    positions = layout_manager.hierarchical_layout(G)
                elif layout_method == 'kamada_kawai':
                    positions = layout_manager.kamada_kawai_layout(G)
                else:
                    positions = layout_manager.spring_layout(G)

                # Snap to grid if requested
                if 'snap' in snap_to_grid:
                    grid_size = 0.1
                    positions = {
                        node: (
                            round(pos[0] / grid_size) * grid_size,
                            round(pos[1] / grid_size) * grid_size,
                        )
                        for node, pos in positions.items()
                    }

                # Convert to serializable format
                layout_data = {node: list(pos) for node, pos in positions.items()}

                return layout_data

            except Exception as e:
                logging.error(f'Error applying layout: {e}')
                logging.error(traceback.format_exc())
                return None

        @self.app.callback(
            Output('network-graph', 'figure'),
            [Input('layout-store', 'data'), Input('network-store', 'data')],
        )
        def update_network_graph(
            layout_data: Optional[Dict], network_data: Optional[Dict]
        ) -> go.Figure:
            """Update network visualization."""
            if not network_data or not layout_data:
                # Return empty figure with instructions
                fig = go.Figure()
                fig.add_annotation(
                    text='Upload data and compute network to visualize',
                    xref='paper',
                    yref='paper',
                    x=0.5,
                    y=0.5,
                    showarrow=False,
                    font={'size': 16, 'color': 'gray'},
                )
                fig.update_layout(
                    xaxis={'visible': False},
                    yaxis={'visible': False},
                    plot_bgcolor='white',
                )
                return fig

            try:
                # Reconstruct network
                G = nx.Graph()
                for node in network_data['nodes']:
                    G.add_node(
                        node['id'],
                        sequence=node['sequence'],
                        frequency=node['frequency'],
                        is_median=node['is_median'],
                        samples=node.get('samples', []),
                    )
                for edge in network_data['edges']:
                    G.add_edge(edge['source'], edge['target'], weight=edge['weight'])

                network = HaplotypeNetwork()
                network.graph = G

                # Convert layout data
                positions = {node: tuple(pos) for node, pos in layout_data.items()}

                # Create interactive plot
                plotter = InteractiveNetworkPlotter()
                fig = plotter.plot_network(network, positions=positions)

                return fig

            except Exception as e:
                fig = go.Figure()
                fig.add_annotation(
                    text=f'Error creating visualization: {str(e)}',
                    xref='paper',
                    yref='paper',
                    x=0.5,
                    y=0.5,
                    showarrow=False,
                    font={'size': 14, 'color': 'red'},
                )
                return fig

        @self.app.callback(
            Output('statistics-display', 'children'), Input('network-store', 'data')
        )
        def update_statistics(network_data: Optional[Dict]) -> html.Div:
            """Update statistics display."""
            if not network_data:
                return html.Div(
                    'Compute a network to see statistics',
                    style={'color': 'gray', 'padding': '20px'},
                )

            try:
                # Reconstruct network
                G = nx.Graph()
                for node in network_data['nodes']:
                    G.add_node(
                        node['id'],
                        sequence=node['sequence'],
                        frequency=node['frequency'],
                        is_median=node['is_median'],
                    )
                for edge in network_data['edges']:
                    G.add_edge(edge['source'], edge['target'], weight=edge['weight'])

                network = HaplotypeNetwork()
                network.graph = G

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
                                html.Li(
                                    f'Network Diameter: {network_metrics.diameter}'
                                ),
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
                        self._format_central_haplotypes(central_haps),
                    ]
                )

            except Exception as e:
                return html.Div(
                    f'Error calculating statistics: {str(e)}', style={'color': 'red'}
                )

        @self.app.callback(
            Output('alignment-display', 'children'), Input('alignment-store', 'data')
        )
        def update_alignment_display(alignment_data: Optional[Dict]) -> str:
            """Update alignment viewer."""
            if not alignment_data:
                return 'Upload data to view alignment'

            try:
                # Format alignment for display
                sequences = alignment_data['sequences'][:50]  # Limit to first 50
                max_id_len = max(len(seq['id']) for seq in sequences)

                lines = []
                for seq in sequences:
                    # Format: ID (padded) + sequence
                    lines.append(f'{seq["id"]:<{max_id_len}}  {seq["data"]}')

                if len(alignment_data['sequences']) > 50:
                    lines.append(
                        f'\n... ({len(alignment_data["sequences"]) - 50} more sequences)'
                    )

                return '\n'.join(lines)

            except Exception as e:
                return f'Error displaying alignment: {str(e)}'

        @self.app.callback(
            Output('download-data', 'data'),
            Input('export-button', 'n_clicks'),
            [
                State('network-store', 'data'),
                State('export-format', 'value'),
                State('network-graph', 'figure'),
            ],
            prevent_initial_call=True,
        )
        def export_network(
            n_clicks: int, network_data: Dict, export_format: str, figure: Dict
        ) -> Dict:
            """Export network in selected format."""
            if not network_data:
                raise PreventUpdate

            try:
                # Reconstruct network
                G = nx.Graph()
                for node in network_data['nodes']:
                    G.add_node(
                        node['id'],
                        sequence=node['sequence'],
                        frequency=node['frequency'],
                        is_median=node['is_median'],
                    )
                for edge in network_data['edges']:
                    G.add_edge(edge['source'], edge['target'], weight=edge['weight'])

                network = HaplotypeNetwork()
                network.graph = G

                # Export based on format
                if export_format == 'graphml':
                    with tempfile.NamedTemporaryFile(
                        mode='w', suffix='.graphml', delete=False
                    ) as tmp:
                        exporter = GraphMLExporter(tmp.name)
                        exporter.export(network)
                        with open(tmp.name, 'r') as f:
                            content = f.read()
                    return {
                        'content': content,
                        'filename': 'network.graphml',
                        'type': 'text/xml',
                    }

                elif export_format == 'gml':
                    content = '\n'.join(nx.generate_gml(G))
                    return {
                        'content': content,
                        'filename': 'network.gml',
                        'type': 'text/plain',
                    }

                elif export_format == 'json':
                    with tempfile.NamedTemporaryFile(
                        mode='w', suffix='.json', delete=False
                    ) as tmp:
                        exporter = JSONExporter(tmp.name)
                        exporter.export(network)
                        with open(tmp.name, 'r') as f:
                            content = f.read()
                    return {
                        'content': content,
                        'filename': 'network.json',
                        'type': 'application/json',
                    }

                elif export_format in ['png', 'svg']:
                    # Export figure as image
                    fig = go.Figure(figure)
                    img_bytes = fig.to_image(format=export_format)
                    return {
                        'content': base64.b64encode(img_bytes).decode(),
                        'filename': f'network.{export_format}',
                        'base64': True,
                    }

            except Exception as e:
                logging.error(f'Error exporting: {e}')
                logging.error(traceback.format_exc())
                raise PreventUpdate from None

    def _format_central_haplotypes(self, central: Dict) -> html.Div:
        """Format central haplotypes for display."""
        try:
            return html.Ul(
                [
                    html.Li(
                        f'Degree Centrality: {central["degree_centrality"]} '
                        f'(degree: {central["degree"]})'
                    ),
                    html.Li(
                        f'Betweenness Centrality: {central["betweenness_centrality"]}'
                    ),
                    html.Li(f'Closeness Centrality: {central["closeness_centrality"]}'),
                ]
            )
        except Exception:
            return html.Div('Unable to identify central haplotypes')

    def run(self) -> None:
        """Run the Dash application."""
        # Use run() for Dash 2.0+, which replaced run_server()
        self.app.run(debug=self.debug, port=self.port)


def main(debug: bool = False, port: int = 8050) -> None:
    """
    Launch the PyPopART GUI application.

    Parameters
    ----------
    debug : bool, default=False
        Enable debug mode
    port : int, default=8050
        Port number for web server
    """
    app = PyPopARTApp(debug=debug, port=port)
    print(f'Starting PyPopART GUI on http://localhost:{port}')
    app.run()


if __name__ == '__main__':
    main(debug=True)

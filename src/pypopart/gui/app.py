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
from pypopart.io.metadata import MetadataReader, extract_coordinates
from pypopart.io.network_export import GraphMLExporter, JSONExporter
from pypopart.layout.algorithms import GeographicLayout, LayoutManager
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
        debug :
            bool, default=False.
            Enable debug mode for development.
        port :
            int, default=8050.
            Port number for the web server.
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
        self.app.layout = html.Div(
            [
                # Header
                html.Div(
                    html.H1(
                        'PyPopART: Haplotype Network Analysis',
                        className='text-center mb-4',
                    ),
                    style={'padding': '20px', 'backgroundColor': 'white'},
                ),
                # Main resizable container
                html.Div(
                    [
                        # Left panel - Controls (resizable sidebar)
                        html.Div(
                            [
                                self._create_upload_card(),
                                html.Br(),
                                self._create_algorithm_card(),
                                html.Br(),
                                self._create_layout_card(),
                                html.Br(),
                                self._create_export_card(),
                            ],
                            style={
                                'minWidth': '250px',
                                'width': '300px',
                                'maxWidth': '600px',
                                'height': '90vh',
                                'overflowY': 'auto',
                                'padding': '20px',
                                'backgroundColor': '#f8f9fa',
                                'resize': 'horizontal',
                                'overflow': 'auto',
                            },
                        ),
                        # Right panel - Visualization
                        html.Div(
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
                                            self._create_haplotype_summary_tab(),
                                            label='Haplotype Summary',
                                        ),
                                        dbc.Tab(
                                            self._create_alignment_tab(),
                                            label='Alignment',
                                        ),
                                    ]
                                )
                            ],
                            style={
                                'flex': '1',
                                'padding': '20px',
                                'overflow': 'auto',
                            },
                        ),
                    ],
                    style={
                        'display': 'flex',
                        'height': '90vh',
                        'overflow': 'hidden',
                    },
                ),
                # Hidden stores for data
                dcc.Store(id='alignment-store'),
                dcc.Store(id='metadata-store'),
                dcc.Store(id='network-store'),
                dcc.Store(id='layout-store'),
                dcc.Store(id='computation-status'),
                dcc.Store(id='geographic-mode', data=False),
            ]
        )

    def _create_upload_card(self) -> dbc.Card:
        """Create file upload card."""
        return dbc.Card(
            [
                dbc.CardHeader(
                    html.H5('1. Upload Data', className='mb-0'),
                ),
                dbc.CardBody(
                    [
                        dbc.Label('Sequence File', className='fw-bold'),
                        html.Small(
                            'Upload aligned sequences in FASTA, NEXUS, or PHYLIP format',
                            className='text-muted d-block mb-2',
                        ),
                        dcc.Upload(
                            id='upload-data',
                            children=dbc.Button(
                                '📁 Select Sequence File',
                                color='primary',
                                className='w-100',
                            ),
                            multiple=False,
                        ),
                        html.Div(id='upload-status', className='mt-2'),
                        html.Hr(),
                        dbc.Label('Metadata File (Optional)', className='fw-bold'),
                        html.Small(
                            'CSV file with population, location, or trait data',
                            className='text-muted d-block mb-2',
                        ),
                        dcc.Upload(
                            id='upload-metadata',
                            children=dbc.Button(
                                '📊 Select Metadata File',
                                color='secondary',
                                outline=True,
                                className='w-100',
                            ),
                            multiple=False,
                        ),
                        html.Div(id='metadata-status', className='mt-2'),
                        html.Div(
                            id='metadata-template-section',
                            children=[
                                html.Hr(),
                                dbc.Button(
                                    '⬇️ Download Metadata Template',
                                    id='download-template-button',
                                    color='info',
                                    outline=True,
                                    size='sm',
                                    className='w-100',
                                    disabled=True,
                                ),
                                dcc.Download(id='download-template'),
                                html.Small(
                                    'Get a CSV template pre-filled with your sequence IDs',
                                    className='text-muted d-block mt-1',
                                ),
                            ],
                        ),
                    ]
                ),
            ]
        )

    def _create_algorithm_card(self) -> dbc.Card:
        """Create algorithm selection and parameter card."""
        return dbc.Card(
            [
                dbc.CardHeader(
                    html.H5('2. Configure Algorithm', className='mb-0'),
                ),
                dbc.CardBody(
                    [
                        dbc.Label('Network Algorithm', className='fw-bold'),
                        html.Small(
                            'Choose the method for constructing the haplotype network',
                            className='text-muted d-block mb-2',
                        ),
                        dcc.Dropdown(
                            id='algorithm-select',
                            options=[
                                {
                                    'label': 'MST - Minimum Spanning Tree',
                                    'value': 'mst',
                                },
                                {
                                    'label': 'MSN - Minimum Spanning Network',
                                    'value': 'msn',
                                },
                                {
                                    'label': 'TCS - Statistical Parsimony',
                                    'value': 'tcs',
                                },
                                {
                                    'label': 'MJN - Median-Joining Network',
                                    'value': 'mjn',
                                },
                            ],
                            value='msn',
                            style={'whiteSpace': 'nowrap'},
                        ),
                        html.Br(),
                        html.Div(id='algorithm-parameters'),
                        html.Br(),
                        dbc.Button(
                            '⚡ Compute Network',
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
                dbc.CardHeader(
                    html.H5('3. Layout Options', className='mb-0'),
                ),
                dbc.CardBody(
                    [
                        dbc.Label('Layout Algorithm', className='fw-bold'),
                        html.Small(
                            'Choose how to position nodes in the visualization',
                            className='text-muted d-block mb-2',
                        ),
                        dcc.Dropdown(
                            id='layout-select',
                            options=[
                                {
                                    'label': 'Spring (Force-directed)',
                                    'value': 'spring',
                                },
                                {'label': 'Circular', 'value': 'circular'},
                                {'label': 'Radial', 'value': 'radial'},
                                {
                                    'label': 'Hierarchical',
                                    'value': 'hierarchical',
                                },
                                {
                                    'label': 'Kamada-Kawai',
                                    'value': 'kamada_kawai',
                                },
                                {
                                    'label': 'Geographic (requires coordinates)',
                                    'value': 'geographic',
                                },
                            ],
                            value='spring',
                            style={'whiteSpace': 'nowrap'},
                        ),
                        html.Br(),
                        html.Div(
                            id='geographic-options',
                            children=[
                                dbc.Label('Map Projection'),
                                dcc.Dropdown(
                                    id='map-projection',
                                    options=[
                                        {'label': 'Mercator', 'value': 'mercator'},
                                        {
                                            'label': 'PlateCarree',
                                            'value': 'platecarree',
                                        },
                                        {
                                            'label': 'Orthographic',
                                            'value': 'orthographic',
                                        },
                                    ],
                                    value='mercator',
                                ),
                                html.Br(),
                                dbc.Label('Zoom Level'),
                                dcc.Slider(
                                    id='map-zoom',
                                    min=1,
                                    max=10,
                                    step=1,
                                    value=2,
                                    marks={i: str(i) for i in range(1, 11)},
                                ),
                            ],
                            style={'display': 'none'},
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
                            '🎨 Apply Layout',
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
                dbc.CardHeader(
                    html.H5('4. Export', className='mb-0'),
                ),
                dbc.CardBody(
                    [
                        dbc.Label('Export Format', className='fw-bold'),
                        html.Small(
                            'Save your network for further analysis or publication',
                            className='text-muted d-block mb-2',
                        ),
                        dcc.Dropdown(
                            id='export-format',
                            options=[
                                {'label': 'GraphML (Cytoscape/Gephi)', 'value': 'graphml'},
                                {'label': 'GML (Graph Format)', 'value': 'gml'},
                                {'label': 'JSON', 'value': 'json'},
                                {'label': 'PNG Image', 'value': 'png'},
                                {'label': 'SVG Image', 'value': 'svg'},
                            ],
                            value='graphml',
                            style={'whiteSpace': 'nowrap'},
                        ),
                        html.Br(),
                        dbc.Button(
                            '💾 Download',
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
                            config={
                                'displayModeBar': True,
                                'displaylogo': False,
                                'editable': True,
                                'edits': {'shapePosition': True},
                            },
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
                        'fontFamily': 'monospace',
                        'whiteSpace': 'pre',
                    },
                )
            ]
        )

    def _create_haplotype_summary_tab(self) -> html.Div:
        """Create haplotype summary tab showing H number to sequence name mapping."""
        return html.Div(
            [
                html.Div(
                    id='haplotype-summary-display',
                    style={'padding': '20px', 'height': '85vh', 'overflow-y': 'auto'},
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
                Output('download-template-button', 'disabled'),
            ],
            Input('upload-data', 'contents'),
            State('upload-data', 'filename'),
        )
        def handle_file_upload(
            contents: Optional[str], filename: Optional[str]
        ) -> Tuple[html.Div, Optional[Dict], bool, bool]:
            """Handle file upload and parse alignment."""
            if contents is None:
                return html.Div(), None, True, True

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
                            [
                                html.Strong('❌ Unsupported file format'),
                                html.Br(),
                                f'File: {filename}',
                                html.Br(),
                                'Please use FASTA (.fasta, .fa), NEXUS (.nex, .nexus), or PHYLIP (.phy, .phylip) format.',
                            ],
                            color='danger',
                        ),
                        None,
                        True,
                        True,
                    )

                # Validate alignment
                if len(alignment) == 0:
                    return (
                        dbc.Alert(
                            [
                                html.Strong('⚠️ Empty alignment'),
                                html.Br(),
                                'The file contains no sequences. Please check your input file.',
                            ],
                            color='warning',
                        ),
                        None,
                        True,
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
                        html.Strong('✅ Success! '),
                        f'Loaded {len(alignment)} sequences '
                        f'of length {alignment.length} bp',
                    ],
                    color='success',
                )

                # Enable both compute button and template download button
                return status, alignment_data, False, False

            except Exception as e:
                self.logger.error(f'Error parsing file: {e}')
                return (
                    dbc.Alert(
                        [
                            html.Strong('❌ Error parsing file'),
                            html.Br(),
                            f'Error: {str(e)}',
                            html.Br(),
                            'Please check that your file is properly formatted.',
                        ],
                        color='danger',
                    ),
                    None,
                    True,
                    True,
                )

        @self.app.callback(
            [Output('metadata-status', 'children'), Output('metadata-store', 'data')],
            Input('upload-metadata', 'contents'),
            State('upload-metadata', 'filename'),
        )
        def handle_metadata_upload(
            contents: Optional[str], filename: Optional[str]
        ) -> Tuple[html.Div, Optional[Dict]]:
            """Handle metadata file upload and parse coordinates."""
            if contents is None:
                return html.Div(), None

            try:
                content_type, content_string = contents.split(',')
                decoded = base64.b64decode(content_string)

                # Parse CSV metadata
                if not (filename.endswith('.csv') or filename.endswith('.txt')):
                    return (
                        dbc.Alert(
                            [
                                html.Strong('❌ Invalid file type'),
                                html.Br(),
                                'Metadata must be a CSV (.csv) or text (.txt) file.',
                            ],
                            color='danger',
                        ),
                        None,
                    )

                with tempfile.NamedTemporaryFile(
                    mode='wb', suffix='.csv', delete=False
                ) as tmp:
                    tmp.write(decoded)
                    tmp.flush()
                    reader = MetadataReader(tmp.name, validate=False)
                    metadata_dict = reader.read_metadata()

                # Extract coordinates where available
                coordinates = {}
                for seq_id, meta in metadata_dict.items():
                    try:
                        coords = extract_coordinates(
                            meta, lat_column='latitude', lon_column='longitude'
                        )
                        if coords:
                            coordinates[seq_id] = coords
                    except (ValueError, KeyError):
                        pass

                metadata_data = {
                    'raw': metadata_dict,
                    'coordinates': coordinates,
                }

                # Build status message
                status_parts = [html.Strong('✅ Success! ')]
                status_parts.append(
                    f'Loaded metadata for {len(metadata_dict)} sequences.'
                )

                if coordinates:
                    status_parts.append(html.Br())
                    status_parts.append(
                        f'📍 Found geographic coordinates for {len(coordinates)} sequences.'
                    )
                else:
                    status_parts.append(html.Br())
                    status_parts.append(
                        '💡 Tip: Add latitude/longitude columns for geographic visualization.'
                    )

                status = dbc.Alert(status_parts, color='success')

                return status, metadata_data

            except Exception as e:
                self.logger.error(f'Error parsing metadata: {e}')
                return (
                    dbc.Alert(
                        [
                            html.Strong('❌ Error parsing metadata'),
                            html.Br(),
                            f'Error: {str(e)}',
                            html.Br(),
                            'Please check your CSV file format.',
                        ],
                        color='danger',
                    ),
                    None,
                )

        @self.app.callback(
            [
                Output('geographic-options', 'style'),
                Output('geographic-mode', 'data'),
            ],
            Input('layout-select', 'value'),
        )
        def toggle_geographic_options(layout: str) -> Tuple[Dict, bool]:
            """Show/hide geographic options based on layout selection."""
            if layout == 'geographic':
                return {'display': 'block'}, True
            else:
                return {'display': 'none'}, False

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
                                'median_vector', False
                            ),
                            'sample_ids': network.graph.nodes[node].get('sample_ids', []),
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

                # Count median/inferred nodes
                n_medians = sum(
                    1 for node in network.graph.nodes()
                    if network.graph.nodes[node].get('is_median', False)
                )

                feedback_parts = [html.Strong('✅ Network computed! ')]
                feedback_parts.append(
                    f'{len(network.graph.nodes)} haplotypes, '
                    f'{len(network.graph.edges)} connections'
                )

                if n_medians > 0:
                    feedback_parts.append(html.Br())
                    feedback_parts.append(f'🔵 {n_medians} inferred median nodes')

                feedback = dbc.Alert(feedback_parts, color='success')

                return network_data, feedback, False, False

            except Exception as e:
                self.logger.error(f'Error computing network: {e}')
                self.logger.error(traceback.format_exc())
                return (
                    None,
                    dbc.Alert(
                        [
                            html.Strong('❌ Error computing network'),
                            html.Br(),
                            f'Error: {str(e)}',
                            html.Br(),
                            'Please try a different algorithm or check your data.',
                        ],
                        color='danger',
                    ),
                    True,
                    True,
                )

        @self.app.callback(
            Output('layout-store', 'data'),
            [Input('apply-layout-button', 'n_clicks'), Input('network-store', 'data')],
            [
                State('layout-select', 'value'),
                State('snap-to-grid', 'value'),
                State('metadata-store', 'data'),
                State('map-projection', 'value'),
            ],
            prevent_initial_call=False,
        )
        def apply_layout(
            n_clicks: Optional[int],
            network_data: Optional[Dict],
            layout_method: str,
            snap_to_grid: List[str],
            metadata_data: Optional[Dict],
            projection: str,
        ) -> Optional[Dict]:
            """Apply layout algorithm to network."""
            if not network_data:
                return None

            try:
                # Reconstruct network
                network = HaplotypeNetwork.from_serialized(network_data)

                # Apply layout
                if layout_method == 'geographic':
                    # Geographic layout requires metadata with coordinates
                    if not metadata_data or not metadata_data.get('coordinates'):
                        # Fall back to spring layout if no coordinates
                        layout_manager = LayoutManager(network)
                        positions = layout_manager.compute_layout('spring')
                    else:
                        geo_layout = GeographicLayout(network)
                        positions = geo_layout.compute(
                            coordinates=metadata_data['coordinates'],
                            projection=projection or 'mercator',
                        )
                else:
                    layout_manager = LayoutManager(network)
                    positions = layout_manager.compute_layout(layout_method)

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
            [
                Input('layout-store', 'data'),
                Input('network-store', 'data'),
                Input('geographic-mode', 'data'),
            ],
            State('metadata-store', 'data'),
        )
        def update_network_graph(
            layout_data: Optional[Dict],
            network_data: Optional[Dict],
            geographic_mode: bool,
            metadata_data: Optional[Dict],
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
                network = HaplotypeNetwork.from_serialized(network_data)

                # Convert layout data
                positions = {node: tuple(pos) for node, pos in layout_data.items()}

                # Create interactive plot
                plotter = InteractiveNetworkPlotter(network)
                fig = plotter.plot(layout=positions)

                # Add geographic context if in geographic mode
                if (
                    geographic_mode
                    and metadata_data
                    and metadata_data.get('coordinates')
                ):
                    # Add annotation to indicate geographic mode
                    fig.add_annotation(
                        text='Geographic Mode',
                        xref='paper',
                        yref='paper',
                        x=0.02,
                        y=0.98,
                        showarrow=False,
                        font={'size': 12, 'color': 'blue'},
                        bgcolor='rgba(255,255,255,0.8)',
                        bordercolor='blue',
                        borderwidth=1,
                        borderpad=4,
                    )
                    # Update axis labels for geographic coordinates
                    fig.update_xaxes(title_text='Longitude')
                    fig.update_yaxes(title_text='Latitude')

                    # Add base map layer (simple grid with light background)
                    # This creates a simple map-like appearance
                    fig.update_layout(
                        plot_bgcolor='#e6f2ff',  # Light blue for ocean
                        xaxis={
                            'showgrid': True,
                            'gridcolor': 'lightgray',
                            'title': 'Longitude',
                            'zeroline': False,
                        },
                        yaxis={
                            'showgrid': True,
                            'gridcolor': 'lightgray',
                            'title': 'Latitude',
                            'zeroline': False,
                        },
                    )

                return fig

            except Exception as e:
                self.logger.error(f'Error creating visualization: {e}')
                self.logger.error(traceback.format_exc())
                fig = go.Figure()
                fig.add_annotation(
                    text=f'Error creating visualization: {str(e)}<br><br>See console for full traceback',
                    xref='paper',
                    yref='paper',
                    x=0.5,
                    y=0.5,
                    showarrow=False,
                    font={'size': 14, 'color': 'red'},
                )
                return fig

        @self.app.callback(
            Output('layout-store', 'data', allow_duplicate=True),
            Input('network-graph', 'relayoutData'),
            State('layout-store', 'data'),
            prevent_initial_call=True,
        )
        def update_node_positions(relayout_data: Optional[Dict], current_layout: Optional[Dict]) -> Optional[Dict]:
            """Update node positions when user drags nodes."""
            if not relayout_data or not current_layout:
                raise PreventUpdate

            # Check if this is a node drag event (contains x and y coordinates for specific traces)
            updated_layout = current_layout.copy()
            try:
                # Plotly relayout data contains keys like 'xaxis.range[0]', etc for zoom/pan
                # For node dragging, we would need custom implementation
                # For now, this maintains the current layout
                pass
            except Exception as e:
                self.logger.error(f'Error updating node positions: {e}')
                raise PreventUpdate from e

            return updated_layout

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
                self.logger.error(f'Error calculating statistics: {e}')
                self.logger.error(traceback.format_exc())
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
                self.logger.error(f'Error displaying alignment: {e}')
                self.logger.error(traceback.format_exc())
                return f'Error displaying alignment: {str(e)}'

        @self.app.callback(
            Output('haplotype-summary-display', 'children'),
            [Input('network-store', 'data'), Input('alignment-store', 'data')],
        )
        def update_haplotype_summary(
            network_data: Optional[Dict], alignment_data: Optional[Dict]
        ) -> html.Div:
            """Update haplotype summary showing H number to sequence name mapping."""
            if not network_data or not alignment_data:
                return html.Div('Compute network to view haplotype summary')

            try:
                # Reconstruct network
                network = HaplotypeNetwork.from_serialized(network_data)

                # Create mapping of H numbers to sequence names
                haplotype_mapping = []
                for i, node_id in enumerate(sorted(network.graph.nodes()), start=1):
                    node_data = network.graph.nodes[node_id]
                    sample_ids = node_data.get('sample_ids', [])
                    is_median = node_data.get('median_vector', False)
                    frequency = node_data.get('frequency', len(sample_ids))

                    h_label = f'H{i}'

                    # Determine if this is an inferred haplotype
                    if is_median or len(sample_ids) == 0:
                        haplotype_type = '🔵 Inferred'
                        sample_display = 'None (inferred ancestral haplotype)'
                    else:
                        haplotype_type = '🟢 Observed'
                        sample_display = ', '.join(sample_ids) if sample_ids else 'Unknown'

                    haplotype_mapping.append({
                        'h_label': h_label,
                        'node_id': node_id,
                        'type': haplotype_type,
                        'frequency': frequency,
                        'samples': sample_display,
                    })

                # Create table
                table_header = [
                    html.Thead(
                        html.Tr([
                            html.Th('H Number'),
                            html.Th('Type'),
                            html.Th('Frequency'),
                            html.Th('Sample IDs'),
                        ])
                    )
                ]

                table_rows = [
                    html.Tr([
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
                    ])
                    for hap in haplotype_mapping
                ]

                table_body = [html.Tbody(table_rows)]

                # Count statistics
                n_observed = sum(1 for h in haplotype_mapping if '🟢' in h['type'])
                n_inferred = sum(1 for h in haplotype_mapping if '🔵' in h['type'])

                return html.Div([
                    html.H4('Haplotype Summary'),
                    html.P([
                        f'Total haplotypes: {len(haplotype_mapping)} ',
                        f'(🟢 {n_observed} observed, 🔵 {n_inferred} inferred)',
                    ]),
                    html.Hr(),
                    dbc.Table(
                        table_header + table_body,
                        bordered=True,
                        hover=True,
                        responsive=True,
                        striped=True,
                        style={'fontSize': '14px'},
                    ),
                ])

            except Exception as e:
                self.logger.error(f'Error creating haplotype summary: {e}')
                self.logger.error(traceback.format_exc())
                return html.Div(
                    [
                        html.H5('Error creating haplotype summary', style={'color': 'red'}),
                        html.P(str(e)),
                        html.Details([
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
                        ]),
                    ],
                    style={'color': 'red', 'padding': '20px'},
                )

        @self.app.callback(
            Output('download-template', 'data'),
            Input('download-template-button', 'n_clicks'),
            State('alignment-store', 'data'),
            prevent_initial_call=True,
        )
        def download_metadata_template(
            n_clicks: int, alignment_data: Optional[Dict]
        ) -> Dict:
            """Generate and download metadata template CSV."""
            if not alignment_data:
                raise PreventUpdate

            try:
                # Generate CSV template with sequence IDs
                sequence_ids = [seq['id'] for seq in alignment_data['sequences']]

                # Create CSV content with headers
                csv_lines = ['sequence_id,population,latitude,longitude,color,notes']

                # Add a row for each sequence with empty fields
                for seq_id in sequence_ids:
                    csv_lines.append(f'{seq_id},,,,,')

                csv_content = '\n'.join(csv_lines)

                return {
                    'content': csv_content,
                    'filename': 'metadata_template.csv',
                    'type': 'text/csv',
                }

            except Exception as e:
                logging.error(f'Error generating template: {e}')
                raise PreventUpdate from None

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
                network = HaplotypeNetwork.from_serialized(network_data)
                G = network._graph

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
    debug :
        bool, default=False.
        Enable debug mode.
    port :
        int, default=8050.
        Port number for web server.
    """
    import argparse

    parser = argparse.ArgumentParser(
        description='PyPopART - Haplotype Network Analysis GUI',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  pypopart-gui                    # Start GUI on default port 8050
  pypopart-gui --port 8080        # Start GUI on custom port
  pypopart-gui --debug            # Start GUI in debug mode

Once started, open your browser to http://localhost:8050
Press Ctrl+C to stop the server.
        """,
    )
    parser.add_argument(
        '--debug',
        action='store_true',
        help='Enable debug mode for development',
    )
    parser.add_argument(
        '--port',
        type=int,
        default=8050,
        help='Port number for web server (default: 8050)',
    )

    # Parse args only if running from command line
    import sys

    if len(sys.argv) > 1:
        args = parser.parse_args()
        debug = args.debug
        port = args.port

    print('=' * 60)
    print('PyPopART GUI - Haplotype Network Analysis')
    print('=' * 60)
    print(f'\n🚀 Starting web server on http://localhost:{port}')
    if debug:
        print('⚠️  Debug mode enabled')
    print('\n📖 Quick Start:')
    print('   1. Upload your sequence alignment (FASTA, NEXUS, or PHYLIP)')
    print('   2. Optionally upload metadata (CSV with population/location data)')
    print('   3. Choose a network algorithm and click "Compute Network"')
    print('   4. Customize the layout and export your results')
    print('\n⚠️  To stop the server, press Ctrl+C')
    print('=' * 60)
    print()

    app = PyPopARTApp(debug=debug, port=port)
    app.run()


if __name__ == '__main__':
    main(debug=True)

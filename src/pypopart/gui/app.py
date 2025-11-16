"""
Dash-based GUI for PyPopART haplotype network analysis.

This module provides a web-based graphical user interface for PyPopART,
allowing users to upload sequence data, configure network algorithms,
visualize results, and export outputs.

Features:
- Interactive network visualization with manual node adjustment
- Haplotype summary tab showing H number to sequence mapping
- Geographic layout mode with base map display
- Support for multiple network algorithms (MST, MSN, TCS, MJN)
"""

import base64
import logging
import tempfile
import traceback
from typing import Dict, List, Optional, Tuple

import dash
import dash_cytoscape as cyto
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
from pypopart.visualization.cytoscape_plot import (
    InteractiveCytoscapePlotter,
    create_cytoscape_network,
)


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
                                            self._create_metadata_tab(),
                                            label='Metadata',
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
                dcc.Store(id='manual-edit-flag', data=False),
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
                            value='hierarchical',
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
                        dbc.Label('Node Spacing', className='fw-bold'),
                        html.Small(
                            'Adjust the spacing between nodes',
                            className='text-muted d-block mb-2',
                        ),
                        dcc.Slider(
                            id='spacing-slider',
                            min=0.5,
                            max=3.0,
                            step=0.1,
                            value=1.0,
                            marks={0.5: '0.5x', 1.0: '1.0x', 2.0: '2.0x', 3.0: '3.0x'},
                            tooltip={"placement": "bottom", "always_visible": False},
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
                # Search bar
                html.Div(
                    [
                        html.Label('Search Haplotype:', style={'marginRight': '10px', 'fontWeight': 'bold'}),
                        dcc.Dropdown(
                            id='haplotype-search',
                            placeholder='Select H number(s) to highlight...',
                            style={'width': '300px', 'display': 'inline-block'},
                            clearable=True,
                            multi=True,
                        ),
                        html.Div(id='search-feedback', style={'display': 'inline-block', 'marginLeft': '10px', 'color': 'red'}),
                    ],
                    style={
                        'position': 'absolute',
                        'bottom': '10px',
                        'right': '10px',
                        'background': 'white',
                        'padding': '10px',
                        'border': '1px solid #ccc',
                        'borderRadius': '5px',
                        'zIndex': 1000,
                        'display': 'flex',
                        'alignItems': 'center',
                    },
                ),
                # Legend display
                html.Div(
                    id='network-legend',
                    style={
                        'position': 'absolute',
                        'top': '10px',
                        'right': '10px',
                        'background': 'white',
                        'padding': '10px',
                        'border': '1px solid #ccc',
                        'borderRadius': '5px',
                        'zIndex': 1000,
                        'maxWidth': '200px',
                    },
                ),
                # Tooltip display on hover
                html.Div(
                    id='node-tooltip',
                    style={
                        'position': 'absolute',
                        'display': 'none',
                        'background': 'rgba(0, 0, 0, 0.8)',
                        'color': 'white',
                        'padding': '10px',
                        'borderRadius': '5px',
                        'zIndex': 2000,
                        'pointerEvents': 'none',
                        'maxWidth': '300px',
                        'fontSize': '12px',
                    },
                ),
                dcc.Loading(
                    id='loading-network',
                    type='default',
                    children=[
                        cyto.Cytoscape(
                            id='network-graph',
                            layout={'name': 'preset'},
                            style={'width': '100%', 'height': '85vh'},
                            elements=[],
                            stylesheet=[],
                            minZoom=0.1,
                            maxZoom=5,
                            wheelSensitivity=0.2,
                            zoom=1,
                            autoungrabify=False,
                        )
                    ],
                )
            ],
            style={'position': 'relative'}
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
                    [
                        dbc.Button(
                            '⬇️ Download Summary CSV',
                            id='download-haplotype-csv-button',
                            color='primary',
                            size='sm',
                            style={'marginRight': '10px'},
                        ),
                        dcc.Download(id='download-haplotype-csv'),
                        dbc.Button(
                            '⬇️ Download Label Template',
                            id='download-h-number-template-button',
                            color='info',
                            outline=True,
                            size='sm',
                            style={'marginRight': '10px'},
                        ),
                        dcc.Download(id='download-h-number-template'),
                        dcc.Upload(
                            id='upload-h-number-mapping',
                            children=dbc.Button(
                                '⬆️ Upload Label Mapping',
                                color='warning',
                                outline=True,
                                size='sm',
                            ),
                            style={'display': 'inline-block', 'marginRight': '10px'},
                        ),
                        html.Div(id='h-number-feedback', style={'display': 'inline-block', 'marginLeft': '10px'}),
                    ],
                    style={'padding': '20px 20px 10px 20px'},
                ),
                html.Div(
                    id='haplotype-summary-display',
                    style={'padding': '0 20px 20px 20px', 'height': '75vh', 'overflow-y': 'auto'},
                ),
                dcc.Store(id='h-number-mapping-store'),
            ]
        )

    def _create_metadata_tab(self) -> html.Div:
        """Create metadata tab showing imported metadata and alignment IDs."""
        return html.Div(
            [
                html.Div(
                    id='metadata-warnings',
                    style={'padding': '20px 20px 10px 20px'},
                ),
                html.Div(
                    id='metadata-display',
                    style={'padding': '0 20px 20px 20px', 'height': '75vh', 'overflow-y': 'auto'},
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
                            coordinates[seq_id] = {'lat': coords[0], 'lon': coords[1]}
                    except (ValueError, KeyError):
                        pass

                # Extract population labels from metadata
                populations = {}
                population_colors = {}
                for seq_id, meta in metadata_dict.items():
                    if 'population' in meta and meta['population']:
                        populations[seq_id] = meta['population']

                # Extract color mappings if provided in metadata
                colors = {}
                for seq_id, meta in metadata_dict.items():
                    if 'color' in meta and meta['color']:
                        colors[seq_id] = meta['color']
                        # Also track population colors
                        if 'population' in meta and meta['population']:
                            pop = meta['population']
                            if pop not in population_colors:
                                population_colors[pop] = meta['color']

                # If population labels provided but no colors, generate colors automatically
                if populations and not population_colors:
                    import colorsys
                    unique_pops = sorted(set(populations.values()))
                    n = len(unique_pops)
                    for i, pop in enumerate(unique_pops):
                        # Generate evenly spaced hues for distinct colors
                        hue = i / n
                        saturation = 0.7
                        value = 0.9
                        r, g, b = colorsys.hsv_to_rgb(hue, saturation, value)
                        hex_color = '#{:02x}{:02x}{:02x}'.format(
                            int(r * 255), int(g * 255), int(b * 255)
                        )
                        population_colors[pop] = hex_color
                    self.logger.info(f'Auto-generated colors for {len(unique_pops)} populations')

                metadata_data = {
                    'raw': metadata_dict,
                    'coordinates': coordinates,
                    'populations': populations,
                    'population_colors': population_colors if population_colors else None,
                    'sequence_ids': list(metadata_dict.keys()),
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
            [Input('apply-layout-button', 'n_clicks'), Input('network-store', 'data'), Input('spacing-slider', 'value')],
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

                # Apply spacing factor to expand/contract the layout
                if spacing_factor and spacing_factor != 1.0:
                    positions = {node: (pos[0] * spacing_factor, pos[1] * spacing_factor)
                                for node, pos in positions.items()}

                # Convert to serializable format
                layout_data = {node: list(pos) for node, pos in positions.items()}

                return layout_data

            except Exception as e:
                logging.error(f'Error applying layout: {e}')
                logging.error(traceback.format_exc())
                return None

        @self.app.callback(
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
            """Update network visualization with Cytoscape."""
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

                # Create Cytoscape elements and stylesheet
                elements, stylesheet = create_cytoscape_network(
                    network,
                    layout=positions,
                    population_colors=population_colors,
                    show_labels=True,
                    show_edge_labels=True,
                    node_labels=node_labels,
                )

                # Add geographic styling if in geographic mode
                if (
                    geographic_mode
                    and metadata_data
                    and metadata_data.get('coordinates')
                ):
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
                self.logger.error(f'Error creating visualization: {e}')
                self.logger.error(traceback.format_exc())
                error_msg = html.Div(
                    [
                        html.Strong('Error creating visualization'),
                        html.Br(),
                        str(e),
                    ],
                    style={'color': 'red'},
                )
                return [], [], error_msg

        @self.app.callback(
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
            """Update node positions when user drags nodes in Cytoscape."""
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
                self.logger.error(f'Error updating node positions: {e}')
                raise PreventUpdate from e

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
        def update_alignment_display(alignment_data: Optional[Dict]):
            """Update alignment viewer with colored nucleotides for polymorphic sites."""
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
                                        'color': 'white' if base.upper() not in ['-', 'N'] else 'black',
                                        'padding': '0 1px',
                                    },
                                )
                            )
                        else:
                            # Keep invariant positions as plain text
                            seq_spans.append(html.Span(base, style={'color': 'black'}))

                    # Combine ID and sequence
                    rows.append(html.Div([id_span] + seq_spans, style={'whiteSpace': 'pre'}))

                # Add note if sequences were truncated
                if len(alignment_data['sequences']) > 50:
                    rows.append(
                        html.Div(
                            f'\n... ({len(alignment_data["sequences"]) - 50} more sequences)',
                            style={'color': 'gray', 'fontStyle': 'italic', 'marginTop': '10px'},
                        )
                    )

                return html.Div(rows)

            except Exception as e:
                self.logger.error(f'Error displaying alignment: {e}')
                self.logger.error(traceback.format_exc())
                return f'Error displaying alignment: {str(e)}'

        @self.app.callback(
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
            """Update haplotype summary showing H number to sequence name mapping."""
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
                            populations_display = ', '.join(sorted(populations)) if populations else ''
                        else:
                            populations_display = ''

                    haplotype_mapping.append({
                        'h_label': h_label,
                        'node_id': node_id,
                        'type': haplotype_type,
                        'frequency': frequency,
                        'samples': sample_display,
                        'populations': populations_display,
                    })

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
                csv_lines = ['id,population,latitude,longitude,color,notes']

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

        # New callbacks for enhanced features

        @self.app.callback(
            [
                Output('haplotype-search', 'options'),
                Output('network-graph', 'stylesheet', allow_duplicate=True),
            ],
            [Input('network-store', 'data'), Input('haplotype-search', 'value')],
            [State('network-graph', 'stylesheet'), State('h-number-mapping-store', 'data')],
            prevent_initial_call=True,
        )
        def update_search_and_highlight(
            network_data: Optional[Dict],
            selected_h_list: Optional[List[str]],
            current_stylesheet: List[Dict],
            h_number_mapping: Optional[Dict],
        ) -> Tuple[List[Dict], List[Dict]]:
            """Update search dropdown options and highlight selected nodes."""
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

                # Remove ALL existing highlight styles (for any node)
                base_stylesheet = [s for s in current_stylesheet
                                  if not (s.get('selector', '').startswith('node[id = "')
                                         and 'border-color' in s.get('style', {})
                                         and s.get('style', {}).get('border-color') == '#FF0000')]

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
                            base_stylesheet.append({
                                'selector': f'node[id = "{node_id}"]',
                                'style': {
                                    'border-width': '5px',
                                    'border-color': '#FF0000',
                                    'border-style': 'solid',
                                }
                            })

                    return h_numbers, base_stylesheet

                # No selection - return clean stylesheet
                return h_numbers, base_stylesheet

            except Exception as e:
                self.logger.error(f'Error updating search: {e}')
                return [], current_stylesheet or []

        @self.app.callback(
            Output('download-haplotype-csv', 'data'),
            Input('download-haplotype-csv-button', 'n_clicks'),
            [
                State('network-store', 'data'),
                State('alignment-store', 'data'),
                State('metadata-store', 'data'),
                State('h-number-mapping-store', 'data'),
            ],
            prevent_initial_call=True,
        )
        def download_haplotype_csv(
            n_clicks: Optional[int],
            network_data: Optional[Dict],
            alignment_data: Optional[Dict],
            metadata_data: Optional[Dict],
            h_number_mapping: Optional[Dict],
        ) -> Optional[Dict]:
            """Download haplotype summary as CSV."""
            if not network_data or not alignment_data:
                raise PreventUpdate

            try:
                import csv
                import io

                # Reconstruct network
                network = HaplotypeNetwork.from_serialized(network_data)

                # Build CSV content
                output = io.StringIO()
                writer = csv.writer(output)

                # Write header
                headers = ['H_Number', 'Type', 'Frequency', 'Sample_IDs']
                if metadata_data and metadata_data.get('populations'):
                    headers.append('Populations')
                writer.writerow(headers)

                # Write data rows
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

                    if is_median or len(sample_ids) == 0:
                        haplotype_type = 'Inferred'
                        sample_display = 'None'
                    else:
                        haplotype_type = 'Observed'
                        sample_display = '; '.join(sample_ids) if sample_ids else 'Unknown'

                    row = [h_label, haplotype_type, frequency, sample_display]

                    # Add populations if metadata available
                    if metadata_data and metadata_data.get('populations'):
                        populations = set()
                        for sid in sample_ids:
                            if sid in metadata_data['populations']:
                                populations.add(metadata_data['populations'][sid])
                        pop_display = '; '.join(sorted(populations)) if populations else ''
                        row.append(pop_display)

                    writer.writerow(row)

                return {
                    'content': output.getvalue(),
                    'filename': 'haplotype_summary.csv',
                }

            except Exception as e:
                self.logger.error(f'Error generating CSV: {e}')
                raise PreventUpdate from e

        @self.app.callback(
            Output('metadata-display', 'children'),
            Output('metadata-warnings', 'children'),
            [Input('alignment-store', 'data'), Input('metadata-store', 'data')],
        )
        def update_metadata_tab(
            alignment_data: Optional[Dict],
            metadata_data: Optional[Dict],
        ) -> Tuple[html.Div, html.Div]:
            """Display metadata with all sequence IDs."""
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
                            'latitude': metadata_data.get('coordinates', {}).get(sid, {}).get('lat', ''),
                            'longitude': metadata_data.get('coordinates', {}).get(sid, {}).get('lon', ''),
                        }

                # Union of all IDs
                all_ids = alignment_ids.union(metadata_ids)

                # Check for duplicates in alignment
                alignment_id_list = [seq['id'] for seq in alignment_data['sequences']]
                alignment_duplicates = [sid for sid in set(alignment_id_list) if alignment_id_list.count(sid) > 1]

                # Check for mismatches
                only_in_alignment = alignment_ids - metadata_ids
                only_in_metadata = metadata_ids - alignment_ids

                # Build warnings
                warnings = []
                if alignment_duplicates:
                    warnings.append(dbc.Alert(
                        f'⚠️ Duplicate IDs found in alignment: {", ".join(alignment_duplicates)}',
                        color='warning',
                    ))
                if only_in_alignment and metadata_data:
                    warnings.append(dbc.Alert(
                        f'⚠️ {len(only_in_alignment)} IDs only in alignment (not in metadata)',
                        color='info',
                    ))
                if only_in_metadata:
                    warnings.append(dbc.Alert(
                        f'⚠️ {len(only_in_metadata)} IDs only in metadata (not in alignment)',
                        color='info',
                    ))

                # Get population colors if available
                population_colors = metadata_data.get('population_colors', {}) if metadata_data else {}

                # Build table
                table_header = [
                    html.Thead(
                        html.Tr([
                            html.Th('Sequence ID'),
                            html.Th('In Alignment'),
                            html.Th('In Metadata'),
                            html.Th('Population'),
                            html.Th('Color'),
                            html.Th('Latitude'),
                            html.Th('Longitude'),
                        ])
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
                        color_display = html.Div([
                            html.Span('●', style={'color': color_hex, 'fontSize': '16px', 'marginRight': '5px'}),
                            html.Span(color_hex, style={'fontSize': '12px'})
                        ])

                    table_rows.append(
                        html.Tr([
                            html.Td(sid),
                            html.Td(in_alignment, style={'textAlign': 'center'}),
                            html.Td(in_metadata, style={'textAlign': 'center'}),
                            html.Td(pop),
                            html.Td(color_display),
                            html.Td(meta.get('latitude', '')),
                            html.Td(meta.get('longitude', '')),
                        ])
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
                self.logger.error(f'Error creating metadata display: {e}')
                self.logger.error(traceback.format_exc())
                return html.Div(f'Error: {str(e)}'), html.Div()

        @self.app.callback(
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
            """Update tooltip content based on hovered node."""
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
                    content = html.Div([
                        html.Strong(h_label or 'Unknown'),
                        html.Br(),
                        html.Span('Inferred median vector'),
                    ])
                else:
                    content = html.Div([
                        html.Strong(h_label or 'Unknown'),
                        html.Br(),
                        html.Span(f'Sequences ({len(sample_ids)}):'),
                        html.Br(),
                        html.Span(', '.join(sample_ids[:10]) + ('...' if len(sample_ids) > 10 else '')),
                    ])

                return content

            except Exception as e:
                self.logger.error(f'Error showing tooltip: {e}')
                return html.Div()

        # Use clientside callback for tooltip positioning
        # This gets the actual rendered position from Cytoscape
        self.app.clientside_callback(
            """
            function(hoverData, edgeHoverData) {
                // Hide tooltip if hovering over edge or no node data
                if ((edgeHoverData && !hoverData) || !hoverData) {
                    return {display: 'none'};
                }

                try {
                    // Get Cytoscape instance
                    const cy = document.getElementById('network-graph')._cyreg.cy;
                    if (!cy) {
                        return {display: 'none'};
                    }

                    // Get the node
                    const nodeId = hoverData.id;
                    const node = cy.getElementById(nodeId);

                    if (!node || node.length === 0) {
                        return {display: 'none'};
                    }

                    // Get rendered position (screen coordinates)
                    const renderedPos = node.renderedPosition();

                    // Position tooltip with offset from node
                    return {
                        display: 'block',
                        position: 'absolute',
                        left: (renderedPos.x + 15) + 'px',
                        top: (renderedPos.y - 40) + 'px',
                        background: 'rgba(0, 0, 0, 0.8)',
                        color: 'white',
                        padding: '10px',
                        borderRadius: '5px',
                        zIndex: 2000,
                        pointerEvents: 'none',
                        maxWidth: '300px',
                        fontSize: '12px',
                        whiteSpace: 'normal',
                    };
                } catch (e) {
                    console.log('Error positioning tooltip:', e);
                    return {display: 'none'};
                }
            }
            """,
            Output('node-tooltip', 'style'),
            [Input('network-graph', 'mouseoverNodeData'), Input('network-graph', 'mouseoverEdgeData')],
        )

        @self.app.callback(
            Output('download-h-number-template', 'data'),
            Input('download-h-number-template-button', 'n_clicks'),
            [State('network-store', 'data'), State('h-number-mapping-store', 'data')],
            prevent_initial_call=True,
        )
        def download_h_number_template(
            n_clicks: Optional[int],
            network_data: Optional[Dict],
            mapping_data: Optional[Dict],
        ) -> Optional[Dict]:
            """Download H number label mapping template as CSV."""
            if not network_data:
                raise PreventUpdate

            try:
                import csv
                import io

                # Reconstruct network
                network = HaplotypeNetwork.from_serialized(network_data)

                # Build CSV content
                output = io.StringIO()
                writer = csv.writer(output)

                # Write header
                writer.writerow(['current_h_number', 'new_label'])

                # Write data rows with current H numbers
                for i, node_id in enumerate(sorted(network.graph.nodes()), start=1):
                    current_label = f'H{i}'
                    # If there's already a custom mapping, use it
                    if mapping_data and node_id in mapping_data:
                        new_label = mapping_data[node_id]
                    else:
                        new_label = current_label
                    writer.writerow([current_label, new_label])

                return {
                    'content': output.getvalue(),
                    'filename': 'h_number_mapping_template.csv',
                }

            except Exception as e:
                self.logger.error(f'Error generating H number template: {e}')
                raise PreventUpdate from e

        @self.app.callback(
            [
                Output('h-number-mapping-store', 'data'),
                Output('h-number-feedback', 'children'),
                Output('network-graph', 'elements', allow_duplicate=True),
            ],
            Input('upload-h-number-mapping', 'contents'),
            [
                State('upload-h-number-mapping', 'filename'),
                State('network-store', 'data'),
                State('layout-store', 'data'),
                State('metadata-store', 'data'),
            ],
            prevent_initial_call=True,
        )
        def upload_h_number_mapping(
            contents: Optional[str],
            filename: Optional[str],
            network_data: Optional[Dict],
            layout_data: Optional[Dict],
            metadata_data: Optional[Dict],
        ) -> Tuple[Optional[Dict], html.Div, List[Dict]]:
            """Process uploaded H number mapping CSV and update graph."""
            if not contents or not network_data:
                raise PreventUpdate

            try:
                import csv
                import io

                # Decode uploaded file
                content_type, content_string = contents.split(',')
                decoded = base64.b64decode(content_string).decode('utf-8')

                # Parse CSV
                csv_reader = csv.DictReader(io.StringIO(decoded))

                # Validate required columns
                if csv_reader.fieldnames is None or set(csv_reader.fieldnames) != {'current_h_number', 'new_label'}:
                    return None, dbc.Alert(
                        [
                            html.Strong('❌ Invalid CSV Format'),
                            html.Br(),
                            'CSV must have exactly two columns: "current_h_number" and "new_label"',
                        ],
                        color='danger',
                        dismissable=True,
                    ), dash.no_update

                # Reconstruct network to get node IDs
                network = HaplotypeNetwork.from_serialized(network_data)
                node_ids = sorted(network.graph.nodes())

                # Build mapping from current H numbers to node IDs
                h_to_node = {}
                for i, node_id in enumerate(node_ids, start=1):
                    h_to_node[f'H{i}'] = node_id

                # Parse the uploaded mapping
                new_mapping = {}
                new_labels_list = []
                errors = []

                for row_num, row in enumerate(csv_reader, start=2):
                    current_h = row.get('current_h_number', '').strip()
                    new_label = row.get('new_label', '').strip()

                    if not current_h:
                        errors.append(f'Row {row_num}: Missing current_h_number')
                        continue

                    if not new_label:
                        errors.append(f'Row {row_num}: Missing new_label for {current_h}')
                        continue

                    if current_h not in h_to_node:
                        errors.append(f'Row {row_num}: Unknown H number "{current_h}"')
                        continue

                    node_id = h_to_node[current_h]
                    new_mapping[node_id] = new_label
                    new_labels_list.append(new_label)

                # Check for duplicate new labels
                seen_labels = {}
                for node_id, label in new_mapping.items():
                    if label in seen_labels:
                        errors.append(
                            f'Duplicate label "{label}" for {seen_labels[label]} and node {node_id}'
                        )
                    else:
                        seen_labels[label] = node_id

                if errors:
                    error_msg = html.Div([
                        html.Strong('❌ Validation Errors:'),
                        html.Ul([html.Li(err) for err in errors[:10]]),
                        html.P(f'({len(errors)} total errors)') if len(errors) > 10 else None,
                    ])
                    return None, dbc.Alert(error_msg, color='danger', dismissable=True), dash.no_update

                # If validation passed, update the graph with new labels
                positions = {node: tuple(pos) for node, pos in layout_data.items()}

                # Extract population colors from metadata if available
                population_colors = None
                if metadata_data and metadata_data.get('populations'):
                    population_colors = metadata_data.get('population_colors', {})

                # Create Cytoscape elements with custom labels
                elements, _ = create_cytoscape_network(
                    network,
                    layout=positions,
                    population_colors=population_colors,
                    show_labels=True,
                    show_edge_labels=True,
                    node_labels=new_mapping,
                )

                success_msg = dbc.Alert(
                    [
                        html.Strong('✅ Success!'),
                        html.Br(),
                        f'Updated {len(new_mapping)} H number labels from "{filename}"',
                    ],
                    color='success',
                    dismissable=True,
                    duration=4000,
                )

                return new_mapping, success_msg, elements

            except Exception as e:
                self.logger.error(f'Error processing H number mapping: {e}')
                self.logger.error(traceback.format_exc())
                return None, dbc.Alert(
                    [
                        html.Strong('❌ Error processing mapping'),
                        html.Br(),
                        str(e),
                    ],
                    color='danger',
                    dismissable=True,
                ), dash.no_update

        # Clientside callback to auto-fit network when elements are updated
        # Only fit when not in manual edit mode
        self.app.clientside_callback(
            """
            function(elements, manualEditFlag) {
                if (!elements || elements.length === 0) {
                    return window.dash_clientside.no_update;
                }
                // Only auto-fit if not manually editing
                if (!manualEditFlag) {
                    // Trigger a fit after elements are loaded
                    setTimeout(function() {
                        try {
                            const cy = document.getElementById('network-graph')._cyreg.cy;
                            if (cy) {
                                cy.fit(null, 50);  // Fit with 50px padding
                                cy.center();
                            }
                        } catch (e) {
                            console.log('Could not auto-fit network:', e);
                        }
                    }, 100);
                }
                return elements.length;
            }
            """,
            Output('network-graph', 'zoom'),
            [Input('network-graph', 'elements'), Input('manual-edit-flag', 'data')],
        )

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

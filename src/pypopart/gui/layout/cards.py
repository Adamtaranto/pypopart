"""Layout builders for the PyPopART GUI: page skeleton, cards and tabs."""

from dash import dcc, html
import dash_bootstrap_components as dbc
import dash_cytoscape as cyto


def build_layout(app) -> None:
    """Set up the application layout with all components."""
    app.layout = html.Div(
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
                            create_upload_card(),
                            html.Br(),
                            create_algorithm_card(),
                            html.Br(),
                            create_layout_card(),
                            html.Br(),
                            create_export_card(),
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
                                        create_network_tab(),
                                        label='Network',
                                    ),
                                    dbc.Tab(
                                        create_statistics_tab(),
                                        label='Statistics',
                                    ),
                                    dbc.Tab(
                                        create_haplotype_summary_tab(),
                                        label='Haplotype Summary',
                                    ),
                                    dbc.Tab(
                                        create_metadata_tab(),
                                        label='Metadata',
                                    ),
                                    dbc.Tab(
                                        create_alignment_tab(),
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
            # Store to trigger window resize handling
            dcc.Store(id='window-size-store'),
        ]
    )


def create_upload_card() -> dbc.Card:
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


def create_algorithm_card() -> dbc.Card:
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
                            {
                                'label': 'PN - Parsimony Network',
                                'value': 'pn',
                            },
                            {
                                'label': 'TSW - Tight Span Walker (Parsimony Network)',
                                'value': 'tsw',
                            },
                        ],
                        value='mst',
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


def create_layout_card() -> dbc.Card:
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
                                'label': 'Hierarchical (Fast)',
                                'value': 'hierarchical',
                            },
                            {
                                'label': 'Spring (Force-directed)',
                                'value': 'spring',
                            },
                            {
                                'label': 'Spring - Proportional Edge Length',
                                'value': 'spring_proportional',
                            },
                            {
                                'label': 'Spectral (Fast, Large networks)',
                                'value': 'spectral',
                            },
                            {'label': 'Circular', 'value': 'circular'},
                            {'label': 'Radial', 'value': 'radial'},
                            {
                                'label': 'Kamada-Kawai (High quality, slow)',
                                'value': 'kamada_kawai',
                            },
                            {
                                'label': 'Kamada-Kawai - Proportional Edge Length',
                                'value': 'kamada_kawai_proportional',
                            },
                            #                                {
                            #'label': 'Geographic (requires coordinates)',
                            #'value': 'geographic',
                            #                                },
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
                        value=2.0,
                        marks={0.5: '0.5x', 1.0: '1.0x', 2.0: '2.0x', 3.0: '3.0x'},
                        tooltip={'placement': 'bottom', 'always_visible': False},
                    ),
                    html.Br(),
                    dbc.Label('Node Size', className='fw-bold'),
                    html.Small(
                        'Adjust the size of nodes',
                        className='text-muted d-block mb-2',
                    ),
                    dcc.Slider(
                        id='node-size-slider',
                        min=10,
                        max=100,
                        step=5,
                        value=40,
                        marks={10: '10', 40: '40', 70: '70', 100: '100'},
                        tooltip={'placement': 'bottom', 'always_visible': False},
                    ),
                    html.Br(),
                    dbc.Label('Edge Width', className='fw-bold'),
                    html.Small(
                        'Adjust the width of edges',
                        className='text-muted d-block mb-2',
                    ),
                    dcc.Slider(
                        id='edge-width-slider',
                        min=1,
                        max=10,
                        step=0.5,
                        value=3,
                        marks={1: '1', 3: '3', 6: '6', 10: '10'},
                        tooltip={'placement': 'bottom', 'always_visible': False},
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


def create_export_card() -> dbc.Card:
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
                            {
                                'label': 'GraphML (Cytoscape/Gephi)',
                                'value': 'graphml',
                            },
                            {'label': 'GML (Graph Format)', 'value': 'gml'},
                            {'label': 'JSON', 'value': 'json'},
                            {'label': 'PNG Image', 'value': 'png'},
                            {'label': 'SVG Image', 'value': 'svg'},
                        ],
                        value='svg',
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


def create_network_tab() -> html.Div:
    """Create network visualization tab."""
    return html.Div(
        [
            # Search bar
            html.Div(
                [
                    html.Label(
                        'Search Haplotype:',
                        style={'marginRight': '10px', 'fontWeight': 'bold'},
                    ),
                    dcc.Dropdown(
                        id='haplotype-search',
                        placeholder='Select H number(s) to highlight...',
                        style={'width': '300px', 'display': 'inline-block'},
                        clearable=True,
                        multi=True,
                    ),
                    html.Div(
                        id='search-feedback',
                        style={
                            'display': 'inline-block',
                            'marginLeft': '10px',
                            'color': 'red',
                        },
                    ),
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
            ),
        ],
        style={'position': 'relative'},
    )


def create_statistics_tab() -> html.Div:
    """Create statistics display tab."""
    return html.Div(
        [
            dcc.Loading(
                html.Div(
                    id='statistics-display',
                    style={'padding': '20px', 'height': '85vh', 'overflow-y': 'auto'},
                )
            )
        ]
    )


def create_alignment_tab() -> html.Div:
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


def create_haplotype_summary_tab() -> html.Div:
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
                    html.Div(
                        id='h-number-feedback',
                        style={'display': 'inline-block', 'marginLeft': '10px'},
                    ),
                ],
                style={'padding': '20px 20px 10px 20px'},
            ),
            html.Div(
                id='haplotype-summary-display',
                style={
                    'padding': '0 20px 20px 20px',
                    'height': '75vh',
                    'overflow-y': 'auto',
                },
            ),
            dcc.Store(id='h-number-mapping-store'),
        ]
    )


def create_metadata_tab() -> html.Div:
    """Create metadata tab showing imported metadata and alignment IDs."""
    return html.Div(
        [
            html.Div(
                id='metadata-warnings',
                style={'padding': '20px 20px 10px 20px'},
            ),
            html.Div(
                id='metadata-display',
                style={
                    'padding': '0 20px 20px 20px',
                    'height': '75vh',
                    'overflow-y': 'auto',
                },
            ),
        ]
    )

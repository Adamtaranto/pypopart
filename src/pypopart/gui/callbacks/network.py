"""Algorithm parameter controls and network computation."""

import traceback
from typing import Dict, List, Optional, Tuple

import dash
from dash import Input, Output, State, dcc, html
from dash.exceptions import PreventUpdate
import dash_bootstrap_components as dbc

from pypopart.core.alignment import Alignment
from pypopart.gui.serialization import network_to_store


def register(app, logger) -> None:
    """
    Register network callbacks on the Dash app.

    Parameters
    ----------
    app : dash.Dash
        The Dash application.
    logger : logging.Logger
        Application logger.
    """

    @app.callback(
        Output('algorithm-parameters', 'children'),
        Input('algorithm-select', 'value'),
    )
    def update_algorithm_parameters(algorithm: str) -> html.Div:
        """
        Update parameter controls based on selected algorithm.

        Parameters
        ----------
        algorithm : str
            Selected algorithm name.

        Returns
        -------
        html.Div
            The parameter controls for the selected algorithm.
        """
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
        elif algorithm == 'pn':
            return html.Div(
                [
                    dbc.Label('Number of Trees'),
                    dcc.Slider(
                        id={'type': 'algorithm-param', 'name': 'n_trees'},
                        min=10,
                        max=500,
                        step=10,
                        value=100,
                        marks={i: str(i) for i in range(0, 501, 100)},
                    ),
                    html.Small(
                        'Number of random parsimony trees to sample',
                        className='text-muted',
                    ),
                ]
            )
        elif algorithm == 'tsw':
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
                    html.Br(),
                    html.Small(
                        'Builds parsimony network using tight span computation',
                        className='text-muted',
                    ),
                ]
            )
        return html.Div()

    @app.callback(
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
        """
        Compute haplotype network using selected algorithm.

        Parameters
        ----------
        n_clicks : int
            Button click count from Dash.
        alignment_data : Dict
            Serialized alignment from the alignment store.
        algorithm : str
            Selected algorithm name.
        param_values : List
            Algorithm parameter values from the controls.

        Returns
        -------
        Tuple[Optional[Dict], html.Div, bool, bool]
            The computed network for its store, status feedback, and the
            disabled state of the layout and export controls.
        """
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

            # Select and configure algorithm via the shared registry.
            # The single pattern-matched parameter means different
            # things per algorithm (distance for mst/msn/tsw, etc.).
            from pypopart.algorithms import build

            algo_kwargs = {}
            if algorithm in ('mst', 'msn', 'tsw'):
                algo_kwargs['distance_method'] = param_value or 'hamming'
            elif algorithm == 'tcs':
                # None = no limit (PopART parity); int caps connections
                algo_kwargs['connection_limit'] = param_value
            elif algorithm == 'mjn':
                algo_kwargs['epsilon'] = param_value or 0
            elif algorithm == 'pn':
                algo_kwargs['n_trees'] = param_value or 20
            algo = build(algorithm, **algo_kwargs)

            # Build network
            network = algo.build_network(alignment)

            # Convert to serializable format
            network_data = network_to_store(network)

            # Count median/inferred nodes (graph nodes use 'median_vector')
            n_medians = sum(
                1
                for node in network.graph.nodes()
                if network.graph.nodes[node].get('median_vector', False)
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
            logger.error(f'Error computing network: {e}')
            logger.error(traceback.format_exc())
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

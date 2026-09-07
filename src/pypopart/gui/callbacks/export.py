"""Network and CSV export downloads."""

import logging
from pathlib import Path
import tempfile
import traceback
from typing import Dict, Optional, Tuple

import dash
from dash import Input, Output, State
from dash.exceptions import PreventUpdate

from pypopart.core.graph import HaplotypeNetwork
from pypopart.io.network_export import GMLExporter, GraphMLExporter, JSONExporter


def register(app, logger) -> None:
    """
    Register export callbacks on the Dash app.

    Parameters
    ----------
    app : dash.Dash
        The Dash application.
    logger : logging.Logger
        Application logger.
    """

    @app.callback(
        [
            Output('download-data', 'data'),
            Output('network-graph', 'generateImage', allow_duplicate=True),
        ],
        Input('export-button', 'n_clicks'),
        [
            State('network-store', 'data'),
            State('export-format', 'value'),
        ],
        prevent_initial_call=True,
    )
    def export_network(
        n_clicks: int, network_data: Dict, export_format: str
    ) -> Tuple[Dict, Dict]:
        """Export network in selected format."""
        if not network_data:
            raise PreventUpdate

        try:
            # Reconstruct network
            network = HaplotypeNetwork.from_serialized(network_data)

            # Export based on format
            if export_format == 'graphml':
                with tempfile.TemporaryDirectory() as tmpdir:
                    path = Path(tmpdir) / 'network.graphml'
                    GraphMLExporter(path).export(network)
                    content = path.read_text()
                return (
                    {
                        'content': content,
                        'filename': 'network.graphml',
                        'type': 'text/xml',
                    },
                    dash.no_update,
                )

            elif export_format == 'gml':
                with tempfile.TemporaryDirectory() as tmpdir:
                    path = Path(tmpdir) / 'network.gml'
                    GMLExporter(path).export(network)
                    content = path.read_text()
                return (
                    {
                        'content': content,
                        'filename': 'network.gml',
                        'type': 'text/plain',
                    },
                    dash.no_update,
                )

            elif export_format == 'json':
                with tempfile.TemporaryDirectory() as tmpdir:
                    path = Path(tmpdir) / 'network.json'
                    JSONExporter(path).export(network)
                    content = path.read_text()
                return (
                    {
                        'content': content,
                        'filename': 'network.json',
                        'type': 'application/json',
                    },
                    dash.no_update,
                )

            elif export_format in ['png', 'svg']:
                # Use Cytoscape's built-in image generation feature
                # This triggers client-side export with automatic download
                image_config = {
                    'type': export_format,
                    'action': 'download',
                    'filename': f'network.{export_format}',
                    'options': {
                        'output': 'base64uri',
                        'bg': 'white',
                        'full': True,
                    },
                }
                return (dash.no_update, image_config)

        except Exception as e:
            logging.error(f'Error exporting: {e}')
            logging.error(traceback.format_exc())
            raise PreventUpdate from None

    # New callbacks for enhanced features

    @app.callback(
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
            logger.error(f'Error generating CSV: {e}')
            raise PreventUpdate from e

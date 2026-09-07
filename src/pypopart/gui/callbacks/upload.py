"""Upload handling: alignment/metadata/H-number files and templates."""

import base64
import logging
import traceback
from typing import Dict, List, Optional, Tuple

import dash
from dash import Input, Output, State, html
from dash.exceptions import PreventUpdate
import dash_bootstrap_components as dbc

from pypopart.core.graph import HaplotypeNetwork
from pypopart.io import FastaReader, NexusReader, PhylipReader
from pypopart.io.metadata import MetadataReader, extract_coordinates
from pypopart.visualization.cytoscape_plot import (
    create_cytoscape_network,
)


def _detect_reader(text: str, filename: str):
    """
    Pick a sequence reader from file content, falling back to extension.

    Parameters
    ----------
    text : str
        Decoded upload content.
    filename : str
        Original file name (extension fallback).

    Returns
    -------
    type or None
        Reader class with a from_string constructor, or None when the
        format cannot be determined.
    """
    stripped = text.lstrip()
    if stripped.startswith('#NEXUS'):
        return NexusReader
    if stripped.startswith('>'):
        return FastaReader
    first_line = stripped.splitlines()[0].split() if stripped else []
    if len(first_line) == 2 and all(tok.isdigit() for tok in first_line):
        return PhylipReader

    name = (filename or '').lower()
    if name.endswith(('.fasta', '.fa', '.fna')):
        return FastaReader
    if name.endswith(('.nex', '.nexus')):
        return NexusReader
    if name.endswith(('.phy', '.phylip')):
        return PhylipReader
    return None


def register(app, logger) -> None:
    """
    Register upload callbacks on the Dash app.

    Parameters
    ----------
    app : dash.Dash
        The Dash application.
    logger : logging.Logger
        Application logger.
    """

    @app.callback(
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
        """
        Handle file upload and parse alignment.

        Parameters
        ----------
        contents : str, optional
            Base64-encoded upload payload from Dash.
        filename : str, optional
            Name of the uploaded file.

        Returns
        -------
        Tuple[html.Div, Optional[Dict], bool, bool]
            Upload feedback, the parsed alignment for its store, and the
            disabled state of the compute and template buttons.
        """
        if contents is None:
            return html.Div(), None, True, True

        try:
            content_type, content_string = contents.split(',')
            decoded = base64.b64decode(content_string)

            # Parse in memory: detect format by content, fall back to
            # the file extension (no temporary files)
            text = decoded.decode('utf-8', errors='replace')
            reader_cls = _detect_reader(text, filename)
            if reader_cls is not None:
                alignment = reader_cls.from_string(text).read_alignment()
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
            logger.error(f'Error parsing file: {e}')
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

    @app.callback(
        [Output('metadata-status', 'children'), Output('metadata-store', 'data')],
        Input('upload-metadata', 'contents'),
        State('upload-metadata', 'filename'),
    )
    def handle_metadata_upload(
        contents: Optional[str], filename: Optional[str]
    ) -> Tuple[html.Div, Optional[Dict]]:
        """
        Handle metadata file upload and parse coordinates.

        Parameters
        ----------
        contents : str, optional
            Base64-encoded upload payload from Dash.
        filename : str, optional
            Name of the uploaded file.

        Returns
        -------
        Tuple[html.Div, Optional[Dict]]
            Upload feedback and the parsed metadata for its store.
        """
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

            reader = MetadataReader.from_string(
                decoded.decode('utf-8', errors='replace'), validate=False
            )
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
                logger.info(f'Auto-generated colors for {len(unique_pops)} populations')

            metadata_data = {
                'raw': metadata_dict,
                'coordinates': coordinates,
                'populations': populations,
                'population_colors': population_colors if population_colors else None,
                'sequence_ids': list(metadata_dict.keys()),
            }

            # Build status message
            status_parts = [html.Strong('✅ Success! ')]
            status_parts.append(f'Loaded metadata for {len(metadata_dict)} sequences.')

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
            logger.error(f'Error parsing metadata: {e}')
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

    @app.callback(
        Output('download-template', 'data'),
        Input('download-template-button', 'n_clicks'),
        State('alignment-store', 'data'),
        prevent_initial_call=True,
    )
    def download_metadata_template(
        n_clicks: int, alignment_data: Optional[Dict]
    ) -> Dict:
        """
        Generate and download metadata template CSV.

        Parameters
        ----------
        n_clicks : int
            Button click count from Dash.
        alignment_data : Dict, optional
            Serialized alignment from the alignment store.

        Returns
        -------
        Dict
            Download payload for the metadata template.
        """
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

    @app.callback(
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
        """
        Download H number label mapping template as CSV.

        Parameters
        ----------
        n_clicks : int, optional
            Button click count from Dash.
        network_data : Dict, optional
            Serialized network from the network store.
        mapping_data : Dict, optional
            Uploaded haplotype label mapping.

        Returns
        -------
        Optional[Dict]
            Download payload for the haplotype label template.
        """
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
            logger.error(f'Error generating H number template: {e}')
            raise PreventUpdate from e

    @app.callback(
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
        """
        Process uploaded H number mapping CSV and update graph.

        Parameters
        ----------
        contents : str, optional
            Base64-encoded upload payload from Dash.
        filename : str, optional
            Name of the uploaded file.
        network_data : Dict, optional
            Serialized network from the network store.
        layout_data : Dict, optional
            Node positions from the layout store.
        metadata_data : Dict, optional
            Serialized metadata from the metadata store.

        Returns
        -------
        Tuple[Optional[Dict], html.Div, List[Dict]]
            The parsed label mapping, upload feedback, and the refreshed
            Cytoscape elements.
        """
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
            if csv_reader.fieldnames is None or set(csv_reader.fieldnames) != {
                'current_h_number',
                'new_label',
            }:
                return (
                    None,
                    dbc.Alert(
                        [
                            html.Strong('❌ Invalid CSV Format'),
                            html.Br(),
                            'CSV must have exactly two columns: "current_h_number" and "new_label"',
                        ],
                        color='danger',
                        dismissable=True,
                    ),
                    dash.no_update,
                )

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
                error_msg = html.Div(
                    [
                        html.Strong('❌ Validation Errors:'),
                        html.Ul([html.Li(err) for err in errors[:10]]),
                        html.P(f'({len(errors)} total errors)')
                        if len(errors) > 10
                        else None,
                    ]
                )
                return (
                    None,
                    dbc.Alert(error_msg, color='danger', dismissable=True),
                    dash.no_update,
                )

            # If validation passed, update the graph with new labels
            positions = {node: tuple(pos) for node, pos in layout_data.items()}

            # Extract population colors and mapping from metadata if available
            population_colors = None
            population_mapping = None
            if metadata_data and metadata_data.get('populations'):
                population_colors = metadata_data.get('population_colors', {})
                population_mapping = metadata_data['populations']

            # Create Cytoscape elements with custom labels
            elements, _ = create_cytoscape_network(
                network,
                layout=positions,
                population_colors=population_colors,
                population_mapping=population_mapping,
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
            logger.error(f'Error processing H number mapping: {e}')
            logger.error(traceback.format_exc())
            return (
                None,
                dbc.Alert(
                    [
                        html.Strong('❌ Error processing mapping'),
                        html.Br(),
                        str(e),
                    ],
                    color='danger',
                    dismissable=True,
                ),
                dash.no_update,
            )

    # Clientside callback to auto-fit network when elements are updated
    # Only fit when not in manual edit mode
    app.clientside_callback(
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

    # Clientside callback to handle window resize and adjust network layout
    app.clientside_callback(
        """
        function(elements) {
            // Set up window resize listener once
            if (!window.cytoscapeResizeSetup) {
                window.cytoscapeResizeSetup = true;

                let resizeTimeout;
                window.addEventListener('resize', function() {
                    // Debounce the resize event
                    clearTimeout(resizeTimeout);
                    resizeTimeout = setTimeout(function() {
                        try {
                            const cy = document.getElementById('network-graph')._cyreg.cy;
                            if (cy) {
                                // Resize the cytoscape instance to fit new container size
                                cy.resize();
                                // Re-fit to maintain visibility
                                cy.fit(null, 50);
                            }
                        } catch (e) {
                            console.log('Could not resize network:', e);
                        }
                    }, 250); // Wait 250ms after resize stops
                });
            }
            return window.dash_clientside.no_update;
        }
        """,
        Output('window-size-store', 'data'),
        Input('network-graph', 'elements'),
        prevent_initial_call=True,
    )

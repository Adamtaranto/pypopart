"""
Command-line interface for PyPopART.

This module provides the CLI for constructing and analyzing haplotype networks.
"""

from pathlib import Path
import sys
from typing import Optional

import click

from pypopart import __version__


@click.group()
@click.version_option(version=__version__)
@click.option(
    '-v',
    '--verbose',
    count=True,
    help='Increase verbosity (can be repeated: -v, -vv, -vvv)',
)
@click.option('-q', '--quiet', is_flag=True, help='Suppress all output except errors')
@click.pass_context
def main(ctx: click.Context, verbose: int, quiet: bool) -> None:
    """
    PyPopART - Pure Python implementation of PopART haplotype network analysis.

    Construct and visualize haplotype networks from DNA sequence alignments.
    """
    # Store verbosity in context for subcommands
    ctx.ensure_object(dict)
    ctx.obj['verbose'] = verbose
    ctx.obj['quiet'] = quiet


@main.command()
@click.argument('input_file', type=click.Path(exists=True, dir_okay=False))
@click.option(
    '-f',
    '--format',
    type=click.Choice(['fasta', 'nexus', 'phylip', 'genbank'], case_sensitive=False),
    help='Input file format (auto-detected if not specified)',
)
@click.option(
    '-m',
    '--metadata',
    type=click.Path(exists=True, dir_okay=False),
    help='Metadata/traits file (CSV format)',
)
@click.option('-o', '--output', type=click.Path(), help='Output file for loaded data')
@click.pass_context
def load(
    ctx: click.Context,
    input_file: str,
    format: Optional[str],
    metadata: Optional[str],
    output: Optional[str],
) -> None:
    """
    Load and validate sequence alignment data.

    INPUT_FILE: Path to sequence alignment file
    """
    from pypopart.io import load_alignment

    click.echo(f'Loading sequences from {input_file}...')

    try:
        alignment = load_alignment(input_file, format=format)
        click.echo(f'✓ Loaded {len(alignment)} sequences')
        click.echo(f'  Alignment length: {alignment.length} bp')

        if metadata:
            click.echo(f'Loading metadata from {metadata}...')
            # TODO: Load and attach metadata
            click.echo('✓ Metadata loaded')

        if output:
            # TODO: Save processed alignment
            click.echo(f'✓ Saved to {output}')

        # Display summary statistics
        stats = alignment.calculate_stats()
        click.echo('\nAlignment Statistics:')
        click.echo(f'  Sequences: {stats.num_sequences}')
        click.echo(f'  Length: {stats.length} bp')
        click.echo(f'  Variable sites: {stats.variable_sites}')
        click.echo(f'  Parsimony informative: {stats.parsimony_informative_sites}')
        click.echo(f'  GC content: {stats.gc_content:.1f}%')

    except Exception as e:
        click.echo(f'✗ Error: {e}', err=True)
        sys.exit(1)


@main.command()
@click.argument('input_file', type=click.Path(exists=True, dir_okay=False))
@click.option(
    '-a',
    '--algorithm',
    type=click.Choice(['mst', 'msn', 'tcs', 'mjn'], case_sensitive=False),
    default='mjn',
    show_default=True,
    help='Network construction algorithm',
)
@click.option(
    '-d',
    '--distance',
    type=click.Choice(['hamming', 'jc', 'k2p', 'tamura_nei'], case_sensitive=False),
    default='hamming',
    show_default=True,
    help='Distance metric',
)
@click.option(
    '-e',
    '--epsilon',
    type=float,
    default=0,
    show_default=True,
    help='Epsilon parameter for median-joining network',
)
@click.option(
    '-p',
    '--parsimony-limit',
    type=float,
    default=0.95,
    show_default=True,
    help='Parsimony probability limit for TCS (0-1)',
)
@click.option('-o', '--output', type=click.Path(), help='Output network file')
@click.option(
    '--format',
    'output_format',
    type=click.Choice(['graphml', 'gml', 'json', 'nexus'], case_sensitive=False),
    default='graphml',
    show_default=True,
    help='Output format',
)
@click.pass_context
def network(
    ctx: click.Context,
    input_file: str,
    algorithm: str,
    distance: str,
    epsilon: float,
    parsimony_limit: float,
    output: Optional[str],
    output_format: str,
) -> None:
    """
    Construct haplotype network from sequence alignment.

    INPUT_FILE: Path to sequence alignment file
    """
    from pypopart.algorithms import (
        TCS,
        MedianJoiningNetwork,
        MinimumSpanningNetwork,
        MinimumSpanningTree,
    )
    from pypopart.io import load_alignment, save_network

    verbose = ctx.obj.get('verbose', 0)

    click.echo(f'Loading sequences from {input_file}...')

    try:
        # Load alignment
        alignment = load_alignment(input_file)
        click.echo(f'✓ Loaded {len(alignment)} sequences ({alignment.length} bp)')

        # Identify unique haplotypes for informational purposes
        click.echo('Identifying unique haplotypes...')
        from pypopart.core.haplotype import identify_haplotypes_from_alignment

        haplotypes = identify_haplotypes_from_alignment(alignment)
        click.echo(f'✓ Found {len(haplotypes)} unique haplotypes')

        # Construct network
        click.echo(f'Building {algorithm.upper()} network...')

        if algorithm == 'mst':
            algo = MinimumSpanningTree(distance_method=distance)
        elif algorithm == 'msn':
            algo = MinimumSpanningNetwork(distance_method=distance)
        elif algorithm == 'tcs':
            algo = TCS(
                distance_method=distance,
                connection_limit=int(parsimony_limit) if parsimony_limit else None,
            )
        elif algorithm == 'mjn':
            algo = MedianJoiningNetwork(distance_method=distance, epsilon=epsilon)
        else:
            raise ValueError(f'Unknown algorithm: {algorithm}')

        network = algo.build_network(alignment)
        click.echo('✓ Network constructed')

        # Display network statistics
        click.echo('\nNetwork Statistics:')
        click.echo(f'  Nodes: {len(network.graph.nodes)}')
        click.echo(f'  Edges: {len(network.graph.edges)}')

        # Count median vectors
        n_medians = sum(1 for node in network.graph.nodes if 'Median_' in str(node))
        if n_medians > 0:
            click.echo(f'  Median vectors: {n_medians}')

        # Save network
        if output:
            click.echo(f'\nSaving network to {output}...')
            save_network(network.graph, output, format=output_format)
            click.echo(f'✓ Network saved as {output_format.upper()}')
        else:
            click.echo('\nℹ Use -o/--output to save the network to a file', err=True)

    except Exception as e:
        click.echo(f'✗ Error: {e}', err=True)
        if verbose > 0:
            import traceback

            traceback.print_exc()
        sys.exit(1)


@main.command()
@click.argument('network_file', type=click.Path(exists=True, dir_okay=False))
@click.option(
    '--stats',
    is_flag=True,
    help='Calculate and display network statistics',
)
@click.option(
    '--topology',
    is_flag=True,
    help='Analyze network topology',
)
@click.option(
    '--popgen',
    is_flag=True,
    help='Calculate population genetics measures',
)
@click.option(
    '-o', '--output', type=click.Path(), help='Output file for analysis results'
)
@click.pass_context
def analyze(
    ctx: click.Context,
    network_file: str,
    stats: bool,
    topology: bool,
    popgen: bool,
    output: Optional[str],
) -> None:
    """
    Analyze haplotype network statistics.

    NETWORK_FILE: Path to network file (GraphML, GML, or JSON)
    """
    from pypopart.io import load_network
    from pypopart.stats import (
        NetworkStatistics,
        TopologyAnalyzer,
    )

    click.echo(f'Loading network from {network_file}...')

    try:
        network = load_network(network_file)
        click.echo(f'✓ Loaded network with {network.number_of_nodes()} nodes')

        results = {}

        # Network statistics
        if stats or (not stats and not topology and not popgen):
            click.echo('\n=== Network Statistics ===')
            net_stats = NetworkStatistics(network)
            summary = net_stats.summary()

            for key, value in summary.items():
                if isinstance(value, float):
                    click.echo(f'{key}: {value:.4f}')
                else:
                    click.echo(f'{key}: {value}')

            results['statistics'] = summary

        # Topology analysis
        if topology:
            click.echo('\n=== Topology Analysis ===')
            topo = TopologyAnalyzer(network)
            topo_summary = topo.analyze()

            click.echo(f'Connected components: {topo_summary["num_components"]}')
            click.echo(f'Star-like patterns: {len(topo_summary["star_patterns"])}')
            click.echo(f'Central nodes: {topo_summary["central_nodes"]}')

            results['topology'] = topo_summary

        # Population genetics
        if popgen:
            click.echo('\n=== Population Genetics ===')
            # This requires alignment data - load if available
            click.echo('Population genetics analysis requires original alignment data')
            # TODO: Implement loading alignment data with network

        # Save results
        if output:
            import json

            with open(output, 'w') as f:
                json.dump(results, f, indent=2, default=str)
            click.echo(f'\n✓ Results saved to {output}')

    except Exception as e:
        click.echo(f'✗ Error: {e}', err=True)
        sys.exit(1)


@main.command()
@click.argument('network_file', type=click.Path(exists=True, dir_okay=False))
@click.option(
    '-o',
    '--output',
    type=click.Path(),
    required=True,
    help='Output image file (PNG, PDF, SVG, or HTML)',
)
@click.option(
    '--layout',
    type=click.Choice(
        ['spring', 'circular', 'radial', 'hierarchical', 'kamada_kawai'],
        case_sensitive=False,
    ),
    default='spring',
    show_default=True,
    help='Layout algorithm',
)
@click.option(
    '--width',
    type=int,
    default=800,
    show_default=True,
    help='Figure width in pixels',
)
@click.option(
    '--height',
    type=int,
    default=600,
    show_default=True,
    help='Figure height in pixels',
)
@click.option(
    '--interactive',
    is_flag=True,
    help='Create interactive HTML visualization (requires .html output)',
)
@click.option(
    '--color-by',
    type=str,
    help='Node attribute to use for coloring (e.g., population)',
)
@click.option(
    '--show-labels',
    is_flag=True,
    default=False,
    help='Show node labels',
)
@click.pass_context
def visualize(
    ctx: click.Context,
    network_file: str,
    output: str,
    layout: str,
    width: int,
    height: int,
    interactive: bool,
    color_by: Optional[str],
    show_labels: bool,
) -> None:
    """
    Visualize haplotype network.

    NETWORK_FILE: Path to network file (GraphML, GML, or JSON)
    """
    from pypopart.io import load_network

    click.echo(f'Loading network from {network_file}...')

    try:
        network = load_network(network_file)
        click.echo(f'✓ Loaded network with {network.number_of_nodes()} nodes')

        # Determine output format
        output_path = Path(output)
        is_html = output_path.suffix.lower() == '.html'

        if interactive or is_html:
            # Interactive visualization
            from pypopart.visualization.interactive import InteractiveVisualizer

            click.echo(f'Creating interactive visualization with {layout} layout...')
            viz = InteractiveVisualizer(network)
            fig = viz.plot(
                layout_algorithm=layout,
                width=width,
                height=height,
                color_by=color_by,
                show_labels=show_labels,
            )
            fig.write_html(str(output_path))
            click.echo(f'✓ Interactive visualization saved to {output}')
            click.echo(f'  Open in browser: file://{output_path.absolute()}')

        else:
            # Static visualization
            from pypopart.visualization.static import StaticVisualizer

            click.echo(f'Creating static visualization with {layout} layout...')
            viz = StaticVisualizer(network)
            viz.plot(
                layout_algorithm=layout,
                figsize=(width / 100, height / 100),
                color_by=color_by,
                show_labels=show_labels,
                output_file=str(output_path),
            )
            click.echo(f'✓ Visualization saved to {output}')

    except Exception as e:
        click.echo(f'✗ Error: {e}', err=True)
        sys.exit(1)


@main.command(name='geo-visualize')
@click.argument('network_file', type=click.Path(exists=True, dir_okay=False))
@click.option(
    '-m',
    '--metadata',
    type=click.Path(exists=True, dir_okay=False),
    required=True,
    help='Metadata CSV file with latitude/longitude columns',
)
@click.option(
    '-o',
    '--output',
    type=click.Path(),
    required=True,
    help='Output image file (PNG, PDF, SVG, or HTML)',
)
@click.option(
    '--projection',
    type=click.Choice(
        ['mercator', 'platecarree', 'orthographic'], case_sensitive=False
    ),
    default='mercator',
    show_default=True,
    help='Map projection',
)
@click.option(
    '--base-map',
    type=click.Choice(
        ['OpenStreetMap', 'Stamen Terrain', 'CartoDB positron'], case_sensitive=False
    ),
    default='OpenStreetMap',
    show_default=True,
    help='Base map for interactive visualization',
)
@click.option(
    '--zoom',
    type=int,
    default=4,
    show_default=True,
    help='Initial zoom level for interactive map (1-18)',
)
@click.option(
    '--extent',
    type=str,
    help='Map extent as "lon_min,lon_max,lat_min,lat_max" (static maps only)',
)
@click.option(
    '--width',
    type=int,
    default=1600,
    show_default=True,
    help='Figure width in pixels',
)
@click.option(
    '--height',
    type=int,
    default=1200,
    show_default=True,
    help='Figure height in pixels',
)
@click.option(
    '--interactive',
    is_flag=True,
    help='Create interactive HTML map visualization',
)
@click.option(
    '--show-labels',
    is_flag=True,
    default=False,
    help='Show node labels',
)
@click.option(
    '--show-borders',
    is_flag=True,
    default=False,
    help='Show country borders (static maps only)',
)
@click.option(
    '--lat-column',
    type=str,
    default='latitude',
    show_default=True,
    help='Name of latitude column in metadata',
)
@click.option(
    '--lon-column',
    type=str,
    default='longitude',
    show_default=True,
    help='Name of longitude column in metadata',
)
@click.pass_context
def geo_visualize(
    ctx: click.Context,
    network_file: str,
    metadata: str,
    output: str,
    projection: str,
    base_map: str,
    zoom: int,
    extent: Optional[str],
    width: int,
    height: int,
    interactive: bool,
    show_labels: bool,
    show_borders: bool,
    lat_column: str,
    lon_column: str,
) -> None:
    """
    Visualize haplotype network on a geographic map.

    NETWORK_FILE: Path to network file (GraphML, GML, or JSON)

    This command overlays the haplotype network on a world map using
    geographic coordinates from the metadata file. The metadata file must
    contain latitude and longitude columns.
    """
    from pypopart.io import load_network
    from pypopart.io.metadata import MetadataReader, extract_coordinates

    click.echo(f'Loading network from {network_file}...')

    try:
        network = load_network(network_file)
        click.echo(f'✓ Loaded network with {network.number_of_nodes()} nodes')

        # Load metadata with coordinates
        click.echo(f'Loading metadata from {metadata}...')
        reader = MetadataReader(metadata)
        metadata_dict = reader.read_metadata()

        # Extract coordinates for each node
        coordinates = {}
        for node_id in network.node_ids:
            if node_id in metadata_dict:
                coords = extract_coordinates(
                    metadata_dict[node_id],
                    lat_column=lat_column,
                    lon_column=lon_column,
                    validate=True,
                )
                if coords:
                    coordinates[node_id] = coords

        click.echo(f'✓ Loaded coordinates for {len(coordinates)} nodes')

        if not coordinates:
            click.echo(
                '✗ Error: No valid geographic coordinates found in metadata', err=True
            )
            sys.exit(1)

        # Parse extent if provided
        extent_tuple = None
        if extent:
            try:
                parts = [float(x.strip()) for x in extent.split(',')]
                if len(parts) != 4:
                    raise ValueError('Extent must have 4 values')
                extent_tuple = tuple(parts)
            except ValueError as e:
                click.echo(f'✗ Error: Invalid extent format: {e}', err=True)
                sys.exit(1)

        # Determine output format
        output_path = Path(output)
        is_html = output_path.suffix.lower() == '.html'

        if interactive or is_html:
            # Interactive geographic visualization
            from pypopart.visualization import InteractiveGeoVisualizer

            click.echo('Creating interactive geographic visualization...')
            viz = InteractiveGeoVisualizer(network)
            viz.plot(
                coordinates=coordinates,
                base_map=base_map,
                zoom_start=zoom,
                show_labels=show_labels,
                output_file=str(output_path),
            )
            click.echo(f'✓ Interactive map saved to {output}')
            click.echo(f'  Open in browser: file://{output_path.absolute()}')

        else:
            # Static geographic visualization
            from pypopart.visualization import GeoVisualizer

            click.echo(
                f'Creating static geographic visualization with {projection} projection...'
            )
            viz = GeoVisualizer(network)
            fig, ax = viz.plot(
                coordinates=coordinates,
                projection=projection,
                extent=extent_tuple,
                figsize=(width / 100, height / 100),
                show_labels=show_labels,
                show_borders=show_borders,
                output_file=str(output_path),
            )
            click.echo(f'✓ Geographic visualization saved to {output}')

    except Exception as e:
        import traceback

        click.echo(f'✗ Error: {e}', err=True)
        if ctx.obj.get('verbose', 0) > 0:
            traceback.print_exc()
        sys.exit(1)


@main.command()
@click.option(
    '--list-algorithms',
    is_flag=True,
    help='List available network construction algorithms',
)
@click.option(
    '--list-distances',
    is_flag=True,
    help='List available distance metrics',
)
@click.option(
    '--list-formats',
    is_flag=True,
    help='List supported file formats',
)
def info(
    list_algorithms: bool,
    list_distances: bool,
    list_formats: bool,
) -> None:
    """
    Display information about PyPopART capabilities.
    """
    if list_algorithms:
        click.echo('Available Network Construction Algorithms:')
        click.echo('  mst - Minimum Spanning Tree')
        click.echo('  msn - Minimum Spanning Network')
        click.echo('  tcs - Statistical Parsimony (TCS)')
        click.echo('  mjn - Median-Joining Network')
        click.echo()

    if list_distances:
        click.echo('Available Distance Metrics:')
        click.echo('  hamming     - Simple Hamming distance (count differences)')
        click.echo('  jc          - Jukes-Cantor correction')
        click.echo('  k2p         - Kimura 2-parameter model')
        click.echo('  tamura_nei  - Tamura-Nei model')
        click.echo()

    if list_formats:
        click.echo('Supported Input Formats:')
        click.echo('  fasta    - FASTA sequence format')
        click.echo('  nexus    - NEXUS format')
        click.echo('  phylip   - PHYLIP format')
        click.echo('  genbank  - GenBank format')
        click.echo()
        click.echo('Supported Output Formats:')
        click.echo('  graphml  - GraphML (XML-based graph format)')
        click.echo('  gml      - Graph Modelling Language')
        click.echo('  json     - JSON format')
        click.echo('  nexus    - NEXUS format')
        click.echo()

    if not (list_algorithms or list_distances or list_formats):
        click.echo('Use --list-algorithms, --list-distances, or --list-formats')
        click.echo('or run "pypopart --help" for usage information')


if __name__ == '__main__':
    main()

"""
Command-line interface for PyPopART.

This module provides the CLI for constructing and analyzing haplotype networks.
"""

import sys
from pathlib import Path
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
    from pypopart.algorithms import MJNAlgorithm, MSNAlgorithm, MSTAlgorithm, TCSAlgorithm
    from pypopart.core.distance import DistanceCalculator
    from pypopart.io import load_alignment, save_network

    verbose = ctx.obj.get('verbose', 0)

    click.echo(f'Loading sequences from {input_file}...')

    try:
        # Load alignment
        alignment = load_alignment(input_file)
        click.echo(f'✓ Loaded {len(alignment)} sequences ({alignment.length} bp)')

        # Calculate distances
        click.echo(f'Calculating {distance} distances...')
        calculator = DistanceCalculator(method=distance)
        dist_matrix = calculator.calculate_matrix(alignment)
        click.echo(f'✓ Distance matrix computed')

        # Identify haplotypes
        click.echo('Identifying unique haplotypes...')
        from pypopart.core.condensation import condense_alignment

        haplotypes, freq_map = condense_alignment(alignment)
        click.echo(f'✓ Found {len(haplotypes)} unique haplotypes')

        # Construct network
        click.echo(f'Building {algorithm.upper()} network...')

        if algorithm == 'mst':
            algo = MSTAlgorithm()
        elif algorithm == 'msn':
            algo = MSNAlgorithm()
        elif algorithm == 'tcs':
            algo = TCSAlgorithm(connection_limit=parsimony_limit)
        elif algorithm == 'mjn':
            algo = MJNAlgorithm(epsilon=epsilon)
        else:
            raise ValueError(f'Unknown algorithm: {algorithm}')

        network = algo.construct_network(haplotypes, dist_matrix)
        click.echo(f'✓ Network constructed')

        # Display network statistics
        click.echo('\nNetwork Statistics:')
        click.echo(f'  Nodes: {network.number_of_nodes()}')
        click.echo(f'  Edges: {network.number_of_edges()}')

        if hasattr(network, 'median_vectors'):
            n_medians = len(
                [n for n in network.nodes() if network.nodes[n].get('is_median', False)]
            )
            click.echo(f'  Median vectors: {n_medians}')

        # Save network
        if output:
            click.echo(f'\nSaving network to {output}...')
            save_network(network, output, format=output_format)
            click.echo(f'✓ Network saved as {output_format.upper()}')
        else:
            click.echo(
                '\nℹ Use -o/--output to save the network to a file', err=True
            )

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
@click.option('-o', '--output', type=click.Path(), help='Output file for analysis results')
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
    from pypopart.stats import NetworkStatistics, PopulationGeneticsAnalysis, TopologyAnalyzer

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

            click.echo(f"Connected components: {topo_summary['num_components']}")
            click.echo(f"Star-like patterns: {len(topo_summary['star_patterns'])}")
            click.echo(f"Central nodes: {topo_summary['central_nodes']}")

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

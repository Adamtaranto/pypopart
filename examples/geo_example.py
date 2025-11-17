"""
Example: Geographic visualization of haplotype networks.

This script demonstrates how to create geographic visualizations
of haplotype networks using PyPopART.
"""

from pathlib import Path

from pypopart.io import load_alignment
from pypopart.io.metadata import MetadataReader, extract_coordinates
from pypopart.algorithms import MinimumSpanningNetwork
from pypopart.core.distance import DistanceCalculator
from pypopart.core.haplotype import identify_haplotypes_from_alignment
from pypopart.visualization import GeoVisualizer, InteractiveGeoVisualizer


def main():
    """Run geographic visualization example."""
    # Setup paths
    data_dir = Path(__file__).parent / 'geo_data'
    sequences_file = data_dir / 'sample_sequences.fasta'
    metadata_file = data_dir / 'sample_metadata.csv'
    output_dir = Path(__file__).parent / 'output'
    output_dir.mkdir(exist_ok=True)

    print('=' * 70)
    print('Geographic Haplotype Network Visualization Example')
    print('=' * 70)
    print()

    # 1. Load sequences
    print('Step 1: Loading sequences...')
    alignment = load_alignment(sequences_file)
    print(f'  ✓ Loaded {len(alignment)} sequences')
    print()

    # 2. Load metadata with coordinates
    print('Step 2: Loading metadata with geographic coordinates...')
    reader = MetadataReader(metadata_file)
    metadata_dict = reader.read_metadata()
    print(f'  ✓ Loaded metadata for {len(metadata_dict)} samples')
    print()

    # 3. Identify haplotypes
    print('Step 3: Identifying unique haplotypes...')
    haplotypes = identify_haplotypes_from_alignment(alignment)
    print(f'  ✓ Found {len(haplotypes)} unique haplotypes')
    print()

    # 4. Calculate distance matrix
    print('Step 4: Calculating genetic distances...')
    calculator = DistanceCalculator(method='hamming')
    dist_matrix = calculator.calculate_matrix(alignment)
    print(f'  ✓ Calculated {dist_matrix.shape[0]}x{dist_matrix.shape[1]} distance matrix')
    print()

    # 5. Construct network
    print('Step 5: Constructing Minimum Spanning Network...')
    msn = MinimumSpanningNetwork()
    network = msn.construct_network(haplotypes, dist_matrix)
    print(f'  ✓ Network has {network.num_nodes} nodes and {network.num_edges} edges')
    print()

    # 6. Extract coordinates and attach to network
    print('Step 6: Extracting geographic coordinates...')
    coordinates = {}
    for hap_id, hap in zip(network.node_ids, haplotypes):
        # Get all sample IDs for this haplotype
        for sample_id in hap.sample_ids:
            if sample_id in metadata_dict:
                coords = extract_coordinates(
                    metadata_dict[sample_id],
                    lat_column='latitude',
                    lon_column='longitude',
                    validate=True,
                )
                if coords:
                    coordinates[hap_id] = coords
                    # Also attach to node metadata
                    network._graph.nodes[hap_id]['latitude'] = str(coords[0])
                    network._graph.nodes[hap_id]['longitude'] = str(coords[1])
                    break  # Use first sample's coordinates

    print(f'  ✓ Extracted coordinates for {len(coordinates)} haplotypes')
    print()

    # 7. Create static geographic visualization
    print('Step 7: Creating static geographic visualizations...')

    # Mercator projection
    print('  Creating Mercator projection map...')
    viz_mercator = GeoVisualizer(network)
    fig, ax = viz_mercator.plot(
        coordinates=coordinates,
        projection='mercator',
        show_labels=True,
        show_borders=True,
        title='Haplotype Network - Mercator Projection',
        figsize=(16, 12),
        output_file=str(output_dir / 'geo_network_mercator.png'),
    )
    print(f'    ✓ Saved to {output_dir / "geo_network_mercator.png"}')

    # PlateCarree projection
    print('  Creating PlateCarree (equirectangular) map...')
    viz_plate = GeoVisualizer(network)
    fig, ax = viz_plate.plot(
        coordinates=coordinates,
        projection='platecarree',
        show_labels=True,
        show_borders=True,
        title='Haplotype Network - PlateCarree Projection',
        figsize=(16, 12),
        output_file=str(output_dir / 'geo_network_platecarree.png'),
    )
    print(f'    ✓ Saved to {output_dir / "geo_network_platecarree.png"}')
    print()

    # 8. Create interactive geographic visualization
    print('Step 8: Creating interactive geographic visualization...')
    viz_interactive = InteractiveGeoVisualizer(network)
    map_obj = viz_interactive.plot(
        coordinates=coordinates,
        base_map='OpenStreetMap',
        zoom_start=2,
        show_labels=True,
        output_file=str(output_dir / 'geo_network_interactive.html'),
    )
    print(f'  ✓ Saved interactive map to {output_dir / "geo_network_interactive.html"}')
    print()

    print('=' * 70)
    print('Geographic visualization complete!')
    print('=' * 70)
    print()
    print('Output files:')
    print(f'  - Static (Mercator): {output_dir / "geo_network_mercator.png"}')
    print(f'  - Static (PlateCarree): {output_dir / "geo_network_platecarree.png"}')
    print(f'  - Interactive: {output_dir / "geo_network_interactive.html"}')
    print()
    print('Open the HTML file in a web browser to explore the interactive map!')


if __name__ == '__main__':
    main()

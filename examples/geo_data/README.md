# Geographic Visualization Example

This directory contains sample data for demonstrating PyPopART's geographic visualization capabilities.

## Files

- `sample_sequences.fasta`: Sample DNA sequences from different geographic locations
- `sample_metadata.csv`: Metadata file with geographic coordinates (latitude/longitude)

## Sample Metadata Format

The metadata CSV file must contain the following columns:
- `id`: Sequence identifier (must match FASTA headers)
- `latitude`: Latitude in decimal degrees (-90 to 90)
- `longitude`: Longitude in decimal degrees (-180 to 180)

Optional columns:
- `population`: Population or group name
- `location`: Human-readable location name
- Any other metadata fields

## Usage

### Using the Python API

See `../geo_example.py` for a complete example:

```python
from pypopart.io import load_alignment
from pypopart.io.metadata import MetadataReader, extract_coordinates
from pypopart.algorithms import MinimumSpanningNetwork
from pypopart.visualization import GeoVisualizer, InteractiveGeoVisualizer

# Load data
alignment = load_alignment('sample_sequences.fasta')
reader = MetadataReader('sample_metadata.csv')
metadata_dict = reader.read_metadata()

# Build network (see example for details)
# ...

# Extract coordinates
coordinates = {}
for node_id in network.node_ids:
    if node_id in metadata_dict:
        coords = extract_coordinates(metadata_dict[node_id])
        if coords:
            coordinates[node_id] = coords

# Create static geographic visualization
viz = GeoVisualizer(network)
fig, ax = viz.plot(
    coordinates=coordinates,
    projection='mercator',
    show_labels=True,
    show_borders=True,
    output_file='geo_network.png',
)

# Create interactive map
viz_interactive = InteractiveGeoVisualizer(network)
map_obj = viz_interactive.plot(
    coordinates=coordinates,
    base_map='OpenStreetMap',
    zoom_start=2,
    output_file='geo_network.html',
)
```

### Using the CLI

First, construct a network from the sequences:

```bash
pypopart network sample_sequences.fasta -o network.graphml
```

Then create a geographic visualization:

```bash
# Static map (PNG/PDF/SVG)
pypopart geo-visualize network.graphml \
    -m sample_metadata.csv \
    -o geo_network.png \
    --projection mercator \
    --show-labels \
    --show-borders

# Interactive map (HTML)
pypopart geo-visualize network.graphml \
    -m sample_metadata.csv \
    -o geo_network.html \
    --interactive \
    --base-map OpenStreetMap \
    --zoom 2
```

## Projections

PyPopART supports multiple map projections:

- **mercator**: Web Mercator projection (preserves angles, distorts area at poles)
- **platecarree**: Equirectangular projection (simple lat/lon mapping)
- **orthographic**: Orthographic projection (3D globe view)

## Interactive Map Base Layers

For interactive HTML maps, you can choose from several base map styles:

- **OpenStreetMap**: Standard OpenStreetMap tiles
- **Stamen Terrain**: Terrain-focused map with hill shading
- **CartoDB positron**: Clean, minimal base map

## Handling Multiple Locations

When a single haplotype is found in multiple geographic locations, PyPopART offers several strategies:

1. **Centroid** (default): Place the node at the geographic centroid of all locations
2. **First**: Use the first location encountered
3. **Jitter**: Add small random offsets to distinguish overlapping nodes

Configure this in the Python API using the `handle_multiple_locations` parameter.

## Example Output

Running the example script creates:
- Static maps showing the haplotype network overlaid on a world map
- Interactive HTML map with clickable nodes showing haplotype details
- Network edges showing genetic relationships between geographically distributed haplotypes

## Tips

- Use higher zoom levels (8-12) for regional studies
- Use lower zoom levels (1-4) for global datasets
- Add borders and coastlines for better geographic context
- Use the PlateCarree projection for polar regions
- Use Mercator for equatorial regions

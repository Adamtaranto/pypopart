# Geographic Visualization Features

## Overview

PyPopART now supports geographic visualization of haplotype networks, allowing you to overlay networks on world maps using latitude/longitude coordinates. This feature is useful for studying the geographic distribution of genetic variation and understanding patterns of dispersal and population structure.

## Features

### 1. Geographic Coordinate Support

- Parse latitude/longitude coordinates from CSV metadata files
- Validate coordinate ranges (latitude: -90 to 90, longitude: -180 to 180)
- Support multiple coordinate formats (decimal degrees, with/without degree symbols)
- Handle missing coordinates gracefully

### 2. Geographic Layout Algorithm

The `GeographicLayout` class positions network nodes based on real geographic coordinates:

- **Multiple Projections**: Mercator, PlateCarree (equirectangular), Orthographic (3D globe)
- **Coordinate Projection**: Transform lat/lon to appropriate map coordinates
- **Missing Data Handling**: Nodes without coordinates are placed using default layouts
- **Jitter Support**: Add small random offsets to distinguish overlapping nodes

### 3. Static Map Visualization

The `GeoVisualizer` class creates publication-quality static maps using matplotlib + cartopy:

- Customizable map projections
- Configurable map extent (zoom to specific regions)
- Optional coastlines, land features, and country borders
- Support for PNG, PDF, and SVG output formats
- Full control over node sizes, colors, and labels
- Edge rendering showing genetic relationships

### 4. Interactive Map Visualization

The `InteractiveGeoVisualizer` class creates interactive HTML maps using folium:

- Multiple base map options (OpenStreetMap, Stamen Terrain, CartoDB positron)
- Adjustable zoom levels
- Clickable nodes with popup information
- Interactive panning and zooming
- Edges showing genetic relationships
- Export to standalone HTML files

### 5. CLI Integration

New `geo-visualize` command provides easy access to geographic visualization:

```bash
# Static map
pypopart geo-visualize network.graphml \
    -m metadata.csv \
    -o geo_network.png \
    --projection mercator \
    --show-labels \
    --show-borders

# Interactive map
pypopart geo-visualize network.graphml \
    -m metadata.csv \
    -o geo_network.html \
    --interactive \
    --base-map OpenStreetMap \
    --zoom 4
```

## Use Cases

1. **Phylogeography**: Study the geographic distribution of haplotypes and their evolutionary relationships
2. **Population Structure**: Visualize genetic structure across geographic space
3. **Migration Patterns**: Identify patterns of gene flow and dispersal
4. **Sampling Design**: Plan future sampling efforts based on geographic gaps
5. **Presentation**: Create compelling visualizations for publications and presentations

## Requirements

Geographic visualization requires two additional dependencies:

- **cartopy**: For static map rendering with matplotlib
- **folium**: For interactive HTML map generation

These are included in the main PyPopART dependencies and will be installed automatically.

## Examples

See the `examples/geo_data/` directory for:
- Sample sequences with geographic distribution
- Sample metadata with coordinates
- Complete working example script (`examples/geo_example.py`)
- Detailed README with usage instructions

## Map Projections

### Mercator
- Preserves angles (conformal)
- Distorts area, especially near poles
- Good for equatorial and mid-latitude regions
- Standard web mapping projection

### PlateCarree (Equirectangular)
- Simple lat/lon mapping (no transformation)
- Equal spacing of parallels and meridians
- Good for global views and polar regions
- Simple but not conformal or equal-area

### Orthographic
- 3D globe perspective
- One hemisphere visible at a time
- Good for presentations and publications
- Natural appearance of Earth from space

## Handling Multiple Locations

When a single haplotype is found in multiple geographic locations, several strategies are available:

1. **Centroid** (default): Place node at geographic center of all locations
2. **First**: Use the first location encountered in the data
3. **Jitter**: Add small random offsets to distinguish overlapping nodes

Configure this in the Python API with the `handle_multiple_locations` parameter.

## Tips for Best Results

1. **Coordinate Quality**: Ensure accurate lat/lon coordinates in your metadata
2. **Projection Choice**: Choose projection based on your study region:
   - Global studies: PlateCarree or Orthographic
   - Regional studies: Mercator
   - Polar studies: PlateCarree
3. **Zoom Level**: Adjust based on geographic extent of your data
4. **Base Maps**: Choose map style based on your needs:
   - Scientific: CartoDB positron (clean, minimal)
   - Presentations: OpenStreetMap (detailed)
   - Topography: Stamen Terrain (terrain features)
5. **Node Colors**: Use population colors to show structure
6. **Labels**: Enable for small networks, disable for large ones

## Python API Example

```python
from pypopart.io import load_alignment
from pypopart.io.metadata import MetadataReader, extract_coordinates
from pypopart.algorithms import MinimumSpanningNetwork
from pypopart.visualization import GeoVisualizer, InteractiveGeoVisualizer

# Load data
alignment = load_alignment('sequences.fasta')
reader = MetadataReader('metadata.csv')
metadata_dict = reader.read_metadata()

# Build network
# ... (construct network as usual)

# Extract coordinates
coordinates = {}
for node_id in network.node_ids:
    if node_id in metadata_dict:
        coords = extract_coordinates(metadata_dict[node_id])
        if coords:
            coordinates[node_id] = coords

# Create static map
viz = GeoVisualizer(network)
fig, ax = viz.plot(
    coordinates=coordinates,
    projection='mercator',
    show_labels=True,
    output_file='geo_network.png'
)

# Create interactive map
viz_interactive = InteractiveGeoVisualizer(network)
map_obj = viz_interactive.plot(
    coordinates=coordinates,
    base_map='OpenStreetMap',
    zoom_start=2,
    output_file='geo_network.html'
)
```

## Contributing

Suggestions for improvements to geographic visualization features are welcome! Please open an issue or submit a pull request.

## Citation

If you use the geographic visualization features in your research, please cite:

```
Taranto, A. (2024). PyPopART: Pure Python implementation of haplotype network analysis.
GitHub repository: https://github.com/adamtaranto/pypopart
```

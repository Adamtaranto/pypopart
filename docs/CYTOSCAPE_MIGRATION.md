# Cytoscape Network Visualization Migration

## Overview

The PyPopART Dash application now uses Dash Cytoscape for network visualization, replacing the previous Plotly-based graphs. This migration provides improved interactivity, manual node repositioning, and better support for complex network features.

## Key Features

### 1. Manual Node Repositioning

**What it does:** Nodes can be dragged and repositioned manually using the mouse.

**How to use:**
- Click and hold on any node in the network
- Drag it to a new position
- Release to place the node
- The new position is automatically saved

**Snap to Grid:** Enable the "Snap to Grid" option in the Layout Options panel to align nodes to a grid when released.

### 2. Population Visualization

#### Single Population Nodes
Nodes containing samples from a single population are displayed with the population's assigned color.

#### Mixed Population Nodes
Nodes containing samples from multiple populations are indicated by:
- A **gold double border** to highlight the mixed composition
- The base color represents the dominant (most frequent) population
- Hover over the node to see the exact population breakdown

### 3. Legend Display

The legend (top right corner) shows:
- **Population colors:** Circle markers with population names
- **Mixed populations indicator:** Gold circle (◉) indicating nodes with samples from multiple populations
- **Median vectors:** Gray square (■) for inferred median vector nodes

### 4. Edge Labels

Edges display the number of mutations between connected haplotypes as numeric labels.

### 5. Interactive Controls

**Zoom:** Use mouse wheel or pinch gesture to zoom in/out

**Pan:** Click and drag on the background to move the entire network

**Select:** Click on nodes to select them (highlighted with red border)

**Hover:** Mouse over nodes or edges to see detailed information

## Technical Details

### Cytoscape Elements Structure

Each node element contains:
- `id`: Unique haplotype identifier (e.g., "H1", "H2")
- `label`: Display label
- `size`: Visual size based on frequency
- `color`: Node color (population-based or default)
- `is_median`: Boolean indicating if this is a median vector
- `has_pie`: Boolean indicating mixed population composition
- `pie_data`: Array of population frequency data (for mixed nodes)

Each edge element contains:
- `id`: Unique edge identifier
- `source`: Source node ID
- `target`: Target node ID
- `weight`: Genetic distance (mutation count)
- `label`: Distance label for display

### Layout Algorithms

Available layout algorithms:
- **Spring (Force-directed):** Physically-simulated layout
- **Circular:** Nodes arranged in a circle
- **Radial:** Radial tree layout
- **Hierarchical:** Hierarchical tree structure
- **Kamada-Kawai:** Force-directed with minimum energy
- **Geographic:** Position based on geographic coordinates (requires metadata)

### Color Generation

When population colors are not specified in metadata:
- Colors are automatically generated using HSV color space
- Evenly distributed hues for maximum distinction
- High saturation and value for vivid colors

## Migration from Plotly

### What Changed

1. **Component Type:** `dcc.Graph` → `cyto.Cytoscape`
2. **Data Format:** Plotly figure traces → Cytoscape elements
3. **Callbacks:** Updated to handle Cytoscape events
4. **Node Positioning:** Now directly editable by users

### What Stayed the Same

- All network algorithms (MST, MSN, TCS, MJN)
- Layout computation methods
- Metadata support
- Export functionality
- Statistics calculations
- Haplotype summary

## API Changes

### InteractiveCytoscapePlotter

New class for creating Cytoscape visualizations:

```python
from pypopart.visualization.cytoscape_plot import (
    InteractiveCytoscapePlotter,
    create_cytoscape_network
)

# Create plotter
plotter = InteractiveCytoscapePlotter(network)

# Generate elements and stylesheet
elements = plotter.create_elements(
    layout=positions,
    population_colors={'PopA': '#FF0000', 'PopB': '#0000FF'},
    show_labels=True,
    show_edge_labels=True,
)

stylesheet = plotter.create_stylesheet(
    population_colors=population_colors,
)

# Or use convenience function
elements, stylesheet = create_cytoscape_network(
    network,
    layout=positions,
    population_colors=population_colors,
)
```

### Convenience Functions

```python
# Auto-generate population colors
colors = plotter.generate_population_colors(['PopA', 'PopB', 'PopC'])

# Create pie chart styles
pie_styles = plotter.create_pie_stylesheet(population_colors)
```

## Backward Compatibility

The original Plotly-based `InteractiveNetworkPlotter` is still available and unchanged:

```python
from pypopart.visualization.interactive_plot import InteractiveNetworkPlotter

# Still works as before
plotter = InteractiveNetworkPlotter(network)
fig = plotter.plot(layout=positions)
```

## Testing

Comprehensive test suite in `tests/unit/test_cytoscape_plot.py`:
- Element creation and structure
- Node sizing and coloring
- Population data handling
- Edge labels
- Stylesheet generation
- Empty and edge cases

Run tests:
```bash
pytest tests/unit/test_cytoscape_plot.py -v
```

## Future Enhancements

Potential improvements for future versions:
- True pie chart rendering using custom Cytoscape.js extensions
- Animation support for layout transitions
- Advanced selection and filtering tools
- Export to Cytoscape.js format
- Custom node shapes for different haplotype types

# PyPopART GUI Usage Guide

## Overview

PyPopART provides a web-based graphical user interface (GUI) built with Dash for interactive haplotype network analysis. The GUI offers an intuitive workflow for:

1. Uploading sequence alignments
2. Computing haplotype networks
3. Visualizing and analyzing results
4. Exporting networks and figures

## Installation

Ensure you have PyPopART installed with GUI dependencies:

```bash
pip install pypopart
# Or install in development mode:
pip install -e .
```

The GUI requires the following additional dependencies (automatically installed):

- dash
- dash-bootstrap-components

## Launching the GUI

### From Python

```python
from pypopart.gui import main

# Launch with default settings (port 8050, debug=False)
main()

# Or customize:
main(debug=True, port=8080)
```

### From Command Line

```bash
# Using the pypopart command (if CLI is implemented)
pypopart gui --port 8050

# Or directly with Python
python -m pypopart.gui.app
```

The GUI will start a local web server. Open your browser to `http://localhost:8050` (or your specified port).

## Using the GUI

### 1. Upload Data

- Click "Select File" in the Upload Data section
- Supported formats:
  - FASTA (.fasta, .fa)
  - NEXUS (.nex, .nexus)
  - PHYLIP (.phy, .phylip)
- After upload, you'll see a success message with sequence count and alignment length

### 2. Configure Algorithm

Choose from four network construction algorithms:

#### MST (Minimum Spanning Tree)

- **Parameters**: Distance metric (Hamming, Jukes-Cantor, Kimura 2-parameter)
- **Best for**: Simple, tree-like relationships

#### MSN (Minimum Spanning Network)

- **Parameters**: Distance metric
- **Best for**: Adding alternative connections at the same distance

#### TCS (Statistical Parsimony)

- **Parameters**: Connection limit (1-20 mutations)
- **Best for**: Recent divergence, intraspecific analysis
- **Default**: 10 mutations

#### MJN (Median-Joining Network)

- **Parameters**: Epsilon value (0 = automatic)
- **Best for**: Complex reticulate relationships, infers median vectors

After selecting an algorithm and parameters, click **Compute Network**.

### 3. Layout Options

Choose how to arrange nodes in the network visualization:

- **Spring (Force-Directed)**: Physics-based layout, good for general use
- **Circular**: Nodes arranged in a circle
- **Radial**: Concentric rings from center
- **Hierarchical**: Tree-like hierarchical arrangement
- **Kamada-Kawai**: Energy-minimization layout

**Snap to Grid**: Check this option to align nodes to a grid for cleaner appearance

Click **Apply Layout** to recompute node positions.

### 4. Viewing Results

The GUI provides three tabs:

#### Network Tab

- **Interactive visualization** of the haplotype network
- Node size: Proportional to haplotype frequency
- Node color: By population (if metadata available)
- Edge thickness: By distance
- **Interactions**:
  - Zoom: Scroll or pinch
  - Pan: Click and drag
  - Hover: View node details
  - Click legend items to show/hide populations

#### Statistics Tab

Displays comprehensive network metrics:

- **Basic Metrics**: Nodes, edges, diameter, clustering, reticulation
- **Diversity Metrics**: Haplotype diversity, Shannon index
- **Central Haplotypes**: Degree, betweenness, and closeness centrality

#### Alignment Tab

- View the uploaded sequence alignment
- Shows first 50 sequences in monospace font
- ID and sequence data aligned for easy comparison

### 5. Export

Export your results in various formats:

- **GraphML**: Network format, opens in Cytoscape, Gephi
- **GML**: Graph Markup Language
- **JSON**: JavaScript Object Notation, for web applications
- **PNG**: Raster image, good for presentations
- **SVG**: Vector image, scalable for publications

Click **Download** to save the file.

## Tips and Best Practices

### Data Preparation

- Ensure sequences are aligned before uploading
- All sequences should have the same length
- Include metadata in NEXUS format for population coloring

### Algorithm Selection

- Start with **MSN** for most datasets (good balance)
- Use **TCS** for closely related sequences (e.g., same species)
- Use **MJN** for complex evolutionary scenarios
- **MST** is fastest but simplest

### Layout Optimization

- Try different layout algorithms to find the clearest visualization
- **Spring layout** works well for most networks
- **Radial layout** highlights star-like patterns
- Enable **Snap to Grid** for publication-ready figures

### Performance

- Large datasets (>100 sequences) may take time to compute
- Progress indicators show when computation is running
- Network visualization is interactive even with many nodes

### Troubleshooting

**File Upload Fails**

- Check file format and extension
- Ensure sequences are properly formatted
- Try a smaller test dataset first

**Network Computation Fails**

- Check that all sequences have the same length
- Reduce connection limit for TCS if very large
- Try a different algorithm

**Visualization Issues**

- Refresh the browser page
- Try a different layout algorithm
- Check browser console for errors (F12)

**Export Doesn't Work**

- Ensure network has been computed
- Check browser's download settings
- Try a different export format

## Advanced Features

### Customization

The GUI can be customized by modifying `src/pypopart/gui/app.py`:

- Adjust default parameters
- Change color schemes
- Add custom layout algorithms
- Modify statistics displayed

### Integration

The GUI components can be used programmatically:

```python
from pypopart.gui.app import PyPopARTApp

app = PyPopARTApp(debug=True, port=8050)
# Access app.app for the Dash application object
# Customize before running
app.run()
```

## Support

For issues, questions, or feature requests:

- GitHub Issues: <https://github.com/adamtaranto/pypopart/issues>
- Documentation: <https://github.com/adamtaranto/pypopart>

## Citation

If you use PyPopART in your research, please cite:

[Citation information to be added]

# PyPopART

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Python Version](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)

**PyPopART** is a pure Python implementation of PopART (Population Analysis with Reticulate Trees) for constructing and visualizing haplotype networks from DNA sequence data.

## Features

- **Multiple Network Algorithms**: MST, MSN, TCS (Statistical Parsimony), and Median-Joining
- **Distance Metrics**: Hamming, Jukes-Cantor, Kimura 2-parameter, Tamura-Nei
- **Comprehensive Analysis**: Network statistics, topology analysis, population genetics measures
- **Rich Visualization**: Static (matplotlib) and interactive (Plotly) network plots
- **Flexible I/O**: Support for FASTA, NEXUS, PHYLIP, GenBank formats
- **Command-Line Interface**: Easy-to-use CLI for all operations
- **Python API**: Programmatic access for custom workflows

## Installation

### From PyPI (when released)

```bash
pip install pypopart
```

### From Source

```bash
git clone https://github.com/adamtaranto/pypopart.git
cd pypopart
pip install -e ".[dev]"
```

### Requirements

- Python 3.9 or higher
- Dependencies: biopython, click, matplotlib, networkx, numpy, pandas, plotly, scipy, scikit-learn, numba

## Quick Start

### Command-Line Interface

#### 1. Load and Validate Sequence Data

```bash
pypopart load sequences.fasta
```

Output:

```text
Loading sequences from sequences.fasta...
✓ Loaded 50 sequences
  Alignment length: 500 bp

Alignment Statistics:
  Sequences: 50
  Length: 500 bp
  Variable sites: 25
  Parsimony informative: 15
  GC content: 48.5%
```

#### 2. Construct a Haplotype Network

```bash
# Median-Joining Network (default)
pypopart network sequences.fasta -o network.graphml

# Statistical Parsimony (TCS)
pypopart network sequences.fasta -a tcs -o network.graphml

# With custom distance metric
pypopart network sequences.fasta -a mjn -d k2p -o network.graphml
```

#### 3. Analyze Network Statistics

```bash
pypopart analyze network.graphml --stats
```

Output:

```text
Loading network from network.graphml...
✓ Loaded network with 12 nodes

=== Network Statistics ===
nodes: 12
edges: 15
diameter: 4
avg_degree: 2.5
clustering_coefficient: 0.3214
reticulation_index: 0.25
```

#### 4. Visualize the Network

```bash
# Static plot (PNG/PDF/SVG)
pypopart visualize network.graphml -o network.png --layout spring

# Interactive HTML
pypopart visualize network.graphml -o network.html --interactive
```

#### 5. Get Information

```bash
# List available algorithms
pypopart info --list-algorithms

# List distance metrics
pypopart info --list-distances

# List supported formats
pypopart info --list-formats
```

### Python API

```python
from pypopart.io import load_alignment
from pypopart.core.distance import DistanceCalculator
from pypopart.core.condensation import condense_alignment
from pypopart.algorithms import MJNAlgorithm
from pypopart.visualization import StaticVisualizer

# Load sequences
alignment = load_alignment('sequences.fasta')

# Calculate distances
calculator = DistanceCalculator(method='k2p')
dist_matrix = calculator.calculate_matrix(alignment)

# Identify unique haplotypes
haplotypes, freq_map = condense_alignment(alignment)

# Construct Median-Joining Network
mjn = MJNAlgorithm(epsilon=0)
network = mjn.construct_network(haplotypes, dist_matrix)

# Visualize
viz = StaticVisualizer(network)
viz.plot(layout_algorithm='spring', output_file='network.png')
```

## Network Construction Algorithms

### Minimum Spanning Tree (MST)

Creates a tree connecting all haplotypes with minimum total distance.

```bash
pypopart network sequences.fasta -a mst -o network.graphml
```

**Use when**: You want the simplest possible network structure.

### Minimum Spanning Network (MSN)

Extends MST by adding alternative connections at equal distance.

```bash
pypopart network sequences.fasta -a msn -o network.graphml
```

**Use when**: You want to show alternative evolutionary pathways at the same genetic distance.

### TCS (Statistical Parsimony)

Connects haplotypes within a parsimony probability limit (default 95%).

```bash
pypopart network sequences.fasta -a tcs -p 0.95 -o network.graphml
```

**Use when**: You want statistically justified connections based on parsimony.

### Median-Joining Network (MJN)

Infers ancestral/median sequences and creates a reticulate network.

```bash
pypopart network sequences.fasta -a mjn -e 0 -o network.graphml
```

**Use when**: You want to infer ancestral haplotypes and show complex evolutionary relationships. The epsilon parameter controls network complexity (0 = maximum simplification).

## Distance Metrics

- **hamming**: Simple count of differences (fastest)
- **jc**: Jukes-Cantor correction for multiple substitutions
- **k2p**: Kimura 2-parameter (transitions vs transversions)
- **tamura_nei**: Accounts for GC content and transition/transversion bias

```bash
pypopart network sequences.fasta -d k2p -o network.graphml
```

## File Format Support

### Input Formats

- **FASTA** (`.fasta`, `.fa`, `.fna`)
- **NEXUS** (`.nexus`, `.nex`) - including traits/metadata
- **PHYLIP** (`.phy`, `.phylip`)
- **GenBank** (`.gb`, `.gbk`)

### Output Formats

- **GraphML** (`.graphml`) - Recommended, preserves all attributes
- **GML** (`.gml`)
- **JSON** (`.json`)
- **NEXUS** (`.nexus`, `.nex`)

### Visualization Formats

- **PNG** (`.png`) - Raster image
- **PDF** (`.pdf`) - Vector, publication-ready
- **SVG** (`.svg`) - Vector, web-friendly
- **HTML** (`.html`) - Interactive

## Advanced Usage

### Working with Metadata

```bash
# Load alignment with population metadata
pypopart load sequences.fasta -m metadata.csv

# Visualize colored by population
pypopart visualize network.graphml -o network.png --color-by population
```

Metadata CSV format:

```csv
id,population,location
Hap1,PopA,Site1
Hap2,PopA,Site1
Hap3,PopB,Site2
```

### Network Analysis

```bash
# Comprehensive statistics
pypopart analyze network.graphml --stats --topology --popgen -o results.json
```

### Topology Analysis

```bash
pypopart analyze network.graphml --topology
```

Identifies:

- Connected components
- Star-like patterns
- Central/hub nodes
- Potential ancestral nodes

### Custom Visualization

```bash
# Circular layout with labels
pypopart visualize network.graphml -o network.pdf \
    --layout circular --show-labels --width 1200 --height 1200

# Radial layout, interactive
pypopart visualize network.graphml -o network.html \
    --layout radial --interactive
```

### Geographic Visualization

PyPopART supports overlaying haplotype networks on geographic maps using latitude/longitude coordinates.

```bash
# Create network from sequences
pypopart network sequences.fasta -o network.graphml

# Create static geographic visualization
pypopart geo-visualize network.graphml \
    -m metadata.csv \
    -o geo_network.png \
    --projection mercator \
    --show-labels \
    --show-borders

# Create interactive geographic map
pypopart geo-visualize network.graphml \
    -m metadata.csv \
    -o geo_network.html \
    --interactive \
    --base-map OpenStreetMap \
    --zoom 4
```

Geographic metadata CSV format:

```csv
id,population,location,latitude,longitude
Hap1,PopA,New York,40.7128,-74.0060
Hap2,PopB,London,51.5074,-0.1278
Hap3,PopC,Tokyo,35.6762,139.6503
```

**Supported projections:**

- `mercator` - Web Mercator (preserves angles)
- `platecarree` - Equirectangular (simple lat/lon)
- `orthographic` - 3D globe view

**Interactive map base layers:**

- `OpenStreetMap` - Standard map tiles
- `Stamen Terrain` - Terrain with hill shading
- `CartoDB positron` - Clean, minimal style

**Python API:**

```python
from pypopart.visualization import GeoVisualizer, InteractiveGeoVisualizer

# Static map
viz = GeoVisualizer(network)
fig, ax = viz.plot(
    coordinates=coordinates,
    projection='mercator',
    show_labels=True,
    show_borders=True,
    output_file='geo_network.png'
)

# Interactive map
viz_interactive = InteractiveGeoVisualizer(network)
map_obj = viz_interactive.plot(
    coordinates=coordinates,
    base_map='OpenStreetMap',
    zoom_start=2,
    output_file='geo_network.html'
)
```

## Examples and Tutorials

Example data and Jupyter notebooks can be found in the `examples/` directory:

- `01_basic_workflow.ipynb` - Complete workflow from sequences to network
- `02_algorithm_comparison.ipynb` - Comparing different network algorithms
- `03_visualization_options.ipynb` - Customizing network plots
- `geo_example.py` - Geographic visualization with real-world coordinates
- `geo_data/` - Sample data with geographic metadata
- `04_population_genetics.ipynb` - Population genetics analysis
- `05_real_world_example.ipynb` - Case study with real data

## Documentation

Full documentation is available at [https://pypopart.readthedocs.io](https://pypopart.readthedocs.io) (coming soon)

Topics covered:

- Installation and setup
- Detailed API reference
- Algorithm descriptions and parameters
- Visualization customization
- Population genetics measures
- File format specifications
- Troubleshooting guide
- Contributing guidelines

## Citation

If you use PyPopART in your research, please cite:

```text
Taranto, A. (2024). PyPopART: Pure Python implementation of haplotype network analysis.
GitHub repository: https://github.com/adamtaranto/pypopart
```

## License

PyPopART is licensed under the GNU General Public License v3.0 or later. See [LICENSE](LICENSE) for details.

## Contributing

Contributions are welcome! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

### Development Setup

```bash
git clone https://github.com/adamtaranto/pypopart.git
cd pypopart
pip install -e ".[dev]"
pre-commit install
```

### Running Tests

```bash
pytest
```

### Code Quality

```bash
# Linting
ruff check src/

# Formatting
ruff format src/

# Type checking
mypy src/
```

## Related Projects

- [PopART](http://popart.otago.ac.nz/) - Original PopART software (Leigh & Bryant, 2015)
- [pegas](https://cran.r-project.org/package=pegas) - R package for population genetics
- [popart-networks](https://github.com/jessicawleigh/popart-current) - PopART source code

## Acknowledgments

PyPopART is inspired by the original PopART software developed by Jessica Leigh and David Bryant.

## Contact

- Author: Adam Taranto
- GitHub: [@adamtaranto](https://github.com/adamtaranto)
- Issues: [GitHub Issues](https://github.com/adamtaranto/pypopart/issues)

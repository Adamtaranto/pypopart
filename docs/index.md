# PyPopART Documentation

Welcome to **PyPopART** - a pure Python implementation of PopART (Population Analysis with Reticulate Trees) for constructing and visualizing haplotype networks from DNA sequence data.

## Overview

PyPopART provides a complete toolkit for:

- **Network Construction**: Six algorithms (MST, MSN, TCS, MJN, PN, TSW)
- **Distance Calculation**: Multiple evolutionary models
- **Network Analysis**: Comprehensive statistics and topology analysis
- **Visualization**: Static, interactive, and web-based network plots
- **Population Genetics**: Diversity measures, FST, Tajima's D, and more
- **Dual Interface**: Command-line tools and web-based GUI

## Key Features

### 🌳 Multiple Algorithms

- **MST** (Minimum Spanning Tree): Simplest tree-based network
- **MSN** (Minimum Spanning Network): Shows alternative connections
- **TCS** (Statistical Parsimony): Statistically justified connections
- **MJN** (Median-Joining): Infers ancestral haplotypes
- **PN** (Parsimony Network): Consensus from multiple trees
- **TSW** (Tight Span Walker): Metric-preserving network construction

### 📏 Distance Metrics

- Hamming distance
- Jukes-Cantor correction
- Kimura 2-parameter
- Tamura-Nei model

### 📊 Comprehensive Analysis

- Network statistics (diameter, clustering, centrality)
- Topology analysis (star patterns, hubs, bridges)
- Population genetics measures (Tajima's D, Fu's Fs, FST)
- Diversity metrics (nucleotide, haplotype, Shannon)

### 🎨 Rich Visualization

- Static plots with matplotlib (PNG, PDF, SVG)
- Interactive plots with Plotly (HTML)
- Multiple layout algorithms
- Customizable colors, sizes, and labels

### 📁 Flexible I/O

- Input: FASTA, NEXUS, PHYLIP, GenBank
- Output: GraphML, GML, JSON, NEXUS
- Metadata support (populations, traits, locations)

## Entry Points

PyPopART offers two interfaces for different workflows:

### Command-Line Interface (CLI)

For automation, scripting, and batch processing:

```bash
# Get help
pypopart --help

# List available algorithms
pypopart info --list-algorithms

# Construct a median-joining network
pypopart network sequences.fasta -o network.graphml

# Visualize the network
pypopart visualize network.graphml -o network.png
```

### Web-based GUI

For interactive analysis and exploration:

```bash
# Start the GUI application
pypopart-gui

# Opens web interface at http://localhost:8050
# or specify custom port:
pypopart-gui --port 8080
```

The GUI provides:
- Drag-and-drop file upload
- Interactive network visualization (zoom, pan, drag nodes)
- Real-time algorithm parameter adjustment
- Population-based coloring with pie charts
- Multiple layout algorithms
- Network statistics and haplotype summaries
- Export to various formats

## Quick Examples

### Command Line Workflow

```bash
# Load and validate sequences
pypopart load sequences.fasta

# Construct network with TCS algorithm
pypopart network sequences.fasta -a tcs -o network.graphml

# Analyze network statistics
pypopart analyze network.graphml --stats

# Create visualization
pypopart visualize network.graphml -o network.png
```

### Python API

```python
from pypopart.io import load_alignment
from pypopart.algorithms import MJNAlgorithm
from pypopart.core.distance import DistanceCalculator
from pypopart.core.condensation import condense_alignment
from pypopart.visualization import StaticVisualizer

# Load data
alignment = load_alignment('sequences.fasta')

# Calculate distances
calc = DistanceCalculator(method='k2p')
distances = calc.calculate_matrix(alignment)

# Build network
haplotypes, _ = condense_alignment(alignment)
mjn = MJNAlgorithm(epsilon=0)
network = mjn.construct_network(haplotypes, distances)

# Visualize
viz = StaticVisualizer(network)
viz.plot(output_file='network.png')
```

## Getting Started

- **[Installation](installation.md)**: Install PyPopART
- **[Quick Start](quickstart.md)**: Your first haplotype network
- **[Basic Concepts](concepts.md)**: Understanding haplotype networks

## User Guide

- **[CLI Guide](guide/cli.md)**: Command-line interface
- **[Python API](guide/api.md)**: Programmatic usage
- **[Algorithms](guide/algorithms.md)**: Choosing the right algorithm
- **[Visualization](guide/visualization.md)**: Creating beautiful plots

## Tutorials

- **[Basic Workflow](tutorials/basic_workflow.md)**: Complete example
- **[Algorithm Comparison](tutorials/algorithm_comparison.md)**: Compare all algorithms
- **[Visualization Options](tutorials/visualization.md)**: Customize your plots
- **[Population Genetics](tutorials/popgen.md)**: Genetic diversity analysis

## API Reference

Detailed documentation for all classes and functions:

- **[Core Classes](api/core/sequence.md)**: Sequence, Alignment, Haplotype
- **[Algorithms](api/algorithms/mst.md)**: Network construction
- **[Visualization](api/visualization/static.md)**: Plotting functions
- **[Statistics](api/stats/statistics.md)**: Analysis tools

## About

PyPopART is developed and maintained by Adam Taranto. It is inspired by the original [PopART software](http://popart.otago.ac.nz/) by Jessica Leigh and David Bryant.

**License**: GNU General Public License v3.0 or later

**Citation**: If you use PyPopART in your research, please cite:

```text
Taranto, A. (2024). PyPopART: Pure Python implementation of haplotype network analysis.
GitHub repository: https://github.com/adamtaranto/pypopart
```

## Need Help?

- **[FAQ](faq.md)**: Frequently asked questions
- **[GitHub Issues](https://github.com/adamtaranto/pypopart/issues)**: Report bugs or request features
- **[Contributing](contributing.md)**: Contribute to PyPopART

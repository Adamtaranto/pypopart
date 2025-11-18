# Command-Line Interface

PyPopART provides a comprehensive command-line interface for scripting and automation.

## Basic Usage

```bash
pypopart [OPTIONS] COMMAND [ARGS]...
```

## Main Commands

### `network` - Build Haplotype Networks

Create a haplotype network from sequence data:

```bash
pypopart network input.fasta -a MST -o output
```

**Options:**

- `-i, --input PATH`: Input sequence file (required)
- `-a, --algorithm ALGORITHM`: Network algorithm (MST, MSN, TCS, MJN, PN, TSW)
- `-d, --distance METRIC`: Distance metric (hamming, jukes-cantor, k2p, tamura-nei)
- `-o, --output PATH`: Output prefix for files
- `-f, --format FORMAT`: Output format (nexus, gml, graphml, json)
- `--epsilon FLOAT`: TCS connection limit (default: 0.95)
- `--plot`: Generate visualization plot
- `--interactive`: Create interactive HTML plot

**Examples:**

```bash
# Basic MST network
pypopart network sequences.fasta -a MST -o my_network

# TCS network with custom epsilon
pypopart network sequences.fasta -a TCS --epsilon 0.99 -o tcs_network

# MJN network with interactive plot
pypopart network sequences.fasta -a MJN --plot --interactive -o mjn_network
```

### `distance` - Calculate Distance Matrix

Compute pairwise genetic distances:

```bash
pypopart distance input.fasta -m k2p -o distances.csv
```

**Options:**

- `-i, --input PATH`: Input sequence file (required)
- `-m, --metric METRIC`: Distance metric (hamming, jukes-cantor, k2p, tamura-nei)
- `-o, --output PATH`: Output file path
- `-f, --format FORMAT`: Output format (csv, tsv, json)

### `stats` - Network Statistics

Calculate network statistics:

```bash
pypopart stats network.gml -o statistics.txt
```

**Options:**

- `-i, --input PATH`: Input network file (required)
- `-o, --output PATH`: Output file path
- `--topology`: Include topology metrics
- `--diversity`: Include diversity indices
- `--popgen`: Include population genetics statistics

### `plot` - Visualize Networks

Create network visualizations:

```bash
pypopart plot network.gml -o figure.png
```

**Options:**

- `-i, --input PATH`: Input network file (required)
- `-o, --output PATH`: Output figure path
- `-l, --layout LAYOUT`: Layout algorithm (spring, circular, kamada-kawai)
- `--width INTEGER`: Figure width in pixels
- `--height INTEGER`: Figure height in pixels
- `--node-size FLOAT`: Node size multiplier
- `--edge-width FLOAT`: Edge width
- `--color-by ATTRIBUTE`: Color nodes by metadata attribute
- `--interactive`: Create interactive HTML plot

## Global Options

- `--version`: Show version and exit
- `--help`: Show help message
- `-v, --verbose`: Enable verbose output
- `--quiet`: Suppress non-error output

## Input Formats

PyPopART supports multiple sequence formats:

- **FASTA** (`.fasta`, `.fa`, `.fna`)
- **NEXUS** (`.nex`, `.nexus`)
- **PHYLIP** (`.phy`, `.phylip`)
- **GenBank** (`.gb`, `.genbank`)

## Output Formats

Networks can be exported in various formats:

- **NEXUS** - Compatible with PopART and other tools
- **GML** - Graph Modeling Language
- **GraphML** - XML-based graph format
- **JSON** - JavaScript Object Notation

## Working with Metadata

Include metadata in NEXUS format:

```nexus
#NEXUS
BEGIN TAXA;
    DIMENSIONS NTAX=4;
    TAXLABELS Seq1 Seq2 Seq3 Seq4;
END;

BEGIN CHARACTERS;
    DIMENSIONS NCHAR=20;
    FORMAT DATATYPE=DNA;
    MATRIX
        Seq1 ATCGATCGATCGATCGATCG
        Seq2 ATCGATCGATCGATCGATCG
        Seq3 ATCGATCGATCGATTGATCG
        Seq4 ATCGATCGATCGATTGATCG
    ;
END;

BEGIN TRAITS;
    DIMENSIONS NTRAITS=2;
    FORMAT LABELS=YES SEPARATOR=,;
    TRAITLABELS Population Location;
    MATRIX
        Seq1 PopA Site1
        Seq2 PopA Site1
        Seq3 PopB Site2
        Seq4 PopB Site2
    ;
END;
```

Use metadata for coloring:

```bash
pypopart plot network.gml --color-by Population -o colored_network.png
```

## Scripting Examples

### Batch Processing

Process multiple files:

```bash
for file in *.fasta; do
    base=$(basename "$file" .fasta)
    pypopart network "$file" -a MST -o "networks/${base}"
done
```

### Pipeline Example

Complete analysis pipeline:

```bash
# 1. Build network
pypopart network input.fasta -a MJN -o mjn_network -f gml

# 2. Calculate statistics
pypopart stats mjn_network.gml --topology --popgen -o statistics.txt

# 3. Create visualizations
pypopart plot mjn_network.gml -o static_plot.png
pypopart plot mjn_network.gml --interactive -o interactive_plot.html
```

### Comparing Algorithms

Compare different algorithms on the same data:

```bash
algorithms=(MST MSN TCS MJN)
for alg in "${algorithms[@]}"; do
    pypopart network sequences.fasta -a "$alg" -o "comparison/${alg}"
done
```

## Tips and Best Practices

1. **Choose the right algorithm**: Start with MST for exploration, use TCS for within-species data, MJN for complex scenarios

2. **Distance metrics**: Use Hamming for very similar sequences, K2P or Tamura-Nei for more divergent data

3. **Visualization**: Generate both static (for publications) and interactive (for exploration) plots

4. **Large datasets**: Consider using MST or MSN first before trying computationally intensive algorithms like MJN

5. **Metadata**: Include population/trait information for richer analyses and visualizations

## Troubleshooting

### Common Issues

**"No such file or directory"**
- Check file paths are correct
- Use absolute paths if needed

**"Invalid sequence format"**
- Verify file format is supported
- Check for file corruption
- Ensure sequences are aligned

**"Algorithm failed to converge"**
- Try different algorithm
- Check data quality
- Adjust algorithm parameters

## Next Steps

- [Python API Guide](api.md): Use PyPopART programmatically
- [Tutorials](../tutorials/basic_workflow.md): Detailed walkthroughs
- [Algorithm Guide](algorithms.md): Choose the right algorithm

# PyPopART Network Construction Algorithm Usage

This document provides examples of how to use the different network construction algorithms in PyPopART.

## Overview

PyPopART provides four main network construction algorithms:

1. **MST (Minimum Spanning Tree)** - Creates a simple tree connecting all haplotypes with minimum total distance
2. **MSN (Minimum Spanning Network)** - Extends MST by adding alternative connections at equal distances
3. **TCS (Statistical Parsimony)** - Connects haplotypes based on 95% parsimony criterion
4. **MJN (Median-Joining Network)** - Infers median vectors (ancestral sequences) to simplify networks

## Basic Usage

```python
from pypopart.core.sequence import Sequence
from pypopart.core.alignment import Alignment
from pypopart.algorithms import (
    MinimumSpanningTree,
    MinimumSpanningNetwork,
    TCS,
    MedianJoiningNetwork
)

# Create alignment
sequences = [
    Sequence("sample1", "ATCGATCG"),
    Sequence("sample2", "ATCGATCC"),
    Sequence("sample3", "ATCGGTCG"),
    Sequence("sample4", "GTCGATCG")
]
alignment = Alignment(sequences)

# 1. Minimum Spanning Tree
mst = MinimumSpanningTree(distance_method="hamming", algorithm="prim")
network_mst = mst.construct_network(alignment)
print(f"MST: {len(network_mst.haplotypes)} haplotypes, {len(network_mst.edges)} edges")

# 2. Minimum Spanning Network
msn = MinimumSpanningNetwork(distance_method="hamming", epsilon=0.0)
network_msn = msn.construct_network(alignment)
print(f"MSN: {len(network_msn.haplotypes)} haplotypes, {len(network_msn.edges)} edges")

# 3. TCS (Statistical Parsimony)
tcs = TCS(distance_method="hamming", confidence=0.95)
network_tcs = tcs.construct_network(alignment)
print(f"TCS: {len(network_tcs.haplotypes)} haplotypes, {len(network_tcs.edges)} edges")

# 4. Median-Joining Network
mjn = MedianJoiningNetwork(distance_method="hamming", simplify=True)
network_mjn = mjn.construct_network(alignment)
print(f"MJN: {len(network_mjn.haplotypes)} haplotypes, {len(network_mjn.edges)} edges")
```

## Algorithm-Specific Parameters

### Minimum Spanning Tree (MST)

```python
mst = MinimumSpanningTree(
    distance_method="hamming",  # Distance metric: hamming, p, jc, k2p, tn
    algorithm="prim"            # MST algorithm: "prim" or "kruskal"
)
```

**When to use:**
- Need simplest representation
- Want guaranteed tree structure (no cycles)
- Have well-separated haplotypes

### Minimum Spanning Network (MSN)

```python
msn = MinimumSpanningNetwork(
    distance_method="hamming",
    epsilon=0.0,              # Tolerance for equal distances
    max_connections=None      # Limit connections per node (None = no limit)
)
```

**When to use:**
- Want to show alternative relationships
- Sequences have equal distances to multiple neighbors
- Need more realistic network than simple tree

### TCS (Statistical Parsimony)

```python
tcs = TCS(
    distance_method="hamming",
    confidence=0.95,          # Parsimony confidence level (0.0-1.0)
    connection_limit=None     # Max distance (auto-calculated if None)
)
```

**When to use:**
- Intraspecific/population-level data
- Want statistically justified connections
- Need to identify separate lineages
- Disconnected networks are acceptable

### Median-Joining Network (MJN)

```python
mjn = MedianJoiningNetwork(
    distance_method="hamming",
    epsilon=0.0,              # Complexity control parameter
    max_median_vectors=None,  # Limit number of inferred nodes
    simplify=True             # Remove unnecessary median vectors
)
```

**When to use:**
- Have missing intermediate haplotypes
- Complex reticulation patterns
- Want to infer ancestral sequences
- Need most parsimonious representation

## Distance Methods

All algorithms support multiple distance calculation methods:

```python
# Hamming distance (count of differences)
algorithm = Algorithm(distance_method="hamming")

# p-distance (proportion of differences)
algorithm = Algorithm(distance_method="p")

# Jukes-Cantor correction
algorithm = Algorithm(distance_method="jc")

# Kimura 2-parameter
algorithm = Algorithm(distance_method="k2p")

# Tamura-Nei
algorithm = Algorithm(distance_method="tn")
```

**Recommendations:**
- **Hamming**: Best for closely related sequences, discrete mutation counts
- **p-distance**: Simple proportion, good for quick analysis
- **Jukes-Cantor**: Corrects for multiple substitutions, assumes equal rates
- **K2P**: Distinguishes transitions vs transversions
- **Tamura-Nei**: Most sophisticated, accounts for base composition and different rates

## Network Analysis

After constructing a network, you can analyze it:

```python
# Check connectivity
print(f"Is connected: {network.is_connected()}")

# Get network statistics
stats = network.calculate_stats()
print(f"Nodes: {stats['num_nodes']}")
print(f"Edges: {stats['num_edges']}")
print(f"Connected: {stats['is_connected']}")

# Get haplotype information
for hap in network.haplotypes:
    print(f"{hap.id}: frequency={hap.frequency}, samples={len(hap.sample_ids)}")
    
# Get edge information
for u, v in network.edges:
    dist = network.get_edge_distance(u, v)
    print(f"{u} -- {v}: distance={dist}")
```

## Comparing Algorithms

```python
algorithms = {
    'MST': MinimumSpanningTree(),
    'MSN': MinimumSpanningNetwork(),
    'TCS': TCS(),
    'MJN': MedianJoiningNetwork()
}

for name, algo in algorithms.items():
    network = algo.construct_network(alignment)
    print(f"\n{name}:")
    print(f"  Haplotypes: {len(network.haplotypes)}")
    print(f"  Edges: {len(network.edges)}")
    print(f"  Connected: {network.is_connected()}")
```

## Best Practices

1. **Start with MST** - Simplest algorithm, good baseline
2. **Try TCS** - Good for population data with statistical justification
3. **Use MSN** - When you need to show alternative relationships
4. **Use MJN** - When you suspect missing intermediates

5. **Choose distance method carefully**:
   - Hamming for closely related sequences
   - Evolutionary models (JC, K2P, TN) for more divergent sequences

6. **Validate results** - Compare networks from different algorithms
7. **Consider biology** - Choose algorithm that matches your evolutionary model

## Performance Tips

- Use Prim's algorithm for MST with dense graphs
- Use Kruskal's algorithm for MST with sparse graphs
- Limit max_median_vectors in MJN for large datasets
- Pre-calculate distance matrix if running multiple algorithms

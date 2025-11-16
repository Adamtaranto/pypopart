# Tight Span Walker (TSW) - Parsimony Network Algorithm

## Overview

The Tight Span Walker (TSW) algorithm constructs haplotype networks using parsimony principles by computing the **tight span** of a distance matrix. The tight span is the smallest metric space that contains all optimal paths between sequences, making it ideal for representing complex evolutionary relationships with reticulation events.

## Algorithm Details

### How It Works

1. **Compute dT Distances (Tree Metric)**
   - For each pair of sequences (i, j), calculate: `dT(i,j) = max over all k of |d(i,k) - d(j,k)|`
   - This represents the minimum distance if sequences were constrained to a tree structure

2. **Build Geodesic Paths**
   - For each pair of haplotypes, construct the geodesic (shortest) path
   - Compare the actual distance with the dT distance
   - If `actual_distance - dT_distance > epsilon`, infer intermediate median vertices

3. **Infer Median Vertices**
   - When needed, create intermediate (ancestral) nodes
   - Medians represent hypothetical ancestral or intermediate sequences

4. **Construct Network**
   - Connect all haplotypes through geodesic paths
   - Include inferred median vertices to maintain metric properties

### Key Features

- **Metric Preservation**: Maintains all distance relationships from original data
- **Reticulate Networks**: Can represent complex relationships
- **Ancestral Inference**: Automatically infers hypothetical ancestral sequences
- **Parsimony-Based**: Uses parsimony principles

## Parameters

- **epsilon** (float, default=1e-6): Tolerance for metric comparisons
- **distance_method** (str, default='hamming'): Distance calculation method

## When to Use TSW

✓ Complex evolutionary relationships with reticulation  
✓ Small to medium datasets (n < 100)  
✓ Accurate metric representation needed  
✓ Ancestral sequence inference desired

## Usage Example

```python
from pypopart.algorithms import TightSpanWalker
from pypopart.io import load_alignment

alignment = load_alignment('sequences.fasta')
tsw = TightSpanWalker(distance_method='hamming')
network = tsw.construct_network(alignment)
```

## References

1. Dress, A. W., & Huson, D. H. (2004). Constructing splits graphs. IEEE/ACM Transactions on Computational Biology and Bioinformatics, 1(3), 109-115.
2. Bryant, D., & Moulton, V. (2004). Neighbor-Net. Molecular Biology and Evolution, 21(2), 255-265.

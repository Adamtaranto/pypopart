# Algorithm Refinements: TSW and PN

## Summary

This document describes the refinements made to the Tight Span Walker (TSW) and Parsimony Network (PN) algorithms to reduce the excess creation of intermediate nodes and improve network accuracy.

## Problem Statement

### TSW Issues
- The original Python implementation was creating an excess of intermediate (median) nodes
- Simplified geodesic path computation didn't match the C++ reference implementation
- Missing pruning logic for unnecessary median vertices

### PN Issues
- The algorithm was automatically subdividing all multi-mutation edges into single-mutation steps
- This created many unnecessary intermediate nodes where edges should be direct
- Deviated from the C++ reference implementation which relies on consensus sampling

## Solutions Implemented

### TSW Refinements

#### 1. Proper Geodesic Computation
The geodesic path algorithm now follows the C++ reference implementation:

- **Bipartite Coloring**: Vertices are classified as "green" (closer to source) or "red" (closer to target) based on their dT distances
- **Delta Calculation**: Computes the splitting parameter delta to determine if an intermediate vertex is needed
- **Vertex Map**: Tracks dT vectors to avoid creating duplicate vertices with the same metric properties
- **Recursive Construction**: Properly creates intermediate vertices only when delta < dT(f,g)

**Key Code Changes:**
```python
# Old approach - overly simplistic
if actual_dist - dt_dist > self.epsilon:
    median_hap = self._create_median_vertex(hap1, hap2)
    # Always created medians when distances didn't match

# New approach - proper geodesic computation
# 1. Compute bipartite coloring
for i in range(n):
    if abs(fi + dt_fg + gi - distance_matrix.matrix[i, idx2]) < self.epsilon:
        green_vertices.add(i)
    elif abs(fi + dt_fg - gi) < self.epsilon:
        red_vertices.add(i)

# 2. Compute delta (splitting parameter)
delta = min([(fi + fj - dij) for i, j in green_pairs]) / 2.0

# 3. Only create vertex if needed
if dt_fg > delta + self.epsilon:
    # Create intermediate vertex and recurse
```

#### 2. Median Vertex Pruning
Added post-processing to remove medians that don't serve a topological purpose:

**Pruning Criteria:**
A median vertex is kept only if it:
1. Connects two or more observed (non-median) vertices (bridges observed nodes), OR
2. Connects to at least one observed vertex and one or more medians (on path between observed), OR  
3. Connects to two or more medians (internal to a path)

**Benefits:**
- Eliminates isolated or non-bridging intermediate nodes
- Produces cleaner, more interpretable networks
- Maintains connectivity while removing redundancy

### PN Refinements

#### 1. Removed Automatic Edge Subdivision
The key change was to stop automatically subdividing multi-mutation edges:

**Old Behavior:**
```python
# Add edge if not already present
if not network.has_edge(id1, id2):
    network.add_edge(id1, id2, distance=distance)
    
    # Automatically subdivide edges > 1 mutation
    if distance > 1:
        self._add_median_vertices_along_edge(
            network, id1, id2, seq1, seq2, int(distance)
        )
```

**New Behavior:**
```python
# Add edge if not already present
# Let the consensus sampling handle network structure
if not network.has_edge(id1, id2):
    network.add_edge(id1, id2, distance=distance)
```

**Rationale:**
- The PN algorithm samples edges from multiple random parsimony trees
- The consensus process (frequency thresholding) naturally determines which edges appear in the final network
- Artificially subdividing edges defeats the purpose of the consensus approach
- The C++ implementation doesn't subdivide edges - it relies on the sampling to create the network structure

#### 2. Removed Median Vertex Creation Method
The `_add_median_vertices_along_edge` method was completely removed as it's no longer needed.

## Validation

### Test Results
All existing unit tests pass:
- TSW: 17/17 tests passing
- PN: 21/21 tests passing

### Example Improvements

#### TSW Example
**Input:** 3 sequences with 4 mutations between seq1 and seq2
- seq1: AAAA
- seq2: TTTT (4 diffs)
- seq3: AATT (2 diffs from each)

**Before:** Potentially created unnecessary medians
**After:** 3 nodes (3 observed, 0 median), 3 edges

#### PN Example
**Input:** 2 sequences with 4 mutations
- seq1: AAAA
- seq2: TTTT

**Before:** Would create 3 intermediate medians to break the edge into 4 single-mutation steps (5 total nodes)
**After:** Direct edge with distance=4.0 (2 nodes, 1 edge)

## Impact

### TSW
✓ Cleaner networks with only necessary medians
✓ Better matches C++ reference implementation behavior
✓ Improved computational efficiency (fewer vertices to process)
✓ More interpretable network structures

### PN
✓ Networks match the intended consensus sampling design
✓ Edges accurately reflect sampled tree topologies
✓ Reduced node count without loss of information
✓ Aligns with C++ reference implementation

## Technical Details

### TSW Implementation Notes
- Extended dT matrix dynamically as new medians are created
- Vertex map uses tuple of dT distances as key for uniqueness
- Bipartite coloring based on path decomposition from Dress & Huson (2004)
- Delta calculation follows the formula: δ = min{(dT(f,i) + dT(f,j) - d(i,j))/2}

### PN Implementation Notes
- Maintains tree sampling and frequency thresholding
- Edge distances calculated directly from sequence differences
- Random seed support ensures reproducibility
- Consensus threshold (min_edge_frequency) controls network complexity

## References

1. Dress, A. W., & Huson, D. H. (2004). Constructing splits graphs. IEEE/ACM Transactions on Computational Biology and Bioinformatics, 1(3), 109-115.

2. Leigh, J. W., & Bryant, D. (2015). PopART: Full-feature software for haplotype network construction. Methods in Ecology and Evolution, 6(9), 1110-1116.

## Future Considerations

### Potential Enhancements
1. **TSW**: Could optimize dT matrix updates to avoid full matrix expansion
2. **TSW**: Could add configurable pruning strategies
3. **PN**: Could expose more tree generation parameters
4. **Both**: Could add visualization highlighting median vs observed nodes

### Performance Notes
- TSW: O(n³) complexity remains, but with fewer vertices created
- PN: Linear in number of edges sampled, no additional overhead from median creation
- Both algorithms suitable for datasets with < 200 sequences

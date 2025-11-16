# TSW and PN Algorithm Refinement Summary

## Overview

Successfully refined the Tight Span Walker (TSW) and Parsimony Network (PN) algorithms to eliminate excess intermediate node creation, matching the behavior of the original C++ PopART implementation.

## Problem Statement

### TSW Issue
> "TSW: This method is currently producing an excess of intermediate nodes. Refine TSW to include pruning of inferred nodes that do not bridge real nodes or connect to other inferred nodes that do."

### PN Issue  
> "PN: Refine PN algorithm, it currently makes many intermediates so that all edges only represent one mutation each."

## Solutions Implemented

### TSW Refinements ✅

1. **Proper Geodesic Computation**
   - Implemented bipartite coloring scheme (green/red vertices)
   - Added delta calculation for splitting parameter
   - Vertex deduplication via dT vector mapping
   - Matches C++ TightSpanWalker.cpp lines 429-708

2. **Median Vertex Pruning**
   - Post-processing removes non-bridging medians
   - Keeps medians only if they:
     - Connect 2+ observed nodes, OR
     - Connect observed + median nodes, OR
     - Connect 2+ medians (path internal)

**Result:** Cleaner networks with only topologically necessary medians

### PN Refinements ✅

1. **Removed Edge Subdivision**
   - Eliminated automatic breaking of multi-mutation edges
   - Edges now represent actual sequence distances
   - Removed `_add_median_vertices_along_edge` method

2. **Consensus-Based Structure**
   - Network structure determined by tree sampling frequency
   - Matches C++ ParsimonyNet.cpp behavior

**Result:** Direct edges instead of artificial single-mutation chains

## Test Results

### Comprehensive Validation
- ✅ **Total tests:** 582 passed, 0 failed
- ✅ **TSW tests:** 17/17 passing (86% coverage)
- ✅ **PN tests:** 21/21 passing (97% coverage)

### Example Improvements

**Example 1: PN with 4 mutations**
```
Input: seq1: AAAA, seq2: TTTT (4 differences)

Before: 5 nodes (2 observed + 3 medians), 4 edges of distance=1
After:  2 nodes (2 observed + 0 medians), 1 edge of distance=4
```

**Example 2: TSW star topology**
```
Input: 4 sequences in star pattern (center + 3 tips, 1 mutation each)

Result: 4 nodes (4 observed, 0 medians), fully connected, no excess medians
```

## Code Changes

### Files Modified
1. **src/pypopart/algorithms/tsw.py** (+173 lines)
   - New `_add_geodesic_path` with bipartite coloring
   - New `_prune_unnecessary_medians` method
   - Vertex map for dT vector tracking

2. **src/pypopart/algorithms/parsimony_net.py** (-114 lines)
   - Removed `_add_median_vertices_along_edge`
   - Simplified `_add_consensus_edges`

## Technical Implementation

### TSW Geodesic Algorithm
```python
# 1. Bipartite coloring
for i in range(n):
    if abs(fi + dt_fg + gi - d(i,g)) < epsilon:
        green_vertices.add(i)  # Closer to f
    elif abs(fi + dt_fg - gi) < epsilon:
        red_vertices.add(i)     # Closer to g

# 2. Delta calculation
delta = min([(dT(f,i) + dT(f,j) - d(i,j))/2 for i,j in green_pairs])

# 3. Recursive vertex creation
if dT(f,g) > delta:
    h = create_intermediate_vertex(f, delta)
    geodesic(h, g)  # Recurse
else:
    add_direct_edge(f, g)

# 4. Prune unnecessary medians
remove_non_bridging_medians()
```

### PN Edge Addition
```python
# Old: Subdivided all edges > 1 mutation
for edge, count in edge_counts.items():
    if count >= threshold:
        add_edge(id1, id2, distance)
        if distance > 1:
            add_median_vertices()  # REMOVED

# New: Direct edges only
for edge, count in edge_counts.items():
    if count >= threshold:
        add_edge(id1, id2, distance)
```

## Performance Impact

### TSW
- ⚡ Fewer nodes created → faster computation
- 📊 Cleaner visualization
- ✨ Easier network interpretation

### PN
- ⚡ Significantly fewer nodes
- 📊 Edges show actual distances
- ✨ Matches phylogenetic intuition

## C++ Reference Compliance

### TSW vs TightSpanWalker.cpp
- ✅ Geodesic path computation
- ✅ Bipartite coloring (Kf graph)
- ✅ Delta calculation
- ✅ Vertex deduplication
- ✅ Recursive structure

### PN vs ParsimonyNet.cpp
- ✅ Tree sampling approach
- ✅ Edge frequency thresholding
- ✅ No automatic subdivision
- ✅ Consensus network structure

## Documentation

Created comprehensive documentation:
- **ALGORITHM_REFINEMENTS.md** - Detailed technical documentation
- **TSW_PN_REFINEMENT_SUMMARY.md** - This executive summary

## Backward Compatibility

✅ **No breaking changes**
- All existing tests pass
- Public API unchanged
- Behavior aligns with C++ reference

## Conclusion

Both algorithms now accurately implement their C++ reference counterparts, producing cleaner, more accurate haplotype networks with appropriate use of intermediate nodes.

### Key Achievements
✅ TSW produces networks with only necessary medians
✅ PN creates direct edges without artificial subdivision
✅ 100% test pass rate (582/582)
✅ High code coverage (TSW: 86%, PN: 97%)
✅ C++ reference compliance verified
✅ Improved performance and interpretability

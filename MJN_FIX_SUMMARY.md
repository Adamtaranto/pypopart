# MJN Algorithm Bug Fix Summary

## Problem
The Median-Joining Network (MJN) algorithm failed with the following error when epsilon > 0:
```
KeyError: 'Edge (Median_52, H4) not found'
```

This error occurred during the network simplification phase when the algorithm attempted to check if a direct edge exists between two neighbors of a degree-2 median vector.

## Root Cause
The `_simplify_network()` method in `src/pypopart/algorithms/mjn.py` had three issues:

1. **Missing edge existence check**: Line 667 called `get_edge_distance(n1, n2)` which raises a KeyError if the edge doesn't exist, instead of first checking with `has_edge()`

2. **Single-pass simplification**: The method only ran once, but should iterate until no more changes (matches C++ behavior)

3. **Incomplete obsolete removal**: Only removed degree-2 medians, but should also remove degree < 2 medians (degree 0 or 1)

## Solution
Fixed all three issues in the `_simplify_network()` method:

```python
# Before (broken):
direct_dist = network.get_edge_distance(n1, n2)  # Raises KeyError if edge doesn't exist
if direct_dist is None:
    network.add_edge(n1, n2, distance=d1 + d2)

# After (fixed):
if not network.has_edge(n1, n2):  # Check existence first
    network.add_edge(n1, n2, distance=d1 + d2)
```

Also added:
- Iterative loop: `while changed:` to match C++ `removeObsoleteVerts()` behavior
- Obsolete median removal: Remove all medians with degree < 2

## Testing
Created comprehensive test suite in `tests/unit/test_mjn_edge_not_found_fix.py`:

- ✅ Test with epsilon values 0, 1, 2, 3, 4, 5
- ✅ Verify all median nodes have degree >= 2 after simplification
- ✅ Verify obsolete medians (degree < 2) are removed
- ✅ Verify degree-2 medians are replaced with direct edges
- ✅ Verify observed haplotypes are never removed
- ✅ All 34 MJN tests pass (29 existing + 5 new)
- ✅ All 81 graph and MJN tests pass
- ✅ No security vulnerabilities (CodeQL passed)

## Comparison with C++ Implementation
Reviewed the original popart C++ implementation (`archive_popart_src/networks/MedJoinNet.cpp`) and confirmed our Python implementation now matches:

### C++ removeObsoleteVerts() behavior:
```cpp
bool MedJoinNet::removeObsoleteVerts() {
    while (changed) {
        for (unsigned i = _nsamples; i < vertexCount(); i++) {
            v = vertex(i);
            if (v->degree() < 2)  // Remove vertices with degree < 2
                obsoleteVerts.push_back(v);
        }
        // Remove obsolete vertices
        // Repeat until no changes
    }
}
```

Our Python implementation now matches this behavior exactly.

## Answer to Original Question
**Question**: "Should isolated median nodes (that don't bridge observed nodes) be preserved in the final network?"

**Answer**: **NO**. According to the C++ reference implementation:
- Median nodes with degree < 2 are considered "obsolete vertices"
- They are removed by `removeObsoleteVerts()` both during the algorithm and in post-processing
- Only median nodes that bridge between other nodes (degree >= 2) are preserved

This makes biological sense: median vectors are inferred ancestral/intermediate haplotypes that should connect observed haplotypes. If they don't connect at least 2 other nodes, they serve no purpose in the network.

## Files Changed
1. `src/pypopart/algorithms/mjn.py`: Fixed `_simplify_network()` method
2. `tests/unit/test_mjn_edge_not_found_fix.py`: Added 5 new comprehensive tests

## Verification
```bash
# All tests pass
python -m pytest tests/unit/test_mjn*.py -v
# 34 passed

# No lint issues
ruff check src/pypopart/algorithms/mjn.py tests/unit/test_mjn_edge_not_found_fix.py
# All checks passed!

# No security issues
# CodeQL: 0 alerts found
```

## Impact
- ✅ MJN algorithm now works correctly with any epsilon value
- ✅ Network simplification matches C++ reference implementation
- ✅ Final networks are properly simplified with only meaningful median vectors
- ✅ No breaking changes to existing functionality

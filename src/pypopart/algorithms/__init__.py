"""
Network construction algorithms for PyPopART.

This module provides various algorithms for constructing haplotype networks
from DNA sequence data.
"""

from .base import NetworkAlgorithm
from .mjn import MedianJoiningNetwork
from .msn import MinimumSpanningNetwork
from .mst import MinimumSpanningTree
from .tcs import TCS

# Convenient aliases
MSTAlgorithm = MinimumSpanningTree
MSNAlgorithm = MinimumSpanningNetwork
TCSAlgorithm = TCS
MJNAlgorithm = MedianJoiningNetwork

__all__ = [
    'NetworkAlgorithm',
    'MinimumSpanningTree',
    'MinimumSpanningNetwork',
    'TCS',
    'MedianJoiningNetwork',
    'MSTAlgorithm',
    'MSNAlgorithm',
    'TCSAlgorithm',
    'MJNAlgorithm',
]

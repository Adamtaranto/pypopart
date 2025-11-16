"""
Network construction algorithms for PyPopART.

This module provides various algorithms for constructing haplotype networks
from DNA sequence data.
"""

from .base import NetworkAlgorithm
from .mjn import MedianJoiningNetwork
from .msn import MinimumSpanningNetwork
from .mst import MinimumSpanningTree
from .parsimony_net import ParsimonyNetwork
from .tcs import TCS
from .tsw import TightSpanWalker

# Convenient aliases
MSTAlgorithm = MinimumSpanningTree
MSNAlgorithm = MinimumSpanningNetwork
TCSAlgorithm = TCS
MJNAlgorithm = MedianJoiningNetwork
PNAlgorithm = ParsimonyNetwork
TSWAlgorithm = TightSpanWalker

__all__ = [
    'NetworkAlgorithm',
    'MinimumSpanningTree',
    'MinimumSpanningNetwork',
    'TCS',
    'MedianJoiningNetwork',
    'ParsimonyNetwork',
    'TightSpanWalker',
    'MSTAlgorithm',
    'MSNAlgorithm',
    'TCSAlgorithm',
    'MJNAlgorithm',
    'PNAlgorithm',
    'TSWAlgorithm',
]

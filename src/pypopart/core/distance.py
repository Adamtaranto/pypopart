"""
Distance calculation and metrics for PyPopART.

Implements various distance metrics for DNA sequences including
Hamming distance and evolutionary models.
"""

from typing import Optional, Callable, List
import numpy as np
from math import log

from .sequence import Sequence
from .alignment import Alignment


def hamming_distance(seq1: Sequence, seq2: Sequence, ignore_gaps: bool = True) -> int:
    """Calculate Hamming distance between two sequences."""
    if len(seq1) != len(seq2):
        raise ValueError(f"Sequences must have same length: {len(seq1)} vs {len(seq2)}")
    
    differences = 0
    for c1, c2 in zip(seq1.data, seq2.data):
        if ignore_gaps and (c1 == '-' or c2 == '-'):
            continue
        if c1 != c2:
            differences += 1
    
    return differences


def p_distance(seq1: Sequence, seq2: Sequence, ignore_gaps: bool = True) -> float:
    """Calculate p-distance (proportion of differing sites)."""
    if len(seq1) != len(seq2):
        raise ValueError(f"Sequences must have same length: {len(seq1)} vs {len(seq2)}")
    
    differences = 0
    compared_sites = 0
    
    for c1, c2 in zip(seq1.data, seq2.data):
        if ignore_gaps and (c1 == '-' or c2 == '-'):
            continue
        compared_sites += 1
        if c1 != c2:
            differences += 1
    
    if compared_sites == 0:
        raise ValueError("No valid sites to compare")
    
    return differences / compared_sites


def jukes_cantor_distance(seq1: Sequence, seq2: Sequence, ignore_gaps: bool = True) -> float:
    """Calculate Jukes-Cantor corrected distance."""
    p = p_distance(seq1, seq2, ignore_gaps)
    
    if p >= 0.75:
        raise ValueError(f"Sequences too divergent for Jukes-Cantor correction (p={p:.3f} >= 0.75)")
    
    distance = -0.75 * log(1 - (4.0 / 3.0) * p)
    return distance


def kimura_2p_distance(seq1: Sequence, seq2: Sequence, ignore_gaps: bool = True) -> float:
    """Calculate Kimura 2-parameter distance."""
    if len(seq1) != len(seq2):
        raise ValueError(f"Sequences must have same length: {len(seq1)} vs {len(seq2)}")
    
    transitions = 0
    transversions = 0
    compared_sites = 0
    
    transition_pairs = {('A', 'G'), ('G', 'A'), ('C', 'T'), ('T', 'C')}
    
    for c1, c2 in zip(seq1.data, seq2.data):
        if ignore_gaps and (c1 == '-' or c2 == '-'):
            continue
        if c1 in 'N?' or c2 in 'N?':
            continue
        
        compared_sites += 1
        
        if c1 != c2:
            if (c1, c2) in transition_pairs:
                transitions += 1
            else:
                transversions += 1
    
    if compared_sites == 0:
        raise ValueError("No valid sites to compare")
    
    P = transitions / compared_sites
    Q = transversions / compared_sites
    
    term1 = 1 - 2*P - Q
    term2 = 1 - 2*Q
    
    if term1 <= 0 or term2 <= 0:
        raise ValueError(f"Sequences too divergent for K2P correction (P={P:.3f}, Q={Q:.3f})")
    
    distance = -0.5 * log(term1 * (term2 ** 0.5))
    return distance


class DistanceMatrix:
    """Store and manage pairwise distance matrix."""
    
    def __init__(self, labels: List[str], matrix: Optional[np.ndarray] = None):
        """Initialize distance matrix."""
        self.labels = labels
        self.n = len(labels)
        self._label_index = {label: i for i, label in enumerate(labels)}
        
        if matrix is not None:
            if matrix.shape != (self.n, self.n):
                raise ValueError(f"Matrix shape {matrix.shape} doesn't match labels {self.n}")
            self.matrix = matrix
        else:
            self.matrix = np.zeros((self.n, self.n))
    
    def get_distance(self, label1: str, label2: str) -> float:
        """Get distance between two sequences by label."""
        i = self._label_index[label1]
        j = self._label_index[label2]
        return self.matrix[i, j]
    
    def set_distance(self, label1: str, label2: str, distance: float) -> None:
        """Set distance between two sequences."""
        i = self._label_index[label1]
        j = self._label_index[label2]
        self.matrix[i, j] = distance
        self.matrix[j, i] = distance
    
    def get_row(self, label: str) -> np.ndarray:
        """Get all distances for a sequence."""
        i = self._label_index[label]
        return self.matrix[i, :]
    
    def get_min_distance(self, exclude_zero: bool = True) -> float:
        """Get minimum distance in matrix."""
        if exclude_zero:
            mask = np.triu(np.ones_like(self.matrix, dtype=bool), k=1)
            return self.matrix[mask].min()
        else:
            return self.matrix.min()
    
    def get_max_distance(self) -> float:
        """Get maximum distance in matrix."""
        return self.matrix.max()
    
    def to_dict(self) -> dict:
        """Convert matrix to dictionary representation."""
        return {'labels': self.labels, 'matrix': self.matrix.tolist()}
    
    @classmethod
    def from_dict(cls, data: dict) -> 'DistanceMatrix':
        """Create distance matrix from dictionary."""
        return cls(labels=data['labels'], matrix=np.array(data['matrix']))
    
    def __str__(self) -> str:
        """String representation."""
        return f"DistanceMatrix({self.n} sequences)"
    
    def __repr__(self) -> str:
        """Detailed representation."""
        return f"DistanceMatrix(n={self.n}, min={self.get_min_distance():.4f}, max={self.get_max_distance():.4f})"


def calculate_distance_matrix(
    alignment: Alignment,
    distance_func: Optional[Callable[[Sequence, Sequence], float]] = None,
    **kwargs
) -> DistanceMatrix:
    """Calculate pairwise distance matrix for alignment."""
    if distance_func is None:
        distance_func = hamming_distance
    
    n = len(alignment)
    labels = alignment.sequence_ids
    matrix = np.zeros((n, n))
    
    sequences = list(alignment)
    for i in range(n):
        for j in range(i + 1, n):
            try:
                dist = distance_func(sequences[i], sequences[j], **kwargs)
                matrix[i, j] = matrix[j, i] = dist
            except ValueError:
                matrix[i, j] = matrix[j, i] = np.inf
    
    return DistanceMatrix(labels, matrix)


def calculate_pairwise_distances(
    alignment: Alignment,
    method: str = "hamming",
    ignore_gaps: bool = True
) -> DistanceMatrix:
    """Calculate pairwise distances using specified method."""
    method = method.lower()
    
    if method == "hamming":
        distance_func = hamming_distance
    elif method == "p":
        distance_func = p_distance
    elif method in ("jc", "jukes-cantor"):
        distance_func = jukes_cantor_distance
    elif method in ("k2p", "kimura"):
        distance_func = kimura_2p_distance
    else:
        raise ValueError(f"Unknown distance method: {method}")
    
    return calculate_distance_matrix(alignment, distance_func=distance_func, ignore_gaps=ignore_gaps)

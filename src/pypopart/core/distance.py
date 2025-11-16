"""
Distance calculation and metrics for PyPopART.

Implements various distance metrics for DNA sequences including
Hamming distance and evolutionary models.
"""

from math import log
from typing import TYPE_CHECKING, Callable, List, Optional, Tuple

import numpy as np

from .alignment import Alignment
from .sequence import Sequence

if TYPE_CHECKING:
    import matplotlib.figure

# Try to import optimized Numba versions
try:
    from .distance_optimized import hamming_distance_optimized

    _NUMBA_AVAILABLE = True
except ImportError:
    _NUMBA_AVAILABLE = False


def hamming_distance(
    seq1: Sequence,
    seq2: Sequence,
    ignore_gaps: bool = True,
    use_numba: bool = True,
) -> int:
    """
    Calculate Hamming distance between two sequences.

    Parameters
    ----------
    seq1 :
        Sequence.
        First sequence.
    seq2 :
        Sequence.
        Second sequence.
    ignore_gaps :
        bool, default=True.
        Whether to ignore gap characters ('-').
    use_numba :
        bool, default=True.
        Use Numba-optimized version if available.

    Returns
    -------
        int        Number of differing positions.

    Raises
    ------
    ValueError
        If sequences have different lengths

    Notes
    -----
    When use_numba=True and Numba is available, uses JIT-compiled
    optimized version for better performance on large datasets.

    N and ? characters are treated as ambiguous and do not count as
    mutations when compared to any base (A, T, G, C) or to each other.
    """
    if len(seq1) != len(seq2):
        raise ValueError(f'Sequences must have same length: {len(seq1)} vs {len(seq2)}')

    # Use Numba optimized version if available and requested
    if use_numba and _NUMBA_AVAILABLE:
        return hamming_distance_optimized(seq1, seq2, ignore_gaps)

    # Fall back to pure Python implementation
    differences = 0
    for c1, c2 in zip(seq1.data, seq2.data):
        if ignore_gaps and (c1 == '-' or c2 == '-'):
            continue
        # Skip positions with N or ? (ambiguous bases)
        if c1 in 'N?' or c2 in 'N?':
            continue
        if c1 != c2:
            differences += 1

    return differences


def p_distance(seq1: Sequence, seq2: Sequence, ignore_gaps: bool = True) -> float:
    """
    Calculate p-distance (proportion of differing sites).

    N and ? characters are treated as ambiguous and do not count as
    mutations when compared to any base (A, T, G, C) or to each other.
    """
    if len(seq1) != len(seq2):
        raise ValueError(f'Sequences must have same length: {len(seq1)} vs {len(seq2)}')

    differences = 0
    compared_sites = 0

    for c1, c2 in zip(seq1.data, seq2.data):
        if ignore_gaps and (c1 == '-' or c2 == '-'):
            continue
        # Skip positions with N or ? (ambiguous bases)
        if c1 in 'N?' or c2 in 'N?':
            continue
        compared_sites += 1
        if c1 != c2:
            differences += 1

    if compared_sites == 0:
        raise ValueError('No valid sites to compare')

    return differences / compared_sites


def jukes_cantor_distance(
    seq1: Sequence, seq2: Sequence, ignore_gaps: bool = True
) -> float:
    """Calculate Jukes-Cantor corrected distance."""
    p = p_distance(seq1, seq2, ignore_gaps)

    if p >= 0.75:
        raise ValueError(
            f'Sequences too divergent for Jukes-Cantor correction (p={p:.3f} >= 0.75)'
        )

    distance = -0.75 * log(1 - (4.0 / 3.0) * p)
    return distance


def kimura_2p_distance(
    seq1: Sequence, seq2: Sequence, ignore_gaps: bool = True
) -> float:
    """Calculate Kimura 2-parameter distance."""
    if len(seq1) != len(seq2):
        raise ValueError(f'Sequences must have same length: {len(seq1)} vs {len(seq2)}')

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
        raise ValueError('No valid sites to compare')

    P = transitions / compared_sites
    Q = transversions / compared_sites

    term1 = 1 - 2 * P - Q
    term2 = 1 - 2 * Q

    if term1 <= 0 or term2 <= 0:
        raise ValueError(
            f'Sequences too divergent for K2P correction (P={P:.3f}, Q={Q:.3f})'
        )

    distance = -0.5 * log(term1 * (term2**0.5))
    return distance


def tamura_nei_distance(
    seq1: Sequence, seq2: Sequence, ignore_gaps: bool = True
) -> float:
    """
    Calculate Tamura-Nei distance.

    The Tamura-Nei model accounts for:
    - Different base frequencies (GC content)
    - Different rates for transitions within purines (A<->G) and pyrimidines (C<->T)
    - Different rate for transversions

    Args:
        seq1: First sequence
        seq2: Second sequence
        ignore_gaps: Whether to ignore gap positions

    Returns
    -------
        Tamura-Nei corrected distance.

    Raises
    ------
        ValueError: If sequences have different lengths or are too divergent

    Reference:
        Tamura K, Nei M (1993) Mol Biol Evol 10(3):512-526
    """
    if len(seq1) != len(seq2):
        raise ValueError(f'Sequences must have same length: {len(seq1)} vs {len(seq2)}')

    # Count base frequencies and differences
    purine_transitions = 0  # A<->G
    pyrimidine_transitions = 0  # C<->T
    transversions = 0
    compared_sites = 0

    base_counts = {'A': 0, 'G': 0, 'C': 0, 'T': 0}

    for c1, c2 in zip(seq1.data, seq2.data):
        if ignore_gaps and (c1 == '-' or c2 == '-'):
            continue
        if c1 in 'N?' or c2 in 'N?':
            continue

        compared_sites += 1

        # Count base frequencies
        if c1 in base_counts:
            base_counts[c1] += 1
        if c2 in base_counts:
            base_counts[c2] += 1

        # Count differences
        if c1 != c2:
            if (c1, c2) in {('A', 'G'), ('G', 'A')}:
                purine_transitions += 1
            elif (c1, c2) in {('C', 'T'), ('T', 'C')}:
                pyrimidine_transitions += 1
            else:
                transversions += 1

    if compared_sites == 0:
        raise ValueError('No valid sites to compare')

    # Calculate base frequencies
    total_bases = sum(base_counts.values())
    if total_bases == 0:
        raise ValueError('No valid bases found')

    freq_A = base_counts['A'] / total_bases
    freq_G = base_counts['G'] / total_bases
    freq_C = base_counts['C'] / total_bases
    freq_T = base_counts['T'] / total_bases

    # Calculate GC content components
    freq_purines = freq_A + freq_G
    freq_pyrimidines = freq_C + freq_T

    # Calculate proportions of differences
    P1 = purine_transitions / compared_sites  # A<->G transitions
    P2 = pyrimidine_transitions / compared_sites  # C<->T transitions
    Q = transversions / compared_sites  # Transversions

    # Calculate GC content and related parameters
    freq_G + freq_C

    # Handle edge cases
    if freq_purines == 0 or freq_pyrimidines == 0:
        # Fall back to Kimura 2-parameter
        P = P1 + P2
        term1 = 1 - 2 * P - Q
        term2 = 1 - 2 * Q
        if term1 <= 0 or term2 <= 0:
            raise ValueError(
                f'Sequences too divergent for TN correction (P={P:.3f}, Q={Q:.3f})'
            )
        return -0.5 * log(term1 * (term2**0.5))

    # Calculate correction terms
    h_R = (
        2 * freq_A * freq_G / freq_purines if freq_purines > 0 else 0
    )  # Purine heterozygosity
    h_Y = (
        2 * freq_C * freq_T / freq_pyrimidines if freq_pyrimidines > 0 else 0
    )  # Pyrimidine heterozygosity

    # Calculate w terms for distance calculation
    # Use small epsilon to avoid division by zero
    epsilon = 1e-10

    if h_R > epsilon:
        w1 = 1 - P1 / h_R - Q / (2 * freq_purines * freq_pyrimidines)
    else:
        # When h_R is very small, use a simplified form
        w1 = (
            1 - P1 - Q / (2 * freq_purines * freq_pyrimidines)
            if freq_purines * freq_pyrimidines > 0
            else 1
        )

    if h_Y > epsilon:
        w2 = 1 - P2 / h_Y - Q / (2 * freq_purines * freq_pyrimidines)
    else:
        # When h_Y is very small, use a simplified form
        w2 = (
            1 - P2 - Q / (2 * freq_purines * freq_pyrimidines)
            if freq_purines * freq_pyrimidines > 0
            else 1
        )

    w3 = (
        1 - Q / (2 * freq_pyrimidines * freq_purines)
        if freq_purines * freq_pyrimidines > 0
        else 1
    )

    # Check for divergence
    if w1 <= 0 or w2 <= 0 or w3 <= 0:
        raise ValueError(
            f'Sequences too divergent for Tamura-Nei correction '
            f'(P1={P1:.3f}, P2={P2:.3f}, Q={Q:.3f})'
        )

    # Calculate Tamura-Nei distance
    if h_R > epsilon and h_Y > epsilon:
        # Full Tamura-Nei formula
        distance = (
            -h_R * log(w1)
            - h_Y * log(w2)
            - (freq_purines * freq_pyrimidines - h_R * h_Y / (h_R + h_Y)) * log(w3)
        )
    else:
        # Simplified form when heterozygosity is low
        distance = -0.5 * log(w1 * w2 * w3)

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
                raise ValueError(
                    f"Matrix shape {matrix.shape} doesn't match labels {self.n}"
                )
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

    def visualize(
        self,
        title: str = 'Distance Matrix',
        cmap: str = 'viridis',
        figsize: Optional[Tuple[int, int]] = None,
        show_values: bool = True,
        save_path: Optional[str] = None,
    ) -> 'matplotlib.figure.Figure':
        """
            Visualize distance matrix as a heatmap.

        Parameters
        ----------
            title :
                Plot title.
            cmap :
                Matplotlib colormap name.
            figsize :
                Figure size (width, height), auto-calculated if None.
            show_values :
                Whether to show distance values in cells.
            save_path :
                Optional path to save figure.

        Returns
        -------
            Matplotlib figure object.

            Raises :
            ImportError :
                If matplotlib is not installed.
        """
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            raise ImportError(
                'Matplotlib is required for visualization. '
                'Install it with: pip install matplotlib'
            ) from None

        # Auto-calculate figure size based on matrix size
        if figsize is None:
            size = max(8, min(20, self.n * 0.5))
            figsize = (size, size)

        fig, ax = plt.subplots(figsize=figsize)

        # Create heatmap
        im = ax.imshow(self.matrix, cmap=cmap, aspect='auto')

        # Set ticks and labels
        ax.set_xticks(np.arange(self.n))
        ax.set_yticks(np.arange(self.n))
        ax.set_xticklabels(self.labels, rotation=45, ha='right')
        ax.set_yticklabels(self.labels)

        # Add colorbar
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label('Distance', rotation=270, labelpad=15)

        # Show values in cells if requested and matrix is not too large
        if show_values and self.n <= 20:
            for i in range(self.n):
                for j in range(self.n):
                    value = self.matrix[i, j]
                    if np.isfinite(value):
                        text_color = (
                            'white' if value > self.matrix.max() * 0.5 else 'black'
                        )
                        ax.text(
                            j,
                            i,
                            f'{value:.2f}',
                            ha='center',
                            va='center',
                            color=text_color,
                            fontsize=8,
                        )

        ax.set_title(title)
        fig.tight_layout()

        if save_path:
            fig.savefig(save_path, dpi=300, bbox_inches='tight')

        return fig

    def to_csv(self, filepath: str, delimiter: str = ',') -> None:
        """
        Export distance matrix to CSV file.

        Parameters
        ----------
        filepath :
            Path to output CSV file.
        delimiter :
            Delimiter character (default: comma).
        """
        with open(filepath, 'w') as f:
            # Write header
            f.write(delimiter.join([''] + self.labels) + '\n')

            # Write data rows
            for i, label in enumerate(self.labels):
                row_data = [label] + [f'{val:.6f}' for val in self.matrix[i]]
                f.write(delimiter.join(row_data) + '\n')

    @classmethod
    def from_csv(cls, filepath: str, delimiter: str = ',') -> 'DistanceMatrix':
        """
            Import distance matrix from CSV file.

        Parameters
        ----------
            filepath :
                Path to CSV file.
            delimiter :
                Delimiter character (default: comma).

        Returns
        -------
            DistanceMatrix object.
        """
        with open(filepath, 'r') as f:
            lines = [line.strip() for line in f if line.strip()]

        # Parse header (column labels)
        header = lines[0].split(delimiter)
        labels = header[1:]  # Skip first empty cell

        # Parse data
        n = len(labels)
        matrix = np.zeros((n, n))

        for i, line in enumerate(lines[1:]):
            parts = line.split(delimiter)
            parts[0]
            values = [float(v) for v in parts[1:]]
            matrix[i] = values

        return cls(labels, matrix)

    def __str__(self) -> str:
        """Return string representation."""
        return f'DistanceMatrix({self.n} sequences)'

    def __repr__(self) -> str:
        """Detailed representation."""
        return f'DistanceMatrix(n={self.n}, min={self.get_min_distance():.4f}, max={self.get_max_distance():.4f})'


def calculate_distance_matrix(
    alignment: Alignment,
    distance_func: Optional[Callable[[Sequence, Sequence], float]] = None,
    **kwargs,
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
    alignment: Alignment, method: str = 'hamming', ignore_gaps: bool = True
) -> DistanceMatrix:
    """Calculate pairwise distances using specified method."""
    method = method.lower()

    if method == 'hamming':
        distance_func = hamming_distance
    elif method == 'p':
        distance_func = p_distance
    elif method in ('jc', 'jukes-cantor'):
        distance_func = jukes_cantor_distance
    elif method in ('k2p', 'kimura'):
        distance_func = kimura_2p_distance
    elif method in ('tn', 'tamura-nei'):
        distance_func = tamura_nei_distance
    else:
        raise ValueError(f'Unknown distance method: {method}')

    return calculate_distance_matrix(
        alignment, distance_func=distance_func, ignore_gaps=ignore_gaps
    )


class DistanceCalculator:
    """
    Convenience class for calculating distance matrices.

    Provides a simple interface for distance calculation with
    different evolutionary models.
    """

    def __init__(self, method: str = 'hamming', ignore_gaps: bool = True):
        """
        Initialize distance calculator.

        Parameters
        ----------
        method :
            str.
            Distance method: 'hamming', 'jc', 'k2p', 'tamura_nei'.
        ignore_gaps :
            bool.
            Whether to ignore gaps in calculations.
        """
        self.method = method.lower()
        self.ignore_gaps = ignore_gaps

        # Map method names to functions
        if self.method == 'hamming':
            self.distance_func = hamming_distance
        elif self.method in ('jc', 'jukes_cantor', 'jukes-cantor'):
            self.distance_func = jukes_cantor_distance
        elif self.method in ('k2p', 'kimura', 'kimura_2p'):
            self.distance_func = kimura_2p_distance
        elif self.method in ('tn', 'tamura_nei', 'tamura-nei'):
            self.distance_func = tamura_nei_distance
        elif self.method == 'p':
            self.distance_func = p_distance
        else:
            raise ValueError(f'Unknown distance method: {method}')

    def calculate(self, seq1: Sequence, seq2: Sequence) -> float:
        """
        Calculate distance between two sequences.

        Parameters
        ----------
        seq1 :
            Sequence.
            First sequence.
        seq2 :
            Sequence.
            Second sequence.

        Returns
        -------
            float            Distance value.
        """
        return self.distance_func(seq1, seq2, ignore_gaps=self.ignore_gaps)

    def calculate_matrix(self, alignment: Alignment) -> np.ndarray:
        """
        Calculate pairwise distance matrix for alignment.

        Parameters
        ----------
        alignment :
            Alignment.
            Sequence alignment.

        Returns
        -------
            np.ndarray            Square distance matrix.
        """
        dist_matrix = calculate_pairwise_distances(
            alignment, method=self.method, ignore_gaps=self.ignore_gaps
        )
        return dist_matrix.matrix

    def __repr__(self) -> str:
        """Return detailed string representation."""
        return f"DistanceCalculator(method='{self.method}', ignore_gaps={self.ignore_gaps})"

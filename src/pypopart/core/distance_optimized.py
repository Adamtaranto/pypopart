"""
Optimized distance calculations using Numba JIT compilation.

This module provides high-performance distance calculations for large datasets
using Numba's just-in-time compilation to native machine code.
"""

from typing import Tuple

import numba
import numpy as np


@numba.jit(nopython=True, cache=True)
def hamming_distance_numba(
    seq1_bytes: np.ndarray, seq2_bytes: np.ndarray, ignore_gaps: bool = True
) -> int:
    """
    Calculate Hamming distance between two sequences using Numba JIT.

    Optimized for performance using compiled native code. Significantly faster
    than pure Python implementation for large sequences or many comparisons.

    Parameters
    ----------
    seq1_bytes : np.ndarray
        First sequence as numpy array of bytes
    seq2_bytes : np.ndarray
        Second sequence as numpy array of bytes
    ignore_gaps : bool, default=True
        Whether to ignore gap characters ('-')

    Returns
    -------
    int
        Number of differing positions

    Notes
    -----
    This function is JIT-compiled and cached for maximum performance.
    First call may be slower due to compilation overhead.
    """
    if len(seq1_bytes) != len(seq2_bytes):
        return -1  # Error indicator

    differences = 0
    gap_byte = ord('-')

    for i in range(len(seq1_bytes)):
        c1 = seq1_bytes[i]
        c2 = seq2_bytes[i]

        if ignore_gaps and (c1 == gap_byte or c2 == gap_byte):
            continue

        if c1 != c2:
            differences += 1

    return differences


@numba.jit(nopython=True, cache=True)
def p_distance_numba(
    seq1_bytes: np.ndarray, seq2_bytes: np.ndarray, ignore_gaps: bool = True
) -> float:
    """
    Calculate p-distance (proportion of differing sites) using Numba JIT.

    Parameters
    ----------
    seq1_bytes : np.ndarray
        First sequence as numpy array of bytes
    seq2_bytes : np.ndarray
        Second sequence as numpy array of bytes
    ignore_gaps : bool, default=True
        Whether to ignore gap characters

    Returns
    -------
    float
        Proportion of differing sites (0.0 to 1.0)
    """
    if len(seq1_bytes) != len(seq2_bytes):
        return -1.0  # Error indicator

    differences = 0
    compared_sites = 0
    gap_byte = ord('-')

    for i in range(len(seq1_bytes)):
        c1 = seq1_bytes[i]
        c2 = seq2_bytes[i]

        if ignore_gaps and (c1 == gap_byte or c2 == gap_byte):
            continue

        compared_sites += 1

        if c1 != c2:
            differences += 1

    if compared_sites == 0:
        return 0.0

    return differences / compared_sites


@numba.jit(nopython=True, cache=True)
def kimura_2p_counts_numba(
    seq1_bytes: np.ndarray, seq2_bytes: np.ndarray, ignore_gaps: bool = True
) -> Tuple[int, int, int]:
    """
    Count transitions and transversions for Kimura 2-parameter distance.

    Parameters
    ----------
    seq1_bytes : np.ndarray
        First sequence as numpy array of bytes
    seq2_bytes : np.ndarray
        Second sequence as numpy array of bytes
    ignore_gaps : bool, default=True
        Whether to ignore gap characters

    Returns
    -------
    tuple of (int, int, int)
        (transitions, transversions, compared_sites)

    Notes
    -----
    Transitions are purine-purine (A<->G) or pyrimidine-pyrimidine (C<->T).
    Transversions are purine-pyrimidine substitutions.
    """
    if len(seq1_bytes) != len(seq2_bytes):
        return (-1, -1, -1)  # Error indicator

    transitions = 0
    transversions = 0
    compared_sites = 0

    # Define byte values for nucleotides
    A = ord('A')
    G = ord('G')
    C = ord('C')
    T = ord('T')
    gap = ord('-')
    N = ord('N')
    question = ord('?')

    for i in range(len(seq1_bytes)):
        c1 = seq1_bytes[i]
        c2 = seq2_bytes[i]

        if ignore_gaps and (c1 == gap or c2 == gap):
            continue

        if c1 == N or c1 == question or c2 == N or c2 == question:
            continue

        compared_sites += 1

        if c1 != c2:
            # Check for transitions
            if (c1 == A and c2 == G) or (c1 == G and c2 == A):
                transitions += 1
            elif (c1 == C and c2 == T) or (c1 == T and c2 == C):
                transitions += 1
            else:
                # All other differences are transversions
                transversions += 1

    return (transitions, transversions, compared_sites)


@numba.jit(nopython=True, cache=True, parallel=True)
def pairwise_hamming_matrix_numba(
    sequences: np.ndarray, ignore_gaps: bool = True
) -> np.ndarray:
    """
    Calculate pairwise Hamming distance matrix for multiple sequences.

    Uses parallel computation for improved performance on multi-core systems.

    Parameters
    ----------
    sequences : np.ndarray
        2D array where each row is a sequence (as bytes)
    ignore_gaps : bool, default=True
        Whether to ignore gap characters

    Returns
    -------
    np.ndarray
        Symmetric distance matrix of shape (n_sequences, n_sequences)

    Notes
    -----
    This function uses Numba's parallel execution for significant speedup
    on multi-core CPUs. The parallel pragma distributes the outer loop
    iterations across available cores.
    """
    n_seqs = sequences.shape[0]
    matrix = np.zeros((n_seqs, n_seqs), dtype=np.int32)

    for i in numba.prange(n_seqs):
        for j in range(i + 1, n_seqs):
            dist = hamming_distance_numba(sequences[i], sequences[j], ignore_gaps)
            matrix[i, j] = dist
            matrix[j, i] = dist

    return matrix


def prepare_sequences_for_numba(sequences: list) -> np.ndarray:
    """
    Convert list of sequence strings to NumPy array suitable for Numba.

    Parameters
    ----------
    sequences : list of str
        List of sequence strings

    Returns
    -------
    np.ndarray
        2D array where each row is a sequence as byte values

    Examples
    --------
    >>> seqs = ["ATCG", "ATCC", "GTCG"]
    >>> seq_array = prepare_sequences_for_numba(seqs)
    >>> seq_array.shape
    (3, 4)
    """
    # Convert strings to byte arrays
    max_len = max(len(s) for s in sequences)

    # Create 2D array
    seq_array = np.zeros((len(sequences), max_len), dtype=np.uint8)

    for i, seq in enumerate(sequences):
        seq_bytes = np.frombuffer(seq.encode('ascii'), dtype=np.uint8)
        seq_array[i, : len(seq_bytes)] = seq_bytes

    return seq_array


# Wrapper functions that convert from Sequence objects to Numba format


def hamming_distance_optimized(seq1, seq2, ignore_gaps: bool = True) -> int:
    """
    Calculate Hamming distance with automatic Numba optimization.

    This is a wrapper around the Numba-optimized function that handles
    Sequence objects and string conversion automatically.

    Parameters
    ----------
    seq1 : Sequence or str
        First sequence
    seq2 : Sequence or str
        Second sequence
    ignore_gaps : bool, default=True
        Whether to ignore gap characters

    Returns
    -------
    int
        Hamming distance

    Notes
    -----
    For single comparisons, the overhead of conversion may outweigh
    the performance benefit. Use this primarily for batch operations
    or repeated calls where JIT compilation is amortized.
    """
    # Convert to strings if Sequence objects
    s1 = seq1.data if hasattr(seq1, 'data') else str(seq1)
    s2 = seq2.data if hasattr(seq2, 'data') else str(seq2)

    # Convert to NumPy byte arrays
    seq1_bytes = np.frombuffer(s1.encode('ascii'), dtype=np.uint8)
    seq2_bytes = np.frombuffer(s2.encode('ascii'), dtype=np.uint8)

    return hamming_distance_numba(seq1_bytes, seq2_bytes, ignore_gaps)


def p_distance_optimized(seq1, seq2, ignore_gaps: bool = True) -> float:
    """
    Calculate p-distance with automatic Numba optimization.

    Parameters
    ----------
    seq1 : Sequence or str
        First sequence
    seq2 : Sequence or str
        Second sequence
    ignore_gaps : bool, default=True
        Whether to ignore gap characters

    Returns
    -------
    float
        Proportion of differing sites
    """
    s1 = seq1.data if hasattr(seq1, 'data') else str(seq1)
    s2 = seq2.data if hasattr(seq2, 'data') else str(seq2)

    seq1_bytes = np.frombuffer(s1.encode('ascii'), dtype=np.uint8)
    seq2_bytes = np.frombuffer(s2.encode('ascii'), dtype=np.uint8)

    return p_distance_numba(seq1_bytes, seq2_bytes, ignore_gaps)


def kimura_2p_counts_optimized(
    seq1, seq2, ignore_gaps: bool = True
) -> Tuple[int, int, int]:
    """
    Count transitions and transversions with Numba optimization.

    Parameters
    ----------
    seq1 : Sequence or str
        First sequence
    seq2 : Sequence or str
        Second sequence
    ignore_gaps : bool, default=True
        Whether to ignore gap characters

    Returns
    -------
    tuple of (int, int, int)
        (transitions, transversions, compared_sites)
    """
    s1 = seq1.data if hasattr(seq1, 'data') else str(seq1)
    s2 = seq2.data if hasattr(seq2, 'data') else str(seq2)

    seq1_bytes = np.frombuffer(s1.encode('ascii'), dtype=np.uint8)
    seq2_bytes = np.frombuffer(s2.encode('ascii'), dtype=np.uint8)

    return kimura_2p_counts_numba(seq1_bytes, seq2_bytes, ignore_gaps)

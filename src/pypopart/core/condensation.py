"""
Haplotype condensation - identifying unique sequences from alignments.

This module provides functions for identifying unique haplotypes from
sequence alignments and calculating their frequencies.
"""

from collections import defaultdict
from typing import Dict, List, Tuple

from .alignment import Alignment
from .haplotype import Haplotype
from .sequence import Sequence


def condense_alignment(
    alignment: Alignment, populations: Dict[str, str] = None
) -> Tuple[List[Haplotype], Dict[str, List[str]]]:
    """
    Identify unique haplotypes from an alignment.

    Condenses identical sequences into unique haplotypes with frequency
    information. Optionally tracks population assignments.

    Parameters
    ----------
    alignment : Alignment
        Sequence alignment to condense
    populations : dict, optional
        Mapping from sequence ID to population name

    Returns
    -------
    haplotypes : list of Haplotype
        List of unique haplotypes
    frequency_map : dict
        Mapping from haplotype ID to list of sample IDs

    Examples
    --------
    >>> from pypopart.io import load_alignment
    >>> alignment = load_alignment('sequences.fasta')
    >>> haplotypes, freq_map = condense_alignment(alignment)
    >>> print(f"Found {len(haplotypes)} unique haplotypes")
    """
    # Group sequences by their data (unique haplotypes)
    haplotype_map: Dict[str, List[str]] = defaultdict(list)
    sequence_dict: Dict[str, Sequence] = {}

    for seq in alignment:
        # Use sequence data as key
        seq_data = seq.data
        haplotype_map[seq_data].append(seq.id)

        # Store one representative sequence for each haplotype
        if seq_data not in sequence_dict:
            sequence_dict[seq_data] = seq

    # Create Haplotype objects
    haplotypes = []
    frequency_map = {}

    for seq_data, sample_ids in haplotype_map.items():
        # Get representative sequence
        seq = sequence_dict[seq_data]

        # Create new sequence for haplotype with unique ID
        hap_id = f'Hap{len(haplotypes) + 1}'
        hap_seq = Sequence(
            id=hap_id,
            data=seq.data,
            metadata=seq.metadata.copy() if seq.metadata else {},
            description=f'Haplotype {len(haplotypes) + 1} ({len(sample_ids)} samples)',
        )

        # Create haplotype
        haplotype = Haplotype(
            sequence=hap_seq, sample_ids=sample_ids, populations=populations
        )

        haplotypes.append(haplotype)
        frequency_map[hap_id] = sample_ids

    return haplotypes, frequency_map

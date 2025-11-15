"""
Haplotype representation and management for PyPopART.

Haplotypes represent unique sequence variants with associated frequency
and population information.
"""

from collections import defaultdict
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Set

from .sequence import Sequence


@dataclass
class HaplotypeFrequency:
    """Frequency information for a haplotype across populations."""

    total: int
    by_population: Dict[str, int]


class Haplotype:
    """
    Represents a unique haplotype (sequence variant) with frequency information.

    A haplotype is a unique DNA sequence that may be shared by multiple
    individuals or samples. This class tracks the sequence, frequency,
    and population assignments.
    """

    def __init__(
        self,
        sequence: Sequence,
        sample_ids: Optional[List[str]] = None,
        populations: Optional[Dict[str, str]] = None,
    ):
        """
        Initialize a haplotype.

        Parameters
        ----------
        sequence :
            The unique sequence for this haplotype.
        sample_ids :
            List of sample IDs sharing this haplotype.
        populations :
            Dictionary mapping sample_id -> population/group.
        """
        self.sequence = sequence
        self._sample_ids: Set[str] = set(sample_ids) if sample_ids else set()
        self._populations: Dict[str, str] = populations if populations else {}

    @property
    def id(self) -> str:
        """Get haplotype ID (same as sequence ID)."""
        return self.sequence.id

    @property
    def data(self) -> str:
        """Get sequence data."""
        return self.sequence.data

    @property
    def frequency(self) -> int:
        """Get total frequency (number of samples with this haplotype)."""
        return len(self._sample_ids)

    @property
    def sample_ids(self) -> List[str]:
        """Get list of sample IDs with this haplotype."""
        return sorted(self._sample_ids)

    def add_sample(self, sample_id: str, population: Optional[str] = None) -> None:
        """
        Add a sample to this haplotype.

        Parameters
        ----------
        sample_id :
            Sample identifier.
        population :
            Optional population/group assignment.
        """
        self._sample_ids.add(sample_id)
        if population:
            self._populations[sample_id] = population

    def remove_sample(self, sample_id: str) -> None:
        """
        Remove a sample from this haplotype.

        Parameters
        ----------
        sample_id :
            Sample identifier.

        Raises :
        KeyError :
            If sample not found.
        """
        if sample_id not in self._sample_ids:
            raise KeyError(f"Sample '{sample_id}' not in haplotype")

        self._sample_ids.discard(sample_id)
        self._populations.pop(sample_id, None)

    def get_population(self, sample_id: str) -> Optional[str]:
        """
            Get population assignment for a sample.

        Parameters
        ----------
            sample_id :
                Sample identifier.

        Returns
        -------
            Population name or None if not assigned.
        """
        return self._populations.get(sample_id)

    def get_populations(self) -> Set[str]:
        """
        Get set of all populations represented in this haplotype.

        Returns
        -------
            Set of population names.
        """
        return set(self._populations.values())

    def get_frequency_by_population(self) -> Dict[str, int]:
        """
        Calculate frequency of this haplotype per population.

        Returns
        -------
            Dictionary mapping population -> count.
        """
        freq_by_pop: Dict[str, int] = defaultdict(int)
        for sample_id in self._sample_ids:
            pop = self._populations.get(sample_id, 'Unassigned')
            freq_by_pop[pop] += 1
        return dict(freq_by_pop)

    def get_frequency_info(self) -> HaplotypeFrequency:
        """
        Get comprehensive frequency information.

        Returns
        -------
            HaplotypeFrequency object with total and per-population frequencies.
        """
        return HaplotypeFrequency(
            total=self.frequency, by_population=self.get_frequency_by_population()
        )

    def __len__(self) -> int:
        """Return sequence length."""
        return len(self.sequence)

    def __eq__(self, other: object) -> bool:
        """Check equality based on sequence data."""
        if not isinstance(other, Haplotype):
            return False
        return self.sequence == other.sequence

    def __hash__(self) -> int:
        """Hash based on sequence data."""
        return hash(self.sequence)

    def __str__(self) -> str:
        """Return string representation."""
        return f'Haplotype {self.id}: {self.frequency} samples, {len(self.get_populations())} populations'

    def __repr__(self) -> str:
        """Detailed representation."""
        return f"Haplotype(id='{self.id}', frequency={self.frequency}, populations={len(self.get_populations())})"

    def to_dict(self) -> Dict[str, Any]:
        """
        Convert haplotype to dictionary representation.

        Returns
        -------
            Dictionary with haplotype data.
        """
        return {
            'id': self.id,
            'sequence': self.sequence.to_dict(),
            'frequency': self.frequency,
            'sample_ids': self.sample_ids,
            'populations': self._populations,
            'frequency_by_population': self.get_frequency_by_population(),
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'Haplotype':
        """
            Create haplotype from dictionary.

        Parameters
        ----------
            data :
                Dictionary with haplotype information.

        Returns
        -------
            New Haplotype object.
        """
        sequence = Sequence.from_dict(data['sequence'])
        return cls(
            sequence=sequence,
            sample_ids=data.get('sample_ids', []),
            populations=data.get('populations', {}),
        )


def identify_haplotypes_from_alignment(
    alignment: 'Alignment', population_map: Optional[Dict[str, str]] = None
) -> List[Haplotype]:
    """
    Identify unique haplotypes from an alignment.

    Groups sequences by unique sequence data (ignoring gaps) and creates
    Haplotype objects with frequency information.

    Args:
        alignment: Multiple sequence alignment
        population_map: Optional dictionary mapping sequence_id -> population

    Returns
    -------
        List of Haplotype objects.
    """
    # Group sequences by unique haplotype
    haplotype_map: Dict[str, List[str]] = {}
    sequence_map: Dict[str, Sequence] = {}

    for seq in alignment:
        # Use ungapped sequence as haplotype key
        ungapped = seq.remove_gaps()
        key = ungapped.data

        if key not in haplotype_map:
            haplotype_map[key] = []
            sequence_map[key] = ungapped

        haplotype_map[key].append(seq.id)

    # Create Haplotype objects
    haplotypes = []
    for i, (key, sample_ids) in enumerate(sorted(haplotype_map.items()), start=1):
        sequence = sequence_map[key]
        # Give haplotype a numbered ID
        haplotype_seq = Sequence(
            id=f'H{i}',
            data=sequence.data,
            metadata=sequence.metadata.copy(),
            description=f'Haplotype {i} (n={len(sample_ids)})',
        )

        # Map samples to populations if provided
        populations = {}
        if population_map:
            for sample_id in sample_ids:
                if sample_id in population_map:
                    populations[sample_id] = population_map[sample_id]

        haplotype = Haplotype(
            sequence=haplotype_seq, sample_ids=sample_ids, populations=populations
        )
        haplotypes.append(haplotype)

    return haplotypes


def calculate_haplotype_diversity(haplotypes: List[Haplotype]) -> Dict[str, float]:
    """
    Calculate haplotype diversity metrics.

    Args:
        haplotypes: List of Haplotype objects

    Returns
    -------
        Dictionary with diversity metrics:.
        - num_haplotypes: Number of unique haplotypes
        - total_samples: Total number of samples
        - haplotype_diversity: Gene diversity (Nei 1987)
        - singleton_count: Number of haplotypes with frequency 1
    """
    total_samples = sum(h.frequency for h in haplotypes)
    num_haplotypes = len(haplotypes)

    if total_samples == 0:
        return {
            'num_haplotypes': 0,
            'total_samples': 0,
            'haplotype_diversity': 0.0,
            'singleton_count': 0,
        }

    # Calculate haplotype diversity (gene diversity)
    # H = (n / (n-1)) * (1 - Σ(pi^2))
    # where pi is the frequency of haplotype i
    freq_sum_squared = sum((h.frequency / total_samples) ** 2 for h in haplotypes)

    if total_samples > 1:
        haplotype_diversity = (total_samples / (total_samples - 1)) * (
            1 - freq_sum_squared
        )
    else:
        haplotype_diversity = 0.0

    # Count singletons (haplotypes appearing only once)
    singleton_count = sum(1 for h in haplotypes if h.frequency == 1)

    return {
        'num_haplotypes': num_haplotypes,
        'total_samples': total_samples,
        'haplotype_diversity': haplotype_diversity,
        'singleton_count': singleton_count,
    }

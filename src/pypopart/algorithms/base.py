"""Base classes for network construction algorithms in PyPopART."""

from abc import ABC, abstractmethod
from typing import Any, Dict, Optional

from ..core.alignment import Alignment
from ..core.distance import DistanceMatrix, calculate_pairwise_distances
from ..core.graph import HaplotypeNetwork


class NetworkAlgorithm(ABC):
    """
    Abstract base class for haplotype network construction algorithms.

    All network construction algorithms should inherit from this class
    and implement the construct_network method.
    """

    def __init__(self, distance_method: str = 'hamming', **kwargs):
        """
        Initialize network algorithm.

        Parameters
        ----------
        distance_method :
            Method for calculating distances (hamming, p, jc, k2p, tn).
        **kwargs :
            Additional algorithm-specific parameters.
        """
        self.distance_method = distance_method
        self.params = kwargs
        self._distance_matrix: Optional[DistanceMatrix] = None

    @abstractmethod
    def construct_network(
        self, alignment: Alignment, distance_matrix: Optional[DistanceMatrix] = None
    ) -> HaplotypeNetwork:
        """
            Construct haplotype network from sequence alignment.

        Parameters
        ----------
            alignment :
                Multiple sequence alignment.
            distance_matrix :
                Optional pre-computed distance matrix.

        Returns
        -------
            Constructed haplotype network.
        """
        pass

    def calculate_distances(self, alignment: Alignment) -> DistanceMatrix:
        """
            Calculate pairwise distances between sequences.

        Parameters
        ----------
            alignment :
                Multiple sequence alignment.

        Returns
        -------
            Distance matrix.
        """
        return calculate_pairwise_distances(
            alignment,
            method=self.distance_method,
            ignore_gaps=self.params.get('ignore_gaps', True),
        )

    def build_network(
        self, alignment: Alignment, distance_matrix: Optional[DistanceMatrix] = None
    ) -> HaplotypeNetwork:
        """
        Build haplotype network from sequence alignment.

        This is an alias for construct_network() to maintain backward compatibility
        with CLI and GUI code.

        Parameters
        ----------
            alignment :
                Multiple sequence alignment.
            distance_matrix :
                Optional pre-computed distance matrix.

        Returns
        -------
            Constructed haplotype network.
        """
        return self.construct_network(alignment, distance_matrix)

    def get_parameters(self) -> Dict[str, Any]:
        """
        Get algorithm parameters.

        Returns
        -------
            Dictionary of parameters.
        """
        return {'distance_method': self.distance_method, **self.params}

    def __str__(self) -> str:
        """Return string representation."""
        return f'{self.__class__.__name__}(distance={self.distance_method})'

    def __repr__(self) -> str:
        """Detailed representation."""
        params_str = ', '.join(f'{k}={v}' for k, v in self.get_parameters().items())
        return f'{self.__class__.__name__}({params_str})'

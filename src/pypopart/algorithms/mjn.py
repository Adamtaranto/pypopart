"""
Median-Joining Network (MJN) algorithm for haplotype network construction.

Implements the method from Bandelt, Forster & Rohl (1999):
"Median-joining networks for inferring intraspecific phylogenies"
Molecular Biology and Evolution 16: 37-48
"""

from typing import Optional, List, Tuple, Set, Dict
import itertools
from ..core.alignment import Alignment
from ..core.graph import HaplotypeNetwork
from ..core.distance import DistanceMatrix, hamming_distance
from ..core.haplotype import Haplotype, identify_haplotypes_from_alignment as identify_haplotypes
from ..core.sequence import Sequence
from .msn import MinimumSpanningNetwork


class MedianJoiningNetwork(MinimumSpanningNetwork):
    """
    Construct haplotype network using Median-Joining algorithm.
    
    The MJN method extends MSN by inferring median vectors (ancestral or
    unsampled haplotypes) that simplify the network structure. It combines
    minimum spanning network principles with median vector inference to
    create a more parsimonious representation of haplotype relationships.
    
    This is particularly useful for datasets with missing intermediate
    haplotypes and complex reticulation patterns.
    """
    
    def __init__(
        self,
        distance_method: str = "hamming",
        epsilon: float = 0.0,
        max_median_vectors: Optional[int] = None,
        simplify: bool = True,
        **kwargs
    ):
        """
        Initialize MJN algorithm.
        
        Args:
            distance_method: Method for calculating distances
            epsilon: Weight parameter for controlling network complexity
            max_median_vectors: Maximum number of median vectors to add
            simplify: Whether to simplify network after median vector addition
            **kwargs: Additional parameters
        """
        super().__init__(distance_method, epsilon=epsilon, **kwargs)
        self.max_median_vectors = max_median_vectors
        self.simplify = simplify
        self._median_counter = 0
    
    def construct_network(
        self,
        alignment: Alignment,
        distance_matrix: Optional[DistanceMatrix] = None
    ) -> HaplotypeNetwork:
        """
        Construct MJN from sequence alignment.
        
        Args:
            alignment: Multiple sequence alignment
            distance_matrix: Optional pre-computed distance matrix
            
        Returns:
            Haplotype network with inferred median vectors
        """
        # Identify unique haplotypes
        haplotypes = identify_haplotypes(alignment)
        
        if len(haplotypes) <= 2:
            # Too few haplotypes for median vector inference
            return super().construct_network(alignment, distance_matrix)
        
        # Calculate or use provided distance matrix
        if distance_matrix is None:
            distance_matrix = self.calculate_distances(alignment)
        
        self._distance_matrix = distance_matrix
        
        # Build initial MSN
        initial_network = super().construct_network(alignment, distance_matrix)
        
        # Infer and add median vectors
        network_with_medians = self._add_median_vectors(
            initial_network,
            alignment.length
        )
        
        # Simplify network if requested
        if self.simplify:
            final_network = self._simplify_network(network_with_medians)
        else:
            final_network = network_with_medians
        
        return final_network
    
    def _add_median_vectors(
        self,
        network: HaplotypeNetwork,
        sequence_length: int
    ) -> HaplotypeNetwork:
        """
        Infer and add median vectors to the network.
        
        Args:
            network: Initial haplotype network
            sequence_length: Length of sequences
            
        Returns:
            Network with median vectors added
        """
        self._median_counter = 0
        median_vectors_added = 0
        
        # Get all triangles (3-cliques) in the network
        triangles = self._find_triplets(network)
        
        for hap1_id, hap2_id, hap3_id in triangles:
            if self.max_median_vectors and median_vectors_added >= self.max_median_vectors:
                break
            
            # Get the three haplotypes
            hap1 = network.get_haplotype(hap1_id)
            hap2 = network.get_haplotype(hap2_id)
            hap3 = network.get_haplotype(hap3_id)
            
            # Calculate median vector
            median_seq = self._calculate_median(
                hap1.sequence,
                hap2.sequence,
                hap3.sequence
            )
            
            if median_seq is None:
                continue
            
            # Check if median already exists in network
            if self._sequence_in_network(median_seq, network):
                continue
            
            # Calculate distances from median to the three haplotypes
            median_hap = Haplotype(
                sequence=median_seq,
                id=f"Median_{self._median_counter}",
                frequency=0
            )
            self._median_counter += 1
            
            dist1 = hamming_distance(median_hap.sequence, hap1.sequence, ignore_gaps=True)
            dist2 = hamming_distance(median_hap.sequence, hap2.sequence, ignore_gaps=True)
            dist3 = hamming_distance(median_hap.sequence, hap3.sequence, ignore_gaps=True)
            
            # Check if median vector simplifies the network
            # (reduces total edge weight)
            original_edges = [
                network.get_edge_distance(hap1_id, hap2_id) or 0,
                network.get_edge_distance(hap1_id, hap3_id) or 0,
                network.get_edge_distance(hap2_id, hap3_id) or 0
            ]
            original_weight = sum(d for d in original_edges if d > 0)
            new_weight = dist1 + dist2 + dist3
            
            if new_weight < original_weight:
                # Add median vector to network
                network.add_haplotype(median_hap)
                
                # Connect median to the three haplotypes
                if dist1 > 0:
                    network.add_edge(median_hap.id, hap1_id, distance=dist1)
                if dist2 > 0:
                    network.add_edge(median_hap.id, hap2_id, distance=dist2)
                if dist3 > 0:
                    network.add_edge(median_hap.id, hap3_id, distance=dist3)
                
                median_vectors_added += 1
        
        return network
    
    def _find_triplets(self, network: HaplotypeNetwork) -> List[Tuple[str, str, str]]:
        """
        Find all triplets (triangles) in the network.
        
        Args:
            network: Haplotype network
            
        Returns:
            List of triplets (id1, id2, id3)
        """
        triplets = []
        hap_ids = list({h.id for h in network.haplotypes})
        
        # Check all combinations of 3 haplotypes
        for hap1, hap2, hap3 in itertools.combinations(hap_ids, 3):
            # Check if they form a triangle (all three edges exist)
            has_12 = network.has_edge(hap1, hap2)
            has_13 = network.has_edge(hap1, hap3)
            has_23 = network.has_edge(hap2, hap3)
            
            if has_12 and has_13 and has_23:
                triplets.append((hap1, hap2, hap3))
        
        return triplets
    
    def _calculate_median(
        self,
        seq1: Sequence,
        seq2: Sequence,
        seq3: Sequence
    ) -> Optional[Sequence]:
        """
        Calculate median sequence of three sequences.
        
        For each position, the median is the most common nucleotide
        among the three sequences. If all three are different, no
        clear median exists for that position.
        
        Args:
            seq1, seq2, seq3: Three Sequence objects
            
        Returns:
            Median Sequence, or None if no clear median exists
        """
        if len(seq1) != len(seq2) or len(seq1) != len(seq3):
            return None
        
        median_data = []
        
        for i in range(len(seq1)):
            c1, c2, c3 = seq1.data[i], seq2.data[i], seq3.data[i]
            
            # Find most common nucleotide
            if c1 == c2:
                median_data.append(c1)
            elif c1 == c3:
                median_data.append(c1)
            elif c2 == c3:
                median_data.append(c2)
            else:
                # All three different - no clear median
                # Use majority rule or first sequence's base
                # For simplicity, we'll say no median exists for this triplet
                return None
        
        median_seq = Sequence(
            id="median_temp",
            data=''.join(median_data)
        )
        
        return median_seq
    
    def _sequence_in_network(
        self,
        sequence: Sequence,
        network: HaplotypeNetwork
    ) -> bool:
        """
        Check if a sequence already exists in the network.
        
        Args:
            sequence: Sequence to check
            network: Haplotype network
            
        Returns:
            True if sequence exists in network
        """
        seq_data = sequence.data
        
        for haplotype in network.haplotypes:
            if haplotype.sequence.data == seq_data:
                return True
        
        return False
    
    def _simplify_network(self, network: HaplotypeNetwork) -> HaplotypeNetwork:
        """
        Simplify network by removing unnecessary median vectors.
        
        A median vector is unnecessary if:
        - It has degree 2 (only connects two nodes)
        - It can be replaced by a direct edge without increasing total weight
        
        Args:
            network: Network with median vectors
            
        Returns:
            Simplified network
        """
        # Find median vectors with degree 2
        to_remove = []
        
        for hap_id in list({h.id for h in network.haplotypes}):
            if not hap_id.startswith("Median_"):
                continue
            
            degree = network.get_degree(hap_id)
            
            if degree == 2:
                # Get the two neighbors
                neighbors = network.get_neighbors(hap_id)
                if len(neighbors) == 2:
                    n1, n2 = neighbors
                    
                    # Get distances
                    d1 = network.get_edge_distance(hap_id, n1) or 0
                    d2 = network.get_edge_distance(hap_id, n2) or 0
                    
                    # Check if direct connection exists
                    direct_dist = network.get_edge_distance(n1, n2)
                    
                    if direct_dist is None:
                        # No direct connection - add it with combined distance
                        network.add_edge(n1, n2, distance=d1 + d2)
                    
                    # Remove the median vector
                    to_remove.append(hap_id)
        
        # Remove marked median vectors
        for hap_id in to_remove:
            network.remove_haplotype(hap_id)
        
        return network
    
    def get_parameters(self) -> dict:
        """Get algorithm parameters."""
        params = super().get_parameters()
        params['max_median_vectors'] = self.max_median_vectors
        params['simplify'] = self.simplify
        return params

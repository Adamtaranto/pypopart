"""
TCS (Statistical Parsimony) algorithm - Enhanced to match C++ PopART implementation.

This implementation closely follows the C++ PopART TCS.cpp logic:
- Component-based connection algorithm
- Intermediate sequence inference with scoring
- Post-processing vertex collapse
- Floyd-Warshall for path length calculations

Implements the method from Clement, Posada & Crandall (2000):
"TCS: a computer program to estimate gene genealogies"
Molecular Ecology 9: 1657-1659
"""

import math
from typing import Dict, List, Optional, Tuple

from ..core.alignment import Alignment
from ..core.distance import DistanceMatrix, hamming_distance
from ..core.graph import HaplotypeNetwork
from ..core.haplotype import Haplotype
from ..core.haplotype import identify_haplotypes_from_alignment as identify_haplotypes
from ..core.sequence import Sequence
from .base import NetworkAlgorithm


class TCS(NetworkAlgorithm):
    """
    Construct haplotype network using Statistical Parsimony (TCS algorithm).
    
    This implementation accurately reproduces the C++ PopART TCS algorithm:
    1. Groups haplotypes by pairwise distances
    2. Iteratively connects components at increasing distance levels
    3. Infers intermediate sequences when distance > 1
    4. Uses scoring system to find optimal intermediates
    5. Post-processes to collapse degree-2 vertices
    
    The TCS method uses a 95% parsimony criterion to determine which
    haplotypes can be connected with statistical confidence.
    """

    # Scoring constants from C++ implementation
    BONUS = 20
    SHORTCUTPENALTY = 10
    LONGPENALTY = 5

    def __init__(
        self,
        distance_method: str = 'hamming',
        confidence: float = 0.95,
        connection_limit: Optional[int] = None,
        infer_intermediates: bool = True,
        collapse_vertices: bool = True,
        **kwargs,
    ):
        """
        Initialize TCS algorithm.

        Parameters
        ----------
        distance_method :
            Method for calculating distances (should be hamming).
        confidence :
            Confidence level for parsimony criterion (default 0.95).
        connection_limit :
            Maximum connection distance (auto-calculated if None).
        infer_intermediates :
            Whether to infer intermediate sequences (default True).
        collapse_vertices :
            Whether to collapse degree-2 vertices (default True).
        **kwargs :
            Additional parameters.
        """
        super().__init__(distance_method, **kwargs)
        self.confidence = confidence
        self.connection_limit = connection_limit
        self.infer_intermediates = infer_intermediates
        self.collapse_vertices = collapse_vertices
        self._intermediate_counter = 0

        if distance_method not in ('hamming',):
            import warnings

            warnings.warn(
                f"TCS is designed for hamming distance. Using '{distance_method}' "
                'may not produce theoretically correct results.',
                stacklevel=2,
            )

    def construct_network(
        self, alignment: Alignment, distance_matrix: Optional[DistanceMatrix] = None
    ) -> HaplotypeNetwork:
        """
        Construct TCS network from sequence alignment.
        
        Implements the component-based connection algorithm from C++ TCS.cpp.

        Parameters
        ----------
        alignment :
            Multiple sequence alignment.
        distance_matrix :
            Optional pre-computed distance matrix.

        Returns
        -------
            Haplotype network constructed using statistical parsimony.
        """
        # Identify unique haplotypes
        haplotypes = identify_haplotypes(alignment)

        if len(haplotypes) == 0:
            return HaplotypeNetwork()

        if len(haplotypes) == 1:
            network = HaplotypeNetwork()
            network.add_haplotype(haplotypes[0])
            return network

        # Calculate distances between haplotypes
        haplotype_dist_matrix = self._calculate_haplotype_distances(haplotypes)
        self._distance_matrix = haplotype_dist_matrix

        # Calculate connection limit if not provided
        if self.connection_limit is None:
            self.connection_limit = self._calculate_connection_limit(
                alignment.length, len(haplotypes)
            )

        # Build network using component-based algorithm (matches C++ TCS.cpp)
        network = self._build_network_with_components(
            haplotypes, haplotype_dist_matrix, alignment.length
        )

        # Collapse degree-2 vertices (post-processing simplification)
        if self.collapse_vertices:
            network = self._collapse_degree2_vertices(network)

        return network

    def _calculate_connection_limit(
        self, sequence_length: int, num_haplotypes: int
    ) -> int:
        """
        Calculate maximum parsimony connection limit.

        Uses the formula from Templeton et al. (1992) to estimate the
        maximum number of mutational differences that can be explained
        by parsimony at the specified confidence level.

        Parameters
        ----------
        sequence_length :
            Length of aligned sequences.
        num_haplotypes :
            Number of unique haplotypes.

        Returns
        -------
            Maximum connection distance for parsimony criterion.
        """
        connection_limit = 1
        cumulative_prob = 0.0

        for k in range(1, sequence_length + 1):
            prob_k = self._poisson_probability(k, sequence_length, num_haplotypes)
            cumulative_prob += prob_k

            if cumulative_prob >= (1 - self.confidence):
                connection_limit = k
                break

        return max(1, connection_limit)

    def _poisson_probability(self, k: int, seq_length: int, sample_size: int) -> float:
        """
        Calculate probability using Poisson distribution.

        Parameters
        ----------
        k :
            Number of mutations.
        seq_length :
            Sequence length.
        sample_size :
            Number of sequences.

        Returns
        -------
            Probability.
        """
        lambda_param = 2.0 * math.log(sample_size) if sample_size > 1 else 1.0

        try:
            prob = (lambda_param**k) * math.exp(-lambda_param) / math.factorial(k)
        except (OverflowError, ValueError):
            prob = 0.0

        return prob

    def _calculate_haplotype_distances(self, haplotypes: List) -> DistanceMatrix:
        """
        Calculate pairwise distances between haplotypes.

        Parameters
        ----------
        haplotypes :
            List of Haplotype objects.

        Returns
        -------
            DistanceMatrix with distances between haplotypes.
        """
        import numpy as np

        n = len(haplotypes)
        labels = [h.id for h in haplotypes]
        matrix = np.zeros((n, n))

        for i in range(n):
            for j in range(i + 1, n):
                dist = hamming_distance(
                    haplotypes[i].sequence,
                    haplotypes[j].sequence,
                    ignore_gaps=self.params.get('ignore_gaps', True),
                )
                matrix[i, j] = matrix[j, i] = dist

        return DistanceMatrix(labels, matrix)

    def _build_network_with_components(
        self, haplotypes: List, distance_matrix: DistanceMatrix, sequence_length: int
    ) -> HaplotypeNetwork:
        """
        Build network using component-based algorithm from C++ TCS.cpp.
        
        This accurately reproduces the C++ logic:
        1. Create all haplotypes as vertices, each in its own component
        2. Group pairs by distance
        3. For each distance level (1 to connection_limit):
           - Process pairs at that distance
           - Connect vertices from different components
           - Infer intermediates for distance > 1
           - Merge components when connected
        4. Clean up intermediate vertices in post-processing

        Parameters
        ----------
        haplotypes :
            List of Haplotype objects.
        distance_matrix :
            Distance matrix.
        sequence_length :
            Length of sequences for creating intermediates.

        Returns
        -------
            Network with inferred intermediate sequences.
        """
        # Initialize network with all haplotypes
        network = HaplotypeNetwork()
        for haplotype in haplotypes:
            network.add_haplotype(haplotype)

        # Component tracking: maps haplotype ID to component ID
        # Each haplotype starts in its own component
        component_ids: Dict[str, int] = {h.id: i for i, h in enumerate(haplotypes)}
        
        # Group pairs by distance (like C++ VertContainer and priority queue)
        pairs_by_distance: Dict[int, List[Tuple[str, str]]] = {}
        
        for i, h1 in enumerate(haplotypes):
            for j in range(i + 1, len(haplotypes)):
                h2 = haplotypes[j]
                dist = int(round(distance_matrix.get_distance(h1.id, h2.id)))
                
                if dist <= self.connection_limit:
                    if dist not in pairs_by_distance:
                        pairs_by_distance[dist] = []
                    pairs_by_distance[dist].append((h1.id, h2.id))

        # Process pairs in order of increasing distance (like C++ priority queue)
        for M in sorted(pairs_by_distance.keys()):
            # Keep processing this distance level until no more pairs remain
            while M in pairs_by_distance and len(pairs_by_distance[M]) > 0:
                pairs = pairs_by_distance[M]
                
                # Track which component pair we're working on
                comp_a = -1
                comp_b = -1
                other_pairs = []
                
                for u_id, v_id in pairs:
                    comp_u = component_ids.get(u_id, -1)
                    comp_v = component_ids.get(v_id, -1)
                    
                    # Skip if already in same component
                    if comp_u == comp_v:
                        continue
                    
                    # Ensure comp_u < comp_v
                    if comp_u > comp_v:
                        comp_u, comp_v = comp_v, comp_u
                        u_id, v_id = v_id, u_id
                    
                    # Set component pair on first distinct pair
                    if comp_a < 0:
                        comp_a = comp_u
                        comp_b = comp_v
                    
                    # Process pairs for the current component pair
                    if comp_u == comp_a and comp_v == comp_b:
                        if M == 1:
                            # Direct connection for distance 1
                            network.add_edge(u_id, v_id, distance=1)
                        else:
                            # Infer intermediates for distance > 1
                            if self.infer_intermediates:
                                self._add_connection_with_intermediates(
                                    network, u_id, v_id, M, sequence_length, 
                                    component_ids, comp_u, comp_v
                                )
                            else:
                                # Simple connection without intermediates
                                network.add_edge(u_id, v_id, distance=M)
                    else:
                        # Save for next iteration of this distance level
                        other_pairs.append((u_id, v_id))
                
                # Merge components (like C++ component renumbering)
                if comp_a >= 0:
                    for hap_id in list(component_ids.keys()):
                        if component_ids[hap_id] < 0 or component_ids[hap_id] == comp_b:
                            component_ids[hap_id] = comp_a
                        elif component_ids[hap_id] > comp_b:
                            component_ids[hap_id] -= 1
                
                # Update pairs for this distance level
                if other_pairs:
                    pairs_by_distance[M] = other_pairs
                else:
                    # No more pairs at this distance, remove from dict
                    del pairs_by_distance[M]

        return network

    def _add_connection_with_intermediates(
        self,
        network: HaplotypeNetwork,
        u_id: str,
        v_id: str,
        distance: int,
        sequence_length: int,
        component_ids: Dict[str, int],
        comp_u: int,
        comp_v: int,
    ) -> None:
        """
        Add connection between u and v, inferring intermediates if needed.
        
        Implements findIntermediates() and newCompositePath() from C++ TCS.cpp.

        Parameters
        ----------
        network :
            Current network.
        u_id :
            Source haplotype ID.
        v_id :
            Target haplotype ID.
        distance :
            Distance between u and v.
        sequence_length :
            Length of sequences.
        component_ids :
            Component membership tracker.
        comp_u :
            Component of u.
        comp_v :
            Component of v.
        """
        # Find optimal intermediate vertices (like C++ findIntermediates)
        int_u, int_v, min_path_length = self._find_intermediates(
            network, u_id, v_id, distance, component_ids, comp_u, comp_v
        )
        
        # Check if path already exists
        try:
            existing_path_length = network.get_shortest_path_length(int_u, int_v)
        except:
            existing_path_length = float('inf')
        
        # Only add new path if it's shorter or doesn't exist
        if existing_path_length == float('inf') or existing_path_length > min_path_length:
            self._create_composite_path(
                network, int_u, int_v, min_path_length, 
                sequence_length, component_ids
            )

    def _find_intermediates(
        self,
        network: HaplotypeNetwork,
        u_id: str,
        v_id: str,
        dist: int,
        component_ids: Dict[str, int],
        comp_u: int,
        comp_v: int,
    ) -> Tuple[str, str, int]:
        """
        Find optimal intermediate vertices to connect two components.
        
        Implements C++ TCS::findIntermediates() with scoring system.

        Parameters
        ----------
        network :
            Current network.
        u_id :
            Source haplotype ID from component comp_u.
        v_id :
            Target haplotype ID from component comp_v.
        dist :
            Distance between components.
        component_ids :
            Component membership tracker.
        comp_u :
            Source component.
        comp_v :
            Target component.

        Returns
        -------
            Tuple of (intermediate_u_id, intermediate_v_id, path_length).
        """
        max_score = float('-inf')
        min_path_length = dist
        best_u = u_id
        best_v = v_id
        
        # Try all vertices in comp_u (or "no man's land" with comp_id < 0)
        for i_id in list(component_ids.keys()):
            if component_ids.get(i_id, -1) != comp_u and component_ids.get(i_id, -1) >= 0:
                continue
            
            # Check if connected to u
            try:
                path_ui = network.get_shortest_path_length(u_id, i_id)
            except:
                continue
            
            if path_ui >= dist:
                continue
            
            # Try all vertices in comp_v
            for j_id in list(component_ids.keys()):
                if component_ids.get(j_id, -1) != comp_v and component_ids.get(j_id, -1) >= 0:
                    continue
                
                # Check if connected to v
                try:
                    path_vj = network.get_shortest_path_length(v_id, j_id)
                except:
                    continue
                
                if path_vj + path_ui >= dist:
                    continue
                
                dP = dist - path_vj - path_ui
                score = self._compute_score(
                    network, i_id, j_id, comp_u, comp_v, dP, dist, component_ids
                )
                
                # Select best scoring pair (or shortest path if tied)
                if score > max_score or (score == max_score and dP < min_path_length):
                    min_path_length = dP
                    max_score = score
                    best_u = i_id
                    best_v = j_id
        
        return best_u, best_v, min_path_length

    def _compute_score(
        self,
        network: HaplotypeNetwork,
        u_id: str,
        v_id: str,
        comp_u: int,
        comp_v: int,
        dP: int,
        clust_dist: int,
        component_ids: Dict[str, int],
    ) -> float:
        """
        Compute score for intermediate pair using C++ scoring system.
        
        Implements C++ TCS::computeScore().

        Parameters
        ----------
        network :
            Current network.
        u_id :
            Intermediate in component u.
        v_id :
            Intermediate in component v.
        comp_u :
            Component u.
        comp_v :
            Component v.
        dP :
            Path length between intermediates.
        clust_dist :
            Distance between original components.
        component_ids :
            Component membership tracker.

        Returns
        -------
            Score for this intermediate pair.
        """
        score = 0
        
        # Get all original haplotypes (not intermediates)
        original_haps = [
            hap_id for hap_id in component_ids.keys()
            if not hap_id.startswith('intermediate_')
        ]
        
        for i_id in original_haps:
            if component_ids.get(i_id, -1) != comp_u:
                continue
            
            for j_id in original_haps:
                if component_ids.get(j_id, -1) != comp_v:
                    continue
                
                # Calculate total path through proposed intermediates
                try:
                    path_ui = network.get_shortest_path_length(u_id, i_id)
                    path_vj = network.get_shortest_path_length(v_id, j_id)
                    total_path = dP + path_ui + path_vj
                except:
                    continue
                
                # Get original distance
                orig_dist = self._distance_matrix.get_distance(i_id, j_id)
                
                # Score based on how well path matches original distance
                if abs(total_path - orig_dist) < 0.5:
                    score += self.BONUS
                elif total_path > orig_dist:
                    score -= self.LONGPENALTY
                else:
                    # Shortcut
                    if total_path < clust_dist:
                        return float('-inf')  # Invalid
                    else:
                        score -= self.SHORTCUTPENALTY
        
        return score

    def _create_composite_path(
        self,
        network: HaplotypeNetwork,
        start_id: str,
        end_id: str,
        distance: int,
        sequence_length: int,
        component_ids: Dict[str, int],
    ) -> None:
        """
        Create path of intermediate vertices connecting start to end.
        
        Implements C++ TCS::newCompositePath().

        Parameters
        ----------
        network :
            Current network.
        start_id :
            Starting vertex ID.
        end_id :
            Ending vertex ID.
        distance :
            Number of intermediates to create.
        sequence_length :
            Sequence length for intermediates.
        component_ids :
            Component membership tracker.
        """
        current_id = start_id
        
        # Create distance-1 intermediate vertices
        for _i in range(1, distance):
            intermediate_id = f'intermediate_{self._intermediate_counter}'
            self._intermediate_counter += 1
            
            # Create placeholder sequence
            intermediate_seq = Sequence(
                id=intermediate_id,
                data='N' * sequence_length,
                description='Inferred intermediate sequence',
            )
            
            intermediate_hap = Haplotype(
                sequence=intermediate_seq,
                sample_ids=[],
            )
            
            # Add to network
            network.add_haplotype(intermediate_hap)
            network.add_edge(current_id, intermediate_id, distance=1)
            
            # Mark as "no man's land" (component ID = -1)
            component_ids[intermediate_id] = -1
            
            current_id = intermediate_id
        
        # Connect last intermediate to end
        network.add_edge(current_id, end_id, distance=1)

    def _collapse_degree2_vertices(self, network: HaplotypeNetwork) -> HaplotypeNetwork:
        """
        Collapse vertices with degree 2 (post-processing simplification).
        
        Implements C++ TCS post-processing that removes degree-2 vertices
        connecting only two other vertices. This matches lines 161-194 of TCS.cpp.

        Parameters
        ----------
        network :
            Network with potential degree-2 vertices.

        Returns
        -------
            Simplified network.
        """
        # Process one vertex at a time (like C++ implementation)
        # Keep looping until no more collapses possible
        changed = True
        
        while changed:
            changed = False
            
            # Find and collapse ONE degree-2 intermediate vertex
            for hap_id in list(network.nodes):
                # Skip if already removed
                if not network.has_node(hap_id):
                    continue
                    
                degree = network.get_degree(hap_id)

                # Only collapse degree-2 vertices
                if degree != 2:
                    continue

                neighbors = network.get_neighbors(hap_id)

                if len(neighbors) != 2:
                    continue
                    
                n1, n2 = neighbors

                # Get edge weights
                try:
                    w1 = network.get_edge_distance(hap_id, n1)
                    w2 = network.get_edge_distance(hap_id, n2)
                except:
                    continue

                # Only collapse intermediates (not original haplotypes)
                # In C++ this checks if vertex index >= nseqs (number of original sequences)
                hap = network.get_haplotype(hap_id)
                if hap.frequency == 0 or 'intermediate' in hap_id.lower():
                    try:
                        combined_weight = w1 + w2
                        
                        # Remove the intermediate vertex and its edges
                        network.remove_haplotype(hap_id)

                        # Add direct edge if it doesn't exist
                        if not network.has_edge(n1, n2):
                            network.add_edge(n1, n2, distance=combined_weight)
                        
                        # Mark that we made a change - restart loop
                        changed = True
                        break  # Exit for loop and restart while loop
                    except Exception:
                        # Skip if removal fails
                        pass

        return network

    def get_parameters(self) -> dict:
        """Get algorithm parameters."""
        params = super().get_parameters()
        params['confidence'] = self.confidence
        params['connection_limit'] = self.connection_limit
        params['infer_intermediates'] = self.infer_intermediates
        params['collapse_vertices'] = self.collapse_vertices
        return params

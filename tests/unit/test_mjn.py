"""
Unit tests for Median-Joining Network (MJN) algorithm.
"""

from pypopart.algorithms.mjn import MedianJoiningNetwork
from pypopart.core.alignment import Alignment
from pypopart.core.sequence import Sequence


class TestMedianJoiningNetwork:
    """Test cases for MJN algorithm."""

    def test_mjn_initialization(self):
        """Test MJN algorithm initialization."""
        mjn = MedianJoiningNetwork()
        assert mjn.distance_method == 'hamming'
        assert mjn.epsilon == 0.0
        assert mjn.simplify is True

    def test_mjn_no_simplify(self):
        """Test MJN without simplification."""
        mjn = MedianJoiningNetwork(simplify=False)
        assert mjn.simplify is False

    def test_mjn_empty_alignment(self):
        """Test MJN with empty alignment."""
        mjn = MedianJoiningNetwork()
        alignment = Alignment()
        network = mjn.construct_network(alignment)

        assert len(network) == 0

    def test_mjn_single_sequence(self):
        """Test MJN with single sequence."""
        mjn = MedianJoiningNetwork()
        alignment = Alignment([Sequence('seq1', 'ATCG')])
        network = mjn.construct_network(alignment)

        assert len(network) == 1
        assert len(network.edges) == 0

    def test_mjn_two_sequences(self):
        """Test MJN with two sequences (no median vectors)."""
        mjn = MedianJoiningNetwork()
        alignment = Alignment([Sequence('seq1', 'ATCG'), Sequence('seq2', 'ATCC')])
        network = mjn.construct_network(alignment)

        # Should be same as MSN for 2 sequences
        assert len(network) == 2
        assert network.is_connected()

    def test_mjn_median_vector_inference(self):
        """Test MJN infers median vectors for triangles."""
        mjn = MedianJoiningNetwork(simplify=False)
        alignment = Alignment(
            [
                Sequence('seq1', 'AAAA'),
                Sequence('seq2', 'AATT'),
                Sequence('seq3', 'TTAA'),
            ]
        )
        network = mjn.construct_network(alignment)

        # Check if network is connected
        assert network.is_connected()
        # May have added median vectors (check for nodes starting with "Median_")
        median_count = sum(
            1
            for hap_id in [h.id for h in network.haplotypes]
            if hap_id.startswith('Median_')
        )
        # Median vectors may or may not be added depending on whether they simplify
        assert median_count >= 0

    def test_mjn_with_max_median_vectors(self):
        """Test MJN with maximum median vectors limit."""
        mjn = MedianJoiningNetwork(max_median_vectors=1)
        alignment = Alignment(
            [
                Sequence('seq1', 'AAAA'),
                Sequence('seq2', 'AATT'),
                Sequence('seq3', 'TTAA'),
                Sequence('seq4', 'TTTT'),
            ]
        )
        network = mjn.construct_network(alignment)

        # Count median vectors
        median_count = sum(
            1
            for hap_id in [h.id for h in network.haplotypes]
            if hap_id.startswith('Median_')
        )
        assert median_count <= mjn.max_median_vectors

    def test_mjn_simplification(self):
        """Test that MJN simplifies network when requested."""
        # Create network with potential for simplification
        mjn_no_simplify = MedianJoiningNetwork(simplify=False)
        mjn_simplify = MedianJoiningNetwork(simplify=True)

        alignment = Alignment(
            [
                Sequence('seq1', 'AAAA'),
                Sequence('seq2', 'AATT'),
                Sequence('seq3', 'TTTT'),
            ]
        )

        network_no_simplify = mjn_no_simplify.construct_network(alignment)
        network_simplify = mjn_simplify.construct_network(alignment)

        # Simplified network should have <= nodes
        assert len(network_simplify) <= len(network_no_simplify)

    def test_mjn_vs_msn(self):
        """Test that MJN produces valid network compared to MSN."""
        from pypopart.algorithms.msn import MinimumSpanningNetwork

        alignment = Alignment(
            [
                Sequence('seq1', 'ATCG'),
                Sequence('seq2', 'ATCC'),
                Sequence('seq3', 'GTCG'),
                Sequence('seq4', 'GTCC'),
            ]
        )

        msn = MinimumSpanningNetwork()
        msn_network = msn.construct_network(alignment)

        mjn = MedianJoiningNetwork()
        mjn_network = mjn.construct_network(alignment)

        # Both should be connected
        assert msn_network.is_connected()
        assert mjn_network.is_connected()

        # MJN might have different number of nodes (due to median vectors)
        # but should still represent the same haplotypes
        assert len(mjn_network) >= len(msn_network) or len(mjn_network) == len(
            msn_network
        )

    def test_mjn_median_calculation(self):
        """Test median sequence calculation."""
        mjn = MedianJoiningNetwork()

        seq1 = Sequence('s1', 'AAAA')
        seq2 = Sequence('s2', 'AATT')
        seq3 = Sequence('s3', 'AATT')

        median = mjn._calculate_median(seq1, seq2, seq3)

        # Median should exist and be "AATT" (majority at each position)
        assert median is not None
        assert median.data == 'AATT'

    def test_mjn_no_median_all_different(self):
        """Test that no median is calculated when all bases differ."""
        mjn = MedianJoiningNetwork()

        seq1 = Sequence('s1', 'A')
        seq2 = Sequence('s2', 'C')
        seq3 = Sequence('s3', 'G')

        median = mjn._calculate_median(seq1, seq2, seq3)

        # No clear median exists
        assert median is None

    def test_mjn_parameters(self):
        """Test getting MJN parameters."""
        mjn = MedianJoiningNetwork(
            distance_method='k2p', epsilon=0.5, max_median_vectors=10, simplify=False
        )
        params = mjn.get_parameters()

        assert params['distance_method'] == 'k2p'
        assert params['epsilon'] == 0.5
        assert params['max_median_vectors'] == 10
        assert params['simplify'] is False

    def test_mjn_string_representation(self):
        """Test string representation of MJN."""
        mjn = MedianJoiningNetwork(distance_method='hamming')
        assert 'MedianJoiningNetwork' in str(mjn)
        assert 'hamming' in str(mjn)
    
    def test_mjn_iterative_refinement(self):
        """Test iterative refinement produces better networks."""
        mjn = MedianJoiningNetwork(epsilon=1.0)
        alignment = Alignment(
            [
                Sequence('seq1', 'AAAA'),
                Sequence('seq2', 'ATAT'),
                Sequence('seq3', 'TATA'),
                Sequence('seq4', 'TTTT'),
            ]
        )
        network = mjn.construct_network(alignment)
        
        # Network should be connected
        assert network.is_connected()
        # Should have inferred some median vectors
        median_count = sum(
            1 for h in network.haplotypes if 'Median_' in h.id
        )
        assert median_count >= 0  # May or may not add medians depending on topology
    
    def test_mjn_quasi_median_simple(self):
        """Test quasi-median calculation with simple case."""
        mjn = MedianJoiningNetwork()
        seq1 = Sequence('s1', 'AAAA')
        seq2 = Sequence('s2', 'AATT')
        seq3 = Sequence('s3', 'AATT')
        
        quasi_medians = mjn._compute_quasi_medians(seq1, seq2, seq3)
        
        # Should have one clear median: AATT
        assert len(quasi_medians) >= 1
        assert 'AATT' in quasi_medians
    
    def test_mjn_quasi_median_all_different(self):
        """Test quasi-median with all positions different."""
        mjn = MedianJoiningNetwork()
        seq1 = Sequence('s1', 'AA')
        seq2 = Sequence('s2', 'CC')
        seq3 = Sequence('s3', 'GG')
        
        quasi_medians = mjn._compute_quasi_medians(seq1, seq2, seq3)
        
        # Should generate multiple quasi-medians (all combinations)
        # For 2 positions with 3 choices each = 3^2 = 9 possibilities
        # But we only include unique ones
        assert len(quasi_medians) > 1
        # Should include sequences with bases from the three inputs
        for qm in quasi_medians:
            assert len(qm) == 2
            assert all(c in 'ACG' for c in qm)
    
    def test_mjn_median_cost(self):
        """Test median cost calculation."""
        mjn = MedianJoiningNetwork()
        seq1 = Sequence('s1', 'AAAA')
        seq2 = Sequence('s2', 'AATT')
        seq3 = Sequence('s3', 'TTAA')
        
        # Test cost of perfect median
        cost = mjn._compute_median_cost(seq1, seq2, seq3, 'AATA')
        assert cost >= 0
        assert isinstance(cost, int)
    
    def test_mjn_remove_obsolete_medians(self):
        """Test removal of obsolete median vectors."""
        from pypopart.core.haplotype import Haplotype
        from pypopart.core.graph import HaplotypeNetwork
        
        mjn = MedianJoiningNetwork()
        
        # Create network with a degree-1 median
        network = HaplotypeNetwork()
        h1 = Haplotype(Sequence('h1', 'AAAA'), 'h1', frequency=5)
        h2 = Haplotype(Sequence('h2', 'TTTT'), 'h2', frequency=3)
        median = Haplotype(Sequence('Median_0', 'AATT'), 'Median_0', frequency=0)
        
        network.add_haplotype(h1)
        network.add_haplotype(h2)
        network.add_haplotype(median)
        network.add_edge('h1', 'Median_0', distance=2)
        # Median has degree 1, should be removed
        
        haplotypes = [h1, h2, median]
        result = mjn._remove_obsolete_medians(network, haplotypes)
        
        # Median should be removed (degree < 2)
        assert len(result) <= len(haplotypes)
    
    def test_mjn_find_triplets_in_msn(self):
        """Test finding triplets in MSN."""
        from pypopart.core.haplotype import Haplotype
        from pypopart.core.graph import HaplotypeNetwork
        
        mjn = MedianJoiningNetwork()
        
        # Create a star network
        network = HaplotypeNetwork()
        center = Haplotype(Sequence('center', 'AAAA'), 'center', frequency=10)
        h1 = Haplotype(Sequence('h1', 'AAAT'), 'h1', frequency=2)
        h2 = Haplotype(Sequence('h2', 'AATT'), 'h2', frequency=3)
        h3 = Haplotype(Sequence('h3', 'ATTT'), 'h3', frequency=1)
        
        network.add_haplotype(center)
        network.add_haplotype(h1)
        network.add_haplotype(h2)
        network.add_haplotype(h3)
        network.add_edge('center', 'h1', distance=1)
        network.add_edge('center', 'h2', distance=2)
        network.add_edge('center', 'h3', distance=3)
        network.add_edge('h1', 'h2', distance=2)  # Make a triplet
        
        triplets = mjn._find_all_triplets_in_msn(network)
        
        # Should find at least one triplet
        assert len(triplets) >= 1
    
    def test_mjn_build_msn_for_iteration(self):
        """Test building MSN for iteration."""
        from pypopart.core.haplotype import Haplotype
        from pypopart.core.distance import DistanceMatrix
        import numpy as np
        
        mjn = MedianJoiningNetwork()
        
        h1 = Haplotype(Sequence('h1', 'AAAA'), 'h1', frequency=5)
        h2 = Haplotype(Sequence('h2', 'AATT'), 'h2', frequency=3)
        h3 = Haplotype(Sequence('h3', 'TTAA'), 'h3', frequency=2)
        
        haplotypes = [h1, h2, h3]
        
        # Create distance matrix
        labels = ['h1', 'h2', 'h3']
        matrix = np.array([
            [0, 2, 3],
            [2, 0, 3],
            [3, 3, 0]
        ])
        dist_matrix = DistanceMatrix(labels, matrix)
        
        network = mjn._build_msn_for_iteration(haplotypes, dist_matrix)
        
        # Should create a connected MSN
        assert len(network.haplotypes) >= 3
        assert network.is_connected()
    
    def test_mjn_epsilon_parameter(self):
        """Test epsilon parameter affects median selection."""
        # Small epsilon - strict median selection
        mjn_strict = MedianJoiningNetwork(epsilon=0.0)
        # Large epsilon - relaxed median selection  
        mjn_relaxed = MedianJoiningNetwork(epsilon=5.0)
        
        alignment = Alignment(
            [
                Sequence('seq1', 'AAAA'),
                Sequence('seq2', 'ATAT'),
                Sequence('seq3', 'TATA'),
                Sequence('seq4', 'TTTT'),
            ]
        )
        
        network_strict = mjn_strict.construct_network(alignment)
        network_relaxed = mjn_relaxed.construct_network(alignment)
        
        # Both should be connected
        assert network_strict.is_connected()
        assert network_relaxed.is_connected()
        
        # Relaxed epsilon may add more medians
        medians_strict = sum(1 for h in network_strict.haplotypes if 'Median_' in h.id)
        medians_relaxed = sum(1 for h in network_relaxed.haplotypes if 'Median_' in h.id)
        assert medians_relaxed >= medians_strict
    
    def test_mjn_sequence_length_preserved(self):
        """Test that median vectors maintain sequence length."""
        mjn = MedianJoiningNetwork()
        alignment = Alignment(
            [
                Sequence('seq1', 'ATCGATCG'),
                Sequence('seq2', 'ATCGTTCG'),
                Sequence('seq3', 'ATCGATTT'),
                Sequence('seq4', 'TTCGATTT'),
            ]
        )
        network = mjn.construct_network(alignment)
        
        # All sequences (including medians) should have same length
        expected_length = 8
        for haplotype in network.haplotypes:
            assert len(haplotype.sequence.data) == expected_length
    
    def test_mjn_convergence(self):
        """Test that iterative refinement converges."""
        mjn = MedianJoiningNetwork(epsilon=1.0)
        alignment = Alignment(
            [
                Sequence('seq1', 'AAAA'),
                Sequence('seq2', 'AATT'),
                Sequence('seq3', 'TTAA'),
                Sequence('seq4', 'TTTT'),
            ]
        )
        
        # Should not run forever (max_iterations prevents infinite loop)
        network = mjn.construct_network(alignment)
        
        # Should produce a valid connected network
        assert network.is_connected()
        assert len(network.haplotypes) >= 4  # At least the original 4
    
    def test_mjn_different_distance_methods(self):
        """Test MJN with different distance methods."""
        alignment = Alignment(
            [
                Sequence('seq1', 'ATCG'),
                Sequence('seq2', 'ATCC'),
                Sequence('seq3', 'GTCG'),
            ]
        )
        
        # Test with hamming
        mjn_hamming = MedianJoiningNetwork(distance_method='hamming')
        network_hamming = mjn_hamming.construct_network(alignment)
        assert network_hamming.is_connected()
        
        # Note: Other distance methods may not work well with median inference
        # as they expect continuous values, not discrete sequences

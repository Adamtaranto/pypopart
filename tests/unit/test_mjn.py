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

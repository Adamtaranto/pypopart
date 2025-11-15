"""
Unit tests for TCS (Statistical Parsimony) algorithm.
"""

from pypopart.algorithms.tcs import TCS
from pypopart.core.alignment import Alignment
from pypopart.core.sequence import Sequence


class TestTCS:
    """Test cases for TCS algorithm."""

    def test_tcs_initialization(self):
        """Test TCS algorithm initialization."""
        tcs = TCS()
        assert tcs.distance_method == 'hamming'
        assert tcs.confidence == 0.95
        assert tcs.connection_limit is None

    def test_tcs_custom_confidence(self):
        """Test TCS with custom confidence level."""
        tcs = TCS(confidence=0.99)
        assert tcs.confidence == 0.99

    def test_tcs_custom_connection_limit(self):
        """Test TCS with custom connection limit."""
        tcs = TCS(connection_limit=5)
        assert tcs.connection_limit == 5

    def test_tcs_empty_alignment(self):
        """Test TCS with empty alignment."""
        tcs = TCS()
        alignment = Alignment()
        network = tcs.construct_network(alignment)

        assert len(network) == 0

    def test_tcs_single_sequence(self):
        """Test TCS with single sequence."""
        tcs = TCS()
        alignment = Alignment([Sequence('seq1', 'ATCG')])
        network = tcs.construct_network(alignment)

        assert len(network) == 1
        assert len(network.edges) == 0

    def test_tcs_two_sequences(self):
        """Test TCS with two sequences."""
        tcs = TCS()
        alignment = Alignment(
            [
                Sequence('seq1', 'ATCG'),
                Sequence('seq2', 'ATCC'),  # 1 difference
            ]
        )
        network = tcs.construct_network(alignment)

        assert len(network) == 2
        # Should be connected if within connection limit
        if tcs.connection_limit and tcs.connection_limit >= 1:
            assert network.is_connected()

    def test_tcs_respects_connection_limit(self):
        """Test that TCS respects connection limit."""
        tcs = TCS(connection_limit=1)
        alignment = Alignment(
            [
                Sequence('seq1', 'ATCG'),
                Sequence('seq2', 'ATCC'),  # 1 diff - should connect
                Sequence('seq3', 'GGGG'),  # 4 diff - should not connect
            ]
        )
        network = tcs.construct_network(alignment)

        assert len(network) == 3
        # seq1 and seq2 should be connected (1 diff)
        # seq3 should not be connected to others (4 diff > limit)
        assert not network.is_connected()

    def test_tcs_connection_limit_calculation(self):
        """Test automatic connection limit calculation."""
        tcs = TCS()
        alignment = Alignment(
            [
                Sequence('seq1', 'ATCGAT'),
                Sequence('seq2', 'ATCGAC'),
                Sequence('seq3', 'ATCGAG'),
            ]
        )
        tcs.construct_network(alignment)

        # Connection limit should have been calculated
        assert tcs.connection_limit is not None
        assert tcs.connection_limit > 0

    def test_tcs_frequency_based_ordering(self):
        """Test that TCS prioritizes more frequent haplotypes."""
        tcs = TCS(connection_limit=1)
        alignment = Alignment(
            [
                Sequence('seq1', 'ATCG'),  # Will be Haplotype_1 (first)
                Sequence('seq2', 'ATCG'),  # Same as seq1
                Sequence('seq3', 'ATCC'),  # 1 diff from seq1/seq2
            ]
        )
        network = tcs.construct_network(alignment)

        # Should have 2 haplotypes (seq1/seq2 collapsed)
        assert len(network) == 2
        # They should be connected
        assert network.is_connected()

    def test_tcs_parameters(self):
        """Test getting TCS parameters."""
        tcs = TCS(distance_method='hamming', confidence=0.99, connection_limit=10)
        params = tcs.get_parameters()

        assert params['distance_method'] == 'hamming'
        assert params['confidence'] == 0.99
        assert params['connection_limit'] == 10

    def test_tcs_string_representation(self):
        """Test string representation of TCS."""
        tcs = TCS(distance_method='hamming')
        assert 'TCS' in str(tcs)
        assert 'hamming' in str(tcs)

    def test_tcs_disconnected_network(self):
        """Test TCS can produce disconnected networks."""
        tcs = TCS(connection_limit=1)
        alignment = Alignment(
            [
                Sequence('seq1', 'AAAA'),
                Sequence('seq2', 'TTTT'),  # 4 differences - too far
                Sequence('seq3', 'GGGG'),  # 4 differences - too far
            ]
        )
        network = tcs.construct_network(alignment)

        assert len(network) == 3
        # Should not be connected (all distances > limit)
        assert not network.is_connected()

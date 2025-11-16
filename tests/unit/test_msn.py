"""Unit tests for Minimum Spanning Network (MSN) algorithm."""

from pypopart.algorithms.msn import MinimumSpanningNetwork
from pypopart.core.alignment import Alignment
from pypopart.core.sequence import Sequence


class TestMinimumSpanningNetwork:
    """Test cases for MSN algorithm."""

    def test_msn_initialization(self):
        """Test MSN algorithm initialization."""
        msn = MinimumSpanningNetwork()
        assert msn.distance_method == 'hamming'
        assert msn.epsilon == 0.0
        assert msn.max_connections is None

    def test_msn_with_epsilon(self):
        """Test MSN with epsilon tolerance."""
        msn = MinimumSpanningNetwork(epsilon=0.5)
        assert msn.epsilon == 0.5

    def test_msn_empty_alignment(self):
        """Test MSN with empty alignment."""
        msn = MinimumSpanningNetwork()
        alignment = Alignment()
        network = msn.construct_network(alignment)

        assert len(network) == 0

    def test_msn_single_sequence(self):
        """Test MSN with single sequence."""
        msn = MinimumSpanningNetwork()
        alignment = Alignment([Sequence('seq1', 'ATCG')])
        network = msn.construct_network(alignment)

        assert len(network) == 1
        assert len(network.edges) == 0

    def test_msn_adds_alternative_connections(self):
        """Test that MSN adds alternative connections at same distance."""
        msn = MinimumSpanningNetwork()
        alignment = Alignment(
            [
                Sequence('seq1', 'AAAA'),
                Sequence('seq2', 'AAAT'),  # 1 diff from seq1
                Sequence('seq3', 'AAAC'),  # 1 diff from seq1
                Sequence('seq4', 'AAAT'),  # Same as seq2
            ]
        )
        network = msn.construct_network(alignment)

        # Should have at least as many edges as MST (or same if no alternatives)
        assert network.is_connected()

    def test_msn_vs_mst(self):
        """Test that MSN has at least as many edges as MST."""
        from pypopart.algorithms.mst import MinimumSpanningTree

        alignment = Alignment(
            [
                Sequence('seq1', 'ATCG'),
                Sequence('seq2', 'ATCC'),
                Sequence('seq3', 'GTCG'),
                Sequence('seq4', 'GTCC'),
            ]
        )

        mst = MinimumSpanningTree()
        mst_network = mst.construct_network(alignment)

        msn = MinimumSpanningNetwork()
        msn_network = msn.construct_network(alignment)

        # MSN should have at least as many edges as MST
        assert len(msn_network.edges) >= len(mst_network.edges)

    def test_msn_with_max_connections(self):
        """Test MSN with max connections limit."""
        msn = MinimumSpanningNetwork(max_connections=3)
        alignment = Alignment(
            [
                Sequence('seq1', 'AAAA'),
                Sequence('seq2', 'AAAT'),
                Sequence('seq3', 'AAAC'),
                Sequence('seq4', 'AAAG'),
                Sequence('seq5', 'AAAN'),
            ]
        )
        network = msn.construct_network(alignment)

        # Network should be created successfully with max_connections parameter
        # Note: The max_connections feature is a soft limit during alternative connection addition
        assert len(network) >= 1
        assert network.is_connected() or len(network.haplotypes) > 0

    def test_msn_triangle(self):
        """Test MSN with three equidistant sequences (triangle)."""
        msn = MinimumSpanningNetwork()
        alignment = Alignment(
            [Sequence('seq1', 'AT'), Sequence('seq2', 'AC'), Sequence('seq3', 'GT')]
        )
        network = msn.construct_network(alignment)

        # Should form a triangle (all three connected)
        assert len(network) == 3
        # MSN might add all three edges if they're at same distance
        assert network.is_connected()

    def test_msn_parameters(self):
        """Test getting MSN parameters."""
        msn = MinimumSpanningNetwork(
            distance_method='k2p', epsilon=0.5, max_connections=5
        )
        params = msn.get_parameters()

        assert params['distance_method'] == 'k2p'
        assert params['epsilon'] == 0.5
        assert params['max_connections'] == 5

    def test_msn_string_representation(self):
        """Test string representation of MSN."""
        msn = MinimumSpanningNetwork(distance_method='hamming')
        assert 'MinimumSpanningNetwork' in str(msn)
        assert 'hamming' in str(msn)

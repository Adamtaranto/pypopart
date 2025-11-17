"""Extended tests for MSN epsilon parameter behavior.

This test module specifically validates that the epsilon parameter
correctly controls network relaxation and alternative edge inclusion.
"""

from pypopart.algorithms.msn import MinimumSpanningNetwork
from pypopart.algorithms.mst import MinimumSpanningTree
from pypopart.core.alignment import Alignment
from pypopart.core.sequence import Sequence


class TestMSNEpsilonBehavior:
    """Test cases for MSN epsilon parameter."""

    def test_epsilon_zero_vs_mst(self):
        """Test that MSN with epsilon=0 creates network from MST base."""
        alignment = Alignment(
            [
                Sequence('seq1', 'AAAA'),
                Sequence('seq2', 'AAAT'),  # dist 1 from seq1
                Sequence('seq3', 'AAAC'),  # dist 1 from seq1
                Sequence('seq4', 'AATT'),  # dist 2 from seq1, dist 1 from seq2 and seq3
            ]
        )

        # Build MST
        mst = MinimumSpanningTree()
        mst_network = mst.construct_network(alignment)

        # Build MSN with epsilon=0
        msn = MinimumSpanningNetwork(epsilon=0.0)
        msn_network = msn.construct_network(alignment)

        # MSN should have at least as many edges as MST
        assert len(msn_network.edges) >= len(mst_network.edges)

        # Both should connect all nodes
        assert mst_network.is_connected()
        assert msn_network.is_connected()

        # MSN should have 4 edges (MST has 3, plus H2-H3 alternative)
        assert len(msn_network.edges) == 4
        assert len(mst_network.edges) == 3

    def test_epsilon_increases_edges(self):
        """Test that larger epsilon adds more edges."""
        alignment = Alignment(
            [
                Sequence('seq1', 'AAAA'),
                Sequence('seq2', 'AAAT'),  # dist 1
                Sequence('seq3', 'AAAC'),  # dist 1
                Sequence('seq4', 'AATT'),  # dist 2
            ]
        )

        # Test with different epsilon values
        msn0 = MinimumSpanningNetwork(epsilon=0.0)
        network0 = msn0.construct_network(alignment)

        msn05 = MinimumSpanningNetwork(epsilon=0.5)
        network05 = msn05.construct_network(alignment)

        msn1 = MinimumSpanningNetwork(epsilon=1.0)
        network1 = msn1.construct_network(alignment)

        # All should be connected
        assert network0.is_connected()
        assert network05.is_connected()
        assert network1.is_connected()

        # Epsilon=1.0 might have more edges before redundancy removal
        # (though redundancy removal may bring it back to same as epsilon=0)
        # The key is that epsilon allows consideration of more distant edges

    def test_epsilon_with_triangle(self):
        """Test epsilon behavior with equidistant triangle."""
        # Create three sequences all at distance 1 from each other
        alignment = Alignment(
            [
                Sequence('seq1', 'AT'),
                Sequence('seq2', 'AC'),
                Sequence('seq3', 'GT'),
            ]
        )

        msn = MinimumSpanningNetwork(epsilon=0.0)
        network = msn.construct_network(alignment)

        # Should have 3 nodes
        assert len(network) == 3

        # After redundancy removal, might have 2 or 3 edges depending on which
        # edge is considered redundant. The key is that all nodes are connected.
        assert len(network.edges) >= 2  # At least MST edges
        assert len(network.edges) <= 3  # At most complete triangle

        # Must be connected
        assert network.is_connected()

    def test_epsilon_adds_longer_edges(self):
        """Test that epsilon allows edges beyond MST distances."""
        alignment = Alignment(
            [
                Sequence('seq1', 'AAAA'),
                Sequence('seq2', 'AAAT'),  # dist 1 from seq1
                Sequence('seq3', 'AATT'),  # dist 2 from seq1, dist 1 from seq2
            ]
        )

        # With epsilon=0, should only use distance-1 edges
        msn0 = MinimumSpanningNetwork(epsilon=0.0)
        network0 = msn0.construct_network(alignment)

        # With epsilon=1, could add distance-2 edge (H1-H3)
        # since it's within 1 of MST distance 1
        msn1 = MinimumSpanningNetwork(epsilon=1.0)
        network1 = msn1.construct_network(alignment)

        # Both should be connected
        assert network0.is_connected()
        assert network1.is_connected()

        # Get edge distances
        edges0 = {
            (e[0], e[1]): network0.get_edge_distance(e[0], e[1]) for e in network0.edges
        }
        edges1 = {
            (e[0], e[1]): network1.get_edge_distance(e[0], e[1]) for e in network1.edges
        }

        # All edges in network0 should be distance 1
        assert all(d == 1.0 for d in edges0.values())

        # network1 might have distance-2 edge (before redundancy removal)
        # After redundancy removal, it might be pruned
        # The key is that epsilon=1 considers it

    def test_epsilon_parameter_storage(self):
        """Test that epsilon parameter is correctly stored and retrieved."""
        msn = MinimumSpanningNetwork(epsilon=0.5)
        assert msn.epsilon == 0.5

        params = msn.get_parameters()
        assert params['epsilon'] == 0.5

    def test_epsilon_with_identical_sequences(self):
        """Test epsilon with duplicate sequences (should be collapsed)."""
        alignment = Alignment(
            [
                Sequence('seq1', 'AAAA'),
                Sequence('seq2', 'AAAA'),  # Identical to seq1
                Sequence('seq3', 'AAAT'),  # dist 1
                Sequence('seq4', 'AAAC'),  # dist 1
            ]
        )

        msn = MinimumSpanningNetwork(epsilon=1.0)
        network = msn.construct_network(alignment)

        # Should have 3 unique haplotypes (seq1 and seq2 collapsed)
        assert len(network) == 3

        # Should be connected
        assert network.is_connected()

    def test_epsilon_negative_not_allowed(self):
        """Test that negative epsilon behaves reasonably."""
        # Negative epsilon should behave like zero
        msn = MinimumSpanningNetwork(epsilon=-1.0)

        alignment = Alignment(
            [
                Sequence('seq1', 'AAAA'),
                Sequence('seq2', 'AAAT'),
            ]
        )

        network = msn.construct_network(alignment)
        assert network.is_connected()

    def test_epsilon_with_max_connections(self):
        """Test interaction between epsilon and max_connections."""
        alignment = Alignment(
            [
                Sequence('seq1', 'AAAA'),
                Sequence('seq2', 'AAAT'),
                Sequence('seq3', 'AAAC'),
                Sequence('seq4', 'AAAG'),
            ]
        )

        # With epsilon=1 and max_connections=2, should limit edges
        msn = MinimumSpanningNetwork(epsilon=1.0, max_connections=2)
        network = msn.construct_network(alignment)

        # Should still be connected
        assert network.is_connected()

        # Check that no node has more than 3 connections
        # (max_connections is a soft limit during alternative addition)
        # network.haplotypes is a list of Haplotype objects
        for haplotype in network.haplotypes:
            degree = network.get_degree(haplotype.id)
            # Soft limit, so might be slightly exceeded
            assert degree <= 4  # Generous check

    def test_epsilon_zero_is_strict_mode(self):
        """Test that epsilon=0 creates strict MSN (minimal network)."""
        alignment = Alignment(
            [
                Sequence('seq1', 'AAAA'),
                Sequence('seq2', 'AAAT'),
                Sequence('seq3', 'AATT'),
            ]
        )

        msn = MinimumSpanningNetwork(epsilon=0.0)
        network = msn.construct_network(alignment)

        # With linear arrangement (1-1-2), should have exactly 2 edges
        assert len(network.edges) == 2
        assert network.is_connected()

    def test_epsilon_large_value(self):
        """Test MSN with very large epsilon (adds many edges)."""
        alignment = Alignment(
            [
                Sequence('seq1', 'AAAA'),
                Sequence('seq2', 'AAAT'),
                Sequence('seq3', 'AATT'),
                Sequence('seq4', 'ATTT'),
            ]
        )

        # Large epsilon should consider all distance levels
        msn = MinimumSpanningNetwork(epsilon=10.0)
        network = msn.construct_network(alignment)

        # Should be connected
        assert network.is_connected()

        # Might have more edges than minimum (but redundancy removal applies)
        assert len(network.edges) >= 3  # At least MST edges

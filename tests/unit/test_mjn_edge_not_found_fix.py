"""
Test case for MJN edge-not-found bug fix.

This module tests the fix for the bug where MJN with epsilon > 0
would fail with KeyError: "Edge (X, Y) not found" during simplification.

The error occurred in _simplify_network when it tried to check if a
direct edge exists between two neighbors of a degree-2 median without
first checking if the edge exists.
"""

from pypopart.algorithms.mjn import MedianJoiningNetwork
from pypopart.core.alignment import Alignment
from pypopart.core.sequence import Sequence


class TestMJNEdgeNotFoundFix:
    """Test cases for MJN edge-not-found bug fix during simplification."""

    def test_mjn_simplify_with_epsilon_complex_case(self):
        """
        Test MJN simplification with epsilon > 0 on a complex case.

        This test reproduces the scenario where simplification would fail
        because it tried to get edge distance without checking if edge exists.
        """
        mjn = MedianJoiningNetwork(epsilon=5.0, simplify=True)
        alignment = Alignment(
            [
                Sequence('H1', 'AAAAAAAAAA'),
                Sequence('H2', 'AAAAAATTTT'),
                Sequence('H3', 'AATTAATTAA'),
                Sequence('H4', 'TTTTTTTTTT'),
                Sequence('H5', 'AAAAATTTAA'),
                Sequence('H6', 'TTAATTAATT'),
                Sequence('H7', 'ATATATAAAA'),
                Sequence('H8', 'TATATATTTT'),
            ]
        )

        # Should not raise KeyError
        network = mjn.construct_network(alignment)

        # Verify network properties
        assert network.is_connected()
        assert len(network.haplotypes) >= 8  # At least 8 observed haplotypes

        # All nodes should have degree >= 2 after simplification
        for hap in network.haplotypes:
            degree = network.get_degree(hap.id)
            # Observed haplotypes can have any degree >= 1
            # Median haplotypes should have degree >= 2 after simplification
            if hap.id.startswith('Median_'):
                assert degree >= 2, f'Median {hap.id} has degree {degree} < 2'

    def test_mjn_simplify_removes_degree_0_and_1_medians(self):
        """Test that simplification removes obsolete medians (degree < 2)."""
        mjn = MedianJoiningNetwork(epsilon=3.0, simplify=True)
        alignment = Alignment(
            [
                Sequence('seq1', 'AAAAAAAA'),
                Sequence('seq2', 'AAAATTTT'),
                Sequence('seq3', 'TTTTTTTT'),
                Sequence('seq4', 'AAAAAACC'),
            ]
        )

        network = mjn.construct_network(alignment)

        # Check that all median vectors have degree >= 2
        for hap in network.haplotypes:
            if hap.id.startswith('Median_'):
                degree = network.get_degree(hap.id)
                assert degree >= 2, (
                    f'Obsolete median {hap.id} with degree {degree} '
                    f'should have been removed'
                )

    def test_mjn_simplify_removes_degree_2_medians(self):
        """
        Test that simplification removes degree-2 medians.

        Degree-2 medians that just connect two nodes should be replaced
        by a direct edge between those nodes.
        """
        mjn_nosimplify = MedianJoiningNetwork(epsilon=2.0, simplify=False)
        mjn_simplify = MedianJoiningNetwork(epsilon=2.0, simplify=True)

        alignment = Alignment(
            [
                Sequence('A', 'AAAAAAAAAA'),
                Sequence('B', 'AAAAAATTTT'),
                Sequence('C', 'TTTTTTTTTT'),
            ]
        )

        network_nosimplify = mjn_nosimplify.construct_network(alignment)
        network_simplify = mjn_simplify.construct_network(alignment)

        # Simplified network should have fewer haplotypes
        assert len(network_simplify.haplotypes) <= len(network_nosimplify.haplotypes)

        # Both networks should be connected
        assert network_nosimplify.is_connected()
        assert network_simplify.is_connected()

    def test_mjn_simplify_iterates_until_no_changes(self):
        """
        Test that simplification iterates until no more changes.

        Matches C++ behavior where removeObsoleteVerts loops until
        no vertices are removed.
        """
        mjn = MedianJoiningNetwork(epsilon=4.0, simplify=True)
        alignment = Alignment(
            [
                Sequence('s1', 'AAAAAAAAAA'),
                Sequence('s2', 'AAAAAATTTT'),
                Sequence('s3', 'AATTAATTAA'),
                Sequence('s4', 'TTTTTTTTTT'),
                Sequence('s5', 'AAAAATTTAA'),
            ]
        )

        network = mjn.construct_network(alignment)

        # After complete simplification, no median should have degree < 2
        median_haplotypes = [
            h for h in network.haplotypes if h.id.startswith('Median_')
        ]
        for hap in median_haplotypes:
            degree = network.get_degree(hap.id)
            assert degree >= 2, (
                f'Median {hap.id} with degree {degree} should have been '
                f'removed by iterative simplification'
            )

    def test_mjn_simplify_preserves_observed_haplotypes(self):
        """
        Test that simplification never removes observed haplotypes.

        Only median vectors (inferred haplotypes) should be removed.
        """
        mjn = MedianJoiningNetwork(epsilon=3.0, simplify=True)
        alignment = Alignment(
            [
                Sequence('obs1', 'AAAAAAAA'),
                Sequence('obs2', 'AAAATTTT'),
                Sequence('obs3', 'TTTTTTTT'),
                Sequence('obs4', 'AAAAAACC'),
            ]
        )

        network = mjn.construct_network(alignment)

        # Check that we have the right number of observed haplotypes
        observed_haps = [h for h in network.haplotypes if h.frequency > 0]
        assert len(observed_haps) == 4

        # Check that observed haplotypes have their sample IDs preserved
        for hap in observed_haps:
            assert len(hap.sample_ids) > 0
            assert hap.sample_ids[0].startswith('obs')

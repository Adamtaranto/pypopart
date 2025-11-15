"""
Test cases for MJN epsilon > 0 bug fix.

This module tests the fix for the bug where MJN with epsilon > 0
would fail with KeyError: "Haplotype 'Median_0' not found".
"""

from pypopart.algorithms.mjn import MedianJoiningNetwork
from pypopart.core.alignment import Alignment
from pypopart.core.sequence import Sequence


class TestMJNEpsilonFix:
    """Test cases for MJN epsilon > 0 bug fix."""

    def test_mjn_epsilon_zero(self):
        """Test MJN with epsilon=0 (baseline)."""
        mjn = MedianJoiningNetwork(epsilon=0.0)
        alignment = Alignment(
            [
                Sequence('seq1', 'AAAA'),
                Sequence('seq2', 'ATAT'),
                Sequence('seq3', 'TATA'),
                Sequence('seq4', 'TTTT'),
            ]
        )
        network = mjn.construct_network(alignment)

        assert network.is_connected()
        assert len(network.haplotypes) >= 4

    def test_mjn_epsilon_one(self):
        """Test MJN with epsilon=1.0 (previously failed with KeyError)."""
        mjn = MedianJoiningNetwork(epsilon=1.0)
        alignment = Alignment(
            [
                Sequence('seq1', 'AAAA'),
                Sequence('seq2', 'ATAT'),
                Sequence('seq3', 'TATA'),
                Sequence('seq4', 'TTTT'),
            ]
        )
        # Should not raise KeyError
        network = mjn.construct_network(alignment)

        assert network.is_connected()
        assert len(network.haplotypes) >= 4

    def test_mjn_epsilon_large(self):
        """Test MJN with large epsilon (more medians expected)."""
        mjn = MedianJoiningNetwork(epsilon=2.0, simplify=False)
        alignment = Alignment(
            [
                Sequence('sample_A', 'AAAAAAAAAA'),
                Sequence('sample_B', 'AAAAAATTTT'),
                Sequence('sample_C', 'AATTAATTAA'),
                Sequence('sample_D', 'TTTTTTTTTT'),
                Sequence('sample_E', 'AAAAATTTAA'),
                Sequence('sample_F', 'TTAATTAATT'),
            ]
        )
        network = mjn.construct_network(alignment)

        # Network should be connected
        assert network.is_connected()

        # Count observed vs inferred haplotypes
        observed = [h for h in network.haplotypes if h.frequency > 0]
        inferred = [h for h in network.haplotypes if h.frequency == 0]

        # Should have 6 observed haplotypes
        assert len(observed) == 6

        # May have some inferred median haplotypes
        assert len(inferred) >= 0

        # Check that observed haplotypes have correct sample IDs
        for hap in observed:
            assert len(hap.sample_ids) > 0
            assert hap.sample_ids[0].startswith('sample_')

        # Check that inferred haplotypes have no sample IDs
        for hap in inferred:
            assert len(hap.sample_ids) == 0
            assert 'Median_' in hap.id

    def test_mjn_median_haplotype_ids_preserved(self):
        """Test that median haplotype IDs are preserved across iterations."""
        mjn = MedianJoiningNetwork(epsilon=2.0, simplify=False)
        alignment = Alignment(
            [
                Sequence('s1', 'AAAAAAAAAA'),
                Sequence('s2', 'AAAAAATTTT'),
                Sequence('s3', 'AATTAATTAA'),
                Sequence('s4', 'TTTTTTTTTT'),
                Sequence('s5', 'AAAAATTTAA'),
                Sequence('s6', 'TTAATTAATT'),
            ]
        )
        network = mjn.construct_network(alignment)

        # Get all haplotype IDs
        hap_ids = [h.id for h in network.haplotypes]

        # Check that median haplotypes have correct ID format
        median_ids = [hid for hid in hap_ids if 'Median_' in hid]
        for mid in median_ids:
            # Should start with "Median_" followed by a number
            assert mid.startswith('Median_')
            # Should be able to extract the number
            num_part = mid.replace('Median_', '')
            assert num_part.isdigit()

    def test_mjn_sample_ids_not_confused_with_haplotype_ids(self):
        """Test that sample IDs are not confused with haplotype IDs."""
        mjn = MedianJoiningNetwork(epsilon=1.0)
        alignment = Alignment(
            [
                Sequence('sample_001', 'AAAA'),
                Sequence('sample_002', 'ATAT'),
                Sequence('sample_003', 'TATA'),
                Sequence('sample_004', 'TTTT'),
            ]
        )
        network = mjn.construct_network(alignment)

        # Check that haplotype IDs follow H pattern
        hap_ids = [h.id for h in network.haplotypes]
        for hid in hap_ids:
            # Should be either HX or Median_X
            assert hid.startswith('H') or hid.startswith('Median_')

        # Check that sample IDs are preserved correctly
        for hap in network.haplotypes:
            if hap.frequency > 0:  # Observed haplotype
                # Sample IDs should be the original sequence IDs
                for sid in hap.sample_ids:
                    assert sid.startswith('sample_')
            else:  # Inferred median
                # Should have no sample IDs
                assert len(hap.sample_ids) == 0

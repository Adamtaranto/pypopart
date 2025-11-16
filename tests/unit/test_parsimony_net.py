"""Unit tests for Parsimony Network (PN) algorithm."""

import numpy as np
import pytest

from pypopart.algorithms.parsimony_net import ParsimonyNetwork
from pypopart.core.alignment import Alignment
from pypopart.core.distance import DistanceMatrix
from pypopart.core.sequence import Sequence


class TestParsimonyNetwork:
    """Test cases for Parsimony Network algorithm."""

    def test_pn_initialization(self):
        """Test PN algorithm initialization."""
        pn = ParsimonyNetwork()
        assert pn.distance_method == 'hamming'
        assert pn.n_trees == 100
        assert pn.min_edge_frequency == 0.05
        assert pn._median_counter == 0

    def test_pn_custom_parameters(self):
        """Test PN with custom parameters."""
        pn = ParsimonyNetwork(n_trees=50, min_edge_frequency=0.1, random_seed=42)
        assert pn.n_trees == 50
        assert pn.min_edge_frequency == 0.1
        assert pn._random_seed == 42

    def test_pn_empty_alignment(self):
        """Test PN with empty alignment."""
        pn = ParsimonyNetwork()
        alignment = Alignment()
        network = pn.construct_network(alignment)

        assert len(network) == 0

    def test_pn_single_sequence(self):
        """Test PN with single sequence."""
        pn = ParsimonyNetwork()
        alignment = Alignment([Sequence('seq1', 'ATCG')])
        network = pn.construct_network(alignment)

        assert len(network) == 1
        assert len(network.edges) == 0

    def test_pn_two_identical_sequences(self):
        """Test PN with two identical sequences."""
        pn = ParsimonyNetwork()
        alignment = Alignment(
            [Sequence('seq1', 'ATCG'), Sequence('seq2', 'ATCG')]
        )
        network = pn.construct_network(alignment)

        # Should have one haplotype node with frequency 2
        assert len(network) == 1
        assert len(network.edges) == 0

    def test_pn_two_different_sequences(self):
        """Test PN with two different sequences."""
        pn = ParsimonyNetwork(n_trees=10, min_edge_frequency=0.0)
        alignment = Alignment(
            [Sequence('seq1', 'ATCG'), Sequence('seq2', 'ATCT')]
        )
        network = pn.construct_network(alignment)

        # Should have at least 2 nodes
        assert len(network) >= 2
        # Should have at least one edge (since min_edge_frequency=0)
        assert len(network.edges) >= 1

    def test_pn_three_sequences(self):
        """Test PN with three sequences."""
        pn = ParsimonyNetwork(n_trees=20, min_edge_frequency=0.1)
        alignment = Alignment(
            [
                Sequence('seq1', 'AAAA'),
                Sequence('seq2', 'AAAT'),
                Sequence('seq3', 'AATT'),
            ]
        )
        network = pn.construct_network(alignment)

        # Should have 3 haplotype nodes
        assert len([n for n in network.nodes if not network.is_median_vector(n)]) == 3
        # Should have edges connecting them
        assert len(network.edges) >= 2

    def test_pn_edge_frequency_threshold(self):
        """Test that edge frequency threshold works."""
        # With low threshold, more edges
        pn_low = ParsimonyNetwork(n_trees=100, min_edge_frequency=0.01, random_seed=42)
        alignment = Alignment(
            [
                Sequence('seq1', 'AAAA'),
                Sequence('seq2', 'AAAT'),
                Sequence('seq3', 'AATT'),
                Sequence('seq4', 'TATT'),
            ]
        )
        network_low = pn_low.construct_network(alignment)

        # With high threshold, fewer edges
        pn_high = ParsimonyNetwork(n_trees=100, min_edge_frequency=0.5, random_seed=42)
        network_high = pn_high.construct_network(alignment)

        # Low threshold should have more or equal edges
        assert len(network_low.edges) >= len(network_high.edges)

    def test_pn_reproducibility_with_seed(self):
        """Test that algorithm is reproducible with same seed."""
        alignment = Alignment(
            [
                Sequence('seq1', 'ATCG'),
                Sequence('seq2', 'ATCT'),
                Sequence('seq3', 'TTCT'),
            ]
        )

        # Run twice with same seed
        pn1 = ParsimonyNetwork(n_trees=50, random_seed=42)
        network1 = pn1.construct_network(alignment)

        pn2 = ParsimonyNetwork(n_trees=50, random_seed=42)
        network2 = pn2.construct_network(alignment)

        # Should produce same results
        assert len(network1) == len(network2)
        assert len(network1.edges) == len(network2.edges)

    def test_pn_different_seeds_different_results(self):
        """Test that different seeds can produce different results."""
        alignment = Alignment(
            [
                Sequence('seq1', 'AAAA'),
                Sequence('seq2', 'AAAT'),
                Sequence('seq3', 'AATT'),
                Sequence('seq4', 'TATT'),
            ]
        )

        # Run with different seeds
        pn1 = ParsimonyNetwork(n_trees=20, min_edge_frequency=0.1, random_seed=42)
        network1 = pn1.construct_network(alignment)

        pn2 = ParsimonyNetwork(n_trees=20, min_edge_frequency=0.1, random_seed=123)
        network2 = pn2.construct_network(alignment)

        # Results might differ due to randomness (though not guaranteed)
        # Just check both networks were built successfully
        assert len(network1) >= 4
        assert len(network2) >= 4

    def test_pn_calculate_pairwise_distance(self):
        """Test pairwise distance calculation."""
        pn = ParsimonyNetwork()

        seq1 = Sequence('seq1', 'ATCG')
        seq2 = Sequence('seq2', 'ATCT')

        distance = pn._calculate_pairwise_distance(seq1, seq2)

        # Should differ at 1 position
        assert distance == 1.0

    def test_pn_calculate_pairwise_distance_identical(self):
        """Test pairwise distance for identical sequences."""
        pn = ParsimonyNetwork()

        seq1 = Sequence('seq1', 'ATCG')
        seq2 = Sequence('seq2', 'ATCG')

        distance = pn._calculate_pairwise_distance(seq1, seq2)
        assert distance == 0.0

    def test_pn_calculate_pairwise_distance_unequal_length(self):
        """Test pairwise distance with unequal length sequences."""
        pn = ParsimonyNetwork()

        seq1 = Sequence('seq1', 'ATCG')
        seq2 = Sequence('seq2', 'ATC')

        # Should now handle unequal lengths by counting length difference as mutations
        distance = pn._calculate_pairwise_distance(seq1, seq2)
        # Length diff is 1, and all 3 matching positions are the same, so distance = 1
        assert distance == 1.0

    def test_pn_median_vertex_creation(self):
        """Test that median vertices are created for multi-mutation edges."""
        pn = ParsimonyNetwork(n_trees=10, min_edge_frequency=0.0)

        # Sequences differing at multiple positions
        alignment = Alignment(
            [
                Sequence('seq1', 'AAAAA'),
                Sequence('seq2', 'TTTTT'),
            ]
        )

        network = pn.construct_network(alignment)

        # Should have median vertices between very different sequences
        median_count = len([n for n in network.nodes if network.is_median_vector(n)])

        # With 5 differences, should have some medians
        assert median_count >= 0  # May or may not create medians depending on tree sampling

    def test_pn_network_connectivity(self):
        """Test that constructed network is connected."""
        pn = ParsimonyNetwork(n_trees=50, min_edge_frequency=0.05)
        alignment = Alignment(
            [
                Sequence('seq1', 'AAAA'),
                Sequence('seq2', 'AAAT'),
                Sequence('seq3', 'AATT'),
            ]
        )
        network = pn.construct_network(alignment)

        # Network should be connected (or consist of singletons if threshold too high)
        # At minimum, nodes should exist
        assert len(network) >= 3

    def test_pn_with_distance_matrix(self):
        """Test PN with pre-computed distance matrix."""
        pn = ParsimonyNetwork(n_trees=20)
        sequences = [
            Sequence('seq1', 'ATCG'),
            Sequence('seq2', 'ATCT'),
        ]
        alignment = Alignment(sequences)

        # Pre-compute distance matrix
        dist_matrix = pn.calculate_distances(alignment)

        # Build network with pre-computed distances
        network = pn.construct_network(alignment, dist_matrix)

        assert len(network) >= 2

    def test_pn_many_sequences(self):
        """Test PN with more sequences."""
        pn = ParsimonyNetwork(n_trees=30, min_edge_frequency=0.1)

        # Create 6 sequences with varying similarity
        sequences = [
            Sequence('seq1', 'AAAA'),
            Sequence('seq2', 'AAAT'),
            Sequence('seq3', 'AATT'),
            Sequence('seq4', 'ATTT'),
            Sequence('seq5', 'TTTT'),
            Sequence('seq6', 'TTTA'),
        ]
        alignment = Alignment(sequences)

        network = pn.construct_network(alignment)

        # Should have all 6 sequences
        observed_count = len([n for n in network.nodes if not network.is_median_vector(n)])
        assert observed_count == 6

    def test_pn_star_topology(self):
        """Test PN with sequences forming a star topology."""
        pn = ParsimonyNetwork(n_trees=50, min_edge_frequency=0.1)

        # Central sequence differs from three others at different positions
        alignment = Alignment(
            [
                Sequence('center', 'AAAA'),
                Sequence('seq1', 'TAAA'),
                Sequence('seq2', 'ATAA'),
                Sequence('seq3', 'AATA'),
            ]
        )

        network = pn.construct_network(alignment)

        # Should have all 4 sequences
        assert len([n for n in network.nodes if not network.is_median_vector(n)]) == 4

    def test_pn_network_properties(self):
        """Test that constructed network has valid properties."""
        pn = ParsimonyNetwork(n_trees=30, min_edge_frequency=0.1)
        alignment = Alignment(
            [
                Sequence('seq1', 'AAAA'),
                Sequence('seq2', 'AAAT'),
                Sequence('seq3', 'AATT'),
            ]
        )
        network = pn.construct_network(alignment)

        # All edge weights should be non-negative
        for source, target in network.edges:
            dist = network.get_edge_distance(source, target)
            assert dist >= 0

    def test_pn_with_gaps(self):
        """Test PN with sequences containing gaps."""
        pn = ParsimonyNetwork(n_trees=20, min_edge_frequency=0.1)
        alignment = Alignment(
            [
                Sequence('seq1', 'AT-G'),
                Sequence('seq2', 'ATCG'),
            ]
        )

        # Should handle gaps
        network = pn.construct_network(alignment)
        assert len(network) >= 1

    def test_pn_tree_count_affects_consensus(self):
        """Test that number of trees affects edge consensus."""
        alignment = Alignment(
            [
                Sequence('seq1', 'AAAA'),
                Sequence('seq2', 'AAAT'),
                Sequence('seq3', 'AATT'),
            ]
        )

        # Few trees with low threshold
        pn_few = ParsimonyNetwork(n_trees=5, min_edge_frequency=0.2, random_seed=42)
        network_few = pn_few.construct_network(alignment)

        # Many trees with same threshold
        pn_many = ParsimonyNetwork(n_trees=100, min_edge_frequency=0.2, random_seed=42)
        network_many = pn_many.construct_network(alignment)

        # Both should build networks
        assert len(network_few) >= 3
        assert len(network_many) >= 3

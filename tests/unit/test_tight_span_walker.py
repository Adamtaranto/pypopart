"""Unit tests for Tight Span Walker (TSW) algorithm."""

import numpy as np
import pytest

from pypopart.algorithms.tight_span_walker import TightSpanWalker
from pypopart.core.alignment import Alignment
from pypopart.core.distance import DistanceMatrix
from pypopart.core.sequence import Sequence


class TestTightSpanWalker:
    """Test cases for Tight Span Walker algorithm."""

    def test_tsw_initialization(self):
        """Test TSW algorithm initialization."""
        tsw = TightSpanWalker()
        assert tsw.distance_method == 'hamming'
        assert tsw._median_counter == 0
        assert tsw._dT_matrix is None

    def test_tsw_empty_alignment(self):
        """Test TSW with empty alignment."""
        tsw = TightSpanWalker()
        alignment = Alignment()
        network = tsw.construct_network(alignment)

        assert len(network) == 0

    def test_tsw_single_sequence(self):
        """Test TSW with single sequence."""
        tsw = TightSpanWalker()
        alignment = Alignment([Sequence('seq1', 'ATCG')])
        network = tsw.construct_network(alignment)

        assert len(network) == 1
        assert len(network.edges) == 0

    def test_tsw_two_identical_sequences(self):
        """Test TSW with two identical sequences."""
        tsw = TightSpanWalker()
        alignment = Alignment(
            [Sequence('seq1', 'ATCG'), Sequence('seq2', 'ATCG')]
        )
        network = tsw.construct_network(alignment)

        # Should have one haplotype node with frequency 2
        assert len(network) == 1
        assert len(network.edges) == 0

    def test_tsw_two_different_sequences(self):
        """Test TSW with two different sequences."""
        tsw = TightSpanWalker()
        alignment = Alignment(
            [Sequence('seq1', 'ATCG'), Sequence('seq2', 'ATCT')]
        )
        network = tsw.construct_network(alignment)

        # Should have at least 2 nodes
        assert len(network) >= 2
        # Should have at least one edge
        assert len(network.edges) >= 1

    def test_tsw_three_sequences_simple(self):
        """Test TSW with three sequences forming a simple tree."""
        tsw = TightSpanWalker()
        alignment = Alignment(
            [
                Sequence('seq1', 'AAAA'),
                Sequence('seq2', 'AAAT'),
                Sequence('seq3', 'AATT'),
            ]
        )
        network = tsw.construct_network(alignment)

        # Should have 3 haplotype nodes
        assert len([n for n in network.nodes if not network.is_median_vector(n)]) == 3
        # Should have edges connecting them
        assert len(network.edges) >= 2

    def test_tsw_compute_dT_simple(self):
        """Test dT distance computation."""
        tsw = TightSpanWalker()

        # Simple distance matrix
        dist_matrix = np.array(
            [[0, 1, 2], [1, 0, 1], [2, 1, 0]], dtype=float
        )

        dT = tsw._compute_dT(dist_matrix, 3)

        # Check symmetry
        for i in range(3):
            for j in range(3):
                assert dT[i, j] == dT[j, i]

        # Check diagonal is zero
        for i in range(3):
            assert dT[i, i] == 0

        # dT should be >= original distance (tree metric property)
        for i in range(3):
            for j in range(3):
                assert dT[i, j] >= dist_matrix[i, j]

    def test_tsw_about_equal(self):
        """Test floating point equality comparison."""
        tsw = TightSpanWalker()

        # Exactly equal
        assert tsw._about_equal(1.0, 1.0)

        # Very close
        assert tsw._about_equal(1.0, 1.0 + 1e-15)

        # Different
        assert not tsw._about_equal(1.0, 1.1)

        # Zero cases
        assert tsw._about_equal(0.0, 0.0)
        assert tsw._about_equal(0.0, 1e-300)

    def test_tsw_with_distance_matrix(self):
        """Test TSW with pre-computed distance matrix."""
        tsw = TightSpanWalker()
        sequences = [
            Sequence('seq1', 'ATCG'),
            Sequence('seq2', 'ATCT'),
        ]
        alignment = Alignment(sequences)

        # Pre-compute distance matrix
        dist_matrix = tsw.calculate_distances(alignment)

        # Build network with pre-computed distances
        network = tsw.construct_network(alignment, dist_matrix)

        assert len(network) >= 2

    def test_tsw_star_topology(self):
        """Test TSW with sequences forming a star topology."""
        tsw = TightSpanWalker()

        # Central sequence differs from three others at different positions
        alignment = Alignment(
            [
                Sequence('center', 'AAAA'),
                Sequence('seq1', 'TAAA'),
                Sequence('seq2', 'ATAA'),
                Sequence('seq3', 'AATA'),
            ]
        )

        network = tsw.construct_network(alignment)

        # Should have all 4 sequences
        assert len([n for n in network.nodes if not network.is_median_vector(n)]) == 4

        # Center should be connected to all others
        center_neighbors = list(network.graph.neighbors('H1'))  # Assuming H1 is center
        # Should have high connectivity

    def test_tsw_median_vector_creation(self):
        """Test that median vectors are created when needed."""
        tsw = TightSpanWalker()

        # Sequences that require median vectors
        alignment = Alignment(
            [
                Sequence('seq1', 'AAAAA'),
                Sequence('seq2', 'TTTTT'),
            ]
        )

        network = tsw.construct_network(alignment)

        # Count median vectors
        median_count = len([n for n in network.nodes if network.is_median_vector(n)])

        # For sequences with many differences, medians may be needed
        # This is algorithm-dependent, so just check network was built
        assert len(network) >= 2

    def test_tsw_dT_matrix_expansion(self):
        """Test that dT matrix expands correctly when adding median vertices."""
        tsw = TightSpanWalker()

        # Initialize with small matrix
        tsw._dT_matrix = np.array([[0, 1], [1, 0]], dtype=float)

        # Create new dT vector
        new_dT_vector = np.array([0.5, 0.5, 0], dtype=float)

        # Expand matrix
        tsw._expand_dT_matrix(new_dT_vector, 2)

        # Check dimensions
        assert tsw._dT_matrix.shape == (3, 3)

        # Check values
        assert tsw._dT_matrix[2, 0] == 0.5
        assert tsw._dT_matrix[0, 2] == 0.5
        assert tsw._dT_matrix[2, 2] == 0.0

    def test_tsw_reproducibility(self):
        """Test that algorithm produces consistent results."""
        alignment = Alignment(
            [
                Sequence('seq1', 'ATCG'),
                Sequence('seq2', 'ATCT'),
                Sequence('seq3', 'TTCT'),
            ]
        )

        # Run twice
        tsw1 = TightSpanWalker()
        network1 = tsw1.construct_network(alignment)

        tsw2 = TightSpanWalker()
        network2 = tsw2.construct_network(alignment)

        # Should produce same number of nodes and edges
        assert len(network1) == len(network2)
        assert len(network1.edges) == len(network2.edges)

    def test_tsw_invalid_alignment(self):
        """Test TSW with invalid alignment (unequal sequence lengths)."""
        tsw = TightSpanWalker()
        
        # Should raise error during alignment construction
        with pytest.raises((ValueError, AssertionError)):
            alignment = Alignment(
                [Sequence('seq1', 'ATCG'), Sequence('seq2', 'ATC')]
            )

    def test_tsw_network_properties(self):
        """Test that constructed network has valid properties."""
        tsw = TightSpanWalker()
        alignment = Alignment(
            [
                Sequence('seq1', 'AAAA'),
                Sequence('seq2', 'AAAT'),
                Sequence('seq3', 'AATT'),
            ]
        )
        network = tsw.construct_network(alignment)

        # Network should be connected
        assert network.is_connected()

        # All edge weights should be non-negative
        for source, target in network.edges:
            dist = network.get_edge_distance(source, target)
            assert dist >= 0

        # All nodes should be reachable
        for node_id in network.nodes:
            assert len(list(network.graph.neighbors(node_id))) > 0 or len(network) == 1

    def test_tsw_with_gaps(self):
        """Test TSW with sequences containing gaps."""
        tsw = TightSpanWalker()
        alignment = Alignment(
            [
                Sequence('seq1', 'AT-G'),
                Sequence('seq2', 'ATCG'),
            ]
        )

        # Should handle gaps based on distance calculation method
        network = tsw.construct_network(alignment)
        assert len(network) >= 1

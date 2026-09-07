"""Tests for figure labelling: node labels and mutation counts."""

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pytest

from pypopart.core.graph import HaplotypeNetwork
from pypopart.core.haplotype import Haplotype
from pypopart.core.sequence import Sequence
from pypopart.visualization.static_plot import StaticNetworkPlotter


@pytest.fixture
def network():
    """Build a two-haplotype network four mutations apart."""
    net = HaplotypeNetwork(name='LabelTest')
    net.add_haplotype(Haplotype(Sequence('Alpha_01', 'ATCG'), sample_ids=['a']))
    net.add_haplotype(Haplotype(Sequence('Beta_01', 'TTCC'), sample_ids=['b']))
    net.add_edge('Alpha_01', 'Beta_01', distance=4)
    return net


def _texts(ax):
    """Collect the rendered text of every label on the axes."""
    return {t.get_text() for t in ax.texts}


class TestNodeLabels:
    """Nodes carry the same labels the interactive view shows."""

    def test_defaults_to_node_ids(self, network):
        """Without a mapping the node's own ID is drawn."""
        _, ax = StaticNetworkPlotter(network).plot(show_mutations=False)
        try:
            assert {'Alpha_01', 'Beta_01'} <= _texts(ax)
        finally:
            plt.close('all')

    def test_uses_supplied_labels(self, network):
        """The GUI passes H numbers so the figure matches the screen."""
        _, ax = StaticNetworkPlotter(network).plot(
            node_labels={'Alpha_01': 'H1', 'Beta_01': 'H2'},
            show_mutations=False,
        )
        try:
            drawn = _texts(ax)
            assert {'H1', 'H2'} <= drawn
            assert 'Alpha_01' not in drawn
        finally:
            plt.close('all')


class TestEdgeMutationCounts:
    """Edges report the real number of mutations."""

    def test_numeral_shows_the_distance(self, network):
        """'weight' defaults to 1.0 for every edge; 'distance' is the count."""
        _, ax = StaticNetworkPlotter(network).plot(
            show_labels=False, show_edge_ticks=False
        )
        try:
            assert '4' in _texts(ax)
            # The old code read weight and labelled every edge '1'.
            assert '1' not in _texts(ax)
        finally:
            plt.close('all')

    def test_ticks_replace_the_numeral(self, network):
        """Below the threshold the count is drawn as strokes, not text."""
        _, ax = StaticNetworkPlotter(network).plot(
            show_labels=False, show_edge_ticks=True, edge_tick_threshold=10
        )
        try:
            assert '4' not in _texts(ax)
            # One line per mutation, on top of the edge itself.
            assert len(ax.lines) == 4
        finally:
            plt.close('all')

    def test_long_edges_fall_back_to_a_numeral(self, network):
        """A comb of 40 strokes is unreadable, so past the threshold: text."""
        network.graph['Alpha_01']['Beta_01']['distance'] = 40
        _, ax = StaticNetworkPlotter(network).plot(
            show_labels=False, show_edge_ticks=True, edge_tick_threshold=10
        )
        try:
            assert '40' in _texts(ax)
            assert not ax.lines
        finally:
            plt.close('all')

    def test_zero_distance_is_not_marked(self, network):
        """Identical haplotypes need no mutation marker."""
        network.graph['Alpha_01']['Beta_01']['distance'] = 0
        _, ax = StaticNetworkPlotter(network).plot(show_labels=False)
        try:
            assert not ax.lines
            assert not _texts(ax) - {'LabelTest'}
        finally:
            plt.close('all')

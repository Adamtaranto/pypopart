"""Tests for the shared palette and population colour assignment."""

from pathlib import Path

import pytest

from pypopart.visualization.style import (
    PALETTE,
    POP_ART_PALETTE,
    generate_population_colors,
)

THEME_CSS = (
    Path(__file__).resolve().parents[2]
    / 'src'
    / 'pypopart'
    / 'gui'
    / 'assets'
    / 'theme.css'
)


class TestPaletteDoesNotDrift:
    """CSS cannot import the Python constants, so pin them together."""

    def test_theme_css_ships(self):
        """The stylesheet has to exist for Dash to auto-serve it."""
        assert THEME_CSS.is_file()

    @pytest.mark.parametrize('color', PALETTE)
    def test_palette_colour_is_in_the_stylesheet(self, color):
        """Every palette colour is declared as a CSS custom property."""
        assert color in THEME_CSS.read_text().lower()

    @pytest.mark.parametrize('color', POP_ART_PALETTE)
    def test_population_colours_are_distinct(self, color):
        """Two populations sharing a colour would be unreadable."""
        assert POP_ART_PALETTE.count(color) == 1


class TestGeneratePopulationColors:
    """Assigning colours to populations."""

    def test_uses_the_palette_while_it_lasts(self):
        """Small studies get the hand-picked colours."""
        colors = generate_population_colors(['PopB', 'PopA'])

        # Sorted, so PopA takes the first palette entry.
        assert colors['PopA'] == POP_ART_PALETTE[0]
        assert colors['PopB'] == POP_ART_PALETTE[1]

    def test_exactly_palette_length_still_uses_palette(self):
        """The boundary case must not tip into the HSV fallback."""
        pops = [f'Pop{i:02d}' for i in range(len(POP_ART_PALETTE))]
        colors = generate_population_colors(pops)

        assert set(colors.values()) == set(POP_ART_PALETTE)

    def test_falls_back_beyond_the_palette(self):
        """Cycling would repeat a colour, so switch to spaced hues."""
        pops = [f'Pop{i:02d}' for i in range(len(POP_ART_PALETTE) + 1)]
        colors = generate_population_colors(pops)

        assert len(set(colors.values())) == len(pops)

    @pytest.mark.parametrize('n', [1, 2, 5, 8, 9, 20, 50])
    def test_colours_are_always_unique(self, n):
        """The guarantee that matters, at every size."""
        colors = generate_population_colors([f'Pop{i:02d}' for i in range(n)])

        assert len(set(colors.values())) == n
        assert all(c.startswith('#') for c in colors.values())

    def test_empty_input(self):
        """No populations, no colours."""
        assert generate_population_colors([]) == {}

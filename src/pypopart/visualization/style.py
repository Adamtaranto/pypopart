"""
Shared node styling for the static, interactive, and Cytoscape plotters.

Single source of truth for the app's palette, for the colour-priority rule
(median colour > explicit node colour > dominant-population colour >
default) and for generating distinct population colours.

The palette constants here are mirrored in ``gui/assets/theme.css``, which
cannot import Python. ``tests/unit/test_style.py`` pins the two together so
they cannot drift apart silently.
"""

import colorsys
from typing import Dict, List, Optional

#: Bauhaus-flavoured pop art palette used across the app chrome, the
#: network canvas and the static exports.
POP_RED = '#d62828'
POP_NAVY = '#003049'
POP_AMBER = '#fcbf49'
POP_BONE = '#eae2b7'
POP_INK = '#1d1d1b'
POP_PAPER = '#ffffff'

#: Every palette colour, in the order the CSS declares them.
PALETTE = (POP_RED, POP_NAVY, POP_AMBER, POP_BONE, POP_INK, POP_PAPER)

#: Colours assigned to populations, in order. Chosen to stay legible as
#: solid node fills and as pie slices, so the near-white bone and the
#: near-black ink of the base palette are not in here.
POP_ART_PALETTE = (
    POP_RED,
    POP_NAVY,
    POP_AMBER,
    '#2a9d8f',  # teal
    '#8338ec',  # violet
    '#f77f00',  # orange
    '#06a77d',  # green
    '#e5989b',  # rose
)

#: Default node colour used by every plotter.
DEFAULT_NODE_COLOR = '#4d8fac'

#: Default colour for inferred median/intermediate vertices.
DEFAULT_MEDIAN_COLOR = '#c9c9c4'


def node_color(
    node: str,
    haplotype,
    is_median: bool,
    node_color_map: Optional[Dict[str, str]] = None,
    population_colors: Optional[Dict[str, str]] = None,
    median_vector_color: str = DEFAULT_MEDIAN_COLOR,
    default: str = DEFAULT_NODE_COLOR,
) -> str:
    """
    Resolve a node's colour with the shared priority rule.

    Parameters
    ----------
    node : str
        Node ID.
    haplotype : Haplotype or None
        The node's haplotype, if any.
    is_median : bool
        Whether the node is an inferred median vector.
    node_color_map : dict, optional
        Explicit per-node colour overrides.
    population_colors : dict, optional
        Colour per population name.
    median_vector_color : str, default=DEFAULT_MEDIAN_COLOR
        Colour for median vectors.
    default : str, default=DEFAULT_NODE_COLOR
        Fallback colour.

    Returns
    -------
    str
        The resolved colour.
    """
    if is_median:
        return median_vector_color
    if node_color_map and node in node_color_map:
        return node_color_map[node]
    if population_colors and haplotype is not None:
        pop_counts = haplotype.get_frequency_by_population()
        if pop_counts:
            dominant = max(pop_counts.items(), key=lambda item: item[1])[0]
            return population_colors.get(dominant, default)
    return default


def generate_population_colors(populations: List[str]) -> Dict[str, str]:
    """
    Assign a distinct colour to each population.

    Uses the hand-picked pop art palette while it lasts and falls back to
    evenly spaced HSV hues beyond it. The palette is not cycled: repeating
    it would hand two populations the same colour, and distinguishable
    colours matter more than staying on-palette.

    Parameters
    ----------
    populations : list of str
        Population names.

    Returns
    -------
    dict
        Mapping of population name to hex colour. Colours are unique for
        any number of populations.
    """
    names = sorted(populations)
    n = len(names)

    if n <= len(POP_ART_PALETTE):
        return {pop: POP_ART_PALETTE[i] for i, pop in enumerate(names)}

    colors = {}
    for i, pop in enumerate(names):
        hue = i / n if n else 0.0
        r, g, b = colorsys.hsv_to_rgb(hue, 0.7, 0.9)
        colors[pop] = '#{:02x}{:02x}{:02x}'.format(
            int(r * 255), int(g * 255), int(b * 255)
        )
    return colors


def apply_pop_art_rcparams() -> Dict[str, object]:
    """
    Apply the app's palette and typography to matplotlib.

    Called by the static plotter so exported figures carry the same look
    as the on-screen network. Fonts fall back to whatever the system has
    when the bundled families are not installed, which matplotlib handles
    on its own.

    Returns
    -------
    dict
        The rcParams that were set, so callers can restore them.
    """
    import matplotlib as mpl

    params = {
        'figure.facecolor': POP_PAPER,
        'axes.facecolor': POP_PAPER,
        'axes.edgecolor': POP_INK,
        'axes.labelcolor': POP_INK,
        'axes.titlecolor': POP_INK,
        'axes.titleweight': 'bold',
        'text.color': POP_INK,
        'xtick.color': POP_INK,
        'ytick.color': POP_INK,
        'font.family': ['Space Grotesk', 'DejaVu Sans', 'sans-serif'],
        'legend.frameon': True,
        'legend.edgecolor': POP_INK,
        'legend.facecolor': POP_PAPER,
    }
    mpl.rcParams.update(params)
    return params

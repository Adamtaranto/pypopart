"""
Shared node styling for the static, interactive, and Cytoscape plotters.

Single source of truth for the colour-priority rule (median colour >
explicit node colour > dominant-population colour > default) and for
generating distinct population colours.
"""

import colorsys
from typing import Dict, List, Optional

#: Default node colour used by every plotter.
DEFAULT_NODE_COLOR = 'lightblue'

#: Default colour for inferred median/intermediate vertices.
DEFAULT_MEDIAN_COLOR = 'lightgray'


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
    median_vector_color : str, default='lightgray'
        Colour for median vectors.
    default : str, default='lightblue'
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
    Generate distinct colours for populations using HSV colour space.

    Parameters
    ----------
    populations : list of str
        Population names.

    Returns
    -------
    dict
        Mapping of population name to hex colour, evenly spaced in hue.
    """
    n = len(populations)
    colors = {}
    for i, pop in enumerate(sorted(populations)):
        hue = i / n if n else 0.0
        r, g, b = colorsys.hsv_to_rgb(hue, 0.7, 0.9)
        colors[pop] = '#{:02x}{:02x}{:02x}'.format(
            int(r * 255), int(g * 255), int(b * 255)
        )
    return colors

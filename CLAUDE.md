# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

PyPopART is a pure Python port of PopART (Population Analysis with Reticulate Trees) for constructing and visualizing haplotype networks from DNA sequence data. License is GPL-3.0-or-later.

## Commands

```bash
# Install for development (editable, with dev tools)
pip install -e ".[dev]"
pre-commit install

# Run all tests (coverage flags come from pyproject.toml)
pytest

# Run a single test file / test
pytest tests/unit/test_mjn.py
pytest tests/unit/test_mjn.py::test_name

# Lint & format (ruff is the source of truth; single quotes, numpy docstrings)
ruff check --fix .
ruff format .

# Run all pre-commit hooks (ruff, numpydoc-validation, yamlfmt, prettier, etc.)
pre-commit run --all-files

# Docs (mkdocs-material; docs deps in the "docs" extra)
mkdocs serve

# Entry points
pypopart --help      # CLI (click group in src/pypopart/cli/main.py)
pypopart-gui         # Dash web GUI at http://localhost:8050 (src/pypopart/gui/app.py)
```

Versioning is hatch-vcs from git tags (`v*`); `src/pypopart/_version.py` is generated — never edit it.

## Architecture

The package (`src/pypopart/`, src layout) is a pipeline: **io → core → algorithms → stats/layout/visualization**, with the CLI and GUI as thin front-ends over the same API.

- **`core/`** — domain model. `Alignment` (alignment.py) wraps sequences; `distance.py` provides `DistanceMatrix` and `calculate_pairwise_distances` (hamming, jc, k2p, tamura_nei; numba-accelerated versions in `distance_optimized.py`); `haplotype.py`/`condensation.py` collapse identical sequences into unique haplotypes with frequency maps; `graph.py` defines `HaplotypeNetwork`, the networkx-based network object everything downstream consumes.
- **`algorithms/`** — one module per network construction method (mst, msn, tcs, mjn, parsimony_net, tsw). All inherit from `NetworkAlgorithm` in `base.py` and implement `construct_network(alignment, distance_matrix) -> HaplotypeNetwork` (`build_network` is a backward-compat alias used by CLI/GUI). New algorithms follow this pattern and get registered in `algorithms/__init__.py`.
- **`io/`** — one module per format (fasta, nexus, phylip, genbank), plus `metadata.py` (CSV population/geo metadata) and `network_export.py` (GraphML/GML/JSON output).
- **`stats/`** — network statistics, topology analysis, and population genetics measures operating on `HaplotypeNetwork`.
- **`layout/`** — node layout algorithms shared by both visualizers.
- **`visualization/`** — `static_plot.py` (matplotlib), `interactive_plot.py` (plotly), `cytoscape_plot.py` (Dash Cytoscape, used by the GUI).
- **`cli/main.py`** — click command group with subcommands `load`, `network`, `analyze`, `visualize`, `info`.
- **`gui/app.py`** — Dash application (upload → configure algorithm → compute → layout → export).

## Conventions

- Ruff enforces single quotes, numpy-style docstrings, and isort-style imports; numpydoc-validation runs in pre-commit (tests/docs excluded), so public functions need compliant numpy docstrings.
- Tests live flat in `tests/unit/`, named `test_*.py` (pre-commit enforces test-first naming).

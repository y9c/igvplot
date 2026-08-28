"""Central visual theme for igvplot figures.

Everything that renders shares one cohesive, modern visual language:
a single palette for the tracks and a set of matplotlib ``rcParams`` applied
when a figure is built via :func:`apply_theme`. Tweak colours here and every
track (and the README gallery) updates consistently.
"""
from __future__ import annotations

import matplotlib as mpl

# --------------------------------------------------------------------------- #
# Palette (modern, colour-blind friendly where possible)
# --------------------------------------------------------------------------- #
# Track fills / lines.
COVERAGE = "#3d7fc4"          # aligned coverage area
COVERAGE_MISMATCH = "#e0492f"  # "variant pile" ticks
STRAND_FORWARD = "#cdd8e8"     # forward-strand reads (pale steel)
STRAND_REVERSE = "#708098"     # reverse-strand reads (slate)
STRAND_UNKNOWN = "#9aa7b8"
INSERTION = "#2c3e50"          # '+' insertion marker
DELETION = "#e0492f"           # red deletion connector
SITE = "#6c5ce7"               # vertical site markers
SASHIMI = "#e63946"            # splice-junction arcs
TAD = "#e63946"                # TAD boundary triangles
SCALE = "#333333"
BED = "#8d99ae"
HIGH_BG = "#ffe08a"
REF_BASE = "#888888"
READ_LINK = "#8d99ae"

# Nucleotide base letters (mismatch / sequence row).
BASE_COLORS = {
    "A": "#2a9d8f",
    "C": "#2f86eb",
    "G": "#f39c12",
    "T": "#e0492f",
    "N": "#9e9e9e",
}

# Base-modification (m6A / m5C / …) colours, keyed on substring match.
BASE_MOD_COLORS = {
    "m6a": "#e63946",
    "m6am": "#e63946",
    "m5c": "#457b9d",
    "5mc": "#457b9d",
    "5mct": "#457b9d",
    "5hmc": "#2a9d8f",
    "m1a": "#9b5de5",
    "m7g": "#f4a261",
    "ac4c": "#e76f51",
    "psi": "#7f8c8d",
    "pseudo": "#7f8c8d",
    "nm": "#00b4d8",
    "default": "#9e9e9e",
}


def apply_theme() -> None:
    """Apply the modern matplotlib defaults used by every rendered figure.

    Called at the start of :meth:`igvplot.GenomeView.render` so a single
    coherent style is used whether you call ``plot_view``, the builder API, or
    the CLI. Idempotent.
    """
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["DejaVu Sans", "Helvetica Neue", "Helvetica", "Arial"],
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
            "axes.edgecolor": "#9aa5b1",
            "axes.labelcolor": "#333333",
            "axes.titlecolor": "#333333",
            "text.color": "#333333",
            "xtick.color": "#666666",
            "ytick.color": "#666666",
            "xtick.labelsize": 9,
            "ytick.labelsize": 8,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "legend.frameon": False,
        }
    )

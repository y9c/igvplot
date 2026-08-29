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
COVERAGE = "#4a7db5"           # aligned coverage area (gradient base colour)
COVERAGE_MISMATCH = "#e5484d"  # "variant pile" ticks
STRAND_FORWARD = "#ccd9ec"     # forward-strand reads (pale steel blue)
STRAND_REVERSE = "#6b84a3"     # reverse-strand reads (slate blue)
STRAND_UNKNOWN = "#9aa7b8"
INSERTION = "#334155"          # '+' insertion marker
DELETION = "#e5484d"           # red deletion connector
SITE = "#94a3b8"               # vertical site markers (subtle slate)
SITE_TEXT = "#64748b"          # site label text
SASHIMI = "#e5484d"            # splice-junction arcs
TAD = "#e5484d"                # TAD boundary triangles
SCALE = "#3d4451"
BED = "#94a3b8"
HIGH_BG = "#ffe08a"
REF_BASE = "#9aa5b1"
READ_LINK = "#8d99ae"
# Gene/feature track (dna_features_viewer bars are recoloured to this).
GENE_FACE = "#5e6fa3"
GENE_EDGE = "#ffffff"
# Chrome: axis lines and secondary text.
SPINE = "#dbe1e8"
TEXT = "#3d4451"
TEXT_SOFT = "#57606a"

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

# Modern muted categorical palette (read-group / group colouring) — replaces
# the saturated matplotlib tableau defaults.
CATEGORICAL = [
    "#4a7db5",  # steel blue
    "#e5484d",  # soft red
    "#2a9d8f",  # teal
    "#f39c12",  # amber
    "#8e7cc3",  # muted violet
    "#e78ac3",  # rose
    "#66a61e",  # green
    "#a0693a",  # sienna
    "#00b4d8",  # cyan
    "#7f8c8d",  # grey
]

# Default feature-type colours for the gene/feature track (dna_features_viewer).
# Keys are lowercased feature types.
FEATURE_COLORS = {
    "gene": "#5e6fa3",      # indigo slate
    "mrna": "#2a9d8f",      # teal
    "exon": "#e07a5f",      # terracotta / coral
    "cds": "#3d5a80",       # deep steel blue
    "5'utr": "#f4a261",     # warm amber
    "3'utr": "#e9c46a",     # soft gold
    "utr": "#f4a261",       # warm amber
    "promoter": "#81b29a",  # sage green
    "enhancer": "#e76f51",  # burnt orange
    "intron": "#9aa5b1",    # cool grey
    "repeat_region": "#b8b8b8",
    "variation": "#e5484d", # soft red
    "misc_feature": "#8d99ae",
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
            "axes.edgecolor": SPINE,
            "axes.linewidth": 0.9,
            "axes.labelcolor": TEXT_SOFT,
            "axes.titlecolor": TEXT,
            "text.color": TEXT,
            "xtick.color": "#8a94a6",
            "ytick.color": "#8a94a6",
            "xtick.labelsize": 9,
            "ytick.labelsize": 8,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "legend.frameon": False,
        }
    )

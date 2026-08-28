"""igvplot: IGV-style genomic visualization in matplotlib.

Stack aligned reads, per-base coverage, gene/feature annotation, splice
junctions, base modifications, Hi-C/TAD and other signal tracks on shared
x-axes, with read-level mismatch / insertion / deletion display.

Two ways to build a figure:

1. One-call convenience::

       view = igvplot.plot_view(
           bam_path="sample.bam", region="chr1:1,000-2,000",
           features="annotation.gb", reference="genome.fa",
           sites={1050: "m6A"}, sashimi=True, out_path="locus.png",
       )

2. Fluent builder (``IGV``/``GenomeView``/``AlignmentView`` are the same class)::

       from igvplot import IGV
       igv = IGV("chr1:1,000-2,000", reference="genome.fa")
       igv.add_reads("rna.bam", color_by="readgroup").add_coverage("rna.bam")
       igv.add_sashimi("rna.bam").add_features("annotation.gb")
       igv.add_base_mods({1050: (1, "m6A")}).add_sites({1050: "m6A"})
       igv.savefig("locus.png")

Track adders include reads (overlay), coverage (strand/overlay), sashimi
(overlay/junctions-BED), features, sequence, sites, base mods, variants, GC,
signal, interaction arcs, Hi-C/TAD/BED/scale and fully custom tracks
(``add_track``). ``igvplot.summary`` returns region QC stats.
"""
from .region import Region
from .reads import (
    BASE_COLORS,
    Read,
    Reference,
    compute_coverage,
    compute_insert_sizes,
    fetch_reads,
    junction_counts,
    open_reference,
    variant_allele_fraction,
)
from .features import load_features
from .plot import build_legend_items, set_font_size
from .theme import apply_theme
from .bigwig import (
    BigWigUnavailableError,
    coverage_from_bedgraph,
    read_bigwig_coverage,
)
from .view import (
    AlignmentView,
    GenomeView,
    IGV,
    insert_size_histogram,
    plot_view,
    summary,
)
from . import plot

__version__ = "0.1.0"

__all__ = [
    "Region",
    "Read",
    "Reference",
    "BASE_COLORS",
    "GenomeView",
    "IGV",
    "AlignmentView",
    "plot_view",
    "insert_size_histogram",
    "summary",
    "fetch_reads",
    "compute_coverage",
    "compute_insert_sizes",
    "junction_counts",
    "variant_allele_fraction",
    "load_features",
    "build_legend_items",
    "set_font_size",
    "apply_theme",
    "open_reference",
    "read_bigwig_coverage",
    "coverage_from_bedgraph",
    "BigWigUnavailableError",
    "plot",
    "__version__",
]

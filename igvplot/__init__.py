"""igvplot: IGV-style genomic visualization in matplotlib.

Stack aligned reads, per-base coverage and gene/feature annotation tracks
(drawn by ``dna_features_viewer``) on shared x-axes, with read-level mismatch,
insertion and deletion display.

Typical usage::

    import igvplot
    view = igvplot.plot_view(
        bam_path="sample.bam",
        region="chr1:1,000-2,000",
        features="annotation.gb",
        reference="genome.fa",
        sites={1050: "m6A site", 1200: "SNP"},
        out_path="locus.png",
    )
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
from .bigwig import (
    BigWigUnavailableError,
    coverage_from_bedgraph,
    read_bigwig_coverage,
)
from .view import AlignmentView, GenomeView, IGV, insert_size_histogram, plot_view
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
    "fetch_reads",
    "compute_coverage",
    "compute_insert_sizes",
    "junction_counts",
    "variant_allele_fraction",
    "load_features",
    "build_legend_items",
    "set_font_size",
    "open_reference",
    "read_bigwig_coverage",
    "coverage_from_bedgraph",
    "BigWigUnavailableError",
    "plot",
    "__version__",
]

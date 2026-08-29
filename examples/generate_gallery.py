#!/usr/bin/env python3
"""Generate the example images used in README.md.

Each call below renders one figure from the synthetic data in ``data/`` and
writes it to ``examples/``. Run from the repo root::

    python examples/generate_gallery.py

The output images are referenced by ``README.md`` as a visual gallery.
"""
from __future__ import annotations

import os
import sys

import matplotlib

matplotlib.use("Agg")  # headless

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, REPO)

import numpy as np  # noqa: E402

import igvplot  # noqa: E402
from igvplot import GenomeView, Region  # noqa: E402

DATA = os.path.join(REPO, "data")
BAM = os.path.join(DATA, "sample.bam")
REF = os.path.join(DATA, "genome.fa")
GB = os.path.join(DATA, "annotation.gb")
SITES = os.path.join(DATA, "sites.bed")
BASEMOD = os.path.join(DATA, "basemods.bed")
VCF = os.path.join(DATA, "variants.vcf")

# m6A GLORI synthetic data (two genes on opposite strands)
M6A_BAM = os.path.join(DATA, "m6a", "m6a_reads.bam")
M6A_REF = os.path.join(DATA, "m6a", "m6a_genome.fa")
M6A_GB = os.path.join(DATA, "m6a", "m6a_annotation.gb")

# m6A sites: forward gene (A, +) and reverse gene (T, -)
M6A_SITES = {
    250: "m6A", 350: "m6A", 550: "m6A", 650: "m6A",   # forward gene (A→G)
    800: "m6A", 850: "m6A", 1020: "m6A", 1070: "m6A",  # reverse gene (T→C)
}

HERO_REGION = "chrTest:6,930-7,200"

# show structural gene features only (drops the tiny point 'variation' features,
# which are still annotatable via the dedicated site-marker track)
STRUCTURAL = {"gene", "CDS", "mRNA", "exon", "transcript", "protein", "repeat_region"}

HIGH = "#ffe0a3"
HIGH2 = "#a8d8c8"


def save(view, name, dpi=120):
    out = os.path.join(HERE, name)
    view.savefig(out, dpi=dpi)
    print(f"saved examples/{name}")


def hero():
    """The default view: sashimi + feature + coverage + reads + sequence + sites."""
    view = igvplot.plot_view(
        bam_path=BAM,
        region=HERO_REGION,
        features=GB,
        reference=REF,
        sites={7000: "C>T", 7010: "A>G", 7020: "Δ3bp", 7030: "+CC"},
        sashimi=True,
        link_mates=True,
        view_as_pairs=True,
        color_by="strand",
        show_sequence=True,
        min_feature_length=3,
        figsize=(15, 9),
    )
    save(view, "gallery_hero.png", dpi=120)
def base_level():
    """Zoomed base-resolution: every read base + colour-coded reference row."""
    view = igvplot.plot_view(
        bam_path=BAM,
        region="chrTest:6,995-7,030",
        reference=REF,
        sites={7000: "C>T", 7010: "A>G", 7020: "Δ3bp"},
        show_all_bases=True,
        show_sequence=True,
        color_by="strand",
        max_reads=70,
        figsize=(15, 7),
    )
    save(view, "gallery_base_level.png", dpi=130)


def color_pair():
    view = igvplot.plot_view(
        bam_path=BAM,
        region=HERO_REGION,
        reference=REF,
        features=GB,
        color_by="pairOrientation",
        group_by="strand",
        link_mates=True,
        view_as_pairs=True,
        show_soft_clips=True,
        min_feature_length=3,
        figsize=(13, 6),
    )
    save(view, "gallery_color_pair.png", dpi=120)


def color_readgroup():
    view = igvplot.plot_view(
        bam_path=BAM,
        region=HERO_REGION,
        reference=REF,
        features=GB,
        color_by="readgroup",
        group_by="readgroup",
        min_feature_length=3,
        figsize=(13, 6),
    )
    save(view, "gallery_color_readgroup.png", dpi=120)


def color_mapq():
    view = igvplot.plot_view(
        bam_path=BAM,
        region=HERO_REGION,
        reference=REF,
        color_by="mapq",
        colormap="Blues",
        group_by="strand",
        figsize=(13, 6),
    )
    save(view, "gallery_color_mapq.png", dpi=120)


def basemod():
    """m6A GLORI overview: strand-coloured reads + per-base A→G / T→C conversion.

    GLORI is a *negative* method: the reagent converts every unmodified A→G (and
    its reverse-strand complement T→C), but an m6A protects its own base so it
    stays A (or T on the complement). The per-base conversion track therefore
    dips to ~0% exactly at the m6A sites — the places that keep the reference
    base. Reads are coloured by strand and the converted bases are shown red.
    """
    view = (
        GenomeView(
            region="chrM6A:150-1150",
            reference=M6A_REF,
            figsize=(17, 9),
            title="m6A GLORI — A→G / T→C conversion (negative method)",
        )
        .add_conversion_fraction(M6A_BAM, reference=M6A_REF, label="A→G / T→C conversion")
        .add_reads(M6A_BAM, reference=M6A_REF, color_by="strand", group_by="strand",
                   max_reads=120, show_all_bases=False, base_fontsize=7.5, weight=6.5,
                   link_mates=True, view_as_pairs=True, sample_seed=7,
                   mismatch_colors={"A>G": "#e63946", "T>C": "#e63946"})  # red = converted
        .add_features(M6A_GB, min_feature_length=3)
        .add_sequence(M6A_REF)
        .add_sites(M6A_SITES)
    )
    save(view, "gallery_basemod.png", dpi=130)


def basemod_zoom():
    """Base-resolution zoom spanning the forward-gene intron.

    Every read base is drawn, so the A→G conversions (red) are readable
    letter-by-letter while the m6A sites (chrM6A:350, 550) keep their A and are
    left grey — the protected signal. Reads that span the intron (CIGAR ``N``)
    are drawn as separate exon blocks so the splice junction is visible.
    """
    view = (
        GenomeView(
            region="chrM6A:350-550",
            reference=M6A_REF,
            figsize=(15, 6),
            title="Zoom — A→G conversions & splice junction (chrM6A:350-550)",
        )
        .add_reads(M6A_BAM, reference=M6A_REF, color_by="strand", group_by="strand",
                   max_reads=45, show_all_bases=True, base_fontsize=8.0, weight=5.5,
                   link_mates=True, view_as_pairs=True, sample_seed=7,
                   mismatch_colors={"A>G": "#e63946", "T>C": "#e63946"})
        .add_features(M6A_GB, min_feature_length=3)
        .add_sequence(M6A_REF)
        .add_sites({350: "m6A", 550: "m6A"})
    )
    save(view, "gallery_basemod_zoom.png", dpi=130)


def overlay():
    """Multi-sample coverage overlay + shaded highlight regions."""
    view = GenomeView(region=HERO_REGION, reference=REF, figsize=(13, 6))
    view.add_coverage_overlay(
        [(BAM, "#4C86C6", "replicate A"), (BAM, "#E8883A", "replicate B")]
    )
    view.add_highlight_regions([(6960, 6990, HIGH, 0.35), (7040, 7070, HIGH2, 0.35)])
    view.add_reads(BAM, reference=REF, color_by="strand", max_reads=80)
    view.add_features(GB, min_feature_length=3)
    save(view, "gallery_overlay.png", dpi=120)


def hic():
    """pyGenomeTracks-style multi-track: scale + Hi-C + TAD + genes + reads."""
    region = Region.from_any(HERO_REGION)

    def make_hic(region, n_bins=48, tads=((6950, 7060), (7060, 7160), (7160, 7230))):
        edges = np.linspace(region.start, region.end, n_bins + 1)
        centers = (edges[:-1] + edges[1:]) / 2

        def tad_of(x):
            for i, (a, b) in enumerate(tads):
                if a <= x < b:
                    return i
            return 0

        mat = np.zeros((n_bins, n_bins))
        rng = np.random.default_rng(3)
        for i in range(n_bins):
            for j in range(n_bins):
                dist = abs(i - j)
                decay = np.exp(-dist / 6.0)
                same = 1.0 if tad_of(centers[i]) == tad_of(centers[j]) else 0.12
                mat[i, j] = decay * (8.0 * same) + rng.uniform(0, 0.4)
        return mat

    gv = GenomeView(region=HERO_REGION, figsize=(15, 9), dpi=120)
    gv.add_scale_bar(window_bp=100)
    gv.add_hic(make_hic(region), cmap="Reds")
    gv.add_tads([7060, 7160])
    gv.add_features(GB, min_feature_length=3)
    gv.add_coverage(BAM, reference=REF)
    gv.add_reads(BAM, reference=REF, color_by="readgroup", group_by="pairOrientation")
    save(gv, "gallery_hic.png", dpi=120)


def compare():
    """Multi-sample comparison: strand coverage + sashimi overlay + reads overlay."""
    view = (
        GenomeView(region=HERO_REGION, reference=REF, figsize=(15, 9))
        .add_coverage_strands(BAM, weight=1.6)
        .add_sashimi_overlay([(BAM, "#e63946", "isoform1"), (BAM, "#457b9d", "isoform2")], weight=1.6)
        .add_reads_overlay(
            [(BAM, "#4c86c6", "replicate A"), (BAM, "#e8883a", "replicate B")],
            max_reads=50,
            weight=3.0,
        )
        .add_features(GB, min_feature_length=3)
    )
    save(view, "gallery_compare.png", dpi=120)


def new_tracks():
    """Signal, GC-content, interaction arcs, variants and a custom track."""
    import numpy as np
    region = Region.from_any(HERO_REGION)
    n = region.length
    gv = (
        GenomeView(region=HERO_REGION, reference=REF, figsize=(15, 10))
        .add_gc(window=60, weight=0.9)
        .add_signal((np.sin(np.linspace(0, 6, n)) * 0.5 + 0.5) * 9, ylabel="motif score", color="#9b5de5", weight=1.0)
        .add_arc([(6930, 7060, 3.0), (6970, 7140, 2.0), (7060, 7200, 1.0)], label="loop", weight=1.8)
        .add_variants(VCF)
        .add_coverage(BAM, reference=REF)
        .add_reads(BAM, reference=REF, color_by="readgroup", group_by="pairOrientation", max_reads=40)
        .add_features(GB, min_feature_length=3)
    )
    save(gv, "gallery_tracks.png", dpi=120)


def epi():
    """Per-base m5C stoichiometry with bisulfite C→T mutations."""
    M5C_BAM = os.path.join(REPO, "data", "m5c", "m5c_reads.bam")
    M5C_REF = os.path.join(REPO, "data", "m5c", "m5c_genome.fa")
    M5C_GB = os.path.join(REPO, "data", "m5c", "m5c_annotation.gb")
    gv = (
        GenomeView(region="chrM5C:120-380", reference=M5C_REF, figsize=(15, 9))
        .add_mod_fraction(
            {150: (1, "m5C", 0.90), 250: (1, "m5C", 0.75), 350: (1, "m5C", 0.60)},
            label="m5C %",
            weight=1.0,
        )
        .add_reads(M5C_BAM, reference=M5C_REF, color_by="strand", max_reads=80)
        .add_features(M5C_GB, min_feature_length=3)
        .add_sites({150: "m5C", 250: "m5C", 350: "m5C"})
    )
    save(gv, "gallery_epigenetics.png", dpi=120)


def variants():
    """A single variant-centred plot with the variant allele fraction in the title."""
    vaf, depth, alt = igvplot.variant_allele_fraction(BAM, "chrTest", 7000, reference=REF)
    view = igvplot.plot_view(
        bam_path=BAM,
        region="chrTest:6,930-7,100",
        reference=REF,
        features=GB,
        sites={7000: "C>T"},
        color_by="strand",
        show_sequence=True,
        show_all_bases=False,
        basemod={7000: (1, "m6A"), 7010: (-1, "m5C")},
        min_feature_length=3,
        figsize=(13, 6),
    )
    view.set_title(f"chrTest:7001 C>T    VAF={vaf:.2f} ({alt}/{depth})")
    save(view, "gallery_variants.png", dpi=120)


def main():
    hero()
    base_level()
    color_pair()
    color_readgroup()
    color_mapq()
    basemod()
    basemod_zoom()
    overlay()
    hic()
    compare()
    new_tracks()
    epi()
    variants()
    print("done.")


if __name__ == "__main__":
    main()

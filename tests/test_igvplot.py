import os
import sys

import matplotlib

matplotlib.use("Agg")  # headless rendering for tests

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
DATA = os.path.join(REPO, "data")

import pytest  # noqa: E402

from igvplot import (  # noqa: E402
    GenomeView,
    Region,
    compute_coverage,
    coverage_from_bedgraph,
    fetch_reads,
    junction_counts,
    load_features,
    plot_view,
)


def _data(fname):
    path = os.path.join(DATA, fname)
    if not os.path.exists(path):
        pytest.skip(f"missing synthetic data: {path} (run make_synthetic_bam.py)")
    return path


@pytest.fixture(autouse=True)
def _reset_global_font():
    """Reset the global font scale after every test so state never leaks."""
    yield
    from igvplot import set_font_size

    set_font_size(8)


# --------------------------------------------------------------------------- #
# Region
# --------------------------------------------------------------------------- #
def test_region_parse_string():
    r = Region.from_string("chr1:1,000-2,000")
    assert r.chrom == "chr1"
    assert r.start == 999  # 0-based inclusive
    assert r.end == 2000  # 0-based exclusive
    assert r.length == 1001


def test_region_coerce():
    r = Region.from_any(("chr1", 100, 200))
    assert (r.chrom, r.start, r.end) == ("chr1", 100, 200)
    assert Region.from_any(r) == r


def test_region_invalid():
    with pytest.raises(ValueError):
        Region("chr1", 500, 100)
    with pytest.raises(ValueError):
        Region.from_string("no-tab-chrom")


# --------------------------------------------------------------------------- #
# Read fetching / mutations on synthetic data
# --------------------------------------------------------------------------- #
def test_fetch_reads_mutations_match_synthetic():
    bam = _data("sample.bam")
    ref = _data("genome.fa")
    region = Region("chrTest", 6900, 7250)
    reads = fetch_reads(bam, region, reference=ref)

    from collections import Counter

    mm, ins, dl = Counter(), Counter(), Counter()
    for r in reads:
        for p in r.mismatches:
            mm[p] += 1
        for p in r.insertions:
            ins[p] += 1
        for a, b in r.deletions:
            dl[a] += 1

    assert mm[7000] == 12  # SNP (12 reads)
    assert mm[7010] == 5  # SNP (5 reads)
    assert ins[7029] == 4  # CC insertion anchored at 7029
    assert dl[7020] == 6  # 3 bp deletion starting at 7020


def test_fetch_reads_has_both_strands():
    bam = _data("sample.bam")
    reads = fetch_reads(bam, Region("chrTest", 6900, 7250), reference=_data("genome.fa"))
    assert any(r.is_reverse for r in reads)
    assert any(not r.is_reverse for r in reads)


def test_compute_coverage_and_variant_pile():
    bam = _data("sample.bam")
    region = Region("chrTest", 6900, 7250)
    depths, mmc = compute_coverage(bam, region, reference=_data("genome.fa"))
    assert depths.shape == (region.length,)
    assert depths.min() > 0  # region is fully tiled
    assert mmc[7000 - region.start] == 12
    assert mmc[7010 - region.start] == 5


def test_print_n_skip_not_deletion():
    # RNA-seq N (intron-skip) ref positions must not be reported as deletions.
    region = Region("chrTest", 6900, 7250)
    reads = fetch_reads(_data("sample.bam"), region, reference=_data("genome.fa"))
    dl = sum(1 for r in reads for a, b in r.deletions if a == 7020)
    assert dl == 6  # only the 6 real 3bp-deletion reads


# --------------------------------------------------------------------------- #
# RNA-seq junctions / sashimi
# --------------------------------------------------------------------------- #
def test_junction_counts_and_mates():
    region = Region("chrTest", 6900, 7250)
    reads = fetch_reads(_data("sample.bam"), region, reference=_data("genome.fa"))
    juncs = junction_counts(reads, region)
    assert juncs[(6970, 7060)] == 25
    assert juncs[(7090, 7140)] == 12
    # mate metadata and read groups are populated
    assert any(r.paired and r.mate_start is not None for r in reads)
    groups = {r.read_group for r in reads if r.read_group}
    assert "isoform1" in groups and "isoform2" in groups
    assert any(r.insert_size > 0 for r in reads)


def test_sashimi_track_renders(tmp_path):
    out = tmp_path / "sashimi.png"
    plot_view(
        bam_path=_data("sample.bam"),
        region="chrTest:6,930-7,200",
        reference=_data("genome.fa"),
        sashimi=True,
        link_mates=True,
        out_path=str(out),
        dpi=80,
        figsize=(10, 6),
    )
    assert out.exists() and out.stat().st_size > 0


# --------------------------------------------------------------------------- #
# Colour modes
# --------------------------------------------------------------------------- #
def test_color_by_modes_render(tmp_path):
    bam = _data("sample.bam")
    region = "chrTest:6,930-7,200"
    for mode, cmap in [("mapq", "plasma"), ("readgroup", "viridis"), ("proper", "viridis"), ("mate", "viridis"), ("none", "viridis")]:
        out = tmp_path / f"color_{mode}.png"
        plot_view(
            bam_path=bam, region=region, color_by=mode, colormap=cmap,
            show_coverage=False, out_path=str(out), dpi=60, figsize=(8, 4),
        )
        assert out.exists() and out.stat().st_size > 0, f"color_by={mode} failed"


def test_color_by_invalid_raises(tmp_path):
    import pytest
    view = GenomeView(region="chrTest:6,930-7,200")
    view.add_reads(_data("sample.bam"), reference=_data("genome.fa"), color_by="bogus")
    with pytest.raises(ValueError):
        view.render()


def test_downsampling_reduces_reads():
    region = Region("chrTest", 0, 20000)
    full = fetch_reads(_data("sample.bam"), region, reference=_data("genome.fa"))
    sampled = fetch_reads(
        _data("sample.bam"), region, reference=_data("genome.fa"),
        sampling_window=2000, max_per_window=10,
    )
    assert len(sampled) < len(full) // 2
    assert len(sampled) > 0


def test_legend_and_highlight_and_mode_render(tmp_path):
    import pytest
    from igvplot import build_legend_items
    region = "chrTest:6,930-7,200"
    bam = _data("sample.bam")
    # legend items present for categorical modes
    reads = fetch_reads(bam, Region.from_any(region), reference=_data("genome.fa"))
    for mode in ("strand", "readgroup", "pairOrientation"):
        assert build_legend_items(mode, reads), f"no legend for {mode}"
    # squished + highlight render
    out = tmp_path / "sq.png"
    plot_view(
        bam_path=bam, region=region, display_mode="squished",
        highlight={"read_00000"}, show_legend=True,
        out_path=str(out), dpi=60, figsize=(8, 4),
    )
    assert out.exists() and out.stat().st_size > 0


def test_legend_none_for_continuous():
    from igvplot import build_legend_items
    reads = fetch_reads(_data("sample.bam"), Region("chrTest", 6900, 7250),
                        reference=_data("genome.fa"))
    assert build_legend_items("mapq", reads) == []


def test_show_all_bases_collects_and_renders(tmp_path):
    # zoomed (small) region -> base-level view collects per-read bases
    region = Region("chrTest", 6990, 7010)
    reads = fetch_reads(_data("sample.bam"), region, reference=_data("genome.fa"),
                        collect_bases=True)
    assert reads and any(r.bases for r in reads)  # full bases collected
    out = tmp_path / "base.png"
    plot_view(bam_path=_data("sample.bam"), region="chrTest:6,990-7,010",
              reference=_data("genome.fa"), out_path=str(out),
              dpi=80, figsize=(9, 6))
    assert out.exists() and out.stat().st_size > 0


def test_full_mode_sort_and_indel_text(tmp_path):
    out = tmp_path / "full.png"
    plot_view(
        bam_path=_data("sample.bam"), region="chrTest:6,990-7,060",
        reference=_data("genome.fa"), display_mode="full", sort_by="mapq",
        show_insertion_text=True, show_deletion_text=True,
        out_path=str(out), dpi=80, figsize=(9, 6),
    )
    assert out.exists() and out.stat().st_size > 0


def test_coverage_overlay_and_highlight_regions(tmp_path):
    gv = GenomeView(region="chrTest:6,950-7,160", figsize=(10, 6))
    gv.add_coverage_overlay(
        [(_data("sample.bam"), "#1f77b4", "A"), (_data("sample.bam"), "#d62728", "B")]
    )
    gv.add_reads(_data("sample.bam"), reference=_data("genome.fa"))
    gv.add_highlight_regions([(6995, 7005), (7015, 7025, "#ef9a9a", 0.4)])
    out = tmp_path / "overlay.png"
    gv.savefig(str(out), dpi=80)
    assert out.exists() and out.stat().st_size > 0


def test_base_mod_track_and_coloring(tmp_path):
    from igvplot import GenomeView
    gv = GenomeView(region="chrTest:6,990-7,050", figsize=(10, 6))
    gv.add_base_mods({7000: (1, "m6A", "#c0392b"), 7010: (-1, "m5C", "#2c6fbb")})
    gv.add_reads(_data("sample.bam"), reference=_data("genome.fa"), color_by="basemod")
    out = tmp_path / "basemod.png"
    gv.savefig(str(out), dpi=80)
    assert out.exists() and out.stat().st_size > 0
    # mutation-site parsing from tuples and 4-tuples
    mods = GenomeView._parse_mod_sites([(7000, 1, "m6A"), (7010, -1, "m5C", "#00ff00")])
    assert len(mods) == 2
    assert mods[7010] == (-1, "m5C", "#00ff00")


def test_hic_scale_tad_bed_tracks(tmp_path):
    import numpy as np
    from igvplot import GenomeView
    region = "chrTest:6,930-7,230"
    mat = np.random.default_rng(1).random((20, 20)) * 8
    bed = tmp_path / "feat.bed"
    bed.write_text("chrTest\t7000\t7005\tpeakA\t800\nchrTest\t7060\t7070\tpeakB\t300\n")
    gv = GenomeView(region=region, figsize=(10, 8))
    gv.add_scale_bar(window_bp=100)
    gv.add_hic(mat)
    gv.add_tads([7010, 7110])
    gv.add_bed_features(str(bed))
    out = tmp_path / "hic.png"
    gv.savefig(str(out), dpi=80)
    assert out.exists() and out.stat().st_size > 0


def test_igv_fluent_api(tmp_path):
    from igvplot import IGV, AlignmentView, GenomeView
    out = tmp_path / "igv.png"
    igv = (
        IGV("chrTest:6,930-7,200", reference=_data("genome.fa"), dpi=80)
        .bam(_data("sample.bam"), color_by="readgroup", group_by="pairOrientation")
        .add_features(_data("annotation.gb"))
        .add_sashimi(_data("sample.bam"))
        .add_sites({7000: "SNP", 7010: "m5C"})
        .add_base_mods({7000: (1, "m6A"), 7010: (-1, "m5C")})
        .add_highlight_regions([(6995, 7005)])
        .savefig(str(out))
    )
    assert out.exists() and out.stat().st_size > 0
    # one class, three names
    assert IGV is GenomeView is AlignmentView


# --------------------------------------------------------------------------- #
# BigWig coverage fallback (bedGraph parser)
# --------------------------------------------------------------------------- #
def test_coverage_from_bedgraph():
    region = Region("chrX", 100, 110)
    text = "\n".join(
        [
            "# comment",
            "track name=x",
            "chrX\t100\t105\t7.5",
            "chrX\t107\t110\t3.0",
            "chr9\t100\t110\t99",  # wrong chrom -> ignored
        ]
    )
    depths = coverage_from_bedgraph(text, region, region.length)
    assert (depths[:5] == 7.5).all()
    assert depths[5] == 0.0 and depths[6] == 0.0
    assert (depths[7:] == 3.0).all()


# --------------------------------------------------------------------------- #
# Feature loading (dna_features_viewer integration)
# --------------------------------------------------------------------------- #
def test_load_features_crops_to_region():
    record = load_features(_data("annotation.gb"), region=Region("chrTest", 7000, 7100))
    assert record.first_index == 7000
    assert record.sequence_length == 100


# --------------------------------------------------------------------------- #
# Plotting
# --------------------------------------------------------------------------- #
def test_plot_view_saves_png(tmp_path):
    out = tmp_path / "locus.png"
    plot_view(
        bam_path=_data("sample.bam"),
        region="chrTest:6,950-7,140",
        features=_data("annotation.gb"),
        reference=_data("genome.fa"),
        sites={7000: "SNP C>A"},
        out_path=str(out),
        dpi=80,
        figsize=(10, 6),
    )
    assert out.exists()
    assert out.stat().st_size > 0


def test_genomeview_builds_figure():
    gv = GenomeView(region="chrTest:6,950-7,140", figsize=(8, 5))
    gv.add_coverage(_data("sample.bam"), reference=_data("genome.fa"))
    gv.add_reads(_data("sample.bam"), reference=_data("genome.fa"), max_reads=30)
    gv.add_features(_data("annotation.gb"))
    gv.add_sites({7000: "site"})
    fig = gv.render()
    assert fig is not None


# --------------------------------------------------------------------------- #
# Global font scaling
# --------------------------------------------------------------------------- #
def test_set_font_size_scales(tmp_path):
    from igvplot import set_font_size
    # the default (8) is a no-op; a larger base must still render cleanly
    out = tmp_path / "font12.png"
    set_font_size(12)
    plot_view(
        bam_path=_data("sample.bam"),
        region="chrTest:6,995-7,040",
        reference=_data("genome.fa"),
        sites={7000: "SNP"},
        out_path=str(out),
        dpi=60,
        figsize=(8, 5),
    )
    assert out.exists() and out.stat().st_size > 0


def test_cli_fontsize_flag(tmp_path):
    from igvplot import cli
    out = tmp_path / "cli.png"
    rc = cli.main(
        [
            _data("sample.bam"),
            "chrTest:6,995-7,040",
            "-o", str(out),
            "-r", _data("genome.fa"),
            "--fontSize", "12",
            "--dpi", "60",
            "--figsize", "8", "5",
        ]
    )
    assert rc == 0 and out.exists() and out.stat().st_size > 0


def test_feature_track_min_length_renders(tmp_path):
    # a clean structural gene track hides the 1bp point-variant features
    out = tmp_path / "feat.png"
    plot_view(
        bam_path=_data("sample.bam"),
        region="chrTest:6,940-7,180",
        features=_data("annotation.gb"),
        min_feature_length=3,
        out_path=str(out),
        dpi=60,
        figsize=(9, 5),
    )
    assert out.exists() and out.stat().st_size > 0


def test_apply_theme_is_safe():
    from igvplot import apply_theme
    apply_theme()  # idempotent, must not raise
    apply_theme()

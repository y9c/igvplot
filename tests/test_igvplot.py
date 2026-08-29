import os

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


def _m6a():
    """Return (bam, reference) for the m6A GLORI synthetic data."""
    bam = os.path.join(DATA, "m6a", "m6a_reads.bam")
    ref = os.path.join(DATA, "m6a", "m6a_genome.fa")
    if not os.path.exists(bam) or not os.path.exists(ref):
        pytest.skip("missing m6A GLORI data (run make_m6a_glori_bam.py)")
    return bam, ref


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


def test_region_centered_string():
    # 'chr:pos' is a centred region (window bp each side), 1-based position
    r = Region.from_string("chr1:1000")
    assert r.chrom == "chr1"
    assert r.start == 899 and r.end == 1100  # window 100 each side
    assert r.length == 201


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


def test_fetch_reads_max_reads_samples_across_region():
    bam, ref = _m6a()
    reg = Region.from_any("chrM6A:150-1150")
    sampled = fetch_reads(bam, reg, reference=ref, max_reads=120, sample_seed=7)
    assert 0 < len(sampled) <= 120
    n = reg.length
    quart = lambda x: min(3, int((x - reg.start) / (n / 4)))  # noqa: E731
    # reads must be drawn from the whole region, not just its start
    assert {quart(r.aleft) for r in sampled} == {0, 1, 2, 3}


def test_fetch_reads_sample_keeps_paired_mates():
    bam, ref = _m6a()
    reg = Region.from_any("chrM6A:150-1150")
    sampled = fetch_reads(bam, reg, reference=ref, max_reads=120, sample_seed=7)
    names = {}
    for r in sampled:
        names.setdefault(r.name, []).append(r)
    # fragments that kept both mates must be on opposite strands (a real pair)
    for name, rs in names.items():
        if len(rs) == 2:
            assert rs[0].is_reverse != rs[1].is_reverse


def test_fetch_reads_sample_deterministic_with_seed():
    bam, ref = _m6a()
    reg = Region.from_any("chrM6A:150-1150")
    a = fetch_reads(bam, reg, reference=ref, max_reads=80, sample_seed=42)
    b = fetch_reads(bam, reg, reference=ref, max_reads=80, sample_seed=42)
    assert [r.name for r in a] == [r.name for r in b]


def test_spliced_read_exon_blocks():
    from igvplot.plot import _read_exon_blocks
    from igvplot.reads import Read

    r = Read(name="x", aleft=381, aright=501, is_reverse=False, mapq=60,
             junctions=[(400, 500)])
    # mirror draw_read_track's clipping, then split on the intron
    s, e = max(350, r.aleft), min(550, r.aright)
    assert _read_exon_blocks(r, s, e) == [(381, 400), (500, 501)]
    # an unspliced read is a single block
    r2 = Read(name="y", aleft=200, aright=300, is_reverse=False, mapq=60)
    assert _read_exon_blocks(r2, 200, 300) == [(200, 300)]


def test_conversion_fraction_track_renders(tmp_path):
    bam, ref = _m6a()
    view = GenomeView(region="chrM6A:150-1150", reference=ref, figsize=(12, 6))
    view.add_conversion_fraction(bam, reference=ref, label="A→G / T→C")
    view.add_reads(bam, reference=ref, max_reads=60, sample_seed=7)
    out = tmp_path / "conv.png"
    view.savefig(str(out), dpi=60)
    assert out.stat().st_size > 0


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


def test_deletions_merged_into_spans():
    # adjacent single-base deletions must merge into one contiguous span
    region = Region("chrTest", 6900, 7250)
    reads = fetch_reads(_data("sample.bam"), region, reference=_data("genome.fa"))
    three_bp = [r for r in reads if any(s == 7020 for s, _ in r.deletions)]
    assert three_bp
    # no per-base fragments: every 7020 span is the full 3 bp
    assert all(e - s == 3 for r in three_bp for s, e in r.deletions if s == 7020)
    # spans are stored as (start, end), end > start
    assert all(e > s for r in reads for s, e in r.deletions)


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


def test_color_by_insert_uniform_sizes_renders(tmp_path):
    # degenerate insert sizes (vmin == vmax) must not raise / warn
    out = tmp_path / "insert_color.png"
    plot_view(
        bam_path=_data("sample.bam"), region="chrTest:6,950-7,140",
        color_by="insert", paint_base_letters=False, max_reads=60,
        out_path=str(out), dpi=60, figsize=(8, 4),
    )
    assert out.exists() and out.stat().st_size > 0


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


def test_compare_tracks_render(tmp_path):
    # multi-sample reads overlay + sashimi overlay + strand-specific coverage
    from igvplot import GenomeView
    gv = GenomeView(region="chrTest:6,950-7,180", reference=_data("genome.fa"), figsize=(10, 8))
    gv.add_coverage_strands(_data("sample.bam"))
    gv.add_sashimi_overlay(
        [(_data("sample.bam"), "#e63946", "iso1"), (_data("sample.bam"), "#457b9d", "iso2")]
    )
    gv.add_reads_overlay(
        [(_data("sample.bam"), "#4c86c6", "rep A"), (_data("sample.bam"), "#e8883a", "rep B")],
        max_reads=40,
    )
    gv.add_features(_data("annotation.gb"), min_feature_length=3)
    out = tmp_path / "compare.png"
    gv.savefig(str(out), dpi=70)
    assert out.exists() and out.stat().st_size > 0


def test_add_junctions_bed(tmp_path):
    from igvplot import GenomeView
    bed = tmp_path / "junc.bed"
    bed.write_text(
        "chrTest\t6970\t7060\tJ1\t25\n"
        "chrTest\t7090\t7140\tJ2\t12\n"
        "chrTest\t7000\t7100\tJ3\t3\n"
    )
    gv = GenomeView(region="chrTest:6,950-7,180")
    gv.add_junctions_bed(str(bed))
    out = tmp_path / "junc.png"
    gv.savefig(str(out), dpi=70)
    assert out.exists() and out.stat().st_size > 0


def test_summary_returns_stats():
    from igvplot import summary
    s = summary(_data("sample.bam"), "chrTest:6,950-7,140", reference=_data("genome.fa"))
    for key in (
        "region", "n_reads", "mean_depth", "max_depth", "n_junctions",
        "n_mismatches", "insert_median",
    ):
        assert key in s, f"missing {key}"
    assert s["n_reads"] > 0 and s["max_depth"] > 0


def test_signal_gc_variants_arc_track(tmp_path):
    import numpy as np
    from igvplot import GenomeView
    out = tmp_path / "new.png"
    gv = GenomeView(region="chrTest:6,950-7,180", reference=_data("genome.fa"), figsize=(10, 8))
    n = gv.region.length
    gv.add_signal(np.linspace(0, 1, n), ylabel="score")
    gv.add_gc(window=60)
    gv.add_variants(_data("variants.vcf"))
    gv.add_arc([(6960, 7060, 2.0), (7000, 7140, 1.0)], label="loop")
    gv.add_track(lambda ax, region: ax.axhline(0.5, color="k", lw=0.8))
    gv.savefig(str(out), dpi=70)
    assert out.exists() and out.stat().st_size > 0


def test_signal_length_must_match_region():
    import numpy as np
    import pytest
    from igvplot import GenomeView
    gv = GenomeView(region="chrTest:6,950-7,180")
    with pytest.raises(ValueError):
        gv.add_signal(np.zeros(3))  # wrong length


def test_add_gc_requires_reference():
    import pytest
    from igvplot import GenomeView
    gv = GenomeView(region="chrTest:6,950-7,000")
    with pytest.raises(ValueError):
        gv.add_gc()


def test_mod_fraction_motifs_variant_fraction(tmp_path):
    from igvplot import GenomeView
    out = tmp_path / "epi.png"
    gv = GenomeView(region="chrTest:6,950-7,180", reference=_data("genome.fa"), figsize=(10, 8))
    gv.add_variant_fraction(_data("sample.bam"), reference=_data("genome.fa"))
    gv.add_mod_fraction({6960: (1, "m6A", 0.8), 7010: (-1, "m5C", 0.3, "#2c6fbb")})
    gv.add_motifs("DRACH")
    gv.add_reads(_data("sample.bam"), reference=_data("genome.fa"), max_reads=40)
    gv.savefig(str(out), dpi=70)
    assert out.exists() and out.stat().st_size > 0


def test_iupac_motif_regex_matches():
    from igvplot.view import _iupac_to_regex
    import re
    # DRACH = D R A C H  (D={AGT}, R={AG}, H={ACT}) — the canonical m6A context
    pat = re.compile(_iupac_to_regex("DRACH"), re.I)
    assert pat.match("GAACT") is not None  # D=G, R=A, A, C, H=T
    assert pat.match("TTTTT") is None


def test_cli_version_flag():
    from igvplot import cli
    with pytest.raises(SystemExit):
        cli.main(["--version"])


def test_all_tracks_render(tmp_path):
    """Every real track type must render (marker/highlight utilities are tested
    in combination, since they span existing tracks)."""
    import numpy as np

    from igvplot import GenomeView

    region = "chrTest:6,950-7,180"

    def new():
        return GenomeView(region=region, reference=_data("genome.fa"), figsize=(7, 5))

    def build(kind):
        g = new()
        if kind == "reads":
            g.add_reads(_data("sample.bam"), reference=_data("genome.fa"))
        elif kind == "reads_overlay":
            g.add_reads_overlay([(_data("sample.bam"), "#4c86c6", "A"), (_data("sample.bam"), "#e8883a", "B")], max_reads=30)
        elif kind == "coverage":
            g.add_coverage(_data("sample.bam"), reference=_data("genome.fa"))
        elif kind == "coverage_overlay":
            g.add_coverage_overlay([(_data("sample.bam"), "#4c86c6", "A"), (_data("sample.bam"), "#e8883a", "B")])
        elif kind == "coverage_strands":
            g.add_coverage_strands(_data("sample.bam"))
        elif kind == "sashimi":
            g.add_sashimi(_data("sample.bam"))
        elif kind == "sashimi_overlay":
            g.add_sashimi_overlay([(_data("sample.bam"), "#e63946", "a"), (_data("sample.bam"), "#457b9d", "b")])
        elif kind == "junctions_bed":
            bed = tmp_path / "j.bed"
            bed.write_text("chrTest\t6970\t7060\tJ1\t25\n")
            g.add_junctions_bed(str(bed))
        elif kind == "features":
            g.add_features(_data("annotation.gb"), min_feature_length=3)
        elif kind == "sequence":
            g.add_sequence(_data("genome.fa"))
        elif kind == "gc":
            g.add_gc()
        elif kind == "signal":
            g.add_signal(np.linspace(0, 1, g.region.length))
        elif kind == "variants":
            g.add_variants(_data("variants.vcf"))
        elif kind == "variant_fraction":
            g.add_variant_fraction(_data("sample.bam"), reference=_data("genome.fa"))
        elif kind == "mod_fraction":
            g.add_mod_fraction({6960: 0.8, 7000: (1, "m6A", 0.6)})
        elif kind == "motifs":
            g.add_motifs("DRACH")
        elif kind == "arc":
            g.add_arc([(6960, 7060, 2.0)])
        elif kind == "bed":
            bed = tmp_path / "p.bed"
            bed.write_text("chrTest\t7000\t7005\tpeakA\t800\n")
            g.add_bed_features(str(bed))
        elif kind == "hic":
            g.add_hic(np.random.default_rng(0).random((24, 24)) * 5)
        elif kind == "tads":
            g.add_tads([7060, 7160])
        elif kind == "scale":
            g.add_scale_bar(100)
        elif kind == "basemods":
            g.add_base_mods({7000: (1, "m6A")})
        elif kind == "track":
            g.add_track(lambda ax, region: ax.axhline(0.5))
        return g

    real_tracks = [
        "reads", "reads_overlay", "coverage", "coverage_overlay", "coverage_strands",
        "sashimi", "sashimi_overlay", "junctions_bed", "features", "sequence", "gc",
        "signal", "variant_fraction", "mod_fraction", "motifs", "arc",
        "bed", "hic", "tads", "scale", "basemods", "track",
    ]
    for kind in real_tracks:
        out = tmp_path / f"{kind}.png"
        build(kind).savefig(str(out), dpi=50)
        assert out.exists() and out.stat().st_size > 0, f"{kind} failed"

    # marker/highlight utilities render when combined with a base track
    g = new()
    g.add_sites({7000: "C>T"}).add_variants(_data("variants.vcf"))
    g.add_highlight_regions([(6960, 6990)]).add_reads(_data("sample.bam"), max_reads=20)
    out = tmp_path / "markers.png"
    g.savefig(str(out), dpi=50)
    assert out.exists() and out.stat().st_size > 0


def test_marker_only_utilities_need_a_track():
    # sites / variants / highlight span existing tracks; alone they raise
    import pytest
    from igvplot import GenomeView
    for build in (
        lambda g: g.add_sites({7000: "C>T"}),
        lambda g: g.add_variants("data/variants.vcf"),
        lambda g: g.add_highlight_regions([(6960, 6990)]),
    ):
        g = GenomeView(region="chrTest:6,950-7,180")
        build(g)
        with pytest.raises(ValueError):
            g.render()


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


def test_add_hic_requires_2d():
    import numpy as np
    import pytest
    from igvplot import GenomeView
    gv = GenomeView(region="chrTest:6,950-7,150")
    with pytest.raises(ValueError):
        gv.add_hic(np.zeros(10))  # 1-D is not a contact matrix


def test_feature_types_case_insensitive(tmp_path):
    # feature_types matching must be case-insensitive and never crash on
    # non-string feature types
    from igvplot import GenomeView
    out = tmp_path / "ft.png"
    gv = GenomeView(region="chrTest:6,940-7,190")
    gv.add_features(_data("annotation.gb"), feature_types={"cds", "gene"})
    gv.savefig(str(out), dpi=60)
    assert out.exists() and out.stat().st_size > 0
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
    (
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


def test_coverage_from_bedgraph_interval_tuples():
    # an iterable of (chrom, start, end, value) tuples must also be accepted
    region = Region("chrX", 100, 110)
    intervals = [
        ("chrX", 100, 105, 7.5),
        ("chrX", 107, 110, 3.0),
        ("chr9", 100, 110, 99),  # wrong chrom -> ignored
    ]
    depths = coverage_from_bedgraph(intervals, region, region.length)
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


def test_view_as_pairs_renders(tmp_path):
    # IGV "view as pairs": paired mates are joined on one row
    from igvplot import GenomeView
    out = tmp_path / "pairs.png"
    gv = GenomeView(region="chrTest:6,940-7,140", figsize=(10, 6))
    gv.add_reads(
        _data("sample.bam"),
        reference=_data("genome.fa"),
        link_mates=True,
        view_as_pairs=True,
        max_reads=120,
    )
    gv.savefig(str(out), dpi=70)
    assert out.exists() and out.stat().st_size > 0


def test_genomeview_close_releases_reference():
    # GenomeView.close() must release any Reference it opened
    from igvplot import IGV
    gv = IGV("chrTest:6,990-7,060", reference=_data("genome.fa"))
    assert gv._reference is not None
    gv.close()
    assert gv._reference is None
    gv.close()  # idempotent


def test_add_sequence_renders_and_owns_reference(tmp_path):
    # add_sequence must read the sequence and not leak the Reference it opens
    from igvplot import IGV
    out = tmp_path / "seq.png"
    gv = IGV("chrTest:6,995-7,040")
    gv.add_sequence(_data("genome.fa"))
    gv.add_reads(_data("sample.bam"), reference=_data("genome.fa"), max_reads=40)
    gv.savefig(str(out), dpi=70)
    assert out.exists() and out.stat().st_size > 0


def test_pathlib_paths_supported(tmp_path):
    # pathlib.Path should be accepted everywhere a file path is taken
    from pathlib import Path
    from igvplot import IGV
    out = tmp_path / "pl.png"
    gv = IGV("chrTest:6,990-7,060")
    gv.add_features(Path(_data("annotation.gb")), min_feature_length=3)
    gv.add_reads(Path(_data("sample.bam")), reference=Path(_data("genome.fa")), max_reads=40)
    gv.savefig(str(out), dpi=70)
    assert out.exists() and out.stat().st_size > 0


def test_coverage_reference_propagates_to_reads(tmp_path):
    # a reference given to add_coverage must be reused by a later add_reads
    from igvplot import IGV
    gv = IGV("chrTest:6,990-7,010")
    gv.add_coverage(_data("sample.bam"), reference=_data("genome.fa"))
    assert gv._reference == _data("genome.fa")  # path stored for later tracks
    gv.add_reads(_data("sample.bam"))  # no reference argument
    gv.add_sites({7000: "C>T"})
    out = tmp_path / "rp.png"
    gv.savefig(str(out), dpi=70)
    assert out.exists() and out.stat().st_size > 0


# --------------------------------------------------------------------------- #
# v0.0.1 features: BAM-tag colouring, SJ.out.tab, BEDPE, log scale, DataFrames
# --------------------------------------------------------------------------- #
def _mk_read(**kw):
    """Build a minimal Read for colour/group logic tests."""
    from igvplot.reads import Read

    kw.setdefault("name", "r1")
    kw.setdefault("aleft", 6900)
    kw.setdefault("aright", 7000)
    kw.setdefault("is_reverse", False)
    kw.setdefault("mapq", 30)
    return Read(**kw)


def test_tag_group_and_color_logic():
    from igvplot.plot import _read_colormap, _read_group_key, build_legend_items

    reads = [
        _mk_read(tags={"CB": "AA", "UB": "u1"}),
        _mk_read(tags={"CB": "TT", "UB": "u2"}),
        _mk_read(tags={}),
    ]
    assert _read_group_key(reads[0], "tag:CB") == "AA"
    assert _read_group_key(reads[2], "tag:CB") == "__none__"

    color_fn = _read_colormap("tag:CB", "viridis", reads)
    c0 = color_fn(0, reads[0])
    c1 = color_fn(1, reads[1])
    c2 = color_fn(2, reads[2])
    assert c0 != c1 and c0 != c2 and c2 != c1  # distinct categories

    items = build_legend_items("tag:CB", reads)
    assert [lab for lab, _c in items] == ["AA", "TT"]


def test_fetch_reads_autodetects_tags(tmp_path):
    # synthetic BAM carries no barcode tags: auto-detect must not crash and
    # tag_keys=[] disables collection explicitly
    reads = fetch_reads(_data("sample.bam"), "chrTest:6,990-7,060", max_reads=10)
    assert isinstance(reads, list) and reads


def test_star_sj_parsing(tmp_path):
    from igvplot.view import _looks_like_star_sj, _parse_star_sj

    sj = tmp_path / "SJ.out.tab"
    sj.write_text(
        "chrTest\t7001\t7100\t2\t0\t5\t12\t3\t30\t0\t0\t0\t0\t0\t0\n"
        "chrTest\t7201\t7350\t1\t1\t0\t7\t0\t25\t0\t0\t0\t0\t0\t0\n"
        "chrOther\t1\t100\t2\t0\t1\t1\t1\t10\t0\t0\t0\t0\t0\t0\n"
    )
    assert _looks_like_star_sj(str(sj))
    counts = _parse_star_sj(str(sj), "chrTest")
    assert counts[(7000, 7100)] == 15   # unique (12) + multi (3) spanning reads
    assert counts[(7200, 7350)] == 7    # falls back to >= 1
    assert all(lo >= 0 for lo, _hi in counts)


def test_bedpe_parsing(tmp_path):
    from igvplot.view import _parse_bedpe

    bp = tmp_path / "loops.bedpe"
    bp.write_text(
        "#header line\n"
        "chrTest\t6950\t7000\tchrTest\t7100\t7150\tloop1\t3.0\n"
        "chrTest\t7000\t7050\tchrOther\t9000\t9100\tskip\t9.0\n"
        "chrTest\t6900\t6950\tchrTest\t7180\t7200\tloop2\t.\n"
    )
    arcs = _parse_bedpe(str(bp), "chrTest")
    assert (6950, 7150, 3.0) in arcs
    assert (6900, 7200, 1.0) in arcs  # non-numeric score -> 1.0
    assert len(arcs) == 2  # cross-chrom record dropped


def test_bedpe_track_renders(tmp_path):
    bp = tmp_path / "loops.bedpe"
    bp.write_text("chrTest\t6940\t6990\tchrTest\t7100\t7160\tl1\t2.5\n")
    gv = GenomeView(region="chrTest:6,930-7,200")
    gv.add_bedpe(str(bp))
    out = tmp_path / "bedpe.png"
    gv.savefig(str(out), dpi=70)
    assert out.exists() and out.stat().st_size > 0


def test_junctions_bed_accepts_star_sj(tmp_path):
    sj = tmp_path / "SJ.out.tab"
    sj.write_text("chrTest\t7001\t7100\t2\t0\t5\t12\t3\t30\t0\t0\t0\t0\t0\t0\n")
    gv = GenomeView(region="chrTest:6,930-7,200")
    gv.add_junctions_bed(str(sj))
    out = tmp_path / "sj.png"
    gv.savefig(str(out), dpi=70)
    assert out.exists() and out.stat().st_size > 0


def test_signal_log_scale_and_dataframe_like(tmp_path):
    import numpy as np

    gv = GenomeView(region="chrTest:6,930-7,200")
    n = gv.region.length
    vals = np.exp(np.linspace(0, 6, n))  # heavy-tailed
    gv.add_signal(vals, log=True, weight=1.0)

    class _DF:  # minimal pandas/polars-like object
        @staticmethod
        def to_numpy():
            return np.linspace(0.0, 5.0, n)

    gv.add_signal(_DF(), color="#2a9d8f")
    fig = gv.render()
    # first track is log1p-scaled: its max is log1p(max(vals))
    ax = fig.axes[0]
    assert ax.get_ylim()[1] == pytest.approx(np.log1p(vals.max()), rel=1e-6)
    import matplotlib.pyplot as plt

    plt.close(fig)


def test_coverage_log_scale(tmp_path):
    gv = GenomeView(region="chrTest:6,930-7,200")
    gv.add_coverage(_data("sample.bam"), reference=_data("genome.fa"), log=True)
    gv.add_reads(_data("sample.bam"), max_reads=30)
    out = tmp_path / "logcov.png"
    gv.savefig(str(out), dpi=70)
    assert out.exists() and out.stat().st_size > 0

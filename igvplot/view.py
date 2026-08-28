"""High-level IGV-style view that stacks feature, coverage and read tracks on
shared x-axes and renders them to a matplotlib figure.

The gene/feature track is drawn by ``dna_features_viewer`` on the same
coordinate axis as the pysam-derived reads/coverage, so genomic loci line up
exactly between tracks.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple, Union

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.gridspec import GridSpec
from matplotlib.figure import Figure
from dna_features_viewer import GraphicRecord

from .bigwig import read_bigwig_coverage
from .features import load_features
from .plot import (
    _fs,
    build_legend_items,
    draw_base_mod_track,
    draw_coverage_track,
    draw_interaction_arc,
    draw_read_track,
    draw_sashimi_overlay_track,
    draw_sashimi_track,
    draw_sequence_track,
    draw_sites,
)
from .reads import (
    Reference,
    Read,
    compute_coverage,
    fetch_reads,
    junction_counts,
)
from .region import Region
from .theme import BED, COVERAGE, HIGH_BG, SASHIMI, TAD, apply_theme

__all__ = ["GenomeView", "plot_view"]


def _xy(ax) -> Tuple[float, float]:
    """Return (xmin, xmax) as floats for an axis."""
    x0, x1 = ax.get_xlim()
    return float(x0), float(x1)


def _human_bp(n: float) -> str:
    """Format a base-pair length human-readable ('1.5 kb', '300 bp')."""
    n = float(n)
    if n >= 1_000_000:
        return f"{n / 1_000_000:.1f} Mb"
    if n >= 1000:
        return f"{n / 1000:.1f} kb"
    return f"{int(n)} bp"


def _parse_bed_features(bed_path: str) -> list:
    """Parse BED3+ -> list of (chrom, start, end, name, score)."""

    features = []
    with open(bed_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith(("track", "browser")):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            try:
                s, e = int(parts[1]), int(parts[2])
            except ValueError:
                continue
            name = parts[3] if len(parts) > 3 else ""
            score = None
            if len(parts) > 4:
                try:
                    score = float(parts[4])
                except ValueError:
                    score = None
            features.append((parts[0], s, e, name, score))
    return features


def _parse_variants(source: str, chrom_filter: Optional[str] = None) -> dict:
    """Parse a VCF or BED of variants -> ``{0-based pos: "REF>ALT"}``."""
    sites: Dict[int, str] = {}
    lower = str(source).lower()
    is_vcf = lower.endswith(".vcf") or lower.endswith(".vcf.gz")
    with open(source) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            if is_vcf:
                # VCF: POS is 1-based, need REF+ALT
                if len(parts) < 5:
                    continue
                try:
                    chrom, pos1, ref, alt = parts[0], int(parts[1]), parts[3], parts[4].split(",")[0]
                except (ValueError, IndexError):
                    continue
            else:
                # BED3+ : start is 0-based
                try:
                    chrom, pos1, ref, alt = parts[0], int(parts[1]) + 1, parts[3], parts[4]
                except (ValueError, IndexError):
                    continue
            if chrom_filter is not None and chrom != chrom_filter:
                continue
            sites[pos1 - 1] = f"{ref}>{alt}"
    return sites


def _parse_junction_bed(bed_path: str, chrom_filter: Optional[str] = None) -> dict:
    """Parse a junction BED -> ``{(start, end): count}``.

    Columns: ``chrom<TAB>start<TAB>end<TAB>[name]<TAB>[count]``. If the 5th
    column is numeric it is used as the supporting-read count, else 1.
    """
    counts: Dict[Tuple[int, int], int] = {}
    with open(bed_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith(("track", "browser")):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            if chrom_filter is not None and parts[0] != chrom_filter:
                continue
            try:
                s, e = int(parts[1]), int(parts[2])
            except ValueError:
                continue
            if e <= s:
                continue
            count = 1
            if len(parts) > 4:
                try:
                    count = max(int(parts[4]), 1)
                except ValueError:
                    pass
            counts[(s, e)] = counts.get((s, e), 0) + count
    return counts


def _iupac_to_regex(motif: str) -> str:
    """Translate an IUPAC motif to a regex character class (case-insensitive)."""
    classes = {
        "A": "A", "C": "C", "G": "G", "T": "T", "U": "T",
        "R": "AG", "Y": "CT", "S": "GC", "W": "AT", "K": "GT",
        "M": "AC", "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG",
        "N": "ACGT", ".": "ACGT",
    }
    return "".join(f"[{classes.get(ch.upper(), ch)}]" for ch in motif)


def _style_gene_labels(ax) -> None:
    """Restyle dna_features_viewer's inline labels to a clean, modern look.

    Strips the boxy default annotation and applies a subtle rounded translucent
    label with neutral text, so the feature names read clearly without the
    dated grey boxes.
    """
    bbox = dict(
        boxstyle="round,pad=0.14,rounding_size=0.35",
        facecolor="white",
        edgecolor="#d9dde3",
        linewidth=0.6,
        alpha=0.92,
    )
    for t in ax.texts:
        t.set_bbox(bbox)
        t.set_color("#333333")
        t.set_fontsize(_fs(11))


@dataclass
class Track:
    """A single stacked track: a drawing recipe plus a height-weight."""

    weight: float
    kind: str
    draw: object

    def render(self, ax, region) -> None:
        self.draw(ax, region)


class GenomeView:
    """Collect tracks and render them as one stacked, shared-x figure."""

    def __init__(
        self,
        region: Union[Region, str, tuple],
        figsize: Tuple[float, float] = (14, 7),
        dpi: int = 120,
        reference: Optional[str] = None,
        title: Optional[str] = None,
    ):
        self.region = Region.from_any(region)
        self.figsize = figsize
        self.dpi = dpi
        self.tracks: List[Track] = []
        self.sites: Dict[int, str] = {}
        self._reference = Reference(reference) if reference is not None else None
        self._legend_items: List[Tuple[str, str]] = []
        self._show_legend = True
        self._highlight_regions: List[Tuple[float, float, str, float]] = []
        self._base_mods: Dict[int, tuple] = {}
        self._title: Optional[str] = None
        # regions this small default to the base-resolution "show all bases" view
        self.auto_base_threshold = 200
        if title:
            self.set_title(title)

    # ------------------------------------------------------------------ #
    # public track adders
    # ------------------------------------------------------------------ #
    def add_features(
        self,
        source=None,
        weight: float = 1.2,
        plot_kwargs: Optional[dict] = None,
        feature_types: Optional[set] = None,
        min_feature_length: int = 0,
    ) -> "GenomeView":
        """Add a gene/feature track rendered by dna_features_viewer.

        ``source`` is anything accepted by :func:`igvplot.features.load_features`
        (a .gb/.gff path, a Biopython SeqRecord, a GraphicRecord, or a list of
        GraphicFeature).

        ``feature_types`` keeps only features whose ``type`` is in the set
        (e.g. ``{"gene", "CDS", "mRNA", "exon"}``) — handy to hide tiny point
        variants and show a clean structural gene track. ``min_feature_length``
        drops features shorter than this many base pairs (default 0 = keep all).
        """
        region = self.region
        if source is None:
            record = GraphicRecord(first_index=region.start, sequence_length=region.length, features=[])
        else:
            record = load_features(source, region=region)

        feats = list(record.features or [])
        if feature_types is not None:
            accepted = {str(t).lower() for t in feature_types}
            feats = [
                f for f in feats
                if str(getattr(f, "feature_type", "")).lower() in accepted
            ]
        if min_feature_length:
            feats = [
                f for f in feats
                if abs(getattr(f, "end", 0) - getattr(f, "start", 0)) >= min_feature_length
            ]
        record = GraphicRecord(
            first_index=region.start,
            sequence_length=region.length,
            features=feats,
        ) if feats != list(record.features or []) else record

        kwargs = dict(
            with_ruler=False,
            draw_line=True,
            annotate_inline=True,
        )
        if plot_kwargs:
            kwargs.update(plot_kwargs)

        def _draw(ax, region):
            record.plot(ax=ax, **kwargs)
            _style_gene_labels(ax)
            # Clip any ruler ticks left behind; the bottom axis owns the ruler.
            ax.set_xticks([])
            ax.tick_params(axis="x", which="both", length=0)
            ax.set_yticks([])

        self.tracks.append(Track(weight=weight, kind="features", draw=_draw))
        return self

    def add_coverage(
        self,
        bam_path: Optional[str] = None,
        bigwig: Optional[str] = None,
        depths: Optional[np.ndarray] = None,
        reference: Optional[Union[str, Reference]] = None,
        min_mapq: int = 0,
        strand: Optional[str] = None,
        weight: float = 1.0,
        fill_color: str = COVERAGE,
        ylabel: str = "depth",
        ymax: Optional[float] = None,
    ) -> "GenomeView":
        """Add a per-base coverage track.

        Supply exactly one of ``bam_path`` (pileup from an indexed BAM/CRAM),
        ``bigwig`` (a BigWig file, for very large files), or precomputed
        ``depths`` (length == region.length, 0-based). ``ymax`` sets a fixed
        maximum for the y-axis (None = autoscale).
        """
        region = self.region
        if depths is None:
            if bigwig is not None:
                depths = read_bigwig_coverage(bigwig, region)
                mismatches = None
            elif bam_path is not None:
                depths, mismatches = compute_coverage(
                    bam_path, region, reference=reference, min_mapq=min_mapq, strand=strand
                )
            else:
                raise ValueError(
                    "add_coverage needs one of bam_path, bigwig or precomputed depths"
                )
        else:
            mismatches = None

        # propagate any reference (path or Reference) so later add_reads /
        # add_sashimi calls reuse it for mismatch detection
        if reference is not None:
            old = self._reference
            self._reference = reference
            if isinstance(old, Reference) and old is not reference:
                old.close()

        def _draw(ax, region):
            draw_coverage_track(
                ax,
                depths,
                region=region,
                fill_color=fill_color,
                mismatch_counts=mismatches,
                ylabel=ylabel,
                ymax=ymax,
            )

        self.tracks.append(Track(weight=weight, kind="coverage", draw=_draw))
        return self

    # ------------------------------------------------------------------ #
    # legend / highlight controls
    # ------------------------------------------------------------------ #
    def add_legend_items(self, items: List[Tuple[str, str]]) -> "GenomeView":
        self._legend_items.extend(items)
        return self

    def set_legend(self, show: bool) -> "GenomeView":
        self._show_legend = show
        return self

    def set_title(self, title: str) -> "GenomeView":
        self._title = title
        return self

    def add_base_mods(
        self,
        sites,
        weight: float = 0.5,
    ) -> "GenomeView":
        """A base-modification (methylation) track, IGV-style.

        ``sites`` maps a 0-based position to either a label (strand +1),
        ``(strand, label)``, or ``(strand, label, color)`` (or is an iterable of
        those tuples). Reads drawn with ``color_by='basemod'`` are coloured by
        the modification they span. Supports a `--basemod` BED in the CLI.
        """
        mods = self._parse_mod_sites(sites)
        if mods:
            self._base_mods.update(mods)

        def _draw(ax, region):
            draw_base_mod_track(ax, self._base_mods, region=region)

        self.tracks.append(Track(weight=weight, kind="base_mod", draw=_draw))
        return self

    def add_hic(
        self,
        matrix,
        cmap: str = "Reds",
        vmin: float = None,
        vmax: float = None,
        weight: float = 3.0,
        ylabel: str = "Hi-C contact",
        interpolation: str = "nearest",
    ) -> "GenomeView":
        """Add a Hi-C contact-map heatmap track (pyGenomeTracks 'hic' style).

        ``matrix`` is a dense ``(n, n)`` array of contact counts sampled across
        the region (bin ``i`` spans genomic bin ``i``); it is rendered as a
        heatmap aligned to the region's genomic coordinates.
        """
        mat = np.asarray(matrix, dtype=float)
        if mat.ndim != 2:
            raise ValueError(
                f"add_hic expects a 2-D contact matrix, got shape {mat.shape}"
            )

        def _draw(ax, region):
            ax.imshow(
                mat,
                aspect="auto",
                cmap=cmap,
                vmin=vmin,
                vmax=vmax,
                interpolation=interpolation,
                extent=[region.start, region.end, region.end, region.start],
                origin="upper",
                zorder=1,
                clip_on=True,
            )
            ax.set_yticks([])
            ax.set_ylabel(ylabel, fontsize=_fs(11))

        self.tracks.append(Track(weight=weight, kind="hic", draw=_draw))
        return self

    def add_tads(self, boundaries, color: str = TAD, label: str = "TAD") -> "GenomeView":
        """Add a TAD/domain-boundary track: dashed vertical lines at ``boundaries``
        (0-based positions) with a small triangle marker."""
        import matplotlib.patches as mpatches

        def _draw(ax, region):
            for b in boundaries:
                b = float(b)
                if not (region.start <= b < region.end):
                    continue
                ax.axvline(b, color=color, ls="--", lw=1.0, alpha=0.8, zorder=2)
                ax.add_patch(
                    mpatches.Polygon(
                        [[b - 6.0, 0.0], [b + 6.0, 0.0], [b, 1.0]],
                        closed=True,
                        facecolor=color,
                        edgecolor="none",
                        zorder=3,
                    )
                )
            ax.axhline(0, color="#c9d1d9", lw=0.8, zorder=1)  # baseline
            ax.set_xlim(region.start, region.end)
            ax.set_ylim(-0.2, 1.0)
            ax.set_yticks([])
            ax.set_ylabel(label, fontsize=_fs(11))

        self.tracks.append(Track(weight=0.4, kind="tads", draw=_draw))
        return self

    def add_scale_bar(self, window_bp: int = None, label: str = None, weight: float = 0.3) -> "GenomeView":
        """Add a scale-bar track (pyGenomeTracks 'scale' style)."""
        import matplotlib.patches as mpatches

        def _draw(ax, region):
            wb = int(window_bp) if window_bp else max(1, int(region.length * 0.25))
            ax.add_patch(
                mpatches.Rectangle(
                    (region.start, -0.12), wb, 0.24, facecolor="#333333", edgecolor="none"
                )
            )
            ax.text(
                region.start + wb,
                0.20,
                label or _human_bp(wb),
                ha="center",
                va="bottom",
                fontsize=_fs(9.5),
                color="#333333",
            )
            ax.set_xlim(region.start, region.end)
            ax.set_ylim(-0.3, 0.6)
            ax.set_yticks([])
            ax.set_ylabel("scale", fontsize=8)

        self.tracks.append(Track(weight=weight, kind="scale", draw=_draw))
        return self

    def add_bed_features(
        self,
        bed_path: str,
        weight: float = 1.0,
        color: str = BED,
        feature_type: str = "annotation",
    ) -> "GenomeView":
        """Add a BED/narrowPeak annotation track (pyGenomeTracks 'bed' style).
        Reads BED3+ (start/end 0-based, col4 = label); narrowPeak score is used
        to scale peak height."""
        import matplotlib.patches as mpatches

        features = _parse_bed_features(bed_path)

        def _draw(ax, region):
            for chrom, s, e, name, score in features:
                if chrom != region.chrom:
                    continue
                s = max(s, region.start)
                e = min(e, region.end)
                if e <= s:
                    continue
                h = min(1.0, 0.2 + float(score or 0) / 1000.0) if score is not None else 0.6
                ax.add_patch(
                    mpatches.Rectangle(
                        (s, 0), e - s, h, facecolor=color, edgecolor="#333333",
                        lw=0.4, zorder=2,
                    )
                )
                if name:
                    ax.text(s, h + 0.05, name, ha="left", va="bottom", fontsize=_fs(9), color="#333333")
            ax.set_xlim(region.start, region.end)
            ax.set_ylim(0, 1.15)
            ax.set_yticks([])
            ax.set_ylabel(feature_type, fontsize=8)

        self.tracks.append(Track(weight=weight, kind="bed", draw=_draw))
        return self

    @staticmethod
    def _parse_mod_sites(sites) -> Dict[int, tuple]:
        from .plot import base_mod_color

        mods: Dict[int, tuple] = {}
        entries = sites.items() if hasattr(sites, "items") else list(sites)
        for entry in entries:
            if hasattr(sites, "items"):
                pos, val = entry
                if isinstance(val, (tuple, list)):
                    if len(val) >= 3:
                        strand, label, color = val[0], val[1], val[2]
                    else:
                        strand, label = val[0], val[1]
                        color = base_mod_color(label)
                else:
                    strand, label, color = 1, val, base_mod_color(val)
            else:
                # (pos, strand, label[, color])
                if len(entry) == 3:
                    pos, strand, label = entry
                    color = base_mod_color(label)
                elif len(entry) == 4:
                    pos, strand, label, color = entry
                else:
                    raise ValueError(f"Unrecognised base-mod entry: {entry!r}")
            mods[int(pos)] = (int(strand), str(label), color)
        return mods

    def add_highlight_regions(
        self,
        regions,
        color: str = HIGH_BG,
        alpha: float = 0.3,
    ) -> "GenomeView":
        """Shade background ``regions`` across every track.

        ``regions`` is an iterable of ``(start, end)`` (0-based) or
        ``(start, end, color, alpha)`` tuples.
        """
        for item in regions:
            if len(item) == 2:
                start, end = item
                rcolor, ralpha = color, alpha
            elif len(item) == 4:
                start, end, rcolor, ralpha = item
            else:
                raise ValueError("highlight regions must be (start, end) or (start, end, color, alpha)")
            self._highlight_regions.append((float(start), float(end), rcolor, ralpha))
        return self

    def add_coverage_overlay(
        self,
        samples: List[Tuple[Optional[str], str, str]],
        reference: Optional[str] = None,
        min_mapq: int = 0,
        weight: float = 1.2,
        ymax: Optional[float] = None,
        label: str = "coverage",
    ) -> "GenomeView":
        """Overlay coverage from several samples on a single shared axis.

        ``samples`` is a list of ``(bam_path_or_None, color, label)``; each
        entry's coverage is drawn as an area on the same axis with a legend.
        """

        def _draw(ax, region):
            x = np.arange(region.start, region.end, dtype=float)
            series, tops = [], []
            for bam_path, color, slabel in samples:
                if bam_path is None:
                    raise ValueError("add_coverage_overlay needs a bam_path per sample")
                depths, _ = compute_coverage(bam_path, region, reference=reference, min_mapq=min_mapq)
                series.append((depths, color, slabel))
                tops.append(float(depths.max()) if depths.size else 0.0)
            top = ymax if ymax is not None else max(tops, default=1.0)
            for depths, color, slabel in series:
                ax.fill_between(x, depths, step="mid", color=color, alpha=0.32, lw=0, label=slabel)
                ax.step(x, depths, where="mid", color=color, lw=1.1, zorder=3)
            ax.set_ylim(0, max(float(top), 1e-9))
            ax.set_yticks([])
            if len(samples) > 1:
                ax.legend(fontsize=_fs(9), frameon=False, ncol=len(samples), loc="upper right")
            ax.set_ylabel(label, fontsize=_fs(11))

        self.tracks.append(Track(weight=weight, kind="coverage_overlay", draw=_draw))
        return self

    def add_coverage_strands(
        self,
        bam_path: str,
        reference: Optional[Union[str, Reference]] = None,
        min_mapq: int = 0,
        weight: float = 1.4,
        forward_color: str = "#457b9d",
        reverse_color: str = "#e63946",
        ylabel: str = "depth",
        ymax: Optional[float] = None,
    ) -> "GenomeView":
        """Add strand-specific coverage (RNA-seq sense/antisense).

        Forward-strand coverage is drawn above the baseline, reverse-strand
        coverage is reflected below it, so you can read strand-specific
        transcription in one track.
        """

        def _draw(ax, region):
            x = np.arange(region.start, region.end, dtype=float)
            fwd, _ = compute_coverage(bam_path, region, reference=reference,
                                      min_mapq=min_mapq, strand="forward")
            rev, _ = compute_coverage(bam_path, region, reference=reference,
                                      min_mapq=min_mapq, strand="reverse")
            ax.fill_between(x, fwd, step="mid", color=forward_color, alpha=0.5, lw=0, zorder=2)
            ax.step(x, fwd, where="mid", color=forward_color, lw=1.1, zorder=3)
            ax.fill_between(x, -rev, step="mid", color=reverse_color, alpha=0.5, lw=0, zorder=2)
            ax.step(x, -rev, where="mid", color=reverse_color, lw=1.1, zorder=3)
            ax.axhline(0, color="#888888", lw=0.8, zorder=1)
            top = ymax if ymax is not None else max(float(fwd.max()), float(rev.max()), 1.0)
            ax.set_ylim(-top, top)
            ax.set_yticks([])
            ax.set_ylabel(ylabel, fontsize=_fs(11))

        self.tracks.append(Track(weight=weight, kind="coverage_strands", draw=_draw))
        return self

    def add_reads_overlay(
        self,
        samples,
        reference: Optional[Union[str, Reference]] = None,
        min_mapq: int = 0,
        max_reads: Optional[int] = None,
        weight: float = 3.0,
        label: str = "reads",
        paint_base_letters: bool = True,
    ) -> "GenomeView":
        """Overlay raw alignment reads from several BAMs on one shared axis.

        ``samples`` is a list of ``(bam_path, color, label)``. Reads from all
        samples are packed together and each sample is drawn in its own flat
        colour (with a legend), so you can compare alignments across conditions
        or replicates directly.
        """
        import matplotlib.patches as mpatches

        ref = reference or self._reference

        def _draw(ax, region):
            combined: List[Read] = []
            colors = {}
            for bam, color, slabel in samples:
                if bam is None:
                    raise ValueError("add_reads_overlay needs a bam_path per sample")
                for r in fetch_reads(bam, region, reference=ref, min_mapq=min_mapq,
                                     max_reads=max_reads):
                    r._overlay_label = slabel
                    combined.append(r)
                colors[slabel] = color

            def color_of(i, r):
                return colors.get(getattr(r, "_overlay_label", None), "#9e9e9e")

            draw_read_track(
                ax, combined, region=region, color_fn=color_of,
                paint_base_letters=paint_base_letters, display_mode="expanded",
                show_soft_clips=False,
            )
            if len(samples) > 1:
                handles = [mpatches.Patch(facecolor=c, edgecolor="none", label=lab) for _, c, lab in samples]
                ax.legend(handles=handles, fontsize=_fs(9), frameon=False,
                          ncol=min(3, len(samples)), loc="upper right")
            ax.set_ylabel(label, fontsize=_fs(11))

        self.tracks.append(Track(weight=weight, kind="reads_overlay", draw=_draw))
        return self

    def add_sashimi_overlay(
        self,
        samples,
        reference: Optional[Union[str, Reference]] = None,
        min_counts: int = 1,
        weight: float = 1.6,
    ) -> "GenomeView":
        """Compare splice-junction usage across samples (ggsashimi-style).

        ``samples`` is a list of ``(bam_path, color, label)``; each sample's
        junctions are drawn in its own colour in a stacked band on one track.
        """
        ref = reference or self._reference

        def _draw(ax, region):
            data = []
            for bam, color, slabel in samples:
                if bam is None:
                    raise ValueError("add_sashimi_overlay needs a bam_path per sample")
                reads = fetch_reads(bam, region, reference=ref)
                counts = junction_counts(reads, region, min_counts=min_counts)
                data.append((color, slabel, counts))
            draw_sashimi_overlay_track(ax, data, region=region)

        self.tracks.append(Track(weight=weight, kind="sashimi_overlay", draw=_draw))
        return self

    def add_junctions_bed(
        self,
        bed_path: str,
        arc_color: str = SASHIMI,
        weight: float = 1.2,
        label: str = "junction",
    ) -> "GenomeView":
        """Draw precomputed splice junctions from a BED file (no BAM needed).

        The BED is ``chrom<TAB>start<TAB>end<TAB>[name]<TAB>[count]``; arc
        thickness/height scale with the count (default 1).
        """
        region = self.region
        counts = _parse_junction_bed(bed_path, region.chrom)

        def _draw(ax, region):
            draw_sashimi_track(ax, counts, region=region, arc_color=arc_color)
            ax.set_ylabel(label, fontsize=_fs(11))

        self.tracks.append(Track(weight=weight, kind="junction_bed", draw=_draw))
        return self

    def add_signal(
        self,
        values,
        color: str = COVERAGE,
        ylabel: str = "signal",
        ymax: Optional[float] = None,
        weight: float = 1.0,
    ) -> "GenomeView":
        """Add a generic numeric signal track from a 1-D array.

        ``values`` has length ``region.length`` and is drawn as a filled
        stepped area aligned to the region. Useful for any score/profile
        (e.g. a signal array, a per-base score from a model).
        """
        region = self.region
        vals = np.asarray(values, dtype=float)
        if vals.ndim != 1:
            raise ValueError(f"add_signal expects a 1-D array, got {vals.ndim}D")
        if vals.size != region.length:
            raise ValueError(
                f"signal length {vals.size} != region length {region.length}"
            )

        def _draw(ax, region):
            x = np.arange(region.start, region.end, dtype=float)
            ax.fill_between(x, vals, step="mid", color=color, alpha=0.4, lw=0, zorder=1)
            ax.step(x, vals, where="mid", color=color, lw=1.2, zorder=2)
            lo = min(0.0, float(vals.min())) if vals.size else 0.0
            hi = ymax if ymax is not None else (float(vals.max()) if vals.size else 1.0)
            ax.set_ylim(min(lo, hi), max(hi, 1e-9))
            ax.set_yticks([])
            ax.set_ylabel(ylabel, fontsize=_fs(11))

        self.tracks.append(Track(weight=weight, kind="signal", draw=_draw))
        return self

    def add_gc(
        self,
        reference: Optional[Union[str, Reference]] = None,
        window: int = 50,
        step: int = 10,
        color: str = "#2a9d8f",
        ylabel: str = "GC%",
        weight: float = 0.8,
    ) -> "GenomeView":
        """Add a GC-content track computed from the reference sequence.

        GC% is computed in ``window``-bp windows every ``step`` bp across the
        region and drawn as a stepped line. Requires a reference fasta.
        """
        region = self.region
        ref = reference or self._reference
        if ref is None:
            raise ValueError("add_gc requires a reference fasta")
        opened = not isinstance(ref, Reference)
        r = ref if isinstance(ref, Reference) else Reference(ref)
        seq = r.get(region.chrom, region.start, region.end)
        if opened:
            r.close()
        if not seq:
            raise ValueError(
                f"reference has no sequence for {region.chrom!r} in {region!r}"
            )

        def _draw(ax, region):
            xs, gcs = [], []
            for start in range(region.start, region.end, step):
                lo = start - region.start
                hi = lo + window
                if hi > len(seq):
                    break
                chunk = seq[lo:hi].upper()
                if not chunk:
                    continue
                g = chunk.count("G") + chunk.count("C")
                gc = 100.0 * g / len(chunk)
                xs.append(start + window / 2)
                gcs.append(gc)
            if not xs:
                ax.axis("off")
                return
            ax.fill_between(xs, gcs, step="mid", color=color, alpha=0.35, lw=0, zorder=1)
            ax.step(xs, gcs, where="mid", color=color, lw=1.2, zorder=2)
            ax.set_xlim(region.start, region.end)
            ax.set_ylim(0, 100)
            ax.set_yticks([])
            ax.set_ylabel(ylabel, fontsize=_fs(11))

        self.tracks.append(Track(weight=weight, kind="gc", draw=_draw))
        return self

    def add_variants(self, source: str, color: str = "#e0492f") -> "GenomeView":
        """Mark variants from a VCF/BED file as ``REF>ALT`` across every track.

        Parses the file and registers the sites (like :meth:`add_sites`), so the
        variant positions and labels appear as markers across all stacked tracks.
        """
        region = self.region
        sites = _parse_variants(source, region.chrom)
        self.sites.update(sites)
        return self

    def add_arc(
        self,
        pairs,
        color: str = "#e63946",
        weight: float = 2.0,
        label: str = "interaction",
        max_height: float = 1.0,
    ) -> "GenomeView":
        """Add generic interaction/loop arcs between genomic positions.

        ``pairs`` is an iterable of ``(start, end)`` or ``(start, end, strength)``
        where ``strength`` scales arc height (and labels arcs with strength > 1.5).
        Useful for loops, RNA-RNA contacts or any pairwise connection.
        """

        def _draw(ax, region):
            heights = []
            for item in pairs:
                if len(item) == 2:
                    s, e, strength = item[0], item[1], 1.0
                else:
                    s, e, strength = item[0], item[1], float(item[2])
                draw_interaction_arc(ax, s, e, region, color=color,
                                     strength=strength, max_height=max_height)
                heights.append(strength)
            top = max_height * max(heights, default=1.0)
            ax.set_xlim(region.start, region.end)
            ax.set_ylim(-0.03, top * 1.12 + 0.03)
            ax.set_yticks([])
            ax.set_ylabel(label, fontsize=_fs(11))

        self.tracks.append(Track(weight=weight, kind="arc", draw=_draw))
        return self

    def add_track(self, callback, weight: float = 1.0) -> "GenomeView":
        """Add a fully custom track drawn by a user callback.

        ``callback(ax, region)`` is called with the target axis and the
        ``Region``; draw anything with matplotlib on ``ax``.
        """
        def _draw(ax, region):
            callback(ax, region)

        self.tracks.append(Track(weight=weight, kind="custom", draw=_draw))
        return self

    def add_mod_fraction(
        self,
        sites,
        color: str = "#e63946",
        ymax: float = 1.0,
        weight: float = 1.0,
        label: str = "mod%",
    ) -> "GenomeView":
        """Draw a per-site base-modification (stoichiometry) bar track.

        ``sites`` maps a 0-based position to a fraction (0-1), or to
        ``(strand, label, fraction)`` / ``(strand, label, fraction, color)``.
        Each site is drawn as a bar of height = fraction, IGV-modification style.
        """
        import matplotlib.patches as mpatches

        sites = self._parse_mod_fraction(sites)

        def _draw(ax, region):
            for pos, frac in sites.items():
                if not region.overlaps(pos, pos + 1):
                    continue
                frac = max(0.0, min(float(ymax), frac))
                ax.add_patch(
                    mpatches.Rectangle(
                        (pos + 0.15, 0), 0.7, frac,
                        facecolor=color, edgecolor="none", alpha=0.9, zorder=2,
                    )
                )
            ax.set_xlim(region.start, region.end)
            ax.set_ylim(0, float(ymax) * 1.08)
            ax.set_yticks([])
            ax.set_ylabel(label, fontsize=_fs(11))

        self.tracks.append(Track(weight=weight, kind="mod_fraction", draw=_draw))
        return self

    def add_variant_fraction(
        self,
        bam_path: str,
        reference: Optional[Union[str, Reference]] = None,
        min_mapq: int = 0,
        color: str = "#e0492f",
        weight: float = 1.0,
        label: str = "VAF",
    ) -> "GenomeView":
        """Draw a per-base variant-allele-fraction (0-1) track from a BAM.

        For each base the fraction of reads whose base differs from the
        reference (mismatch pileup / depth) is drawn as a stepped area — useful
        for SNP/allelic-imbalance hotspots.
        """
        region = self.region
        ref = reference or self._reference
        depths, mism = compute_coverage(bam_path, region, reference=ref, min_mapq=min_mapq)
        frac = np.zeros(region.length, dtype=float)
        if mism is not None:
            np.divide(mism, depths, out=frac, where=depths > 0)

        def _draw(ax, region):
            x = np.arange(region.start, region.end, dtype=float)
            ax.fill_between(x, frac, step="mid", color=color, alpha=0.85, lw=0, zorder=2)
            ax.step(x, frac, where="mid", color=color, lw=1.1, zorder=3)
            ax.set_xlim(region.start, region.end)
            top = float(frac.max()) if frac.size else 1.0
            ax.set_ylim(0, max(top, 1e-6) * 1.1)
            ax.set_yticks([])
            ax.set_ylabel(label, fontsize=_fs(11))

        self.tracks.append(Track(weight=weight, kind="variant_fraction", draw=_draw))
        return self

    def add_motifs(
        self,
        motif: str,
        reference: Optional[Union[str, Reference]] = None,
        color: str = "#9b5de5",
        weight: float = 0.7,
        label: str = "motif",
        max_hits: int = 500,
    ) -> "GenomeView":
        """Mark occurrences of a motif in the reference (IUPAC-enabled).

        ``motif`` may contain IUPAC codes (e.g. ``DRACH`` for the m6A context);
        each occurrence is drawn as a small bar spanning the motif length.
        """
        import re as _re

        region = self.region
        ref = reference or self._reference
        if ref is None:
            raise ValueError("add_motifs requires a reference fasta")
        opened = not isinstance(ref, Reference)
        r = ref if isinstance(ref, Reference) else Reference(ref)
        seq = r.get(region.chrom, region.start, region.end)
        if opened:
            r.close()
        pattern = _re.compile(_iupac_to_regex(motif), _re.I)
        hits = []
        for m in _re.finditer(pattern, seq):
            ln = len(m.group(0) or motif)
            hits.append((region.start + m.start(), ln))
            if len(hits) >= max_hits:
                break

        def _draw(ax, region):
            import matplotlib.patches as mpatches
            for pos, ln in hits:
                ax.add_patch(
                    mpatches.Rectangle(
                        (pos, 0.05), ln, 0.55,
                        facecolor=color, edgecolor="none", alpha=0.9, zorder=2,
                    )
                )
            ax.set_xlim(region.start, region.end)
            ax.set_ylim(-0.1, 0.7)
            ax.set_yticks([])
            ax.set_ylabel(label, fontsize=_fs(11))

        self.tracks.append(Track(weight=weight, kind="motifs", draw=_draw))
        return self

    @staticmethod
    def _parse_mod_fraction(sites) -> Dict[int, float]:
        """Normalise a mod-fraction spec to ``{0-based pos: fraction}``."""
        out: Dict[int, float] = {}
        for pos, vals in sites.items() if hasattr(sites, "items") else sites:
            if hasattr(sites, "items"):
                if isinstance(vals, (tuple, list)):
                    frac = float(vals[0])
                else:
                    frac = float(vals)
            else:
                # (pos, frac) or (pos, strand, label, frac)
                pos, frac = vals[0], float(vals[-1])
            out[int(pos)] = max(0.0, min(1.0, frac))
        return out

    def add_reads(
        self,
        bam_path: Optional[str] = None,
        reads: Optional[List[Read]] = None,
        reference: Optional[Union[str, Reference]] = None,
        min_mapq: int = 0,
        max_reads: Optional[int] = None,
        keep_duplicates: bool = False,
        keep_secondary: bool = False,
        paint_base_letters: bool = True,
        base_fontsize: float = 11.0,
        color_by: str = "strand",
        colormap: str = "viridis",
        link_mates: bool = False,
        group_by: str = "none",
        show_soft_clips: bool = False,
        display_mode: str = "expanded",
        view_as_pairs: bool = False,
        highlight: Optional[set] = None,
        highlight_color: str = "#e67e22",
        sampling_window: int = 0,
        max_per_window: int = 0,
        show_all_bases: Optional[bool] = None,
        sort_by: str = "start",
        show_insertion_text: bool = False,
        show_deletion_text: bool = False,
        basemod_sites=None,
        sort_base_pos: Optional[int] = None,
        weight: float = 3.5,
    ) -> "GenomeView":
        """Add an IGV-style aligned-reads pileup track.

        ``color_by`` may be 'strand', 'firstOfPairStrand', 'pairOrientation',
        'mapq', 'insert'/'tlen', 'unexpectedPair', 'readgroup', 'proper',
        'mate' or 'none'. ``group_by`` clusters reads vertically by the same
        set of attributes. ``link_mates`` draws connectors between the two
        mates of a paired read when both are present. ``show_soft_clips``
        draws soft-clipped ends as thin bars. ``display_mode`` is 'expanded' or
        'squished'. ``highlight`` names get a coloured outline.
        ``sampling_window``/``max_per_window`` enable IGV-style downsampling.
        ``show_all_bases`` renders every read base as a letter (base
        resolution, IGV "show all bases"); defaults to automatic for small
        regions.
        """
        region = self.region
        if show_all_bases is None:
            show_all_bases = region.length <= self.auto_base_threshold
        collect_bases = bool(show_all_bases)

        if reads is None:
            if bam_path is None:
                raise ValueError("add_reads needs bam_path or a list of reads")
            reads = fetch_reads(
                bam_path,
                region,
                reference=reference or self._reference,
                min_mapq=min_mapq,
                max_reads=max_reads,
                keep_duplicates=keep_duplicates,
                keep_secondary=keep_secondary,
                sampling_window=sampling_window,
                max_per_window=max_per_window,
                collect_bases=collect_bases,
            )

        def _draw(ax, region):
            draw_read_track(
                ax,
                reads,
                region=region,
                paint_base_letters=paint_base_letters,
                base_fontsize=base_fontsize,
                color_by=color_by,
                colormap=colormap,
                link_mates=link_mates,
                group_by=group_by,
                show_soft_clips=show_soft_clips,
                display_mode=display_mode,
                view_as_pairs=view_as_pairs,
                highlight=highlight,
                highlight_color=highlight_color,
                show_all_bases=show_all_bases,
                sort_by=sort_by,
                show_insertion_text=show_insertion_text,
                show_deletion_text=show_deletion_text,
                basemod_sites=basemod_sites or self._base_mods,
                sort_base_pos=sort_base_pos,
            )

        self.tracks.append(Track(weight=weight, kind="reads", draw=_draw))
        if self._show_legend:
            self.add_legend_items(build_legend_items(color_by, reads, group_by, basemod_sites or self._base_mods))
        return self

    def add_sashimi(
        self,
        bam_path: Optional[str] = None,
        reads: Optional[List[Read]] = None,
        min_counts: int = 1,
        arc_color: str = SASHIMI,
        weight: float = 1.2,
    ) -> "GenomeView":
        """Add a splice-junction (sashimi) track above the reads.

        Junctions are counted from CIGAR 'N' (intron-skipped) operations in
        the aligned reads via :func:`igvplot.reads.junction_counts`.
        """
        region = self.region
        if reads is None:
            if bam_path is None:
                raise ValueError("add_sashimi needs bam_path or a list of reads")
            reads = fetch_reads(bam_path, region, reference=self._reference)
        counts = junction_counts(reads, region, min_counts=min_counts)

        def _draw(ax, region):
            draw_sashimi_track(ax, counts, region=region, arc_color=arc_color)

        self.tracks.append(Track(weight=weight, kind="sashimi", draw=_draw))
        return self

    def add_sequence(
        self,
        reference: Union[str, Reference],
        weight: float = 0.5,
    ) -> "GenomeView":
        """Add a reference-base sequence row (requires the reference fasta)."""
        region = self.region
        opened = not isinstance(reference, Reference)
        ref = reference if isinstance(reference, Reference) else Reference(reference)
        seq = ref.get(region.chrom, region.start, region.end)
        if opened:
            ref.close()

        def _draw(ax, region):
            draw_sequence_track(ax, seq, region=region)

        self.tracks.append(Track(weight=weight, kind="sequence", draw=_draw))
        return self

    def bam(
        self,
        bam_path: str,
        reference: Optional[Union[str, Reference]] = None,
        coverage: bool = True,
        reads: bool = True,
        **read_kwargs,
    ) -> "GenomeView":
        """Convenience: add alignment coverage and read pile-up from one BAM."""
        ref = reference or getattr(self, "_reference", None)
        if coverage:
            self.add_coverage(bam_path, reference=ref)
        if reads:
            self.add_reads(bam_path, reference=ref, **read_kwargs)
        return self

    @property
    def figure(self):
        """Build and return the matplotlib Figure."""
        return self.render()

    def add_sites(self, sites: Dict[int, str]) -> "GenomeView":
        """Add vertical site markers across every track: ``{0-based pos: label}``."""
        self.sites.update(sites)
        return self

    def add_blank(self, weight: float = 0.3, label: Optional[str] = None) -> "GenomeView":
        """Add an empty spacer track (vertical gap) of the given height weight."""
        def _draw(ax, region):
            ax.set_xlim(region.start, region.end)
            ax.set_yticks([])
            ax.set_xticks([])
            if label:
                ax.set_ylabel(label, fontsize=_fs(11))

        self.tracks.append(Track(weight=weight, kind="spacer", draw=_draw))
        return self

    # ------------------------------------------------------------------ #
    # rendering
    # ------------------------------------------------------------------ #
    def render(self, fig: Optional[Figure] = None) -> Figure:
        if not self.tracks:
            raise ValueError("GenomeView has no tracks to render")

        apply_theme()
        region = self.region
        n = len(self.tracks)
        fig = fig or plt.figure(figsize=self.figsize, dpi=self.dpi)
        # clamp track weights to a small positive so a 0/negative weight can't
        # produce a zero-height (or degenerate) subplot
        heights = [max(float(t.weight), 1e-3) for t in self.tracks]
        gs = GridSpec(
            n,
            1,
            figure=fig,
            height_ratios=heights,
            hspace=0.08,
            left=0.06,
            right=0.985,
            bottom=0.06,
            top=0.96,
        )

        axes: List = []
        for i, track in enumerate(self.tracks):
            ax = fig.add_subplot(gs[i], sharex=axes[0] if axes else None)
            axes.append(ax)
            track.render(ax, region)

        # Site markers across all tracks.
        if self.sites:
            draw_sites(axes, region, self.sites)

        # Background highlight regions across all tracks.
        if self._highlight_regions:
            for start, end, hcolor, halpha in self._highlight_regions:
                start = max(start, region.start)
                end = min(end, region.end)
                if end <= start:
                    continue
                for ax in axes:
                    ax.axvspan(start, end, color=hcolor, alpha=halpha, zorder=0)

        # Finalize axes: single shared x range, 1-based genomic ruler on bottom.
        bottom = axes[-1]
        for ax in axes[:-1]:
            ax.set_xticks([])
            ax.tick_params(axis="x", length=0)
            ax.spines["top"].set_visible(False)
            ax.spines["bottom"].set_visible(False)
            for sp in ("left", "right"):
                ax.spines[sp].set_visible(True)
            ax.set_xlim(region.start, region.end)
        bottoms = {
            "left": True,
            "right": False,
            "top": False,
            "bottom": True,
        }
        for sp, vis in bottoms.items():
            bottom.spines[sp].set_visible(vis)

        # Ensure dna_features_viewer didn't shrink the shared window.
        for ax in axes:
            ax.set_xlim(region.start, region.end)

        bottom.set_xlim(region.start, region.end)
        bottom.set_xticks(bottom.get_xticks())
        fmt = mticker.FuncFormatter(lambda x, p: f"{int(x + 1):,}")
        bottom.xaxis.set_major_formatter(fmt)
        bottom.set_xlabel(
            f"{region.chrom} position (1-based)", fontsize=_fs(12)
        )

        fig.align_ylabels(axes)
        if self._show_legend and self._legend_items:
            from matplotlib.patches import Patch

            handles = [
                Patch(facecolor=color, edgecolor="grey", label=label)
                for label, color in self._legend_items
            ]
            fig.legend(
                handles=handles,
                loc="upper left",
                bbox_to_anchor=(1.01, 1.0),
                frameon=False,
                fontsize=_fs(9.5),
                title="alignments",
                title_fontsize=_fs(9.5),
            )
        if self._title:
            fig.suptitle(self._title, fontsize=_fs(16))
            fig.subplots_adjust(top=0.93)
        return fig

    def savefig(self, out_path: str, dpi: Optional[int] = None, **kwargs) -> None:
        fig = self.render()
        fig.savefig(out_path, dpi=dpi or self.dpi, bbox_inches="tight", **kwargs)
        plt.close(fig)
        self.close()

    def close(self) -> None:
        """Release any :class:`Reference` held open by this view."""
        if self._reference is not None:
            obj, self._reference = self._reference, None
            close = getattr(obj, "close", None)
            if callable(close):
                close()

    def show(self) -> None:
        self.render()
        plt.show()

    def __repr__(self) -> str:
        return (
            f"{type(self).__name__}(region={self.region!r}, "
            f"tracks={len(self.tracks)}, dpi={self.dpi})"
        )


# Friendly aliases for the single builder class. ``IGV`` and ``AlignmentView``
# are the same object as ``GenomeView``; pick whichever reads best.
IGV = AlignmentView = GenomeView


def plot_view(
    bam_path: Optional[str] = None,
    region=None,
    out_path: Optional[str] = None,
    features=None,
    reference: Optional[str] = None,
    sites: Optional[Dict[int, str]] = None,
    min_mapq: int = 0,
    max_reads: Optional[int] = None,
    paint_base_letters: bool = True,
    show_coverage: bool = True,
    show_sequence: bool = False,
    sashimi: bool = False,
    sashimi_min_counts: int = 1,
    color_by: str = "strand",
    colormap: str = "viridis",
    link_mates: bool = False,
    group_by: str = "none",
    show_soft_clips: bool = False,
    display_mode: str = "expanded",
    view_as_pairs: bool = False,
    highlight: Optional[set] = None,
    show_legend: bool = True,
    sampling_window: int = 0,
    max_per_window: int = 0,
    show_all_bases: Optional[bool] = None,
    sort_by: str = "start",
    show_insertion_text: bool = False,
    show_deletion_text: bool = False,
    highlight_regions=None,
    basemod=None,
    sort_base_pos: Optional[int] = None,
    title: Optional[str] = None,
    bigwig: Optional[str] = None,
    coverage_ymax: Optional[float] = None,
    figsize: Tuple[float, float] = (14, 8),
    dpi: int = 150,
    feature_types: Optional[set] = None,
    min_feature_length: int = 0,
    coverage_strand: bool = False,
    variants=None,
    gc: bool = False,
    **kwargs,
) -> GenomeView:
    """One-call convenience to build and (optionally) save a full view.

    Parameters
    ----------
    bam_path:
        Sorted, indexed BAM/CRAM file (required unless only features are used).
    region:
        Anything acceptable to :meth:`Region.from_any`.
    out_path:
        If given, save the figure here.
    features:
        Annotation source for the gene track (see :func:`load_features`).
    reference:
        Reference fasta path; enables mismatch detection and the sequence row.
    sites:
        ``{0-based position: label}`` vertical markers across all tracks.
    min_mapq / max_reads / paint_base_letters:
        Read-track options.
    show_coverage / show_sequence:
        Toggle the coverage and reference-sequence rows.
    sashimi / sashimi_min_counts:
        Add a splice-junction (sashimi) track and its minimum read support.
    color_by / colormap / link_mates / group_by / show_soft_clips /
    display_mode:
        Read colouring mode, colormap for continuous modes, mate linking,
        vertical read clustering, soft-clip display and display mode
        ('expanded' / 'squished').
    highlight:
        Set of read names to outline in a highlighted colour.
    show_legend:
        Draw an auto-generated alignment legend.
    sampling_window / max_per_window:
        IGV-style downsampling for very deep/large files.
    show_all_bases:
        Render every read base as a letter (base resolution). Defaults to
        automatic for small regions and forces the reference-sequence row.
    sort_by / show_insertion_text / show_deletion_text:
        Read ordering attribute and inline indel size labels.
    highlight_regions:
        Iterable of ``(start, end)`` or ``(start, end, color, alpha)`` to
        shade across all tracks.
    basemod:
        Base-modification (methylation) sites, e.g. ``{pos: "m6A"}`` or
        ``[(pos, strand, label), ...]``. Also enables ``color_by="basemod"``.
    bigwig:
        Optional BigWig file for coverage (avoids pileup on huge BAMs).
    coverage_ymax:
        Optional fixed maximum for the coverage y-axis.
    """
    if region is None:
        raise ValueError("region is required")
    region = Region.from_any(region)
    view = GenomeView(region=region, figsize=figsize, dpi=dpi)
    view.set_legend(show_legend)

    # Automatic base-resolution: at single-base zoom, show all read bases and
    # the reference row so mutations can be read position by position.
    if show_all_bases is None:
        show_all_bases = region.length <= view.auto_base_threshold
    if show_all_bases and reference is not None:
        show_sequence = True

    ok_to_add = False

    if sashimi and bam_path is not None:
        view.add_sashimi(bam_path, min_counts=sashimi_min_counts)
        ok_to_add = True

    if features is not None:
        view.add_features(
            features,
            feature_types=feature_types,
            min_feature_length=min_feature_length,
        )
        ok_to_add = True

    if show_coverage:
        if coverage_strand and bam_path is not None:
            view.add_coverage_strands(bam_path)
            ok_to_add = True
        elif bigwig is not None:
            view.add_coverage(bigwig=bigwig, ymax=coverage_ymax)
            ok_to_add = True
        elif bam_path is not None:
            view.add_coverage(
                bam_path, reference=reference, min_mapq=min_mapq, ymax=coverage_ymax
            )
            ok_to_add = True

    if variants is not None:
        view.add_variants(variants)
        ok_to_add = True

    if gc and reference is not None:
        view.add_gc(reference)
        ok_to_add = True

    if basemod is not None:
        view.add_base_mods(basemod)
        ok_to_add = True

    if bam_path is not None:
        view.add_reads(
            bam_path,
            reference=reference,
            min_mapq=min_mapq,
            max_reads=max_reads,
            paint_base_letters=paint_base_letters,
            color_by=color_by,
            colormap=colormap,
            link_mates=link_mates,
            group_by=group_by,
            show_soft_clips=show_soft_clips,
            display_mode=display_mode,
            view_as_pairs=view_as_pairs,
            highlight=highlight,
            sampling_window=sampling_window,
            max_per_window=max_per_window,
            show_all_bases=show_all_bases,
            sort_by=sort_by,
            show_insertion_text=show_insertion_text,
            show_deletion_text=show_deletion_text,
            sort_base_pos=sort_base_pos,
        )
        ok_to_add = True

    if highlight_regions:
        view.add_highlight_regions(highlight_regions)

    if title:
        view.set_title(title)

    if show_sequence and reference is not None:
        view.add_sequence(reference)

    if sites:
        view.add_sites(sites)

    if not ok_to_add:
        raise ValueError("nothing to plot: provide bam_path and/or features")

    if out_path:
        view.savefig(out_path, dpi=dpi)
    return view


def summary(
    bam_path: str,
    region,
    reference: Optional[Union[str, Reference]] = None,
    min_mapq: int = 0,
) -> dict:
    """Return QC-style statistics for a BAM region.

    Returns a dict with read counts, per-base depth stats, variant/indel and
    junction counts and the insert-size distribution.
    """
    from .reads import (
        compute_coverage,
        compute_insert_sizes,
        fetch_reads,
        junction_counts,
    )

    region = Region.from_any(region)
    reads = fetch_reads(bam_path, region, reference=reference, min_mapq=min_mapq)
    depths, mism = compute_coverage(bam_path, region, reference=reference, min_mapq=min_mapq)
    sizes = compute_insert_sizes(bam_path, region, min_mapq=min_mapq)
    juncs = junction_counts(reads, region)
    n = depths.size
    return {
        "region": region,
        "n_reads": len(reads),
        "n_reverse": sum(1 for r in reads if r.is_reverse),
        "n_mismatches": sum(len(r.mismatches) for r in reads),
        "n_insertions": sum(len(r.insertions) for r in reads),
        "n_deletions": sum(len(r.deletions) for r in reads),
        "n_junctions": sum(juncs.values()),
        "mean_depth": float(depths.mean()) if n else 0.0,
        "max_depth": int(depths.max()) if n else 0,
        "n_variant_positions": int(np.count_nonzero(mism)) if mism is not None else 0,
        "insert_median": float(np.median(sizes)) if sizes.size else 0.0,
        "insert_mean": float(sizes.mean()) if sizes.size else 0.0,
    }

"""Command-line interface for igvplot: ``igvplot BAM REGION [-o OUT] ...``."""
from __future__ import annotations

import argparse
import sys

from . import __version__
from .region import Region


def _parse_sites(bed_path: str) -> dict:
    """Parse a BED3/4 file: 0-based start, optional name -> {pos: label}."""
    sites = {}
    with open(bed_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            try:
                chrom, start, _ = parts[0], int(parts[1]), int(parts[2])
            except (ValueError, IndexError):
                continue
            label = parts[3] if len(parts) > 3 else f"{chrom}:{start}"
            sites[start] = label
    return sites


def _parse_highlight(bed_path: str) -> list:
    """Parse a BED3+ file of regions -> list of (start, end) 0-based intervals."""
    regions = []
    with open(bed_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith(("track", "browser")):
                continue
            parts = line.split("\t")
            try:
                regions.append((int(parts[1]), int(parts[2])))
            except (ValueError, IndexError):
                continue
    return regions


def _parse_variants(vcf_path: str = None, bed_path: str = None) -> list:
    """Parse a VCF or BED of variant sites -> list of (chrom, pos_1based, ref, alt)."""
    variants = []
    if vcf_path:
        with open(vcf_path) as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) < 5:
                    continue
                try:
                    variants.append((parts[0], int(parts[1]), parts[3], parts[4].split(",")[0]))
                except (ValueError, IndexError):
                    continue
    elif bed_path:
        with open(bed_path) as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#") or line.startswith(("track", "browser")):
                    continue
                parts = line.split("\t")
                if len(parts) < 3:
                    continue
                try:
                    start = int(parts[1])
                except ValueError:
                    continue
                ref = parts[3] if len(parts) > 3 else "."
                alt = parts[4] if len(parts) > 4 else "."
                variants.append((parts[0], start + 1, ref, alt))
    return variants


def _parse_basemod(bed_path: str) -> dict:
    """Parse a BED file of base-modification sites -> {pos: (strand, label)}.
    Columns: chrom, start, end, [label], [, strand ('+'/'-')]."""
    mods = {}
    with open(bed_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith(("track", "browser")):
                continue
            parts = line.split("\t")
            try:
                start = int(parts[1])
            except (ValueError, IndexError):
                continue
            label = parts[3] if len(parts) > 3 else "mod"
            strand = 1
            if len(parts) > 5 and parts[5] in ("+", "-"):
                strand = 1 if parts[5] == "+" else -1
            mods[start] = (strand, label)
    return mods


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="igvplot",
        description="IGV-style BAM read pileup, coverage and gene-feature "
        "plots in matplotlib (dna_features_viewer for the gene track).",
    )
    p.add_argument("bam", help="Sorted, indexed BAM/CRAM file")
    p.add_argument("region", nargs="?", default=None,
                   help="Region, e.g. 'chr1:1,000-2,000' (1-based inclusive) or 'chr1:1000' (centered, with --window)")
    p.add_argument("-o", "--out", default=None, help="Output image path (e.g. locus.png)")
    p.add_argument("--vcf", default=None, help="VCF of sites: batch-plot one centered figure per variant")
    p.add_argument("--bed", default=None, help="BED of sites: batch-plot one centered figure per variant")
    p.add_argument("--window", type=int, default=100, help="Flanking window for centered variant plots")
    p.add_argument("--prefix", default="igvplot", help="Output prefix for batch/single plots")
    p.add_argument("--out-format", default=None, help="Output format for batch plots (png/pdf/svg)")
    p.add_argument("--sort-by-variant", action="store_true", help="Sort reads by the variant base at each site")
    p.add_argument("--vaf", action="store_true", help="Annotate the variant allele fraction in the title")
    p.add_argument("--features", default=None, help="Gene/feature annotation (.gb, .gff, .gtf) for the gene track")
    p.add_argument("--reference", "-r", default=None, help="Reference fasta (enables mismatch detection)")
    p.add_argument("--sites", default=None, help="BED file of sites to mark across all tracks")
    p.add_argument("--min-mapq", type=int, default=0, help="Drop reads below this mapping quality")
    p.add_argument("--max-reads", type=int, default=None, help="Cap the number of reads drawn")
    p.add_argument("--no-base-letters", action="store_true", help="Do not paint mismatch base letters")
    p.add_argument("--no-coverage", action="store_true", help="Omit the coverage track")
    p.add_argument("--show-sequence", action="store_true", help="Add a reference sequence row")
    p.add_argument("--bigwig", default=None, help="Use a BigWig file for coverage (avoids pileup on huge BAMs)")
    p.add_argument("--color-by", default="strand", choices=["strand", "firstOfPairStrand", "pairOrientation", "mapq", "insert", "tlen", "unexpectedPair", "readgroup", "proper", "mate", "basemod", "none"], help="How to colour reads")
    p.add_argument("--colormap", default="viridis", help="Matplotlib colormap for continuous colour modes (mapq/insert)")
    p.add_argument("--group-by", default="none", choices=["none", "strand", "readgroup", "proper", "mate", "firstOfPairStrand", "pairOrientation", "insert"], help="Cluster reads vertically by attribute")
    p.add_argument("--link-mates", action="store_true", help="Draw connectors between paired reads")
    p.add_argument("--show-soft-clips", action="store_true", help="Draw soft-clipped ends as thin bars")
    p.add_argument("--sort-by", default="start", choices=["start", "strand", "mapq", "insert", "mate_start", "name"], help="Order reads within the view")
    p.add_argument("--display-mode", default="expanded", choices=["expanded", "squished", "full"], help="Read display mode")
    p.add_argument("--view-as-pairs", action="store_true", help="Place each paired fragment on one row (IGV 'view as pairs')")
    p.add_argument("--show-insertion-text", action="store_true", help="Show inserted-base counts inline")
    p.add_argument("--show-deletion-text", action="store_true", help="Show deletion lengths inline")
    p.add_argument("--highlight-regions", default=None, help="BED file of (start,end) regions to shade across tracks")
    p.add_argument("--basemod", default=None, help="BED file of base-modification sites (col4=label, col6=strand)")
    p.add_argument("--sashimi", action="store_true", help="Add an RNA-seq splice-junction (sashimi) track")
    p.add_argument("--sashimi-min-counts", type=int, default=1, help="Minimum junction support for the sashimi track")
    p.add_argument("--figsize", type=float, nargs=2, default=(14, 8), metavar=("W", "H"), help="Figure size in inches")
    p.add_argument("--dpi", type=int, default=150, help="Output resolution")
    p.add_argument("--font-size", "--fontSize", dest="font_size", type=float, default=None,
                   help="Globally scale all label fonts (points; 8 = default design size)")
    p.add_argument("-V", "--version", action="version", version=f"igvplot {__version__}")
    p.add_argument("--coverage-strand", action="store_true", help="Strand-specific (sense/antisense) coverage")
    p.add_argument("--variants", default=None, help="VCF/BED of variants to mark (REF>ALT)")
    p.add_argument("--gc", action="store_true", help="Add a GC-content track (needs --reference)")
    return p


def _region_from_arg(region_str, window):
    """Accept 'chr:start-end' or a centered 'chr:pos' (with ``window``)."""
    import re

    if region_str is None:
        return None
    m = re.match(r"^(?P<chrom>[A-Za-z0-9_.]+):(?P<pos>\d+)$", region_str.strip().replace(",", ""))
    if m:
        from .region import Region

        return Region.centered(m.group("chrom"), int(m.group("pos")), window)
    from .region import Region

    return Region.from_string(region_str)


def _centered(chrom, pos1, window):

    return Region.centered(chrom, pos1, window)


def sites_from_args(args):
    if not getattr(args, "_sites_cache", None) and args.sites:
        try:
            args._sites_cache = _parse_sites(args.sites)
        except OSError as exc:
            raise SystemExit(f"igvplot: could not read sites file {args.sites!r}: {exc}")
    return getattr(args, "_sites_cache", None)


def _common_kwargs(args, region, out_path, title=None, sites=None, sort_base_pos=None, sort_by=None):
    """Build shared plot_view kwargs for single and batch modes."""

    return dict(
        bam_path=args.bam,
        region=region,
        out_path=out_path,
        features=args.features,
        reference=args.reference,
        sites=sites if sites is not None else sites_from_args(args),
        min_mapq=args.min_mapq,
        max_reads=args.max_reads,
        paint_base_letters=not args.no_base_letters,
        show_coverage=not args.no_coverage,
        show_sequence=args.show_sequence,
        bigwig=args.bigwig,
        sashimi=args.sashimi,
        sashimi_min_counts=args.sashimi_min_counts,
        color_by=args.color_by,
        colormap=args.colormap,
        link_mates=args.link_mates,
        group_by=args.group_by,
        show_soft_clips=args.show_soft_clips,
        display_mode=args.display_mode,
        view_as_pairs=args.view_as_pairs,
        sort_by=sort_by or args.sort_by,
        sort_base_pos=sort_base_pos,
        show_insertion_text=args.show_insertion_text,
        show_deletion_text=args.show_deletion_text,
        highlight_regions=(
            _parse_highlight(args.highlight_regions) if args.highlight_regions else None
        ),
        basemod=basemod_from_args(args),
        coverage_strand=args.coverage_strand,
        variants=args.variants,
        gc=args.gc,
        title=title,
        figsize=tuple(args.figsize),
        dpi=args.dpi,
    )


def basemod_from_args(args):
    if not getattr(args, "_basemod_cache", None) and args.basemod:
        args._basemod_cache = _parse_basemod(args.basemod)
    return getattr(args, "_basemod_cache", None)


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    if getattr(args, "font_size", None):
        from .plot import set_font_size

        set_font_size(args.font_size)
    out_fmt = args.out_format or ("png" if (args.out and "." in args.out) else "png")

    # --- variant / region batch mode -------------------------------------
    if args.vcf or args.bed:
        variants = _parse_variants(args.vcf, args.bed)
        if not variants:
            print("igvplot: no variants parsed from input", file=sys.stderr)
            return 2
        from .reads import variant_allele_fraction
        from .view import plot_view

        for chrom, pos, ref, alt in variants:
            region = _centered(chrom, pos, args.window)
            site_label = f"{ref}> {alt}".replace("> ", ">") if ref and alt else "variant"
            title = f"{chrom}:{pos} {ref}>{alt}" if ref and alt else f"{chrom}:{pos}"
            if args.vaf:
                vaf, depth, altc = variant_allele_fraction(
                    args.bam, chrom, pos - 1, reference=args.reference, min_mapq=args.min_mapq
                )
                title += f"  VAF={vaf:.2f} ({altc}/{depth})"
            out = f"{args.prefix}_{chrom}_{pos}.{out_fmt}"
            plot_view(
                **_common_kwargs(
                    args,
                    region,
                    out,
                    title=title,
                    sites={pos - 1: site_label},
                    sort_base_pos=(pos - 1) if args.sort_by_variant else None,
                    sort_by="base" if args.sort_by_variant else None,
                )
            )
            print(f"wrote {out}")
        return 0

    # --- single region ----------------------------------------------------
    if args.region is None and not args.out:
        # allow --vcf/--bed path where region was omitted
        print("igvplot: a region (or --vcf/--bed) is required", file=sys.stderr)
        build_parser().print_help()
        return 2

    region = _region_from_arg(args.region, args.window) if args.region else args.region
    try:
        kwargs = _common_kwargs(args, region, args.out)
        from .view import plot_view

        plot_view(**kwargs)
    except Exception as exc:  # keep CLI errors user-friendly
        print(f"igvplot error: {exc}", file=sys.stderr)
        return 1

    if args.out:
        print(f"wrote {args.out}")
    else:
        import matplotlib.pyplot as plt

        plt.show()
    return 0


if __name__ == "__main__":
    sys.exit(main())

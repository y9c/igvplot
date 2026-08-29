# Changelog

All notable changes to **igvplot** are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/); this project adheres to
[Semantic Versioning](https://semver.org/).

## [0.0.1] - 2026-08-29

First public release.

### Changed
- Visual refresh (modern pyGenomeTracks-style look) across every track:
  gradient-filled coverage/GC/signal areas, a filled modern gene bar (recoloured
  from dna_features_viewer's defaults and sized to its track), horizontal track
  labels, a slim in-figure legend band (no more floating outside legend),
  subtler dotted site markers with halo'd labels, softer categorical palette
  for read-group colouring, hairline read edges and receding match letters in
  the base-resolution view.

### Fixed
- Gene/feature track no longer renders as a thin bar in a mostly-empty axis:
  dna_features_viewer reserves a tall annotation canvas by default; the axis is
  now clamped to the feature stack and the bar scaled to fill it.

### Added
- Single-cell / long-read colouring: `color_by="tag:CB"` / `group_by="tag:UB"`
  — categorical read colouring and clustering by any auxiliary BAM tag
  (simple scalar tags auto-collected, or pass `tag_keys` explicitly).
- Input formats: junction BED **or STAR `SJ.out.tab`** (auto-detected) for
  `.add_junctions_bed`; BEDPE arcs via the new `.add_bedpe`; pandas/polars
  Series accepted anywhere an array is.
- `log=True` for `.add_coverage` / `.add_signal` (log1p scaling for
  heavy-tailed coverage and bigwig-style signals).
- Hi-C heatmap (`.add_hic`), TAD boundaries (`.add_tads`), scale bar, BED
  annotations, highlight regions.
- Multi-sample comparison tracks: `.add_reads_overlay`, `.add_sashimi_overlay`
  (ggsashimi-style), `.add_coverage_strands` (strand-specific coverage),
  `.add_junctions_bed` (precomputed junctions), `add_coverage(strand=...)`.
- Generic/utility tracks: `.add_signal` (1-D signal), `.add_gc` (GC content),
  `.add_variants` (VCF/BED variants), `.add_arc` (interaction/loop arcs),
  `.add_track` (custom callback).
- Epitranscriptomics tracks: `.add_variant_fraction` (per-base VAF),
  `.add_mod_fraction` (per-site modification stoichiometry),
  `.add_motifs` (IUPAC motif scan, e.g. `DRACH`).
- `igvplot.summary(bam, region, reference)` → region QC statistics.
- `python -m igvplot`, `--version`/`-V`, and CLI flags `--coverage-strand`,
  `--variants`, `--gc`; `docs/api.md` + `Makefile` + `.editorconfig`.
- `draw_interaction_arc` primitive and `add_coverage(strand=...)` pileup filter.

### Fixed
- Reference lifecycle: `add_coverage` propagates a reference to later tracks and
  closes a replaced held `Reference`; `add_sequence` closes the fasta it opens.
- `GenomeView.close()` releases a held reference (called from `savefig`).
- Deletions merged into contiguous spans (correct `show_deletion_text` lengths).
- `add_coverage_overlay` no longer scans the BAM twice per sample.
- `_pack_pairs_into_rows` dedup (pairs processed once), IGV "view as pairs".
- `coverage_from_bedgraph` accepts interval tuples; `pathlib.Path` supported
  everywhere; centred `chr:pos` regions; degenerate insert-colormap guard.
- README converted to an image-first gallery; CI adds a ruff lint job and runs
  tests with warnings-as-errors.

### Changed
- Removed the non-track QC histograms (`insert_size_histogram`, `read_length_histogram`
  and the `--insert-length` / `--read-length` CLI modes); the insert-size statistic
  remains available via `igvplot.summary`.

## [0.1.0] — 2026

### Added
- IGV-style **read pileup** with per-base mismatch letters, `+` insertions and
  red deletion bars, packed into rows and clipped to the region.
- Per-base **coverage** track with the IGV "variant pile" (red ticks where reads
  mismatch the reference).
- **Gene / feature** track rendered by `dna_features_viewer` (GenBank / GFF /
  GTF / `GraphicRecord`), cropped to the region on a shared x-axis.
- **Reference sequence** row and **site markers** with auto-stacked labels.
- **Sashimi** (splice-junction) arcs: thickness ∝ supporting reads, height ∝
  intron length, with read-count labels.
- **Color-by / group-by** matching igv.js (`strand`, `pairOrientation`,
  `readgroup`, `mapq`, `insert`, `proper`, `mate`, `basemod`, …), `link_mates`,
  `show_soft_clips`, `display_mode` (`expanded`/`squished`/`full`) and `sort_by`.
- **Base-resolution** "show all bases" zoom with a colour-coded reference row.
- **Base modifications** (m6A / m5C / …) as strand-aware markers, and
  `color_by="basemod"` to colour reads by the modification they span.
- **Hi-C** contact heatmap, **TAD** boundaries, **scale bar** and **BED/narrowPeak**
  annotation tracks.
- Multi-sample **coverage overlay** and **highlight regions**.
- BigWig coverage via `pyBigWig` or the UCSC `bigWigToBedGraph` fallback.
- `plot_view()` one-call API, a fluent `GenomeView` / `IGV`/`AlignmentView`
  object builder, and a CLI (`single`, variant/`--vcf`/`--bed` batch,
  `--insert-length`).
- Global font scaling (`igvplot.set_font_size` / CLI `--fontSize`).

### Changed
- 2026: README converted to an image-first gallery (see `examples/generate_gallery.py`).

[0.1.0]: https://github.com/yecheng/igvplot/releases/tag/v0.1.0

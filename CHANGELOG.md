# Changelog

All notable changes to **igvplot** are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/); this project adheres to
[Semantic Versioning](https://semver.org/).

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

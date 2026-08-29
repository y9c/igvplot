# igvplot API reference

`igvplot` renders IGV-style genomic tracks in matplotlib on shared x-axes. All
coordinates are **0-based, half-open** internally; human regions (`chr1:1,000-2,000`)
are 1-based inclusive.

## Quick start

```python
import igvplot
from igvplot import IGV, GenomeView, plot_view, summary

# one-call
view = plot_view(bam_path="sample.bam", region="chr1:1,000-2,000",
                 features="annotation.gb", reference="genome.fa",
                 sites={1050: "m6A"}, sashimi=True, out_path="locus.png")

# fluent builder (IGV == GenomeView == AlignmentView)
igv = IGV("chr1:1,000-2,000", reference="genome.fa")   # dpi=120, figsize=(14,7)
igv.add_coverage("rna.bam").add_reads("rna.bam", color_by="readgroup",
                 group_by="pairOrientation", view_as_pairs=True, link_mates=True)
igv.add_sashimi("rna.bam").add_features("annotation.gb", min_feature_length=3)
igv.add_sites({1050: "m6A"}).add_base_mods({1050: (1, "m6A")})
igv.savefig("locus.png")          # or igv.plot()/igv.figure -> Figure; igv.show()
```

## Regions

```python
from igvplot import Region
Region("chr1", 999, 2000)                     # 0-based, half-open
Region.from_string("chr1:1,000-2,000")        # 1-based inclusive
Region.from_any("chr1:1,000-2,000")           # str | (chrom,start,end) | Region
Region.from_string("chr1:5000", window=100)   # centred 'chr:pos' (window each side)
Region.centered("chr1", 5000, window=100)
r.chrom, r.start, r.end, r.length
```


## Track methods (`GenomeView` / `IGV`)

Each returns `self` for chaining; call order = top-to-bottom stack order.

### Alignment
| Method | Notes |
| --- | --- |
| `add_reads(bam_path=None, reads=None, reference=None, min_mapq=0, max_reads=None, keep_duplicates=False, keep_secondary=False, paint_base_letters=True, base_fontsize=11.0, color_by="strand", colormap="viridis", link_mates=False, group_by="none", show_soft_clips=False, display_mode="expanded", view_as_pairs=False, highlight=None, sampling_window=0, max_per_window=0, show_all_bases=None, sort_by="start", show_insertion_text=False, show_deletion_text=False, basemod_sites=None, sort_base_pos=None, mismatch_colors=None, weight=3.5)` | stacked read pileup. `mismatch_colors` colours mismatch letters by substitution type (`"A>G"`) or alternate letter, resolved against `reference` when available |
| `add_reads_overlay(samples, reference=None, min_mapq=0, max_reads=None, weight=3.0, label="reads", paint_base_letters=True)` | `samples=[(bam,color,label),...]`, per-sample colour |
| `bam(bam_path, reference=None, coverage=True, reads=True, **read_kwargs)` | convenience: coverage + reads |

`color_by`: `strand` · `firstOfPairStrand` · `pairOrientation` · `mapq` · `insert`/`tlen` · `unexpectedPair` · `readgroup` · `proper` · `mate` · `basemod` · `none`.
`group_by`: `none` · `strand` · `readgroup` · `proper` · `mate` · `firstOfPairStrand` · `pairOrientation` · `insert`.
`display_mode`: `expanded` · `squished` · `full`. `sort_by`: `start` · `strand` · `mapq` · `insert` · `mate_start` · `name`.

### Coverage
| Method | Notes |
| --- | --- |
| `add_coverage(bam_path=None, bigwig=None, depths=None, reference=None, min_mapq=0, strand=None, weight=1.0, fill_color=..., ylabel="depth", ymax=None)` | per-base depth (+variant pile) |
| `add_coverage_overlay(samples, reference=None, min_mapq=0, weight=1.2, ymax=None, label="coverage")` | overlaid areas `[(bam,color,label),...]` |
| `add_coverage_strands(bam_path, reference=None, min_mapq=0, weight=1.4, forward_color=..., reverse_color=..., ylabel="depth", ymax=None)` | sense/antisense |

### Splicing
| Method | Notes |
| --- | --- |
| `add_sashimi(bam_path=None, reads=None, min_counts=1, arc_color=..., weight=1.2)` | junction arcs from CIGAR `N` |
| `add_sashimi_overlay(samples, reference=None, min_counts=1, weight=1.6)` | comparison, `[(bam,color,label),...]` |
| `add_junctions_bed(bed_path, arc_color=..., weight=1.2, label="junction")` | precomputed junctions (no BAM) |

### Annotation / sequence / signal
| Method | Notes |
| --- | --- |
| `add_features(source=None, weight=1.2, plot_kwargs=None, feature_types=None, min_feature_length=0)` | genes via `dna_features_viewer` |
| `add_sequence(reference, weight=0.5)` | reference bases |
| `add_gc(reference=None, window=50, step=10, color=..., ylabel="GC%", weight=0.8)` | GC content |
| `add_signal(values, color=..., ylabel="signal", ymax=None, weight=1.0)` | numeric array, length = region.length |
| `add_variants(source, color=...)` | mark VCF/BED variants |
| `add_variant_fraction(bam_path, reference=None, min_mapq=0, color=..., weight=1.0, label="VAF")` | per-base allelic fraction |
| `add_mod_fraction(sites, color=..., ymax=1.0, weight=1.0, label="mod%")` | per-site stoichiometry bars |
| `add_motifs(motif, reference=None, color=..., weight=0.7, label="motif", max_hits=500)` | IUPAC motif scan (e.g. `DRACH`) |
| `add_arc(pairs, color=..., weight=2.0, label="interaction", max_height=1.0)` | interaction/loop arcs `(start,end[,strength])` |
| `add_bed_features(bed_path, weight=1.0, color=..., feature_type="annotation")` | BED/narrowPeak |

### Hi-C / domains / markers
| Method | Notes |
| --- | --- |
| `add_hic(matrix, cmap="Reds", vmin=None, vmax=None, weight=3.0, ylabel="Hi-C contact", interpolation="nearest")` | dense `(n,n)` contact map |
| `add_tads(boundaries, color=..., label="TAD")` | domain boundaries |
| `add_scale_bar(window_bp=None, label=None, weight=0.3)` | scale bar |
| `add_sites(sites)` | `{0-based pos: label}` markers across all tracks |
| `add_base_mods(sites, weight=0.5)` | strand-aware mod markers |
| `add_highlight_regions(regions, color=..., alpha=0.3)` | shade `(start,end[,color,alpha])` |
| `add_blank(weight=0.3, label=None)` | spacer |
| `add_track(callback, weight=1.0)` | fully custom `callback(ax, region)` |

### Figure controls
`set_title(text)` · `set_legend(bool)` · `add_legend_items(items)` · `figure` (property) ·
`plot()` / `render(fig=None)` → `Figure` · `savefig(path, dpi=None, **kw)` · `show()` ·
`close()` · `repr`.

## Reading data & stats

```python
from igvplot import fetch_reads, compute_coverage, compute_insert_sizes, variant_allele_fraction, summary
reads   = fetch_reads("sample.bam", region, reference="genome.fa", min_mapq=20, max_reads=500)
depths, mism = compute_coverage("sample.bam", region, reference="genome.fa", strand="forward")
sizes   = compute_insert_sizes("sample.bam", region)
vaf, depth, alt = variant_allele_fraction("sample.bam", "chr1", 1050, reference="genome.fa")
stats   = summary("sample.bam", region, reference="genome.fa")
# stats: region, n_reads, n_reverse, n_mismatches, n_insertions, n_deletions,
#        n_junctions, mean_depth, max_depth, n_variant_positions,
#        insert_median, insert_mean
```

`Read` dataclass fields: `name, aleft, aright, is_reverse, mapq, mismatches, bases,
insertions, deletions, junctions, paired, properly_paired, pairend_first, pairend_second,
mate_is_reverse, mate_chrom, mate_start, insert_size, read_group, clip_left, clip_right`.

## Styling

```python
from igvplot import set_font_size, apply_theme, BASE_COLORS
set_font_size(12)      # scale every label (8 = default design)
apply_theme()          # apply modern matplotlib defaults (auto on render)
# palette lives in igvplot/theme.py: COVERAGE, STRAND_FORWARD, STRAND_REVERSE,
# SITE, SASHIMI, BASE_COLORS, BASE_MOD_COLORS, ...
```

## Command line

```bash
igvplot sample.bam chr1:1,000-2,000 -o locus.png \
    --features annotation.gb -r genome.fa --sites sites.bed \
    --sashimi --link-mates --view-as-pairs --color-by pairOrientation --group-by strand \
    --coverage-strand --variants variants.vcf --gc --show-sequence
igvplot sample.bam --vcf variants.vcf -r genome.fa --window 100 --sort-by-variant --vaf
igvplot --version
```

## Extras
- `pyBigWig`/`bigWigToBedGraph` BigWig coverage (`add_coverage(bigwig=...)`).
- `feature_types` / `min_feature_length` on `add_features` filter the gene track
  (for GenBank, types collapse to `'feature'`, so prefer `min_feature_length`).

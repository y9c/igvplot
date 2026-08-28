# igvplot

IGV-style **BAM read pileup**, **coverage**, **gene features**, **sashimi
junctions** and **epigenetic base modifications** in matplotlib — built on
[pysam](https://github.com/pysam-developers/pysam) and
[dna_features_viewer](https://github.com/Edinburgh-Genome-Foundry/DnaFeaturesViewer).

Reads, coverage and genes are stacked on **shared x-axes** so genomic positions
line up exactly, with **per-base mismatch / insertion / deletion** detail, IGV
colouring, and base-resolution zoom.

---

## The default view

Sashimi arcs → gene arrows → coverage (+ variant pile) → aligned reads → reference sequence.

![default view](examples/gallery_hero.png)

## Base-level zoom

Every read base as a letter, aligned to a colour-coded reference row — read
SNPs, insertions (`+`) and deletions (`−`) position by position.

![base-level zoom](examples/gallery_base_level.png)

## Colouring & grouping (IGV parity)

Colour and cluster reads by attribute via `color_by` / `group_by`.
`link_mates` draws a connector between the two mates; `view_as_pairs` places
each paired fragment on a single joined row (IGV "view as pairs").

| pair orientation + strand | read group | mapping quality |
| --- | --- | --- |
| ![pair orientation](examples/gallery_color_pair.png) | ![read group](examples/gallery_color_readgroup.png) | ![mapq](examples/gallery_color_mapq.png) |

## Epigenetic base modifications

Strand-aware m6A / m5C / … markers in a dedicated track (`add_base_mods`). Reads
can be coloured by the modification they span (`color_by="basemod"`).

![base modifications](examples/gallery_basemod.png)

## Multi-sample comparison

Compare several samples on shared axes — coverage overlay, **comparison sashimi**
(`add_sashimi_overlay`), **strand-specific coverage** (`add_coverage_strands`),
and **read overlay** (`add_reads_overlay`) — then shade regions of interest.

![multi-sample overlay](examples/gallery_overlay.png)

![multi-sample comparison](examples/gallery_compare.png)

## Hi-C multi-track (pyGenomeTracks-style)

Scale bar + contact heatmap + TAD boundaries + genes + coverage + reads —
`GenomeView` track API.

![Hi-C multi-track](examples/gallery_hic.png)

## Variant-centric + allele fraction

One centred plot per variant, with the **variant allele fraction** in the title
and the base change called out at the site (plus any base modifications).

![variant VAF](examples/gallery_variants.png)

---

## Install

```bash
pip install -e .
# optional: GFF/GTF feature support
pip install -e ".[gff]"
```

## Python API

One-liner:

```python
import igvplot

view = igvplot.plot_view(
    bam_path="sample.bam",
    region="chr1:1,000-2,000",        # 1-based inclusive
    features="annotation.gb",
    reference="genome.fa",
    sites={1050: "m6A"},
    sashimi=True,
    link_mates=True,
    out_path="locus.png",             # fig = view.render() to keep it in memory
)
```

Or stack tracks fluently with `GenomeView` (a.k.a. `IGV`):

```python
from igvplot import IGV

igv = (
    IGV("chr1:1,000-2,000", reference="genome.fa", dpi=150)
    .add_sashimi("rna.bam")
    .add_coverage("rna.bam")
    .add_reads("rna.bam", color_by="pairOrientation", group_by="strand",
               link_mates=True, show_soft_clips=True)
    .add_features("annotation.gb")
    .add_sequence("genome.fa")
    .add_sites({1050: "m6A"})
    .add_base_mods({1050: (1, "m6A")})
    .add_highlight_regions([(1050, 1080)])
)
igv.savefig("locus.png")
```

## Command line

```bash
igvplot sample.bam chr1:1,000-2,000 -o locus.png \
    --features annotation.gb -r genome.fa \
    --sites sites.bed --sashimi --link-mates \
    --color-by pairOrientation --group-by strand --show-sequence
```

Variant batch (one figure per VCF/BED site, with VAF):

```bash
igvplot sample.bam --vcf variants.vcf -r genome.fa --window 100 \
    --prefix variant --sort-by-variant --vaf
```

## Regenerate the gallery

```bash
python examples/generate_gallery.py
```

## Reference

| | |
| --- | --- |
| **Region** | `chr1:1,000-2,000` (1-based inclusive), `(chrom, start, end)` tuple, or `Region`. Internal coords are **0-based, half-open**. |
| **Read modes** | `color_by`: `strand` · `pairOrientation` · `readgroup` · `mapq` · `insert` · `proper` · `mate` · `basemod` · `none`. `group_by` cluster reads vertically; `sort_by`: `start` · `strand` · `mapq` · `insert` · `name`. |
| **Display** | `display_mode`: `expanded` · `squished` · `full`. `show_all_bases` for base resolution, `show_soft_clips`, `link_mates`, `view_as_pairs` (IGV "view as pairs"), `highlight={...}`. |
| **Tracks** | `.add_reads` · `.add_reads_overlay` · `.add_coverage` · `.add_coverage_strands` · `.add_coverage_overlay` · `.add_sashimi` · `.add_sashimi_overlay` · `.add_junctions_bed` · `.add_signal` · `.add_gc` · `.add_variants` · `.add_variant_fraction` · `.add_mod_fraction` · `.add_motifs` · `.add_arc` · `.add_track` · `.add_features` · `.add_sequence` · `.add_bed_features` · `.add_hic` · `.add_tads` · `.add_scale_bar` · `.add_sites` · `.add_base_mods` · `.add_highlight_regions` |
| **Stats** | `igvplot.summary(bam, region, reference)` → reads, depth, variant/indel/junction counts, insert-size median. |

![signal / GC / loops / variants / custom](examples/gallery_tracks.png)

![per-base VAF / mod stoichiometry / motifs](examples/gallery_epigenetics.png)

GC, signal, interaction arcs, variants & custom tracks — and the epitranscriptomics
tracks: per-base variant fraction, modification stoichiometry and IUPAC motifs
(`DRACH` for the m6A context).
| **Styling** | `set_font_size(n)` / `--fontSize n` scales every label; `set_title`, `set_legend`, `figsize`, `dpi`. |

MIT — see [LICENSE](LICENSE).

---

## Development & publishing

```bash
pip install -e ".[dev]"        # build, pytest, pytest-cov
python -m pytest -q            # run the test suite
python examples/generate_gallery.py   # regenerate the images above
# or: make test / make lint / make gallery / make build
```

- [docs/api.md](docs/api.md) — full API reference (regions, every track method, stats, CLI).
- [CONTRIBUTING.md](CONTRIBUTING.md) — setup, style, module map.
- [CHANGELOG.md](CHANGELOG.md) — release history.
- CI runs tests on Python 3.9 & 3.11 and builds a wheel on every push
  (`.github/workflows/ci.yml`).

### Releases

```bash
python scripts/bump_version.py --patch   # bump to v0.1.1 (update pyproject + __init__)
python scripts/bump_version.py --patch --tag   # bump, commit, tag v0.1.1
git push && git push --tags             # tag triggers .github/workflows/release.yml
```

Pushing a `v*` tag builds the wheel + sdist, verifies the wheel imports, and
uploads the artifacts to a GitHub Release.

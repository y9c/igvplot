# Contributing

Thanks for helping with **igvplot**! This is a small, focused package: keep
changes consistent with its structure and naming, and always back them with a
test.

## Setup

```bash
git clone <repo-url> igvplot
cd igvplot
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"
```

The package needs no external data to run tests — a small synthetic BAM/fasta
lives in `data/` (regenerate with `python scripts/make_synthetic_bam.py`).

## Tests & gallery

```bash
python -m pytest -q                 # run the suite (uses a headless Agg backend)
python examples/generate_gallery.py # regenerate the README example images
```

## Code style & conventions

- **Coordinates** are always **0-based, half-open** internally (pysam
  convention); human region strings (`chr1:1,000-2,000`) are 1-based inclusive
  and converted on parse.
- **One vocabulary**: the single builder class is `GenomeView` (aliases `IGV`
  and `AlignmentView`); all track adders use the `add_*` prefix, plus the
  composite `.bam()`.
- **Data integrity**: only use real values from primary sources (BAM/FASTA
  reads, index coordinates); never invent numbers.
- Public API lives in `igvplot/__init__.py` — export anything users need there.
- Format with [ruff](https://docs.astral.sh/ruff/) (`line-length = 100`) and
  keep the public surface small and discoverable.

## Where things live

| Area | Module |
| ---- | ------ |
| Region / coordinates | `igvplot/region.py` |
| Reads, coverage, junctions, insert sizes | `igvplot/reads.py` |
| Gene/feature loading | `igvplot/features.py` |
| BigWig / bedGraph | `igvplot/bigwig.py` |
| Low-level drawing primitives | `igvplot/plot.py` |
| Track builder `GenomeView` + `plot_view` | `igvplot/view.py` |
| Command-line interface | `igvplot/cli.py` |
| Synthetic test data generator | `scripts/make_synthetic_bam.py` |

## Before you open a PR

1. `python -m pytest -q` passes.
2. Add a test for any new behaviour.
3. If a feature is user-facing, add it to the README's reference table / gallery.
4. Update `CHANGELOG.md`.

#!/usr/bin/env python3
"""Render a pyGenomeTracks-style Hi-C multi-track figure from synthetic data.

Top-to-bottom: scale bar, Hi-C contact heatmap, TAD boundaries, gene track,
and an alignment/coverage view. Run: python examples/hic_demo.py
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, REPO)

from igvplot import GenomeView, Region  # noqa: E402

DATA = os.path.join(REPO, "data")
REGION = "chrTest:6,930-7,230"


def make_hic_matrix(region, n_bins=40, tads=((6930, 7020), (7020, 7120), (7120, 7230))):
    """Build a block-structured, decaying-contact Hi-C matrix for the region."""
    start, end = region.start, region.end
    edges = np.linspace(start, end, n_bins + 1)
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


def main():
    region = Region.from_any(REGION)
    mat = make_hic_matrix(region)
    boundaries = [7020, 7120]

    gv = GenomeView(region=REGION, figsize=(14, 8), dpi=120)
    gv.add_scale_bar(window_bp=100)
    gv.add_hic(mat, cmap="Reds")
    gv.add_tads(boundaries)
    gv.add_features(os.path.join(DATA, "annotation.gb"))
    gv.add_coverage(os.path.join(DATA, "sample.bam"), reference=os.path.join(DATA, "genome.fa"))
    gv.add_reads(os.path.join(DATA, "sample.bam"), reference=os.path.join(DATA, "genome.fa"),
                 color_by="readgroup", group_by="pairOrientation")
    gv.savefig(os.path.join(HERE, "hic_demo.png"), dpi=120)
    print("saved examples/hic_demo.png")


if __name__ == "__main__":
    main()

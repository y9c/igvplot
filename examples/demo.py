#!/usr/bin/env python3
"""Render the igvplot demo figure from the synthetic data.

Run from the repo root:  python examples/demo.py
Produces: examples/locus.png (reads + coverage + gene features + sites + ref)
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, REPO)

import igvplot  # noqa: E402

DATA = os.path.join(REPO, "data")


def main():
    igvplot.plot_view(
        bam_path=os.path.join(DATA, "sample.bam"),
        region="chrTest:6,930-7,200",
        features=os.path.join(DATA, "annotation.gb"),
        reference=os.path.join(DATA, "genome.fa"),
        sites={
            7000: "SNP C>A",
            7010: "SNP",
            7020: "DEL 3 bp",
            7030: "INS CC",
        },
        sashimi=True,
        link_mates=True,
        color_by="strand",
        show_sequence=True,
        out_path=os.path.join(HERE, "locus.png"),
        dpi=130,
        figsize=(14, 8),
    )
    print("saved examples/locus.png")


if __name__ == "__main__":
    main()

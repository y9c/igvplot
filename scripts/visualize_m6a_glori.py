#!/usr/bin/env python3
"""Visualize m6A GLORI sequencing data.

Creates a comprehensive plot showing:
- Gene structure (reverse strand with 3 exons)
- Spliced reads from both strands
- A→G mutations at m6A sites (GLORI conversion)
- m6A site annotations

Output: data/m6a/m6a_glori_visualization.png
"""
from __future__ import annotations

import os
import sys

import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..'))

from igvplot import plot_view  # noqa: E402

DATA = os.path.join(HERE, "..", "data", "m6a")
OUTPUT = os.path.join(DATA, "m6a_glori_visualization.png")

# View the entire transcript region
REGION = "chrM6A:150-1150"

# m6A sites to highlight
M6A_SITES = [250, 350, 450, 750, 800, 1000, 1050]


def main():
    # Load data
    bam_path = os.path.join(DATA, "m6a_reads.bam")
    annotation_path = os.path.join(DATA, "m6a_annotation.gb")
    reference_path = os.path.join(DATA, "m6a_genome.fa")
    
    # Create visualization
    view = plot_view(
        bam_path=bam_path,
        annotation_path=annotation_path,
        reference_path=reference_path,
        region=REGION,
        title="m6A GLORI Sequencing (A→G mutations, reverse strand gene)",
        color_by="strand",
        show_coverage=True,
        show_junctions=True,
        min_mapq=0,
        max_reads=500,  # Show all reads
        figsize=(16, 10),
    )
    
    # Render to get the figure
    fig = view.render()
    
    # Highlight m6A sites with vertical lines on all axes
    for ax in fig.axes:
        for site in M6A_SITES:
            ax.axvline(x=site, color='purple', alpha=0.5, linestyle='--', linewidth=1)
    
    # Add legend for m6A sites
    fig.axes[0].axvline(x=M6A_SITES[0], color='purple', alpha=0.5, 
                         linestyle='--', linewidth=1, label='m6A sites (A→G)')
    fig.axes[0].legend(loc='upper right', fontsize=9, framealpha=0.8)
    
    # Save
    fig.savefig(OUTPUT, dpi=150, bbox_inches='tight')
    print(f"Visualization saved to: {OUTPUT}")
    
    # Also create a zoomed view of one m6A site
    zoom_output = os.path.join(DATA, "m6a_glori_zoom.png")
    zoom_view = plot_view(
        bam_path=bam_path,
        annotation_path=annotation_path,
        reference_path=reference_path,
        region="chrM6A:230-270",  # Zoom on first m6A site at 250
        title="Zoom on m6A site at position 250 (A→G mutation)",
        color_by="strand",
        show_coverage=True,
        show_junctions=True,
        min_mapq=0,
        max_reads=200,
        figsize=(12, 8),
    )
    
    # Render zoomed view
    zoom_fig = zoom_view.render()
    
    # Add vertical line at m6A site
    for ax in zoom_fig.axes:
        ax.axvline(x=250, color='purple', alpha=0.7, 
                   linestyle='--', linewidth=2, label='m6A site')
    zoom_fig.axes[0].legend(loc='upper right', fontsize=10, framealpha=0.9)
    
    zoom_fig.savefig(zoom_output, dpi=150, bbox_inches='tight')
    print(f"Zoomed visualization saved to: {zoom_output}")


if __name__ == "__main__":
    main()

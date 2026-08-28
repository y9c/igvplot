"""Low-level matplotlib drawing primitives for the IGV-style tracks.

Each function draws one track onto a given axis in **global** reference
coordinates and returns the number of visual rows/lines used (used by the
view to auto-size axes).
"""
from __future__ import annotations

from typing import Dict, List, Optional, Tuple, Union

import numpy as np
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import matplotlib.path as mpath
import matplotlib.patches as mpatches

from .reads import BASE_COLORS, Read, junction_counts
from .region import Region
from .theme import (
    BASE_MOD_COLORS,
    COVERAGE,
    COVERAGE_MISMATCH,
    DELETION,
    INSERTION,
    REF_BASE,
    SASHIMI,
    SITE,
    STRAND_FORWARD,
    STRAND_REVERSE,
    STRAND_UNKNOWN,
)

__all__ = [
    "BASE_COLORS",
    "STRAND_COLORS",
    "COLOR_SCHEMES",
    "pack_reads_into_rows",
    "draw_read_track",
    "draw_coverage_track",
    "draw_sequence_track",
    "draw_sashimi_track",
    "draw_sites",
]

# IGV-style: forward reads light, reverse reads dark (modern palette).
STRAND_COLORS = {
    "forward": STRAND_FORWARD,
    "reverse": STRAND_REVERSE,
    "unknown": STRAND_UNKNOWN,
}
MISMATCH_INS_COLOR = INSERTION
DELETION_COLOR = DELETION

# Per-category colour palettes for discrete colouring (read-group, proper-pair).
COLOR_SCHEMES = {
    "readgroup": list(mcolors.TABLEAU_COLORS)
    + ["#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"],
    "proper": {"proper": "#1f77b4", "improper": "#d62728", "unpaired": "#9e9e9e"},
    "mate": {"connected": "#1f77b4", "unconnected": "#c8c8c8"},
}

# --------------------------------------------------------------------------- #
# Global font sizing
# --------------------------------------------------------------------------- #
# All rendered text is sized relative to this design base, so a single
# ``set_font_size(12)`` (or ``--fontSize 12``) scales *every* track's fonts
# proportionally instead of editing each draw call.
_BASE_DESIGN_FONT = 8.0
_FONT_SCALE = 1.0


def set_font_size(size) -> None:
    """Globally scale every rendered label by ``size`` points.

    ``size`` is the new base (in points) for the currently "8.0"-design text;
    a value of 8 leaves the figure unchanged, 12 makes everything 1.5x larger.
    Affects reads, coverage, features, sashimi, sequence, site and axis labels
    across all tracks. Call before ``render``/``savefig``.
    """
    global _FONT_SCALE
    _FONT_SCALE = float(size) / _BASE_DESIGN_FONT


def _fs(size):
    """Scale a design font size by the global font factor."""
    return size * _FONT_SCALE


def _rgb_hex(color: Tuple[float, float, float]) -> str:
    return "#{:02x}{:02x}{:02x}".format(
        *(int(round(c * 255)) for c in color)
    )


def pack_reads_into_rows(
    reads: List[Read],
    region: Region,
    group_keys: Optional[List] = None,
    sort_keys: Optional[List] = None,
) -> Dict[int, int]:
    """Greedily assign each read to a row so that reads packed in the same row
    do not overlap on the reference. Returns ``{read_index: row_index}``.

    Reads are processed left-to-right (or by ``sort_keys``); each is placed in
    the lowest row whose rightmost occupied position ends at or before the
    read's (clipped) start. If ``group_keys`` is given, reads are processed
    grouped by that key (in order of first appearance) so that reads sharing a
    key cluster together vertically.
    """
    clipped = []
    for r in reads:
        cs = max(region.start, r.aleft)
        ce = min(region.end, r.aright)
        clipped.append((cs, ce))
    valid = [i for i, (s, e) in enumerate(clipped) if e - s >= 1]

    if group_keys is None:
        group_keys = [0] * len(reads)
    if sort_keys is None:
        sort_keys = [clipped[i][0] for i in range(len(reads))]
    labels: List = []
    for g in group_keys:
        if g not in labels:
            labels.append(g)
    rank = {l: i for i, l in enumerate(labels)}

    order = sorted(
        valid,
        key=lambda i: (rank[group_keys[i]], sort_keys[i], clipped[i][0]),
    )
    row_right_end: List[int] = []
    mapping: Dict[int, int] = {}
    for i in order:
        s, e = clipped[i]
        placed = False
        for row, right_end in enumerate(row_right_end):
            if right_end <= s:
                row_right_end[row] = max(right_end, e)
                mapping[i] = row
                placed = True
                break
        if not placed:
            row_right_end.append(e)
            mapping[i] = len(row_right_end) - 1
    return mapping


def _first_pair_strand_key(r: Read) -> str:
    """Return 'forward'/'reverse' for the strand of the first-of-pair mate."""
    if r.pairend_first:
        return "forward" if not r.is_reverse else "reverse"
    return "forward" if not r.mate_is_reverse else "reverse"


def _pair_orientation(r: Read) -> str:
    """Return 'FF'/'FR'/'RF'/'RR'/'unpaired' describing mate orientation."""
    if not r.paired:
        return "unpaired"
    own, mate = r.is_reverse, r.mate_is_reverse
    if not own and not mate:
        return "FF"
    if not own and mate:
        return "FR"
    if own and not mate:
        return "RF"
    return "RR"


_PAIR_ORIENTATION_COLORS = {
    "FF": "#1f77b4",
    "FR": "#d62728",
    "RF": "#ff7f0e",
    "RR": "#9467bd",
    "unpaired": "#9e9e9e",
}

# Default colours for common nucleotide modifications come from the theme.


def base_mod_color(label: str) -> str:
    """Map a modification label/type to a colour (by substring match on known
    modifications, else grey)."""
    text = str(label).lower()
    for key, color in BASE_MOD_COLORS.items():
        if key in text:
            return color
    return BASE_MOD_COLORS["default"]


def _read_sort_value(r: Read, sort_by: str, base_pos: Optional[int] = None):
    """Return a numeric/orderable key for ordering reads within a group."""
    if sort_by == "start" or sort_by == "position":
        return (0, r.aleft)
    if sort_by == "strand":
        return (0 if not r.is_reverse else 1, r.aleft)
    if sort_by == "mapq":
        return (-int(r.mapq), r.aleft)
    if sort_by == "insert":
        return (-int(abs(r.insert_size)), r.aleft)
    if sort_by == "mate_start":
        return (r.mate_start if r.mate_start is not None else -1, r.aleft)
    if sort_by == "name":
        return (r.name or "", r.aleft)
    if sort_by == "base":
        # sort reads by whether they carry the (variant) base at base_pos:
        # alt/mismatch reads first, then the rest.
        carry = 0 if (base_pos is not None and base_pos in r.mismatches) else 1
        return (carry, r.aleft)
    raise ValueError(
        f"Unknown sort_by={sort_by!r}; choose from start, strand, mapq, insert, mate_start, name, base"
    )


def _read_group_key(r: Read, group_by: str):
    if group_by == "none":
        return None
    if group_by == "strand":
        return "forward" if not r.is_reverse else "reverse"
    if group_by == "readgroup":
        return r.read_group or "__none__"
    if group_by == "proper":
        if not r.paired:
            return "unpaired"
        return "proper" if r.properly_paired else "improper"
    if group_by == "mate":
        return r.mate_chrom or "__none__"
    if group_by == "firstOfPairStrand":
        return _first_pair_strand_key(r)
    if group_by == "pairOrientation":
        return _pair_orientation(r)
    if group_by == "insert":
        return "insert"
    raise ValueError(
        f"Unknown group_by={group_by!r}; choose from none, strand, readgroup, "
        "proper, mate, firstOfPairStrand, pairOrientation, insert"
    )


def build_legend_items(
    color_by: str,
    reads: List[Read],
    group_by: str = "none",
    basemod_sites: Optional[Dict[int, tuple]] = None,
) -> List[Tuple[str, str]]:
    """Return ``[(label, color), ...]`` for an auto-generated legend of the
    given colour/group mode. Empty for continuous/trivial modes."""
    if color_by == "none":
        # fall back to grouping if colouring is off
        color_by = group_by
    items: List[Tuple[str, str]] = []
    add = items.append

    if color_by in ("strand", "firstOfPairStrand"):
        add(("forward", STRAND_COLORS["forward"]))
        add(("reverse", STRAND_COLORS["reverse"]))
    elif color_by == "pairOrientation":
        orders = ("FF", "FR", "RF", "RR")
        sampled = {_pair_orientation(r) for r in reads}
        for k in orders:
            if k in sampled:
                add((k, _PAIR_ORIENTATION_COLORS[k]))
    elif color_by == "readgroup":
        groups = sorted({r.read_group for r in reads if r.read_group is not None})
        for i, g in enumerate(groups):
            add((g, COLOR_SCHEMES["readgroup"][i % len(COLOR_SCHEMES["readgroup"])]))
    elif color_by == "proper":
        add(("proper", COLOR_SCHEMES["proper"]["proper"]))
        add(("improper", COLOR_SCHEMES["proper"]["improper"]))
        add(("unpaired", COLOR_SCHEMES["proper"]["unpaired"]))
    elif color_by == "mate":
        add(("connected", COLOR_SCHEMES["mate"]["connected"]))
        add(("unconnected", COLOR_SCHEMES["mate"]["unconnected"]))
    elif color_by in ("insert", "tlen"):
        add(("insert < 10th pct", "#2c6fbb"))
        add(("insert > 90th pct", "#c0392b"))
    elif color_by == "unexpectedPair":
        add(("insert too short", "#2c6fbb"))
        add(("insert too long", "#c0392b"))
        add(("improper pair", "#e67e22"))
    elif color_by == "basemod":
        seen = {}
        for pos, vals in (basemod_sites or {}).items():
            label = vals[1] if isinstance(vals, (tuple, list)) else str(vals)
            color = base_mod_color(label)
            if isinstance(vals, (tuple, list)) and len(vals) > 2:
                color = vals[2]
            seen[label] = color
        for label, color in seen.items():
            add((label, color))
    return items

def _read_colormap(color_by: str, colormap: str, reads: List[Read], basemod_sites=None):
    """Return a ``color(read_index) -> hex`` callable for a colour mode.

    Supported ``color_by`` values: 'strand', 'firstOfPairStrand',
    'pairOrientation', 'mapq', 'insert'/'tlen', 'unexpectedPair', 'readgroup',
    'proper', 'mate', 'none'.
    """
    cmap = mpl.colormaps[colormap]

    # IGV-style insert-size thresholds (blue = too short, red = too long).
    sizes = sorted(abs(r.insert_size) for r in reads if r.insert_size > 0)
    if sizes:
        min_tlen = np.percentile(sizes, 10)
        max_tlen = np.percentile(sizes, 90)
    else:
        min_tlen = max_tlen = 0

    def color_by_strand(r):
        key = "forward" if not r.is_reverse else "reverse"
        return STRAND_COLORS.get(key, STRAND_COLORS["unknown"])

    def color_by_first_pair_strand(r):
        key = _first_pair_strand_key(r)
        return STRAND_COLORS.get(key, STRAND_COLORS["unknown"])

    def color_by_pair_orientation(r):
        return _PAIR_ORIENTATION_COLORS[_pair_orientation(r)]

    def color_by_mapq(r):
        norm = mcolors.Normalize(vmin=0, vmax=60)
        return _rgb_hex(cmap(norm(r.mapq)))

    def color_by_insert(r):
        if r.insert_size and max_tlen:
            if r.insert_size < min_tlen:
                return "#2c6fbb"  # too short
            if r.insert_size > max_tlen:
                return "#c0392b"  # too long
            return _rgb_hex(cmap(mcolors.Normalize(vmin=min_tlen, vmax=max_tlen)(r.insert_size)))
        return STRAND_COLORS["unknown"]

    def color_by_unexpected(r):
        if not r.paired:
            return STRAND_COLORS["unknown"]
        if r.insert_size and r.insert_size < min_tlen:
            return "#2c6fbb"
        if r.insert_size and r.insert_size > max_tlen:
            return "#c0392b"
        if not r.properly_paired:
            return "#e67e22"
        return STRAND_COLORS["unknown"]

    def color_by_readgroup(r):
        groups = sorted({r.read_group for r in reads if r.read_group is not None})
        lookup = {
            g: COLOR_SCHEMES["readgroup"][i % len(COLOR_SCHEMES["readgroup"])]
            for i, g in enumerate(groups)
        }
        if r.read_group is not None:
            return lookup[r.read_group]
        return STRAND_COLORS["unknown"]

    def color_by_proper(r):
        if not r.paired:
            return COLOR_SCHEMES["proper"]["unpaired"]
        return COLOR_SCHEMES["proper"]["proper" if r.properly_paired else "improper"]

    def color_by_mate(r):
        if r.paired and r.mate_start is not None and r.mate_chrom is not None:
            return COLOR_SCHEMES["mate"]["connected"]
        return COLOR_SCHEMES["mate"]["unconnected"]

    def color_by_basemod(r):
        # Colour reads that span a modification site by that modification's
        # colour (IGV "color by base mod"): reads overlapping a mod site are
        # highlighted, unmodified-site reads stay grey.
        mods = basemod_sites or {}
        if r.aleft <= r.aright:
            for pos, vals in mods.items():
                if r.aleft <= pos < r.aright:
                    strand = vals[0] if isinstance(vals, (tuple, list)) else 1
                    color = vals[2] if isinstance(vals, (tuple, list)) and len(vals) > 2 else base_mod_color(vals[1] if isinstance(vals, (tuple, list)) else str(vals))
                    return color
        return STRAND_COLORS["unknown"]

    dispatch = {
        "strand": color_by_strand,
        "firstOfPairStrand": color_by_first_pair_strand,
        "pairOrientation": color_by_pair_orientation,
        "mapq": color_by_mapq,
        "insert": color_by_insert,
        "tlen": color_by_insert,
        "unexpectedPair": color_by_unexpected,
        "readgroup": color_by_readgroup,
        "proper": color_by_proper,
        "mate": color_by_mate,
        "basemod": color_by_basemod,
        "none": lambda r: STRAND_COLORS["unknown"],
    }
    if color_by not in dispatch:
        raise ValueError(
            f"Unknown color_by={color_by!r}; choose from "
            + ", ".join(sorted(dispatch))
        )
    return lambda i, r: dispatch[color_by](r)


def draw_read_track(
    ax,
    reads: List[Read],
    region: Region,
    paint_base_letters: bool = True,
    base_fontsize: float = 11.0,
    strand_colors: Optional[Dict[str, str]] = None,
    base_colors: Optional[Dict[str, str]] = None,
    body_alpha: float = 1.0,
    color_by: str = "strand",
    colormap: str = "viridis",
    link_mates: bool = False,
    mate_color: str = "#444444",
    group_by: str = "none",
    show_soft_clips: bool = False,
    soft_clip_color: str = "#7f8c8d",
    deletion_color: str = DELETION_COLOR,
    insertion_color: str = MISMATCH_INS_COLOR,
    display_mode: str = "expanded",
    highlight: Optional[set] = None,
    highlight_color: str = "#e67e22",
    show_all_bases: bool = False,
    sort_by: str = "start",
    show_insertion_text: bool = False,
    show_deletion_text: bool = False,
    basemod_sites: Optional[Dict[int, tuple]] = None,
    sort_base_pos: Optional[int] = None,
) -> int:
    """Draw a stacked pile of aligned reads (IGV-style).

    The read body is a rectangle coloured by ``color_by`` (strand /
    firstOfPairStrand / pairOrientation / mapq / insert or tlen /
    unexpectedPair / readgroup / proper / mate / none). Per-base mismatches are
    drawn as coloured nucleotide letters, insertions as '+', deletions as a red
    connecting line. Reads can be clustered vertically with ``group_by``, and
    soft-clipped ends can be shown. If ``link_mates``, paired reads that are
    both present are joined by a thin line.

    Returns the number of rows used.
    """
    strand_colors = strand_colors or STRAND_COLORS
    base_colors = base_colors or BASE_COLORS
    color_fn = _read_colormap(color_by, colormap, reads, basemod_sites)
    group_keys = [_read_group_key(r, group_by) for r in reads] if group_by != "none" else None

    squished = display_mode == "squished"
    full = display_mode == "full"
    show_letters = paint_base_letters and not squished
    # Opaque reads that fill their row: no translucent seams, no white gaps
    # between rows, so the pile reads as smooth IGV-style bands.
    row_half = 0.5 if not squished else 0.22
    body_alpha = 0.5 if squished else body_alpha

    group_keys = [_read_group_key(r, group_by) for r in reads] if group_by != "none" else None
    sort_keys = [_read_sort_value(r, sort_by, sort_base_pos) for r in reads]
    mapping = pack_reads_into_rows(reads, region, group_keys, sort_keys)
    if full:
        # one read per row, ordered by group then sort key
        labels = []
        if group_keys is not None:
            for g in group_keys:
                if g not in labels:
                    labels.append(g)
            rank = {l: i for i, l in enumerate(labels)}
        else:
            rank = {i: 0 for i in range(len(reads))}
        order = sorted(
            range(len(reads)),
            key=lambda i: (rank[group_keys[i]] if group_keys is not None else 0, sort_keys[i]),
        )
        mapping = {i: pos for pos, i in enumerate(order)}

    clipped = {}
    for idx, r in enumerate(reads):
        clipped[idx] = (
            max(region.start, r.aleft),
            min(region.end, r.aright),
        )

    used_rows = set()

    for idx, r in enumerate(reads):
        row = mapping.get(idx)
        if row is None:
            continue
        used_rows.add(row)
        s, e = clipped[idx]
        if e - s < 1:
            continue

        body_color = color_fn(idx, r)

        # Soft-clipped ends (thin lighter bars just outside the aligned span).
        if show_soft_clips and (r.clip_left or r.clip_right):
            if r.clip_left:
                cs = max(region.start, s - r.clip_left)
                if e - s and cs < s:
                    ax.add_patch(
                        mpatches.Rectangle(
                            (cs, row - row_half * 0.3),
                            s - cs,
                            row_half * 0.6,
                            facecolor=soft_clip_color,
                            edgecolor="none",
                            alpha=0.8,
                            zorder=1,
                        )
                    )
            if r.clip_right:
                ce = min(region.end, e + r.clip_right)
                if ce > e:
                    ax.add_patch(
                        mpatches.Rectangle(
                            (e, row - row_half * 0.3),
                            ce - e,
                            row_half * 0.6,
                            facecolor=soft_clip_color,
                            edgecolor="none",
                            alpha=0.8,
                            zorder=1,
                        )
                    )

        # Read body.
        ax.add_patch(
            mpatches.Rectangle(
                (s, row - row_half),
                e - s,
                2 * row_half,
                facecolor=body_color,
                edgecolor="none",
                alpha=body_alpha,
                zorder=2,
            )
        )

        # Highlight outline for user-selected reads.
        if highlight and r.name in highlight:
            ax.add_patch(
                mpatches.Rectangle(
                    (s, row - row_half),
                    e - s,
                    2 * row_half,
                    facecolor="none",
                    edgecolor=highlight_color,
                    lw=1.8,
                    zorder=5,
                )
            )

        # Deletions: connecting red line within the read.
        for ds, de in r.deletions:
            ds = max(ds, s)
            de = min(de, e)
            if de - ds < 1:
                continue
            ax.plot([ds, de], [row, row], color=deletion_color, lw=2.0, zorder=3)
            if show_deletion_text:
                ax.text(
                    (ds + de) / 2.0,
                    row + row_half + 0.12,
                    f"{de - ds}bp",
                    ha="center",
                    va="bottom",
                    fontsize=_fs(max(5, base_fontsize - 1.5)),
                    color=deletion_color,
                    zorder=6,
                )

        # Insertions: '+' (optionally with the inserted-base count).
        for ipos, ibases in r.insertions.items():
            if s <= ipos < e and show_letters:
                txt = f"+{len(ibases)}" if show_insertion_text else "+"
                ax.text(
                    ipos + 0.5,
                    row,
                    txt,
                    ha="center",
                    va="center",
                    fontsize=_fs(base_fontsize),
                    color=insertion_color,
                    weight="bold",
                    zorder=6,
                )

        # Base-level "show all bases" view: render every aligned base as a
        # letter and leave '−' for deletions, so the reads read like a
        # multiple sequence alignment against the reference row below.
        if show_all_bases:
            del_pos = {
                p
                for (da, db) in r.deletions
                for p in range(max(da, s), min(db, e))
            }
            for p in range(s, e):
                if p in del_pos:
                    ax.text(
                        p + 0.5,
                        row,
                        "−",
                        ha="center",
                        va="center",
                        fontsize=_fs(base_fontsize),
                        color="#9aa5b1",
                        zorder=6,
                    )
                elif p in r.bases:
                    base = r.bases[p]
                    is_mismatch = p in r.mismatches
                    ax.text(
                        p + 0.5,
                        row,
                        base,
                        ha="center",
                        va="center",
                        fontsize=_fs(base_fontsize),
                        color=base_colors.get(base, "#333333") if is_mismatch else "#aeb8c4",
                        weight="bold" if is_mismatch else "normal",
                        zorder=6,
                    )
        elif show_letters:
            # Mismatch-only view: coloured letters at positions differing from
            # the reference.
            for pos, base in sorted(r.mismatches.items()):
                if s <= pos < e:
                    ax.text(
                        pos + 0.5,
                        row,
                        base,
                        ha="center",
                        va="center",
                        fontsize=_fs(base_fontsize),
                        color=base_colors.get(base, "#333333"),
                        weight="bold",
                        zorder=6,
                    )

    # ---- mate-pair links -------------------------------------------------
    if link_mates:
        by_name: Dict[str, List[int]] = {}
        for idx, r in enumerate(reads):
            by_name.setdefault(r.name, []).append(idx)
        done = set()
        for name, idxs in by_name.items():
            if len(idxs) < 2:
                continue
            idxs = [i for i in idxs if i in mapping]
            if len(idxs) < 2:
                continue
            # connect the two reads' nearest ends
            i1, i2 = idxs[0], idxs[1]
            r1, r2 = reads[i1], reads[i2]
            row1, row2 = mapping[i1], mapping[i2]
            s1, e1 = clipped[i1]
            s2, e2 = clipped[i2]
            if (name,) in done:
                continue
            done.add((name,))
            # choose midpoint for a gentle link
            left_pt = (e1, row1) if e1 < e2 else (e2, row2)
            right_pt = (s1, row1) if s1 > s2 else (s2, row2)
            ax.plot(
                [left_pt[0], right_pt[0]],
                [left_pt[1], right_pt[1]],
                color=mate_color,
                lw=0.8,
                alpha=0.7,
                zorder=1,
            )

    n_used = len(used_rows) if used_rows else 0
    if n_used:
        ax.set_ylim(n_used - 0.5, -0.5)  # row 0 on top
    else:
        ax.set_ylim(-0.5, 0.5)

    ax.set_yticks([])
    ax.set_ylabel(f"{n_used} reads", fontsize=_fs(11))
    return n_used


def draw_coverage_track(
    ax,
    depths: np.ndarray,
    region: Region,
    fill_color: str = COVERAGE,
    mismatch_counts: Optional[np.ndarray] = None,
    mismatch_color: str = COVERAGE_MISMATCH,
    ylabel: str = "depth",
    ymax: Optional[float] = None,
) -> float:
    """Draw per-base coverage as a stepped area plot, plus (optionally) red
    tick bars for positions where reads mismatch the reference ('variant
    pile'). ``ymax`` fixes the top of the y-axis (None = autoscale to the max
    depth). Returns the max depth (for shared y-axis auto-scaling)."""
    n = region.length
    x = np.arange(region.start, region.end, dtype=float)

    ax.fill_between(
        x,
        depths,
        step="mid",
        color=fill_color,
        alpha=0.42,
        edgecolor="none",
        linewidth=0,
        zorder=1,
    )
    # crisp stepped bounding line for a cleaner, more modern look
    ax.step(x, depths, where="mid", color=fill_color, lw=1.3, zorder=2)
    if ymax is None:
        ymax = float(depths.max()) if n else 1.0
    ax.set_ylim(0, max(ymax, 1e-9))

    if mismatch_counts is not None and mismatch_counts.any():
        pos = np.nonzero(mismatch_counts)[0]
        # stagger variant ticks so they are visible above the depth line
        heights = mismatch_counts[pos]
        ax.vlines(
            x[pos] + 0.5,
            [0] * len(pos),
            heights,
            color=mismatch_color,
            lw=1.2,
            zorder=3,
        )

    ax.set_yticks([])
    ax.set_ylabel(ylabel, fontsize=_fs(11))
    return float(depths.max()) if n else 0.0


def draw_sequence_track(
    ax,
    seq: str,
    region: Region,
    base_colors: Optional[Dict[str, str]] = None,
    base_fontsize: float = 11.0,
    min_px_per_char: float = 7.5,
) -> None:
    """Draw the reference nucleotide sequence on its own thin row.

    Letters are only drawn when the axis is wide enough to give each base
    ``min_px_per_char`` pixels (i.e. base resolution); otherwise a plain
    baseline is drawn so letters never overlap (IGV shows the sequence only at
    base resolution).
    """
    base_colors = base_colors or BASE_COLORS
    renderer = ax.figure.canvas.get_renderer()
    try:
        ax_width_px = ax.get_window_extent(renderer).width
    except Exception:
        ax_width_px = float("inf")
    legible = len(seq) * min_px_per_char <= ax_width_px

    if not legible:
        ax.plot([region.start, region.end], [0, 0], color=REF_BASE, lw=1.2, zorder=1)
        ax.set_xlim(region.start, region.end)
        ax.set_ylim(-0.4, 0.4)
        ax.set_yticks([])
        ax.set_ylabel("ref", fontsize=_fs(11))
        return
    for i, base in enumerate(seq.upper()):
        if base not in base_colors:
            base = "N"
        pos = region.start + i
        ax.text(
            pos + 0.5,
            0,
            base,
            ha="center",
            va="center",
            fontsize=_fs(base_fontsize),
            color=base_colors.get(base, "#333333"),
            family="monospace",
            weight="bold",
            zorder=4,
        )
    ax.set_xlim(region.start, region.end)
    ax.set_ylim(-0.6, 0.6)
    ax.set_yticks([])
    ax.set_ylabel("ref", fontsize=_fs(11))


def draw_sashimi_track(
    ax,
    counts: Dict[Tuple[int, int], int],
    region: Region,
    arc_color: str = SASHIMI,
    min_counts: int = 1,
    max_height: float = 1.0,
    label_fontsize: float = 9.5,
    thickness: Tuple[float, float] = (1.0, 4.0),
    jitter: float = 0.04,
) -> float:
    """Draw RNA-seq splice-junction arcs (sashimi plot).

    Junctions are drawn as dome-shaped cubic Bezier arcs (the IGV/ggsashimi
    convention): arc **thickness** scales with the number of supporting reads
    and arc **height** scales with the intron distance. Counts are labelled
    near each arc's peak. A deterministic jitter keeps nearby arcs separable.

    ``counts`` maps ``(skipped_start, skipped_end) -> n_supporting_reads`` (as
    returned by :func:`igvplot.reads.junction_counts`).

    Returns the maximum arc height used.
    """
    items = []
    for (s, e), cnt in sorted(counts.items(), key=lambda kv: (kv[0][0], kv[0][1])):
        s = max(s, region.start)
        e = min(e, region.end)
        if e - s >= 1:
            items.append((s, e, cnt))
    if not items:
        ax.axis("off")
        return 0.0

    max_count = max(cnt for _, _, cnt in items)
    largest_span = max(e - s for s, e, _ in items)
    lw_min, lw_max = thickness

    for i, (s, e, cnt) in enumerate(items):
        span = e - s
        # arc height proportional to intron distance, normalised so the
        # largest junction fills the track (no dead headroom)
        h = max_height * (span / largest_span) if largest_span else max_height
        # deterministic jitter around the height to reduce overdraw
        h *= 1.0 + jitter * ((i % 3) - 1)

        # line thickness proportional to read support
        lw = lw_min + (lw_max - lw_min) * (cnt / max_count)

        verts = [
            (s, 0.0),
            (s, h),
            (e, h),
            (e, 0.0),
        ]
        codes = [
            mpl.path.Path.MOVETO,
            mpl.path.Path.CURVE4,
            mpl.path.Path.CURVE4,
            mpl.path.Path.CURVE4,
        ]
        patch = mpl.patches.PathPatch(
            mpl.path.Path(verts, codes),
            facecolor="none",
            lw=lw,
            edgecolor=arc_color,
            zorder=2,
        )
        ax.add_patch(patch)

        ax.text(
            (s + e) / 2.0,
            h * 0.9,
            str(cnt),
            ha="center",
            va="center",
            fontsize=_fs(label_fontsize),
            color=arc_color,
            weight="bold",
            zorder=3,
        )

    ax.set_ylim(-0.03, max_height * 1.12)
    ax.set_yticks([])
    ax.set_ylabel("junction", fontsize=_fs(11))
    return max_height


def draw_base_mod_track(
    ax,
    mods: Dict[int, tuple],
    region: Region,
    label_fontsize: float = 9.5,
    marker_size: float = 60,
) -> None:
    """Draw a strand-aware base-modification track.

    ``mods`` maps 0-based position -> ``(strand, label, color)``. Each site is
    drawn as a coloured triangle (up for + strand, down for - strand) with its
    label above. Used for m6A / m5C / etc. modification sites (the "base mod"
    track in IGV).
    """
    if not mods:
        ax.axis("off")
        return
    for pos, vals in sorted(mods.items()):
        if not region.overlaps(pos, pos + 1):
            continue
        strand, label, color = vals
        marker = "^" if strand >= 0 else "v"
        ax.scatter([pos + 0.5], [0], marker=marker, c=[color], s=marker_size, zorder=4)
        ax.text(
            pos + 0.5,
            0.55 if strand >= 0 else -0.15,
            label,
            ha="center",
            va="bottom" if strand >= 0 else "top",
            fontsize=_fs(label_fontsize),
            color=color,
            weight="bold",
            zorder=5,
        )
    ax.set_xlim(region.start, region.end)
    ax.set_ylim(-0.5, 1.0)
    ax.set_yticks([])
    ax.set_ylabel("bases", fontsize=_fs(11))


def draw_sites(
    axes,
    region,
    sites: Dict[int, str],
    color: str = SITE,
    label_fontsize: float = 9.5,
) -> None:
    """Draw vertical site markers across all stacked axes, with level-stacking
    of their labels so close sites do not overlap (mirrors dna_features_viewer's
    annotation layout). ``sites`` maps 0-based position -> label."""
    if not sites:
        return
    items = sorted(sites.items())
    ax = axes[0]
    renderer = ax.figure.canvas.get_renderer()
    x0, x1 = ax.get_xlim()
    axwidth_pts = ax.get_window_extent(renderer).width

    # measured width in data units + a padded margin so close sites do not
    # get packed onto the same label row
    pad_pts = label_fontsize * 0.5

    def label_width_data(label: str) -> float:
        t = ax.text(0, 0, label, fontsize=_fs(label_fontsize))
        w_pts = t.get_window_extent(renderer=renderer).width + pad_pts
        t.remove()
        return w_pts * (x1 - x0) / axwidth_pts

    y_top = ax.get_ylim()[1]
    y_bottom = ax.get_ylim()[0]
    y_step = max(0.6, (y_top - y_bottom) / 6)

    # level -> rightmost occupied data x on that label row
    row_rights: List[float] = []

    for pos, label in items:
        cx = pos + 0.5
        if not (region.start <= pos < region.end):
            continue
        half = label_width_data(label) / 2
        left = cx - half
        # pick the first row that is clear of the label to its left
        level = 0
        for i, right in enumerate(row_rights):
            if left > right:
                level = i
                break
            level = i + 1
        if level >= len(row_rights):
            row_rights.append(cx + half)
        else:
            row_rights[level] = max(row_rights[level], cx + half)

        y = y_top + 0.6 * y_step - level * y_step
        for a in axes:
            a.axvline(cx, color=color, ls="--", lw=1.0, zorder=0, alpha=0.7)
        ax.plot(cx, y, marker="o", ms=3.2, color=color, zorder=8, alpha=0.9)
        ax.text(
            cx,
            y,
            label,
            ha="center",
            va="bottom",
            fontsize=_fs(label_fontsize),
            color=color,
            zorder=7,
        )


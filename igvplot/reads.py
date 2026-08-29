"""Fetch aligned reads and coverage from an indexed BAM/CRAM file with
pysam, and extract per-base mutation events (mismatches vs reference,
insertions, deletions) used for IGV-style plotting.

Unlike tools that shell out to ``samtools`` and parse SAM text, this module
works in-process on the pysam objects, so it is faster, robust, and gives
direct access to read flags, CIGAR and per-base events.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Tuple, Union

import numpy as np
import pysam

from os import fspath

from .region import Region
from .theme import BASE_COLORS

__all__ = [
    "Read",
    "Reference",
    "BASE_COLORS",
    "fetch_reads",
    "compute_coverage",
    "junction_counts",
    "compute_insert_sizes",
    "variant_allele_fraction",
    "open_reference",
]

# Canonical IGV-ish DNA base colours (see igvplot.theme for the palette).
_VALID_BASES = frozenset("ACGTN")


@dataclass
class Read:
    """An aligned read with per-base events relative to the reference.

    Coordinates are 0-based, half-open reference coordinates (pysam
    convention). Only events overlapping the requested ``Region`` are kept so
    plotting is bounded by the visible window.
    """

    name: str
    aleft: int  # leftmost aligned reference position (0-based)
    aright: int  # rightmost aligned reference position + 1
    is_reverse: bool
    mapq: int
    # ref_pos -> query base that differs from the reference at that position
    mismatches: Dict[int, str] = field(default_factory=dict)
    # ref_pos -> query base at every aligned reference position (for base-level
    # "show all bases" view). Populated only when collect_bases is requested.
    bases: Dict[int, str] = field(default_factory=dict)
    # anchored ref_pos -> inserted bases (present in read, absent in ref)
    insertions: Dict[int, str] = field(default_factory=dict)
    # list of (start, end) reference ranges deleted within this read
    deletions: List[Tuple[int, int]] = field(default_factory=list)
    # RNA-seq splice junctions: (skipped_start, skipped_end) from CIGAR 'N'
    junctions: List[Tuple[int, int]] = field(default_factory=list)
    # mate / read-group metadata for pairing and colouring (may be unset)
    paired: bool = False
    properly_paired: bool = False
    pairend_first: bool = False  # this read is first-of-pair (read1)
    pairend_second: bool = False  # this read is second-of-pair (read2)
    mate_is_reverse: bool = False  # does the mate map to the reverse strand
    mate_chrom: Optional[str] = None
    mate_start: Optional[int] = None
    insert_size: int = 0
    read_group: Optional[str] = None
    # optional auxiliary BAM tags (e.g. CB/CR/UB barcodes for single-cell or
    # long-read data), populated by fetch_reads(tag_keys=...) or auto-detected
    tags: Dict[str, object] = field(default_factory=dict)
    # soft-clipped bases at each end (5' / 3') of the stored query
    clip_left: int = 0
    clip_right: int = 0

    @property
    def start(self) -> int:
        return self.aleft

    @property
    def end(self) -> int:
        return self.aright

    def __repr__(self) -> str:
        return (
            f"Read(name={self.name!r}, {self.aleft}-{self.aright}, "
            f"{'R' if self.is_reverse else 'F'}, mapq={self.mapq}, "
            f"{len(self.mismatches)} mismatches, {len(self.insertions)} ins, "
            f"{len(self.deletions)} del, {len(self.junctions)} junctions)"
        )


class Reference:
    """Lazy per-region reference-sequence accessor.

    Serves single bases from small cached slices of the reference, which is
    far cheaper than one ``FastaFile.fetch`` call per read. Falls back to
    inert behaviour (base() -> None) when no fasta is available.
    """

    def __init__(self, fasta_path: Optional[str] = None, fasta_obj: Optional[pysam.FastaFile] = None):
        self._fasta = fasta_obj
        if fasta_path is not None:
            self._fasta = pysam.FastaFile(fspath(fasta_path))
        self._cache: Dict[Tuple[str, int, int], str] = {}

    def close(self) -> None:
        if self._fasta is not None:
            try:
                self._fasta.close()
            except Exception:
                pass
        self._fasta = None

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()

    @property
    def available(self) -> bool:
        return self._fasta is not None

    def get(self, chrom: str, start: int, end: int) -> str:
        if self._fasta is None:
            return ""
        key = (chrom, start, end)
        seq = self._cache.get(key)
        if seq is None:
            try:
                seq = self._fasta.fetch(chrom, start, end).upper()
            except (KeyError, ValueError):
                seq = ""
            self._cache[key] = seq
        return seq

    def base(self, chrom: str, pos: int) -> Optional[str]:
        """Reference base at 0-based ``pos`` (uppercase), or None if unknown."""
        for span in ((max(0, pos - 64), pos + 64), (pos, pos + 1)):
            seq = self.get(chrom, span[0], span[1])
            idx = pos - span[0]
            if 0 <= idx < len(seq):
                return seq[idx]
        return None


def open_reference(fasta_path: Optional[str] = None) -> Reference:
    """Return a :class:`Reference` accessor. With no path it is inert (no
    mismatches are computed)."""
    return Reference(fasta_path)


def variant_allele_fraction(
    bam_path: str,
    chrom: str,
    pos: int,
    reference: Optional[Union[str, Reference]] = None,
    min_mapq: int = 0,
) -> Tuple[float, int, int]:
    """Estimate the variant allele fraction at a 0-based ``pos``.

    Returns ``(vaf, depth, alt_count)`` where ``vaf`` is the fraction of reads
    whose base differs from the reference (a proxy for a bi-allelic SNP's ALT
    allele fraction), ``depth`` the total aligned depth and ``alt_count`` the
    number of non-reference reads.
    """
    region = Region(chrom, int(pos), int(pos) + 1)
    depths, mism = compute_coverage(bam_path, region, reference=reference, min_mapq=min_mapq)
    depth = int(depths[0]) if len(depths) else 0
    alt = int(mism[0]) if mism is not None and len(mism) else 0
    if depth == 0:
        return 0.0, 0, 0
    return alt / depth, depth, alt


def compute_insert_sizes(
    bam_path: str,
    region: Union[Region, str, tuple],
    min_mapq: int = 0,
    keep_duplicates: bool = False,
) -> np.ndarray:
    """Collect absolute insert sizes (TLEN) of paired reads overlapping ``region``.

    Returns a numpy array of insert sizes (0 if none). Used by
    :func:`igvplot.summary` for insert-size statistics.
    """
    region = Region.from_any(region)
    sizes: List[int] = []
    with pysam.AlignmentFile(fspath(bam_path), "rb") as bam:
        for seg in bam.fetch(region.chrom, region.start, region.end):
            if not seg.is_paired or seg.is_unmapped or seg.is_supplementary or seg.is_secondary:
                continue
            if not keep_duplicates and seg.is_duplicate:
                continue
            if min_mapq and seg.mapping_quality < min_mapq:
                continue
            tl = seg.template_length
            if tl and abs(tl) > 0:
                # count each fragment once (use the first-of-pair side)
                if not seg.is_read2:
                    sizes.append(abs(int(tl)))
    return np.array(sizes, dtype=np.int64) if sizes else np.zeros(0, dtype=np.int64)


def junction_counts(
    reads: List[Read],
    region: Region,
    min_counts: int = 1,
) -> Dict[Tuple[int, int], int]:
    """Count RNA-seq splice junctions (from CIGAR 'N') across ``reads``.

    Only junctions whose skipped (intron) interval overlaps ``region`` are
    counted. Returns ``{(start, end): n_supporting_reads}`` for junctions with
    at least ``min_counts`` supporting reads.
    """
    counts: Dict[Tuple[int, int], int] = {}
    for r in reads:
        for s, e in r.junctions:
            if not region.overlaps(s, e):
                continue
            key = (max(s, region.start), min(e, region.end))
            counts[key] = counts.get(key, 0) + 1
    return {k: v for k, v in counts.items() if v >= min_counts}


# Tags excluded from auto-detection: alignment bookkeeping that would only
# waste memory per read (users can still request them via tag_keys).
_TAG_BLOCKLIST = {
    "MD", "NM", "SA", "MC", "MQ", "AS", "XS", "NH", "HI", "IH", "ZF",
}


def _simple_tag_value(seg, key: str):
    """Return the tag ``key`` from ``seg`` if it is a simple scalar, else None."""
    try:
        v = seg.get_tag(key)
    except KeyError:
        return None
    if isinstance(v, (bytes, bytearray, list, tuple)):
        return None
    if isinstance(v, (str, int, float)):
        return v
    return None


def _events_for_segment(
    seg,
    reference: Optional[Reference],
    region: Region,
    collect_bases: bool = False,
    tag_keys: Optional[set] = None,
) -> Optional[Read]:
    """Construct a :class:`Read` with per-base events for ``seg`` within ``region``."""
    if seg.is_unmapped:
        return None

    aleft = seg.reference_start
    aright = seg.reference_end if seg.reference_end is not None else aleft

    a_start = seg.query_alignment_start if seg.query_alignment_start is not None else 0
    a_end = seg.query_alignment_end if seg.query_alignment_end is not None else 0
    qseq = seg.query_sequence or ""

    read = Read(
        name=seg.query_name or "?",
        aleft=aleft,
        aright=aright,
        is_reverse=seg.is_reverse,
        mapq=seg.mapping_quality,
    )

    # ---- mate / read-group metadata -------------------------------------
    read.paired = seg.is_paired
    read.properly_paired = seg.is_proper_pair
    read.pairend_first = seg.is_read1
    read.pairend_second = seg.is_read2
    read.mate_is_reverse = bool(seg.mate_is_reverse)
    if seg.is_paired and seg.next_reference_name is not None:
        read.mate_chrom = seg.next_reference_name
        read.mate_start = seg.next_reference_start
    read.insert_size = int(abs(seg.template_length) or 0) if seg.template_length else 0
    try:
        read.read_group = seg.get_tag("RG")
    except (KeyError, ValueError):
        read.read_group = None

    # ---- optional auxiliary tags (barcodes etc.) -------------------------
    if tag_keys:
        for key in tag_keys:
            v = _simple_tag_value(seg, key)
            if v is not None:
                read.tags[key] = v

    # soft-clipped bases at each end
    read.clip_left = seg.query_alignment_start if seg.query_alignment_start else 0
    read.clip_right = (
        (len(seg.query_sequence) - seg.query_alignment_end)
        if seg.query_sequence and seg.query_alignment_end is not None
        else 0
    )

    # ---- RNA-seq splice junctions (CIGAR 'N' skipped reference) ---------
    if seg.cigartuples is not None:
        ref_cursor = aleft
        for op, length in seg.cigartuples:
            if op == 0 or op == 7 or op == 8:  # M / = / X
                ref_cursor += length
            elif op == 2:  # D
                ref_cursor += length
            elif op == 3:  # N (intron / skipped region)
                read.junctions.append((ref_cursor, ref_cursor + length))
                ref_cursor += length
            elif op == 1:  # I
                pass
            elif op in (4, 5, 6):  # S / H / P -- no reference consumption
                pass

    mismatches: Dict[int, str] = {}
    insertions: Dict[int, str] = {}
    deletions: List[Tuple[int, int]] = []
    prev_rpos: Optional[int] = None

    for qpos, rpos in seg.get_aligned_pairs(matches_only=False):
        in_query = qpos is not None and a_start <= qpos < a_end

        if rpos is None:
            # Insertion (query bases with no reference position).
            if (
                in_query
                and prev_rpos is not None
                and region.overlaps(prev_rpos, prev_rpos + 1)
            ):
                base = qseq[qpos].upper() if qpos < len(qseq) else ""
                insertions[prev_rpos] = insertions.get(prev_rpos, "") + base
            continue

        if not region.overlaps(rpos, rpos + 1):
            prev_rpos = rpos
            continue

        if not in_query or qpos >= len(qseq):
            # A reference position with no query base is either a real
            # deletion ('D') or an RNA splice-junction ref-skip ('N'); N
            # positions must NOT be treated as deletions.
            rpos_in_skipped = False
            for jstart, jend in read.junctions:
                if jstart <= rpos < jend:
                    rpos_in_skipped = True
                    break
            if not rpos_in_skipped:
                deletions.append((rpos, rpos + 1))
        else:
            base = qseq[qpos].upper()
            if base in _VALID_BASES:
                if collect_bases:
                    read.bases[rpos] = base
                if reference is not None and reference.available:
                    ref_base = reference.base(region.chrom, rpos)
                    # Reads are stored reference-oriented in BAM (reverse-strand
                    # reads are reverse-complemented on write), so query bases
                    # compare directly against the reference.
                    if ref_base is not None and ref_base != base:
                        mismatches[rpos] = base
        prev_rpos = rpos

    read.mismatches = mismatches
    read.insertions = insertions
    # merge adjacent single-base deletions into contiguous spans so the deletion
    # reads as one event (with the correct length) instead of per-base entries
    merged: List[Tuple[int, int]] = []
    for s, e in sorted(deletions):
        if merged and s == merged[-1][1]:
            merged[-1] = (merged[-1][0], e)
        else:
            merged.append((s, e))
    read.deletions = [
        (max(s, region.start), min(e, region.end)) for s, e in merged if e > s
    ]
    return read


def fetch_reads(
    bam_path: str,
    region: Union[Region, str, tuple],
    reference: Optional[Union[str, Reference]] = None,
    min_mapq: int = 0,
    max_reads: Optional[int] = None,
    keep_duplicates: bool = False,
    keep_secondary: bool = False,
    sampling_window: int = 0,
    max_per_window: int = 0,
    collect_bases: bool = False,
    tag_keys: Optional[Iterable[str]] = None,
) -> List[Read]:
    """Fetch reads overlapping ``region`` from a sorted, indexed BAM/CRAM.

    Parameters
    ----------
    bam_path:
        Path to a sorted, indexed ``.bam`` or ``.cram`` file.
    region:
        Anything acceptable to :meth:`Region.from_any`.
    reference:
        Optional reference: a ``.fa`` path or a :class:`Reference` instance.
        When provided, per-base mismatches against the reference are computed.
    min_mapq:
        Drop reads with mapping quality strictly below this threshold.
    max_reads:
        Stop after collecting this many reads.
    keep_duplicates / keep_secondary:
        By default PCR/optical duplicates and secondary/supplementary
        alignments are skipped (standard IGV behaviour).
    sampling_window / max_per_window:
        IGV-style downsampling: bucket reads by ``reference_start //
        sampling_window`` and keep at most ``max_per_window`` reads per bucket
        (in coordinate order). Useful to bound memory/plot size for very deep
        or very large files.
    collect_bases:
        Also store the query base at every aligned reference position on each
        read (used by the base-resolution "show all bases" view). Set this only
        for small regions, since it increases memory.
    tag_keys:
        Auxiliary BAM tags to collect on every read (e.g. ``["CB", "UB"]`` for
        single-cell barcodes). When None (default) the simple scalar tags
        present in the file are auto-detected from a peek at the first reads,
        excluding alignment bookkeeping tags. Reads fetched via
        :meth:`igvplot.GenomeView.add_reads` can then be coloured or clustered
        with ``color_by="tag:CB"`` / ``group_by="tag:UB"``.
    """
    region = Region.from_any(region)
    ref = reference if isinstance(reference, Reference) else Reference(reference)
    close_ref = not isinstance(reference, Reference)

    reads: List[Read] = []
    seen_bins = {}
    keys: Optional[set] = set(tag_keys) if tag_keys is not None else None
    try:
        with pysam.AlignmentFile(fspath(bam_path), "rb") as bam:
            if keys is None:
                # auto-detect collectable simple tags from a peek at the first
                # reads in the region (barcodes etc., minus bookkeeping tags)
                detected: set = set()
                for probe_seg in bam.fetch(region.chrom, region.start, region.end):
                    for tag in probe_seg.get_tags():
                        key = tag[0]
                        if key not in _TAG_BLOCKLIST and _simple_tag_value(probe_seg, key) is not None:
                            detected.add(key)
                    if len(detected) > 32 or probe_seg.reference_start > region.start + 1000:
                        break
                keys = detected
            for seg in bam.fetch(region.chrom, region.start, region.end):
                if seg.is_unmapped or seg.is_supplementary:
                    continue
                if not keep_secondary and seg.is_secondary:
                    continue
                if not keep_duplicates and seg.is_duplicate:
                    continue
                if seg.is_qcfail:
                    continue
                if min_mapq and seg.mapping_quality < min_mapq:
                    continue
                if sampling_window and max_per_window:
                    b = seg.reference_start // sampling_window
                    n = seen_bins.get(b, 0)
                    if n >= max_per_window:
                        continue
                    seen_bins[b] = n + 1
                read = _events_for_segment(
                    seg, ref, region, collect_bases=collect_bases, tag_keys=keys
                )
                if read is not None:
                    reads.append(read)
                if max_reads is not None and len(reads) >= max_reads:
                    break
    finally:
        if close_ref:
            ref.close()
    return reads


def compute_coverage(
    bam_path: str,
    region: Union[Region, str, tuple],
    reference: Optional[Union[str, Reference]] = None,
    min_mapq: int = 0,
    min_base_quality: int = 0,
    strand: Optional[str] = None,
) -> Tuple[np.ndarray, Optional[np.ndarray]]:
    """Compute per-base sequencing depth over ``region``.

    Parameters
    ----------
    bam_path:
        Sorted, indexed BAM/CRAM file.
    region:
        Anything acceptable to :meth:`Region.from_any`.
    reference:
        Optional reference ``.fa`` path or :class:`Reference`. When provided,
        a second array with the per-base count of *mismatching* reads (the
        "variant pile") is also returned.
    strand:
        If ``'forward'`` or ``'reverse'``, only count reads of that strand.
    min_mapq / min_base_quality:
        Filtering thresholds.

    Returns
    -------
    (depths, mismatch_counts):
        ``depths`` has shape ``(region.length,)``. ``mismatch_counts`` is
        None unless a reference was supplied.
    """
    region = Region.from_any(region)
    n = region.length
    depths = np.zeros(n, dtype=np.int64)
    mismatches = np.zeros(n, dtype=np.int64) if reference is not None else None

    ref = reference if isinstance(reference, Reference) else Reference(reference)
    close_ref = not isinstance(reference, Reference)
    try:
        with pysam.AlignmentFile(fspath(bam_path), "rb") as bam:
            for col in bam.pileup(
                region.chrom,
                region.start,
                region.end,
                min_base_quality=min_base_quality,
                stepper="samtools",
                truncate=True,
            ):
                pos = col.reference_pos
                idx = pos - region.start
                if idx < 0 or idx >= n:
                    continue
                depth = 0
                for pread in col.pileups:
                    seg = pread.alignment
                    if seg.is_unmapped or seg.is_supplementary or seg.is_secondary:
                        continue
                    if seg.is_duplicate:
                        continue
                    if min_mapq and seg.mapping_quality < min_mapq:
                        continue
                    if strand and seg.is_reverse != (strand == "reverse"):
                        continue
                    # Only count reads that actually provide an aligned base
                    # here; deletions ('D') and splice ref-skips ('N') are
                    # not exon coverage (matches IGV behaviour).
                    if pread.is_del or pread.is_refskip:
                        continue
                    depth += 1
                    if mismatches is not None:
                        qpos = pread.query_position
                        if qpos is not None:
                            base = (seg.query_sequence or "")[qpos : qpos + 1].upper()
                            if base in _VALID_BASES:
                                ref_base = ref.base(region.chrom, pos)
                                if ref_base is not None and ref_base != base:
                                    mismatches[idx] += 1
                depths[idx] = depth
    finally:
        if close_ref:
            ref.close()
    return depths, mismatches

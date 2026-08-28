#!/usr/bin/env python3
"""Generate synthetic reference + annotation + BAM for testing/developing igvplot.

Produces (in ``data/``):
    genome.fa           reference contig 'chrTest' (20 kb) + genome.fa.fai
    annotation.gb       GenBank annotation with 4 genes + site features
    sites.bed           BED of the known mutation/indel positions
    sample.bam          reads (incl. known mutations/indels) + sample.bam.bai

The 'feature' region around chrTest:6940-7160 (in geneB, reverse strand)
contains reads with deliberate, verifiable mutation events:
    SNP C->T  at 7000   (12 reads)
    SNP A->G  at 7010   (5 reads)
    DEL 3 bp  at 7020-7023 (6 reads)
    INS "CC"  at 7030   (4 reads)

Run:  python scripts/make_synthetic_bam.py
"""
from __future__ import annotations

import os
import random
from collections import OrderedDict

import numpy as np
import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
CONTIG = "chrTest"
LENGTH = 20000
SEED = 42

READ_LEN = 120
COVERAGE = 50  # nominal depth -> fragment step = READ_LEN // COVERAGE

# ---- reference -----------------------------------------------------------
rng = random.Random(SEED)
seq_chars = "".join(
    rng.choice("ACGT") for _ in range(LENGTH)
)

# Genes (0-based, half-open).
GENES = [
    ("geneA", 1000, 3200, +1),
    ("geneB", 6800, 7600, -1),
    ("geneC", 13000, 15600, +1),
    ("geneD", 17000, 19300, -1),
]

# Known events inside the demo region.
PLOT_REGION = (6940, 7160)
EVENTS = [
    ("snp_7000_CT", 7000, "SNP", "C>T", 12),
    ("snp_7010_AG", 7010, "SNP", "A>G", 5),
    ("del_7020_3bp", 7020, "DEL", "3bp deletion", 6),
    ("ins_7030_CC", 7030, "INS", "CC insertion", 4),
]


def alt_base(b):
    return "A" if b != "A" else "C"


def make_read_query(ref, start, rlen, edits=None):
    """Build a forward-strand query sequence + CIGAR from the reference.

    ``edits`` is a dict with optional keys:
        mismatches: {rel_offset: base}
        insertions: {rel_offset (left boundary): bases}
        deletions:  [(rel_start, rel_end)]
    Returns (query, cigar).
    """
    mismatches = edits.get("mismatches", {}) if edits else {}
    insertions = edits.get("insertions", {}) if edits else {}
    deletions = [(a, b) for (a, b) in (edits.get("deletions", []) if edits else []) if b > a]

    query_parts = []
    cigar = []
    run = 0
    cur = 0
    ins_start = {off for off in insertions}
    del_start = {a for a, _ in deletions}

    while cur < rlen:
        if cur in del_start:
            if run:
                cigar.append((0, run))
                run = 0
            _, b = next((a, b) for a, b in deletions if a == cur)
            cigar.append((2, b - cur))
            cur = b
            continue
        if cur in ins_start:
            if run:
                cigar.append((0, run))
                run = 0
            bases = insertions[cur]
            query_parts.append(bases)
            cigar.append((1, len(bases)))
        base = mismatches.get(cur, ref[start + cur])
        query_parts.append(base)
        run += 1
        cur += 1
    if run:
        cigar.append((0, run))
    return "".join(query_parts), cigar


def revcomp(s):
    return s.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def revcomp(s):
    return s.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def main():
    os.makedirs(DATA, exist_ok=True)
    genome_fa = os.path.join(DATA, "genome.fa")
    annot_gb = os.path.join(DATA, "annotation.gb")
    sites_bed = os.path.join(DATA, "sites.bed")
    bam_path = os.path.join(DATA, "sample.bam")

    # ---- 1. reference fasta + index --------------------------------------
    with open(genome_fa, "w") as fh:
        fh.write(f">{CONTIG}\n{seq_chars}\n")
    pysam.faidx(genome_fa)

    # ---- 2. GenBank annotation ---------------------------------------------
    features = []
    for name, gstart, gend, strand in GENES:
        features.append(
            SeqFeature(
                FeatureLocation(gstart, gend, strand=strand),
                type="CDS",
                qualifiers={"label": name, "gene": name},
            )
        )
    for label, pos, ftype, note, _count in EVENTS:
        features.append(
            SeqFeature(
                FeatureLocation(pos, pos + 1, strand=None),
                type="variation",
                qualifiers={"label": label, "note": note},
            )
        )
    record = SeqRecord(
        Seq(seq_chars),
        id=CONTIG,
        name=CONTIG,
        description="synthetic contig for igvplot demo",
        annotations={"molecule_type": "DNA"},
        features=features,
    )
    with open(annot_gb, "w") as fh:
        SeqIO.write([record], fh, "genbank")

    # ---- 3. sites BED ------------------------------------------------------
    with open(sites_bed, "w") as fh:
        for label, pos, ftype, note, _count in EVENTS:
            fh.write(f"{CONTIG}\t{pos}\t{pos + 1}\t{label}({note})\n")

    # ---- 4. reads ----------------------------------------------------------
    step = max(1, READ_LEN // COVERAGE)  # ~COVERAGE x depth
    records = []
    read_id = 0

    def add_read(read_start, rlen, edits=None, force_forward=False, mapq=60,
                 read_group=None, insert_size=0):
        nonlocal read_id
        is_reverse = (not force_forward) and (read_id % 2 == 1)
        query, cigar = make_read_query(seq_chars, read_start, rlen, edits)
        # BAM convention: reads are stored reference-oriented (reverse-strand
        # reads are reverse-complemented on write). For our synthetic reads we
        # store the reference-oriented sequence directly for both strands and
        # mark orientation with flag 0x10, mirroring real alignment files.
        seg = pysam.AlignedSegment()
        seg.query_name = f"read_{read_id:05d}"
        seg.query_sequence = query
        seg.query_qualities = pysam.qualitystring_to_array("I" * len(query))
        seg.flag = 16 if is_reverse else 0
        seg.reference_id = 0
        seg.reference_start = read_start
        seg.mapping_quality = mapq
        seg.cigar = cigar
        if read_group:
            seg.set_tag("RG", read_group)
        if insert_size:
            seg.template_length = insert_size
        records.append((read_start, seg))
        read_id += 1
        return read_id - 1

    def join_read(exons, force_forward=False, mapq=60, read_group=None):
        """Create a spliced read whose aligned exons are listed in ``exons`` as
        ``(start, end)`` reference intervals; introns are skipped (CIGAR N)."""
        nonlocal read_id
        query_parts, cigar = [], []
        for i, (s, e) in enumerate(exons):
            if i:
                prev_end = exons[i - 1][1]
                gap = s - prev_end
                if gap > 0:
                    cigar.append((3, gap))  # N (intron / skipped region)
            cigar.append((0, e - s))  # M exon
            query_parts.append(seq_chars[s:e])
        is_reverse = (not force_forward) and (read_id % 2 == 1)
        seg = pysam.AlignedSegment()
        seg.query_name = f"read_{read_id:05d}"
        seg.query_sequence = "".join(query_parts)
        seg.query_qualities = pysam.qualitystring_to_array("I" * len(seg.query_sequence))
        seg.flag = 16 if is_reverse else 0
        seg.reference_id = 0
        seg.reference_start = exons[0][0]
        seg.mapping_quality = mapq
        seg.cigar = cigar
        if read_group:
            seg.set_tag("RG", read_group)
        records.append((exons[0][0], seg))
        read_id += 1

    def add_pair(frag_start, mate_start, frag_len, read_group=None, mapq=60):
        """Add a properly-paired read pair with a known insert size."""
        nonlocal read_id
        # read1: forward at frag_start; read2: reverse at mate_start.
        r1 = frag_start
        r2 = mate_start
        seg = pysam.AlignedSegment()
        seg.query_name = f"pair_{read_id:05d}"
        seg.query_sequence = seq_chars[r1:r1 + READ_LEN]
        seg.query_qualities = pysam.qualitystring_to_array("I" * READ_LEN)
        seg.flag = 0x1 | 0x2 | 0x40  # paired, proper, first-of-pair
        seg.reference_id = 0
        seg.reference_start = r1
        seg.mapping_quality = mapq
        seg.cigar = [(0, READ_LEN)]
        seg.next_reference_id = 0
        seg.next_reference_start = r2
        seg.template_length = (r2 + READ_LEN) - r1
        if read_group:
            seg.set_tag("RG", read_group)
        records.append((r1, seg))

        seg2 = pysam.AlignedSegment()
        seg2.query_name = f"pair_{read_id:05d}"
        seg2.query_sequence = seq_chars[r2:r2 + READ_LEN]
        seg2.query_qualities = pysam.qualitystring_to_array("I" * READ_LEN)
        seg2.flag = 0x1 | 0x2 | 0x80 | 0x10  # paired, proper, second-of-pair, reverse
        seg2.reference_id = 0
        seg2.reference_start = r2
        seg2.mapping_quality = mapq
        seg2.cigar = [(0, READ_LEN)]
        seg2.next_reference_id = 0
        seg2.next_reference_start = r1
        seg2.template_length = (r2 + READ_LEN) - r1
        if read_group:
            seg2.set_tag("RG", read_group)
        records.append((r2, seg2))
        read_id += 1

    # baseline tiling fragments (dense => deep coverage)
    for start in range(0, LENGTH - READ_LEN, step):
        # occasional random base error so the view looks "real" -- kept OUTSIDE
        # the demo region so the exact, tested mutation counts stay deterministic
        edits = None
        if rng.random() < 0.004:
            off = rng.randrange(READ_LEN)
            pos = start + off
            if not (PLOT_REGION[0] <= pos < PLOT_REGION[1]):
                edits = {"mismatches": {off: alt_base(seq_chars[pos])}}
        # vary mapping quality so the mapq colour mode shows a real gradient
        add_read(start, READ_LEN, edits, mapq=60 - (start % 20))

    # targeted mutant reads in the demo region
    for label, pos, ftype, _note, count in EVENTS:
        if ftype == "SNP":
            start = pos - 70
            off = pos - start
            for _ in range(count):
                # pick a base that differs from the reference so it is a real SNP
                refb = seq_chars[pos]
                add_read(start, READ_LEN, {"mismatches": {off: alt_base(refb)}}, force_forward=(pos % 2 == 0))
        elif ftype == "DEL":
            start = pos - 50
            off = pos - start
            for _ in range(count):
                add_read(start, READ_LEN, {"deletions": [(off, off + 3)]}, force_forward=(pos % 2 == 0))
        elif ftype == "INS":
            start = pos - 30
            off = pos - start
            for _ in range(count):
                add_read(start, READ_LEN + 2, {"insertions": {off: "CC"}}, force_forward=(pos % 2 == 0))

    # spliced reads (RNA-seq) -> two splice junctions for the sashimi track
    for k in range(25):
        st = 6935 - (k % 5)
        join_read([(st, 6970), (7060, 7080)], mapq=60, read_group="isoform1")
    for k in range(12):
        st = 7050 - (k % 4)
        join_read([(st, 7090), (7140, 7180)], mapq=60, read_group="isoform2")
    # a few low-mapq reads to exercise the mapq colour mode
    for k in range(6):
        add_read(6950 + 5 * k, READ_LEN, mapq=15 + (k % 3) * 10)

    # properly-paired reads for mate linking / insert-size colouring
    for k in range(10):
        add_pair(6880 + 10 * k, 6900 + 10 * k, READ_LEN, read_group="paired", mapq=60)

    # soft-clipped reads to exercise the show_soft_clips path
    for k in range(8):
        start = 6960 + 6 * k
        mlen = 90
        clip = 12 + (k % 3) * 4
        q = (
            "".join(rng.choice("ACGT") for _ in range(clip))
            + seq_chars[start : start + mlen]
            + "".join(rng.choice("ACGT") for _ in range(clip))
        )
        seg = pysam.AlignedSegment()
        seg.query_name = f"soft_{k:02d}"
        seg.query_sequence = q
        seg.query_qualities = pysam.qualitystring_to_array("I" * len(q))
        seg.flag = 16 if (k % 2 == 1) else 0
        seg.reference_id = 0
        seg.reference_start = start
        seg.mapping_quality = 60
        seg.cigar = [(4, clip), (0, mlen), (4, clip)]
        records.append((start, seg))

    # write sorted by coordinate
    records.sort(key=lambda r: r[0])
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": CONTIG, "LN": LENGTH}],
    }
    with pysam.AlignmentFile(bam_path, "wb", header=header) as out:
        for _start, seg in records:
            out.write(seg)
    pysam.index(bam_path)

    print(f"[synthetic] wrote:\n  {genome_fa}\n  {annot_gb}\n  {sites_bed}\n  {bam_path}")
    print(f"[synthetic] contig={CONTIG} length={LENGTH} reads={len(records)}")
    print(f"[synthetic] demo region {CONTIG}:{PLOT_REGION[0]}-{PLOT_REGION[1]}")
    for label, pos, ftype, note, count in EVENTS:
        print(f"  {label:<14} {pos}: {note} ({count} reads)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

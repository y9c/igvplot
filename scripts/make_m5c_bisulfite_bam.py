#!/usr/bin/env python3
"""Generate synthetic bisulfite RNA-seq data for m5C detection.

Bisulfite chemistry: unmethylated C → U (reads as T), methylated C (5mC) → unchanged C.

Produces (in ``data/m5c/``):
    m5c_genome.fa       reference transcript with exons
    m5c_annotation.gb   GenBank with exon/CDS features
    m5c_reads.bam       spliced bisulfite reads + .bai
    m5c_sites.bed       m5C site positions

The transcript has:
    - 2 exons with intron splicing (reads span exon-exon junctions)
    - Multiple C positions: some are m5C (methylated), others unmethylated
    - Reads show T at unmethylated C, C at m5C positions

Run:  python scripts/make_m5c_bisulfite_bam.py
"""
from __future__ import annotations

import os
import random

import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data", "m5c")
os.makedirs(DATA, exist_ok=True)

CONTIG = "chrM5C"
LENGTH = 1000
SEED = 123

READ_LEN = 100
COVERAGE = 30

rng = random.Random(SEED)

# Reference sequence
seq_chars = "".join(rng.choice("ACGT") for _ in range(LENGTH))

# Transcript structure: exon1 (100-400), intron (400-600), exon2 (600-900)
EXON1_START, EXON1_END = 100, 400
INTRON_START, INTRON_END = 400, 600
EXON2_START, EXON2_END = 600, 900

# m5C sites (0-based positions in reference)
M5C_SITES = [150, 250, 350, 650, 750, 850]

# Ensure reference has C at m5C positions
seq_list = list(seq_chars)
for pos in M5C_SITES:
    seq_list[pos] = "C"
seq_chars = "".join(seq_list)


def make_spliced_read(ref, start, end, strand="+"):
    """Build a spliced read that may span exon-exon junction.
    
    Returns (query, cigar, ref_start) for the aligned read.
    """
    query_parts = []
    cigar = []
    
    # Determine which exons this read overlaps
    in_exon1 = start < EXON1_END and end > EXON1_START
    in_exon2 = start < EXON2_END and end > EXON2_START
    spans_intron = start < INTRON_START and end > INTRON_END
    
    if spans_intron:
        # Read spans the intron: align to exon1 part, skip intron, align to exon2 part
        exon1_part = min(end, EXON1_END) - max(start, EXON1_START)
        exon2_part = min(end, EXON2_END) - max(start, EXON2_START)
        
        # Query sequence: concat of both exon parts
        for pos in range(max(start, EXON1_START), min(end, EXON1_END)):
            query_parts.append(ref[pos])
        for pos in range(max(start, EXON2_START), min(end, EXON2_END)):
            query_parts.append(ref[pos])
        
        # CIGAR: match exon1, skip intron, match exon2
        cigar = [
            (0, exon1_part),  # M
            (3, INTRON_END - INTRON_START),  # N (intron skip)
            (0, exon2_part),  # M
        ]
        ref_start = max(start, EXON1_START)
    else:
        # Read within single exon
        actual_start = max(start, EXON1_START if in_exon1 else EXON2_START)
        actual_end = min(end, EXON1_END if in_exon1 else EXON2_END)
        
        for pos in range(actual_start, actual_end):
            query_parts.append(ref[pos])
        
        cigar = [(0, actual_end - actual_start)]
        ref_start = actual_start
    
    query = "".join(query_parts)
    
    # Reverse complement if reverse strand
    if strand == "-":
        query = query.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    
    return query, cigar, ref_start


def apply_bisulfite_conversion(query, ref_positions, m5c_set):
    """Apply bisulfite conversion to query sequence.
    
    C → T (U) except at m5C positions where C is protected.
    """
    result = []
    for i, base in enumerate(query):
        if base == "C" and i < len(ref_positions):
            ref_pos = ref_positions[i]
            if ref_pos in m5c_set:
                result.append("C")  # m5C protected
            else:
                result.append("T")  # bisulfite converted
        else:
            result.append(base)
    return "".join(result)


def get_ref_positions_for_read(cigar, ref_start):
    """Get reference positions consumed by the read."""
    positions = []
    pos = ref_start
    for op, length in cigar:
        if op == 0:  # M
            positions.extend(range(pos, pos + length))
            pos += length
        elif op == 2:  # D
            pos += length
        elif op == 3:  # N (skip)
            pos += length
    return positions


def main():
    genome_fa = os.path.join(DATA, "m5c_genome.fa")
    annot_gb = os.path.join(DATA, "m5c_annotation.gb")
    bam_path = os.path.join(DATA, "m5c_reads.bam")
    sites_bed = os.path.join(DATA, "m5c_sites.bed")

    # 1. Reference FASTA
    with open(genome_fa, "w") as fh:
        fh.write(f">{CONTIG}\n{seq_chars}\n")
    pysam.faidx(genome_fa)

    # 2. GenBank annotation
    features = []
    
    # Gene (full transcript span)
    features.append(SeqFeature(
        FeatureLocation(EXON1_START, EXON2_END, strand=1),
        type="gene",
        qualifiers={"label": "gene1"},
    ))
    
    # Exons
    features.append(SeqFeature(
        FeatureLocation(EXON1_START, EXON1_END, strand=1),
        type="exon",
        qualifiers={"label": "exon1"},
    ))
    features.append(SeqFeature(
        FeatureLocation(EXON2_START, EXON2_END, strand=1),
        type="exon",
        qualifiers={"label": "exon2"},
    ))
    
    # CDS (spans both exons)
    features.append(SeqFeature(
        FeatureLocation(EXON1_START, EXON1_END, strand=1),
        type="CDS",
        qualifiers={"label": "CDS"},
    ))
    features.append(SeqFeature(
        FeatureLocation(EXON2_START, EXON2_END, strand=1),
        type="CDS",
        qualifiers={"label": "CDS"},
    ))
    
    record = SeqRecord(
        Seq(seq_chars),
        id=CONTIG,
        name=CONTIG,
        description="m5C bisulfite RNA transcript",
        annotations={"molecule_type": "RNA"},
        features=features,
    )
    with open(annot_gb, "w") as fh:
        SeqIO.write([record], fh, "genbank")

    # 3. m5C sites BED
    with open(sites_bed, "w") as fh:
        for pos in M5C_SITES:
            fh.write(f"{CONTIG}\t{pos}\t{pos+1}\tm5C\t100\t+\n")

    # 4. BAM file with spliced bisulfite reads
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": CONTIG, "LN": LENGTH}],
    }
    
    m5c_set = set(M5C_SITES)
    reads = []
    
    # Generate reads spanning transcript
    step = READ_LEN // COVERAGE
    for start in range(EXON1_START - 50, EXON2_END + 50, step):
        end = start + READ_LEN
        
        # Only generate reads that overlap exons
        if end <= EXON1_START or start >= EXON2_END:
            continue
        
        # Skip reads entirely in intron
        if start >= INTRON_START and end <= INTRON_END:
            continue
        
        strand = rng.choice(["+", "-"])
        query, cigar, ref_start = make_spliced_read(seq_chars, start, end, strand)
        
        # Get reference positions for bisulfite conversion
        ref_positions = get_ref_positions_for_read(cigar, ref_start)
        
        # Apply bisulfite conversion (C→T except at m5C)
        query = apply_bisulfite_conversion(query, ref_positions, m5c_set)
        
        reads.append({
            "query": query,
            "cigar": cigar,
            "ref_start": ref_start,
            "strand": strand,
        })
    
    # Write BAM
    with pysam.AlignmentFile(bam_path, "wb", header=header) as bam:
        for i, r in enumerate(reads):
            aln = pysam.AlignedSegment(bam.header)
            aln.query_name = f"read_{i:05d}"
            aln.query_sequence = r["query"]
            aln.is_reverse = (r["strand"] == "-")
            aln.mapping_quality = 60
            aln.reference_start = r["ref_start"]
            aln.cigar = r["cigar"]
            aln.query_qualities = pysam.qualitystring_to_array("I" * len(r["query"]))
            bam.write(aln)
    
    pysam.sort("-o", bam_path, bam_path)
    pysam.index(bam_path)

    print(f"[m5C bisulfite] wrote:")
    print(f"  {genome_fa}")
    print(f"  {annot_gb}")
    print(f"  {sites_bed}")
    print(f"  {bam_path}")
    print(f"[m5C bisulfite] contig={CONTIG} length={LENGTH} reads={len(reads)}")
    print(f"[m5C bisulfite] m5C sites: {M5C_SITES}")
    print(f"[m5C bisulfite] exon1: {EXON1_START}-{EXON1_END}, exon2: {EXON2_START}-{EXON2_END}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Generate synthetic GLORI RNA-seq data for m6A detection.

GLORI chemistry: A residues can be modified to A→G in sequencing reads when m6A is present.
Regular A stays as A in reads.

Produces (in ``data/m6a/``):
    m6a_genome.fa       reference sequence with gene on negative strand
    m6a_annotation.gb   GenBank with exon/CDS features (reverse strand)
    m6a_reads.bam       spliced GLORI reads with A→G mutations + .bai
    m6a_sites.bed       m6A site positions

The transcript has:
    - Gene on NEGATIVE strand
    - 3 exons with intron splicing (reads span exon-exon junctions)
    - Multiple A positions: some are m6A (show A→G), others unmethylated (stay A)
    - Reads from BOTH forward and reverse strands

Run:  python scripts/make_m6a_glori_bam.py
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
DATA = os.path.join(HERE, "..", "data", "m6a")
os.makedirs(DATA, exist_ok=True)

CONTIG = "chrM6A"
LENGTH = 1200
SEED = 456

READ_LEN = 120
COVERAGE = 40

rng = random.Random(SEED)

# Reference sequence
seq_chars = "".join(rng.choice("ACGT") for _ in range(LENGTH))

# Gene on NEGATIVE strand: exon1 (200-500), intron1 (500-700), 
# exon2 (700-850), intron2 (850-950), exon3 (950-1100)
EXON1_START, EXON1_END = 200, 500
INTRON1_START, INTRON1_END = 500, 700
EXON2_START, EXON2_END = 700, 850
INTRON2_START, INTRON2_END = 850, 950
EXON3_START, EXON3_END = 950, 1100

GENE_STRAND = -1  # negative strand

# m6A sites (0-based positions in reference) - these will show A→G in reads
M6A_SITES = [250, 350, 450, 750, 800, 1000, 1050]

# Ensure reference has A at m6A positions
seq_list = list(seq_chars)
for pos in M6A_SITES:
    seq_list[pos] = "A"
seq_chars = "".join(seq_list)


def make_spliced_read(ref, start, end, strand="+"):
    """Build a spliced read that may span multiple exon-exon junctions.
    
    Returns (query, cigar, ref_start) for the aligned read.
    """
    query_parts = []
    cigar = []
    
    # Check which exons this read overlaps
    exons = [
        (EXON1_START, EXON1_END),
        (EXON2_START, EXON2_END),
        (EXON3_START, EXON3_END),
    ]
    introns = [
        (INTRON1_START, INTRON1_END),
        (INTRON2_START, INTRON2_END),
    ]
    
    # Find overlapping exons
    overlaps = []
    for exon_start, exon_end in exons:
        if start < exon_end and end > exon_start:
            overlap_start = max(start, exon_start)
            overlap_end = min(end, exon_end)
            overlaps.append((overlap_start, overlap_end))
    
    if not overlaps:
        return None, None, None
    
    # Build query and CIGAR
    ref_start = overlaps[0][0]
    
    for i, (overlap_start, overlap_end) in enumerate(overlaps):
        # Add exon sequence to query
        for pos in range(overlap_start, overlap_end):
            query_parts.append(ref[pos])
        
        # Add match operation
        cigar.append((0, overlap_end - overlap_start))
        
        # If there's a next exon, add intron skip
        if i < len(overlaps) - 1:
            next_start = overlaps[i + 1][0]
            intron_len = next_start - overlap_end
            cigar.append((3, intron_len))  # N operation for intron skip
    
    query = "".join(query_parts)
    
    # Reverse complement if reverse strand
    if strand == "-":
        query = query.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    
    return query, cigar, ref_start


def apply_glori_conversion(query, ref_positions, m6a_set, mutation_rate=0.8):
    """Apply GLORI conversion to query sequence.
    
    A → G at m6A positions (with mutation_rate probability).
    Regular A stays as A.
    """
    result = []
    for i, base in enumerate(query):
        if base == "A" and i < len(ref_positions):
            ref_pos = ref_positions[i]
            if ref_pos in m6a_set and rng.random() < mutation_rate:
                result.append("G")  # GLORI mutation at m6A
            else:
                result.append("A")  # unchanged
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
    genome_fa = os.path.join(DATA, "m6a_genome.fa")
    annot_gb = os.path.join(DATA, "m6a_annotation.gb")
    bam_path = os.path.join(DATA, "m6a_reads.bam")
    sites_bed = os.path.join(DATA, "m6a_sites.bed")

    # 1. Reference FASTA
    with open(genome_fa, "w") as fh:
        fh.write(f">{CONTIG}\n{seq_chars}\n")
    pysam.faidx(genome_fa)

    # 2. GenBank annotation - gene on NEGATIVE strand
    features = []
    
    # Gene (full transcript span, reverse strand)
    features.append(SeqFeature(
        FeatureLocation(EXON1_START, EXON3_END, strand=GENE_STRAND),
        type="gene",
        qualifiers={"label": "gene1", "gene": "GLORI_test"},
    ))
    
    # mRNA (reverse strand)
    features.append(SeqFeature(
        FeatureLocation(EXON1_START, EXON3_END, strand=GENE_STRAND),
        type="mRNA",
        qualifiers={"label": "transcript1", "gene": "GLORI_test"},
    ))
    
    # Exons (reverse strand)
    exons = [
        (EXON1_START, EXON1_END, "exon1"),
        (EXON2_START, EXON2_END, "exon2"),
        (EXON3_START, EXON3_END, "exon3"),
    ]
    for exon_start, exon_end, label in exons:
        features.append(SeqFeature(
            FeatureLocation(exon_start, exon_end, strand=GENE_STRAND),
            type="exon",
            qualifiers={"label": label, "gene": "GLORI_test"},
        ))
    
    # CDS (reverse strand, spans all exons)
    for exon_start, exon_end, _ in exons:
        features.append(SeqFeature(
            FeatureLocation(exon_start, exon_end, strand=GENE_STRAND),
            type="CDS",
            qualifiers={"label": "CDS", "gene": "GLORI_test"},
        ))
    
    record = SeqRecord(
        Seq(seq_chars),
        id=CONTIG,
        name=CONTIG,
        description="m6A GLORI RNA transcript (reverse strand)",
        annotations={"molecule_type": "RNA"},
        features=features,
    )
    with open(annot_gb, "w") as fh:
        SeqIO.write([record], fh, "genbank")

    # 3. m6A sites BED
    with open(sites_bed, "w") as fh:
        for pos in M6A_SITES:
            fh.write(f"{CONTIG}\t{pos}\t{pos+1}\tm6A\t100\t+\n")

    # 4. BAM file with spliced GLORI reads
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": CONTIG, "LN": LENGTH}],
    }
    
    m6a_set = set(M6A_SITES)
    reads = []
    
    # Define introns for checking
    introns = [
        (INTRON1_START, INTRON1_END),
        (INTRON2_START, INTRON2_END),
    ]
    
    # Generate reads spanning transcript region
    step = READ_LEN // COVERAGE
    for start in range(EXON1_START - 100, EXON3_END + 100, step):
        end = start + READ_LEN
        
        # Only generate reads that overlap at least one exon
        overlaps_exon = False
        for exon_start, exon_end, _ in exons:
            if start < exon_end and end > exon_start:
                overlaps_exon = True
                break
        
        if not overlaps_exon:
            continue
        
        # Skip reads entirely in introns
        in_intron_only = True
        for intron_start, intron_end in introns:
            if start >= intron_start and end <= intron_end:
                in_intron_only = True
                break
        else:
            in_intron_only = False
        
        if in_intron_only:
            continue
        
        # Generate reads from BOTH strands
        for strand in ["+", "-"]:
            query, cigar, ref_start = make_spliced_read(seq_chars, start, end, strand)
            
            if query is None:
                continue
            
            # Get reference positions for GLORI conversion
            ref_positions = get_ref_positions_for_read(cigar, ref_start)
            
            # Apply GLORI conversion (A→G at m6A sites)
            query = apply_glori_conversion(query, ref_positions, m6a_set, mutation_rate=0.85)
            
            reads.append({
                "query": query,
                "cigar": cigar,
                "ref_start": ref_start,
                "strand": strand,
            })
    
    # Write BAM
    temp_bam = bam_path + ".unsorted.bam"
    with pysam.AlignmentFile(temp_bam, "wb", header=header) as bam:
        for i, r in enumerate(reads):
            aln = pysam.AlignedSegment(bam.header)
            aln.query_name = f"read_{i:05d}"
            aln.query_sequence = r["query"]
            aln.is_reverse = (r["strand"] == "-")
            aln.mapping_quality = 60
            aln.reference_id = 0  # Set reference ID to the first (and only) chromosome
            aln.reference_start = r["ref_start"]
            aln.cigar = r["cigar"]
            aln.query_qualities = pysam.qualitystring_to_array("I" * len(r["query"]))
            bam.write(aln)
    
    # Sort and index
    pysam.sort("-o", bam_path, temp_bam)
    os.remove(temp_bam)
    pysam.index(bam_path)

    print(f"[m6A GLORI] wrote:")
    print(f"  {genome_fa}")
    print(f"  {annot_gb}")
    print(f"  {sites_bed}")
    print(f"  {bam_path}")
    print(f"[m6A GLORI] contig={CONTIG} length={LENGTH} reads={len(reads)}")
    print("[m6A GLORI] gene strand: NEGATIVE (-)")
    print(f"[m6A GLORI] m6A sites (A→G): {M6A_SITES}")
    print(f"[m6A GLORI] exons: {EXON1_START}-{EXON1_END}, {EXON2_START}-{EXON2_END}, {EXON3_START}-{EXON3_END}")
    print("[m6A GLORI] reads from both +/- strands")


if __name__ == "__main__":
    main()

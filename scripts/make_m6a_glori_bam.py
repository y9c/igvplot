#!/usr/bin/env python3
"""Generate synthetic GLORI RNA-seq data for m6A detection with two genes.

GLORI chemistry:
- Forward strand gene: A residues at m6A sites → A→G in reads
- Reverse strand gene: T residues at m6A sites → T→C in reads (reverse complement)

Produces (in ``data/m6a/``):
    m6a_genome.fa       reference sequence with genes on both strands
    m6a_annotation.gb   GenBank with exon/CDS features (both strands)
    m6a_reads.bam       spliced GLORI reads with strand-specific mutations + .bai
    m6a_sites.bed       m6A site positions

The region chrM6A:200-1100 has:
    - Forward strand gene (exons: 200-400, 500-700) with m6A sites showing A→G
    - Reverse strand gene (exons: 750-900, 1000-1100) with m6A sites showing T→C
    - Paired-end reads from both strands with proper linking

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
INSERT_SIZE = 250
COVERAGE = 100  # Increased from 30 to show more reads

rng = random.Random(SEED)

# Reference sequence
seq_chars = "".join(rng.choice("ACGT") for _ in range(LENGTH))

# Forward strand gene: exons at 200-400, 500-700
FWD_EXON1_START, FWD_EXON1_END = 200, 400
FWD_EXON2_START, FWD_EXON2_END = 500, 700

# Reverse strand gene: exons at 750-900, 1000-1100
REV_EXON1_START, REV_EXON1_END = 750, 900
REV_EXON2_START, REV_EXON2_END = 1000, 1100

# m6A sites for forward gene (A positions, will show A→G)
FWD_M6A_SITES = [250, 350, 550, 650]

# m6A sites for reverse gene (T positions, will show T→C)
REV_M6A_SITES = [800, 850, 1020, 1070]

# Ensure reference has A at forward m6A positions and T at reverse m6A positions
seq_list = list(seq_chars)
for pos in FWD_M6A_SITES:
    seq_list[pos] = "A"
for pos in REV_M6A_SITES:
    seq_list[pos] = "T"
seq_chars = "".join(seq_list)


def make_spliced_read(ref, start, end, exons):
    """Build a spliced read spanning multiple exons.
    
    Returns (query, cigar, ref_start) for the aligned read.
    """
    query_parts = []
    cigar = []
    
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
    return query, cigar, ref_start


def apply_glori_conversion_forward(query, ref_positions, m6a_set, conversion_rate=0.85):
    """Apply GLORI conversion for forward strand.
    
    GLORI is a negative method:
    - Unmodified A → G (reagent converted it)
    - m6A-protected A → stays A (methylation blocked conversion)
    """
    result = []
    for i, base in enumerate(query):
        if base == "A" and i < len(ref_positions):
            ref_pos = ref_positions[i]
            if ref_pos in m6a_set:
                # m6A site - protected, stays as A
                result.append("A")
            else:
                # Unmodified A - convert to G with conversion_rate probability
                if rng.random() < conversion_rate:
                    result.append("G")
                else:
                    result.append("A")
        else:
            result.append(base)
    return "".join(result)


def apply_glori_conversion_reverse(query, ref_positions, m6a_set, conversion_rate=0.85):
    """Apply GLORI conversion for reverse strand.
    
    GLORI is a negative method:
    - Unmodified T → C (reagent converted the complement A)
    - m6A-protected T → stays T (methylation on forward strand blocked conversion)
    """
    result = []
    for i, base in enumerate(query):
        if base == "T" and i < len(ref_positions):
            ref_pos = ref_positions[i]
            if ref_pos in m6a_set:
                # m6A site - protected, stays as T
                result.append("T")
            else:
                # Unmodified T (complement of A) - convert to C
                if rng.random() < conversion_rate:
                    result.append("C")
                else:
                    result.append("T")
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

    # 2. GenBank annotation - genes on both strands
    features = []
    
    # Forward strand gene
    features.append(SeqFeature(
        FeatureLocation(FWD_EXON1_START, FWD_EXON2_END, strand=1),
        type="gene",
        qualifiers={"label": "FWD_gene", "gene": "GLORI_forward"},
    ))
    features.append(SeqFeature(
        FeatureLocation(FWD_EXON1_START, FWD_EXON1_END, strand=1),
        type="exon",
        qualifiers={"label": "FWD_exon1", "gene": "GLORI_forward"},
    ))
    features.append(SeqFeature(
        FeatureLocation(FWD_EXON2_START, FWD_EXON2_END, strand=1),
        type="exon",
        qualifiers={"label": "FWD_exon2", "gene": "GLORI_forward"},
    ))
    
    # Reverse strand gene
    features.append(SeqFeature(
        FeatureLocation(REV_EXON1_START, REV_EXON2_END, strand=-1),
        type="gene",
        qualifiers={"label": "REV_gene", "gene": "GLORI_reverse"},
    ))
    features.append(SeqFeature(
        FeatureLocation(REV_EXON1_START, REV_EXON1_END, strand=-1),
        type="exon",
        qualifiers={"label": "REV_exon1", "gene": "GLORI_reverse"},
    ))
    features.append(SeqFeature(
        FeatureLocation(REV_EXON2_START, REV_EXON2_END, strand=-1),
        type="exon",
        qualifiers={"label": "REV_exon2", "gene": "GLORI_reverse"},
    ))
    
    record = SeqRecord(
        Seq(seq_chars),
        id=CONTIG,
        name=CONTIG,
        description="m6A GLORI RNA with genes on both strands",
        annotations={"molecule_type": "RNA"},
        features=features,
    )
    with open(annot_gb, "w") as fh:
        SeqIO.write([record], fh, "genbank")

    # 3. m6A sites BED (all sites, both strands)
    with open(sites_bed, "w") as fh:
        for pos in FWD_M6A_SITES:
            fh.write(f"{CONTIG}\t{pos}\t{pos+1}\tm6A_fwd\t100\t+\n")
        for pos in REV_M6A_SITES:
            fh.write(f"{CONTIG}\t{pos}\t{pos+1}\tm6A_rev\t100\t-\n")

    # 4. BAM file with spliced GLORI reads from both genes
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": CONTIG, "LN": LENGTH}],
    }
    
    fwd_m6a_set = set(FWD_M6A_SITES)
    rev_m6a_set = set(REV_M6A_SITES)
    reads = []
    frag_id = 0  # mates of one fragment share this query name (for PE join)
    
    # Generate paired-end reads for forward gene
    fwd_exons = [(FWD_EXON1_START, FWD_EXON1_END), (FWD_EXON2_START, FWD_EXON2_END)]
    step = READ_LEN // COVERAGE
    for start in range(FWD_EXON1_START - 50, FWD_EXON2_END + 50, step):
        # Skip reads entirely in intron
        if start >= FWD_EXON1_END and start < FWD_EXON2_START:
            continue
        
        # Read 1 (forward strand)
        query1, cigar1, ref_start1 = make_spliced_read(seq_chars, start, start + READ_LEN, fwd_exons)
        if query1 is None:
            continue
        frag_id += 1
        name = f"frag_{frag_id:05d}"

        # Read 2 (reverse strand, mate)
        mate_start = start + INSERT_SIZE
        query2, cigar2, ref_start2 = make_spliced_read(seq_chars, mate_start, mate_start + READ_LEN, fwd_exons)
        
        # Apply GLORI conversion to forward gene reads (A→G)
        ref_positions1 = get_ref_positions_for_read(cigar1, ref_start1)
        query1 = apply_glori_conversion_forward(query1, ref_positions1, fwd_m6a_set)
        
        reads.append({
            "name": name,
            "query": query1,
            "cigar": cigar1,
            "ref_start": ref_start1,
            "strand": "+",
            "mate_start": mate_start if query2 else None,
            "read_num": 1,
        })
        
        if query2 is not None:
            ref_positions2 = get_ref_positions_for_read(cigar2, ref_start2)
            query2 = apply_glori_conversion_reverse(query2, ref_positions2, fwd_m6a_set)
            
            reads.append({
                "name": name,
                "query": query2,
                "cigar": cigar2,
                "ref_start": ref_start2,
                "strand": "-",
                "mate_start": start,
                "read_num": 2,
            })
    
    # Generate paired-end reads for reverse gene
    rev_exons = [(REV_EXON1_START, REV_EXON1_END), (REV_EXON2_START, REV_EXON2_END)]
    for start in range(REV_EXON1_START - 50, REV_EXON2_END + 50, step):
        # Skip reads entirely in intron
        if start >= REV_EXON1_END and start < REV_EXON2_START:
            continue
        
        # Read 1 (forward strand)
        query1, cigar1, ref_start1 = make_spliced_read(seq_chars, start, start + READ_LEN, rev_exons)
        if query1 is None:
            continue
        frag_id += 1
        name = f"frag_{frag_id:05d}"

        # Read 2 (reverse strand, mate)
        mate_start = start + INSERT_SIZE
        query2, cigar2, ref_start2 = make_spliced_read(seq_chars, mate_start, mate_start + READ_LEN, rev_exons)
        
        # Apply GLORI conversion to reverse gene reads (T→C for reverse complement)
        ref_positions1 = get_ref_positions_for_read(cigar1, ref_start1)
        query1 = apply_glori_conversion_reverse(query1, ref_positions1, rev_m6a_set)
        
        reads.append({
            "name": name,
            "query": query1,
            "cigar": cigar1,
            "ref_start": ref_start1,
            "strand": "+",
            "mate_start": mate_start if query2 else None,
            "read_num": 1,
        })
        
        if query2 is not None:
            ref_positions2 = get_ref_positions_for_read(cigar2, ref_start2)
            query2 = apply_glori_conversion_forward(query2, ref_positions2, rev_m6a_set)
            
            reads.append({
                "name": name,
                "query": query2,
                "cigar": cigar2,
                "ref_start": ref_start2,
                "strand": "-",
                "mate_start": start,
                "read_num": 2,
            })
    
    # Write BAM
    temp_bam = bam_path + ".unsorted.bam"
    with pysam.AlignmentFile(temp_bam, "wb", header=header) as bam:
        for r in reads:
            aln = pysam.AlignedSegment(bam.header)
            aln.query_name = r["name"]
            aln.query_sequence = r["query"]
            aln.is_reverse = (r["strand"] == "-")
            aln.is_read1 = (r["read_num"] == 1)
            aln.is_read2 = (r["read_num"] == 2)
            aln.is_paired = True
            aln.mapping_quality = 60
            aln.reference_id = 0
            aln.reference_start = r["ref_start"]
            aln.cigar = r["cigar"]
            aln.query_qualities = pysam.qualitystring_to_array("I" * len(r["query"]))
            
            if r["mate_start"] is not None:
                aln.next_reference_id = 0
                aln.next_reference_start = r["mate_start"]
                aln.template_length = r["mate_start"] - r["ref_start"]
                aln.mate_is_reverse = (r["strand"] == "+")  # opposite of this read
                aln.is_proper_pair = True
            
            bam.write(aln)
    
    # Sort and index
    pysam.sort("-o", bam_path, temp_bam)
    os.remove(temp_bam)
    pysam.index(bam_path)

    print("[m6A GLORI] wrote:")
    print(f"  {genome_fa}")
    print(f"  {annot_gb}")
    print(f"  {sites_bed}")
    print(f"  {bam_path}")
    print(f"[m6A GLORI] contig={CONTIG} length={LENGTH} reads={len(reads)}")
    print(f"[m6A GLORI] Forward gene m6A sites (A→G): {FWD_M6A_SITES}")
    print(f"[m6A GLORI] Reverse gene m6A sites (T→C): {REV_M6A_SITES}")


if __name__ == "__main__":
    main()

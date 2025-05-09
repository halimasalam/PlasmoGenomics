#!/bin/bash
# Usage: bash fst_pipeline.sh <ref_fasta> <high_bam_list> <low_bam_list> <gff_file> <output_prefix>

set -euo pipefail

# Input arguments
REF="$1"
HIGH_BAMS="$2"
LOW_BAMS="$3"
GFF_FILE="$4"
PREFIX="$5"

# Tools
ANGSD="angsd"
BEDTOOLS="bedtools"

# Extract long contigs (≥50kb) from .fai
awk '$2 >= 50000' "${REF}.fai" | cut -f1 > "${PREFIX}_long_contigs.txt"

# Generate SAF files
$ANGSD -b "$HIGH_BAMS" -ref "$REF" -anc "$REF" -out "${PREFIX}_high" \
  -dosaf 1 -GL 1 -P 4 -minMapQ 30 -minQ 20 -rf "${PREFIX}_long_contigs.txt"

$ANGSD -b "$LOW_BAMS" -ref "$REF" -anc "$REF" -out "${PREFIX}_low" \
  -dosaf 1 -GL 1 -P 4 -minMapQ 30 -minQ 20 -rf "${PREFIX}_long_contigs.txt"

# Estimate site frequency spectrum (SFS)
$ANGSD/misc/realSFS "${PREFIX}_high.saf.idx" "${PREFIX}_low.saf.idx" -P 4 > "${PREFIX}_sfs.sfs"

# Calculate FST
$ANGSD/misc/realSFS fst index "${PREFIX}_high.saf.idx" "${PREFIX}_low.saf.idx" \
  -sfs "${PREFIX}_sfs.sfs" -fstout "${PREFIX}"

$ANGSD/misc/realSFS fst stats "${PREFIX}.fst.idx" > "${PREFIX}_global_fst.txt"
$ANGSD/misc/realSFS fst stats2 "${PREFIX}.fst.idx" -win 50000 -step 10000 > "${PREFIX}_windowed_fst.txt"

# Clean FST output
cut -f2- "${PREFIX}_windowed_fst.txt" > "${PREFIX}_clean_fst.txt"

# Intersect top 1% outliers (generated later in R) with GFF
# Must match formatting of outlier BED from R
$BEDTOOLS intersect -a "${PREFIX}_top1_percent_windows.bed" -b "$GFF_FILE" -wa -wb > "${PREFIX}_outlier_genes.tsv"

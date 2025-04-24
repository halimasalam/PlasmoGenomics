#!/bin/bash

# This script calculates the mutational burden per isolate

# Input VCF file
VCF_FILE="/Users/halimaabdulsalam/Plasmodium-WGS-ARISE/PlasmoGenomics_copy/plasmo_genomics/data/haploid_variants_ann.vcf"
# VCF_FILE="data/haploid_varaints_ann.vcf"
OUTPUT_DIR="analysis_results/mutation_burden"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Query mutation data from VCF
bcftools query -l $VCF_FILE | paste - <(bcftools query -f '[%GT\t]\n' $VCF_FILE | awk '{for(i=1;i<=NF;i++) if($i=="1") count[i]++} END {for(i in count) print count[i]}') | awk '{print $1 " " $2}' > "$OUTPUT_DIR/mutation_burden.txt"

echo "Mutational burden file saved as '$OUTPUT_DIR/mutation_burden.txt'"

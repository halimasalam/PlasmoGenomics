#!/usr/bin/env bash

# PCA Preprocessing Pipeline with PLINK
# Usage: ./run_pca.sh -i high_quality_variants.vcf -o results/pca_output

set -euo pipefail

# Colors
GREEN="\033[0;32m"
NC="\033[0m"

# Usage message
usage() {
  echo -e "Usage: $0 -i <input.vcf> -o <output_dir>"
  exit 1
}

# Parse flags
while getopts ":i:o:" opt; do
  case ${opt} in
    i ) input_vcf=$OPTARG ;;
    o ) output_dir=$OPTARG ;;
    * ) usage ;;
  esac
done

# Validate inputs
if [ -z "${input_vcf:-}" ] || [ -z "${output_dir:-}" ]; then
  usage
fi

# Resolve absolute paths
input_vcf=$(realpath "$input_vcf")
output_dir=$(realpath "$output_dir")

# Create output dir if it doesn't exist
mkdir -p "$output_dir"

# Derive prefix from input filename
prefix=$(basename "$input_vcf" .vcf)

echo -e "${GREEN}[INFO] Input VCF: ${input_vcf}${NC}"
echo -e "${GREEN}[INFO] Output directory: ${output_dir}${NC}"
echo -e "${GREEN}[INFO] Prefix for PLINK files: ${prefix}${NC}"

cd "$output_dir"

# Step 1: Convert VCF to PLINK binary format
plink --vcf "$input_vcf" --snps-only just-acgt --make-bed --out "${prefix}" --allow-extra-chr

# Step 2: Filter SNPs by missingness
plink --bfile "${prefix}" --geno 0.9 --make-bed --out filtered_data --allow-extra-chr

# Step 3: LD pruning
plink --bfile filtered_data --indep-pairwise 50 5 0.8 --out pruned_data --allow-extra-chr

# Step 4: Extract pruned SNPs
plink --bfile filtered_data --extract pruned_data.prune.in --make-bed --out pruned_variants --allow-extra-chr

# Step 5: Run PCA
plink --bfile pruned_variants --pca --out pca_results --allow-extra-chr

# Step 6: Missingness
plink --bfile "${prefix}" --missing --allow-extra-chr

# Step 7: Chromosome-wise PCA (if applicable)
# Extract standard chromosomes
awk '{if ($1 ~ /^NC_/) print $1}' pruned_variants.bim | sort -u > chr_list.txt

while read chr; do
    echo "🔹 Running PCA for $chr"
    plink --bfile pruned_variants --chr "$chr" --pca 10 --out "pca_$chr" --allow-extra-chr \
        || echo "⚠️ Skipped $chr due to missing genotype data."
done < chr_list.txt

echo -e "${GREEN}[DONE] PCA pipeline finished. Results in ${output_dir}${NC}"




#!/bin/bash

# ----------------------------
# Script: rehh.sh
# Description: Prepares .hap and .map files from phased VCF for rehh iHS analysis
# Usage:
#   bash rehh.sh phased_data.vcf.gz /path/to/ref.fna output_dir
# ----------------------------

set -e  # Exit on error

# Input arguments
VCF_GZ="$1"
REF="$2"
OUTDIR="$3"

# Validate inputs
if [[ -z "$VCF_GZ" || -z "$REF" || -z "$OUTDIR" ]]; then
    echo "❗ Usage: bash rehh.sh phased_data.vcf.gz /path/to/ref.fna output_dir"
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p "$OUTDIR"

# Extract all unique chromosomes from the VCF
echo "🔍 Extracting chromosome list from $VCF_GZ..."
CHROM_LIST=$(bcftools query -f '%CHROM\n' "$VCF_GZ" | sort -u)

# Loop over chromosomes and process each one
for CHR in $CHROM_LIST; do
    echo "🔄 Processing chromosome: $CHR"
    
    # Subset the VCF to current chromosome
    CHR_VCF="${OUTDIR}/${CHR}.vcf.gz"
    bcftools view -r "$CHR" "$VCF_GZ" | bgzip -c > "$CHR_VCF"
    tabix -p vcf "$CHR_VCF"

    # Run your vcf2rehh Python script
    python vcf2rehh.py \
        --vcf "$CHR_VCF" \
        --ref "$REF" \
        --out "${OUTDIR}/${CHR}"
done

echo "✅ rehh processing complete! Files saved to: $OUTDIR"
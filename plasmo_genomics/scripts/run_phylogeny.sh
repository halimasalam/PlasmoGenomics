#!/bin/bash

# ------------------------------
# Phylogenetic Analysis Pipeline
# ------------------------------

# # Set input VCF file and output directory
VCF_INPUT="data/haploid_variants.vcf"
OUTPUT_DIR="analysis_results/phylogeny"
VCF2PHYLIP_PATH="scripts/vcf2phylip.py" # Path to vcf2phylip conversion script

# Check that the input VCF file exists
if [ ! -f "$VCF_INPUT" ]; then
  echo "Error: Input VCF file not found at: $VCF_INPUT"
  exit 1
fi

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# ---------------------------------------
# Step 1: Filter VCF using vcftools
echo "Running vcftools filtering..."

# vcftools filters:
# --max-missing 0.3     → Keep variants in ≥30% of samples
# --minQ 30             → Remove variants with QUAL < 30
# --min-meanDP 5        → Keep variants with average DP ≥ 5

vcftools --vcf "$VCF_INPUT" \
  --max-missing 0.3 \
  --minQ 30 \
  --min-meanDP 5 \
  --recode \
  --out "$OUTPUT_DIR/high_quality_variants"

# Check if filtering worked and output file exists
FILTERED_VCF="$OUTPUT_DIR/high_quality_variants.recode.vcf"
if [ ! -f "$FILTERED_VCF" ]; then
  echo "Error: Filtered VCF was not generated. Exiting."
  exit 1
fi

# ---------------------------------------
# Step 2: Convert filtered VCF to PHYLIP
# vcf2phylip converts multi-sample VCF into an alignment for phylogenetics
echo "Converting to PHYLIP format..."
python "$VCF2PHYLIP_PATH" \
       -i "$OUTPUT_DIR/high_quality_variants.recode.vcf" \
       -o "$OUTPUT_DIR/high_quality_variants.phy"

# ---------------------------------------
# Step 3: Build tree using IQ-TREE
echo "Running IQ-TREE..."
cd "$OUTPUT_DIR" || exit

# -m GTR+G   → GTR substitution model with gamma rate heterogeneity
# -bb 1000   → Ultrafast bootstrap with 1000 replicates
# -alrt 1000 → SH-like approximate likelihood ratio test
# -nt AUTO   → Automatically choose optimal number of threads

iqtree \
  -s "high_quality_variants.recode.min4.phy" \
  -m "GTR+G" \
  -bb 1000 \
  -alrt 1000 \
  -nt AUTO                          











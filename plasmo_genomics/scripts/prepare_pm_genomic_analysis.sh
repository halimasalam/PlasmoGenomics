#!/bin/bash

set -euo pipefail

# Define directories and filenames
WORK_DIR="analysis_results/variants_analysis/files/"  # Set this to your working directory
VCF_INPUT="data/haploid_variants.vcf"  # Modify with the correct raw VCF file input path
VCF_OUTPUT="${WORK_DIR}/haploid_annotated_variants.vcf"

# Ensure VCF input file exists
if [ ! -f "$VCF_INPUT" ]; then
    echo "ERROR: Input VCF file does not exist at $VCF_INPUT"
    exit 1
fi

echo "Step 1: Generate High-Quality Annotated VCF (if not provided)"

# Placeholder command for generating the high-quality annotated VCF, adjust as per your workflow
bcftools view -i 'QUAL>=30' $VCF_INPUT -o $VCF_OUTPUT

echo "High-quality annotated VCF generated: $VCF_OUTPUT"


echo "Step 2: Minor Allele Frequency (MAF) Filtering"
# Filter variants with sample coverage (AN >= 5)
bcftools view -i 'AN>=5' $VCF_OUTPUT -o $WORK_DIR/filtered_AN5.vcf

# Keep only variants with MAF between 0.05 and 0.95
bcftools view -i 'AC/AN>0.05 & AC/AN<0.95' $WORK_DIR/filtered_AN5.vcf -o $WORK_DIR/filtered_AN5_MAF05.vcf

echo "MAF filtering complete."


echo "Step 3: Identity-by-Descent (IBD) Calculation"
# Use PLINK to calculate pairwise IBD
plink --bfile pruned_variants_filtered --genome --out $WORK_DIR/ibd_pm --chr-set 52 --allow-extra-chr

# Count pairs with PI_HAT > 0
echo "Pairs with PI_HAT > 0:"
awk '$10 > 0' $WORK_DIR/ibd_pm.genome | wc -l

echo "IBD calculation complete."


echo "Step 4: LD Analysis Prep"
# Create unpruned BED data from raw VCF
plink --vcf $VCF_OUTPUT --snps-only just-acgt --make-bed --out $WORK_DIR/ld_raw_data --allow-extra-chr

# Apply filters for LD (MAF > 0.05, missing rate < 10%)
plink --bfile $WORK_DIR/ld_raw_data --geno 0.9 --maf 0.05 --make-bed --out $WORK_DIR/ld_filtered --allow-extra-chr

# Compute LD (R²) with large window
plink --bfile $WORK_DIR/ld_filtered --r2 --ld-window-kb 1000 --ld-window 999999 --ld-window-r2 0.2 --out $WORK_DIR/ld_results_filtered --allow-extra-chr

echo "LD analysis setup complete."


echo "Step 5: iHS Analysis Prep"

# Download Beagle for phasing (if not already downloaded)
if [ ! -f beagle.28Jun21.220.jar ]; then
    echo "Downloading Beagle..."
    curl -O https://faculty.washington.edu/browning/beagle/beagle.28Jun21.220.jar
else
    echo "Beagle already downloaded."
fi

# Filter contigs with sufficient SNPs (for Pv specifically)
echo "Filtering contigs with sufficient SNPs..."
bcftools query -f '%CHROM\n' $VCF_OUTPUT | sort | uniq -c | sort -nr > contig_counts.txt
awk '$1 >= 5 {print $2}' contig_counts.txt > valid_contigs.txt
awk '{print $0"\t0\t999999999"}' valid_contigs.txt > regions.bed

bcftools view -R regions.bed $VCF_OUTPUT -o $WORK_DIR/high_quality_fv.vcf -Ov

# Phase the VCF using Beagle
echo "Phasing VCF with Beagle..."
java -jar beagle.28Jun21.220.jar gt=$WORK_DIR/high_quality_fv.vcf out=phased_data nthreads=4

# Compress and index phased data
bgzip phased_data.vcf
tabix -p vcf phased_data.vcf.gz

# Run rehh chromosome-wise using rehh.sh
bash rehh.sh

echo "All preprocessing steps complete."

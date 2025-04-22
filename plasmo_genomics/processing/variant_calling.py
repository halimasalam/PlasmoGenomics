# Script for Variant Calling Using FreeBayes
# This script takes BAM files with read groups added, runs variant calling using FreeBayes,
# filters the results based on quality, and outputs the final filtered VCF.
# It checks for the required tools and input files before starting the process.

import subprocess
import sys
from pathlib import Path
from plasmo_genomics.utils.helper import ensure_dir, check_tool_installed

def run_variant_calling(input_dir, output_dir, reference):
    # Convert string paths to Path objects and expand ~ to full user directory
    input_path = Path(input_dir).expanduser()
    output_path = Path(output_dir).expanduser()
    reference = Path(reference).expanduser()

    # Make sure the input directory exists
    if not input_path.exists():
        print(f"Input directory '{input_path}' does not exist.")
        sys.exit(1)

    # Check if the reference genome file is present
    if not reference.exists():
        print(f"Reference genome '{reference}' not found.")
        sys.exit(1)

    # Look for BAM files that already have read groups added
    bams = list(input_path.glob("*_rg.bam"))
    if not bams:
        print(f"No *_rg.bam files found in {input_path}")
        sys.exit(1)

    # Create the output directory if it doesn’t exist
    ensure_dir(output_path)

    # Confirm required tools are installed
    check_tool_installed("freebayes")
    check_tool_installed("bcftools")

    # Set file paths
    raw_vcf = output_path / "variants_raw.vcf"
    output_vcf = output_path / "variants.vcf"  # Final filtered output

    print(f"🔹 Calling variants for all samples...")

    # Run FreeBayes to perform variant calling on all BAMs
    subprocess.run([
        "freebayes",
        "-f", str(reference),
        "-C", "2",
        "--min-base-quality", "20",
        "--min-mapping-quality", "30",
        "--genotype-qualities",
        "--min-coverage", "5",
        "--vcf", str(raw_vcf),
        *map(str, bams)
    ], check=True)

    print("✅ Variant calling completed. Filtering by QUAL > 30...")

    # Filter VCF based on quality score and overwrite the final output
    subprocess.run([
        "bcftools", "filter",
        "-i", "QUAL > 30",
        "-o", str(output_vcf),
        str(raw_vcf)
    ], check=True)

    print(f"✅ Filtering complete. Final VCF saved to {output_vcf}")


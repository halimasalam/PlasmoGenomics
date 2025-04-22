# Script to Add Read Groups to BAM Files Using GATK
# This script processes BAM files that have been marked for duplicates, adds read group information 
# using GATK's AddOrReplaceReadGroups, and then indexes the resulting BAMs with samtools.

import subprocess
from pathlib import Path
import sys
from plasmo_genomics.utils.helper import ensure_dir, check_tool_installed

def add_read_groups(input_dir, output_dir):
    # Convert input/output paths and expand ~ to full paths
    input_path = Path(input_dir).expanduser()
    output_path = Path(output_dir).expanduser()

    # Ensure the input directory exists
    if not input_path.exists():
        print(f"Input directory '{input_path}' does not exist.")
        sys.exit(1)

    # Find BAM files that have been marked for duplicates
    bams = list(input_path.glob("*_markdup.bam"))
    if not bams:
        print(f"No *_markdup.bam files found in {input_path}")
        sys.exit(1)

    # Create the output directory if needed
    ensure_dir(output_path)

    # Check if required tools are available
    check_tool_installed("gatk")
    check_tool_installed("samtools")

    # Loop through each BAM file and add read group information
    for bam in bams:
        # Derive sample name by removing the suffix
        sample = bam.stem.replace("_markdup", "")
        out_bam = output_path / f"{sample}_rg.bam"

        print(f"Adding read groups for {sample}...")

        # Run GATK to add read groups
        subprocess.run([
            "gatk", "AddOrReplaceReadGroups",
            "-I", str(bam),
            "-O", str(out_bam),
            "-RGID", sample,       # Read Group ID
            "-RGLB", "lib1",       # Library name
            "-RGPL", "ILLUMINA",   # Platform (e.g., ILLUMINA)
            "-RGPU", "unit1",      # Platform unit
            "-RGSM", sample        # Sample name
        ], check=True)
        
        # Index the new BAM file
        subprocess.run(["samtools", "index", str(out_bam)], check=True)
    
    print("✅ Read groups added and indexed.")

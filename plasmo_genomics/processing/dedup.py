# Script to Remove PCR Duplicates from BAM Files Using Samtools
# This step takes sorted BAM files, performs name sorting, fixmate correction,
# coordinate sorting, and marks/removes PCR duplicates. Final BAMs are indexed.

import subprocess
from pathlib import Path
import sys
from plasmo_genomics.utils.helper import ensure_dir, check_tool_installed

def run_dedup(input_dir, output_dir):
    # Expand input/output paths
    input_path = Path(input_dir).expanduser()
    output_path = Path(output_dir).expanduser()

    # Make sure the input directory exists
    if not input_path.exists():
        print(f"Input directory '{input_path}' does not exist.")
        sys.exit(1)

    # Find BAM files that were previously sorted
    bams = list(input_path.glob("*_sorted.bam"))
    if not bams:
        print(f"No *_sorted.bam files found in {input_path}")
        sys.exit(1)

    # Ensure output directory exists
    ensure_dir(output_path)

    # Confirm samtools is installed
    check_tool_installed("samtools")

    for bam in bams:
        sample = bam.name.replace("_sorted.bam", "")
        print(f"Processing {sample}...")

        # Temporary intermediate file paths
        nsorted = output_path / f"{sample}_nsorted.bam"
        fixmate = output_path / f"{sample}_fixmate.bam"
        sorted_fixmate = output_path / f"{sample}_fixmate_sorted.bam"
        final_bam = output_path / f"{sample}_markdup.bam"

        # Step 1: Name sort the BAM
        subprocess.run(["samtools", "sort", "-n", "-o", nsorted, bam], check=True)
        # Step 2: Fix mate information
        subprocess.run(["samtools", "fixmate", "-m", nsorted, fixmate], check=True)
        # Step 3: Coordinate sort the fixmate BAM
        subprocess.run(["samtools", "sort", "-o", sorted_fixmate, fixmate], check=True)
        # Step 4: Mark and remove duplicates
        subprocess.run(["samtools", "markdup", "-s", "-r", sorted_fixmate, final_bam], check=True)
        # Step 5: Index the deduplicated BAM
        subprocess.run(["samtools", "index", final_bam], check=True)

        # Clean up intermediate files
        for tmp in [nsorted, fixmate, sorted_fixmate]:
            tmp.unlink(missing_ok=True)

        print(f"Done with {sample}")
    
    print(" Deduplication complete for all samples.")

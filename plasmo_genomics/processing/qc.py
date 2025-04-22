
# Script to Trim Paired-End FASTQ Files Using Trim Galore
# This script scans for paired-end FASTQ files in the input directory, runs Trim Galore 
# to trim adapters and low-quality bases, and saves the cleaned reads in the output directory.

import subprocess
from pathlib import Path
import sys
from plasmo_genomics.utils.helper import ensure_dir, check_tool_installed

def run_trim_galore(input_dir, output_dir):
    # Convert input/output paths and expand ~ to full paths
    input_path = Path(input_dir).expanduser()
    output_path = Path(output_dir).expanduser()
    
    # Make sure the input directory exists
    if not input_path.exists():
        print(f"Input directory '{input_path}' does not exist.")
        sys.exit(1)

    # Look for R1 FASTQ files (paired-end reads)
    fastq_files = list(input_path.glob("*_1.fastq.gz"))
    if not fastq_files:
        print(f"No FASTQ files matching *_1.fastq.gz found in {input_path}")
        sys.exit(1)

    # Create the output directory if needed
    ensure_dir(output_path)

    # Check that Trim Galore is available
    check_tool_installed("trim_galore")
    
    # Loop through all R1 files and find their matching R2 files
    for r1 in fastq_files:
        # Construct the name of the paired R2 file
        r2 = r1.with_name(r1.name.replace("_1.fastq.gz", "_2.fastq.gz"))

        # Skip this pair if the R2 file doesn't exist
        if not r2.exists():
            print(f"Skipping {r1.name} — matching R2 file not found: {r2.name}")
            continue

        print(f"Trimming {r1.name} and {r2.name}")
        # Run Trim Galore in paired-end mode
        subprocess.run([
            "trim_galore", "--paired", str(r1), str(r2), "-o", str(output_path)
        ], check=True)

    print("✅ All samples trimmed successfully!")

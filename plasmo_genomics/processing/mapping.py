# Script to Align Paired-End Reads Using BWA-MEM2 and Sort with Samtools
# This script aligns quality-trimmed FASTQ files to a reference genome using BWA-MEM2, 
# converts the output to BAM format, sorts it, and creates BAM indexes.

import subprocess
from pathlib import Path
import sys
from plasmo_genomics.utils.helper import ensure_dir, check_tool_installed

def run_alignment(input_dir, output_dir, reference, threads=4):
    # Expand and normalize paths
    input_path = Path(input_dir).expanduser()
    output_path = Path(output_dir).expanduser()
    reference = Path(reference).expanduser()

    # Check for input directory and reference genome
    if not input_path.exists():
        print(f"Input directory '{input_path}' does not exist.")
        sys.exit(1)
    if not reference.exists():
        print(f"Reference genome '{reference}' not found.")
        sys.exit(1)

    # Find trimmed R1 FASTQ files (output from Trim Galore)
    fastqs = list(input_path.glob("*_1_val_1.fq.gz"))
    if not fastqs:
        print(f"No *_1_val_1.fq.gz files found in {input_path}")
        sys.exit(1)

    # Make sure output directory exists
    ensure_dir(output_path)

    # Confirm required tools are available
    check_tool_installed("bwa-mem2")
    check_tool_installed("samtools")

    # Loop through each pair of R1 and R2 reads
    for r1 in fastqs:
        r2 = r1.with_name(r1.name.replace("_1_val_1.fq.gz", "_2_val_2.fq.gz"))
        if not r2.exists():
            print(f"Skipping {r1.name}: matching {r2.name} not found")
            continue

        # Derive sample name from file name
        sample = r1.name.replace("_1_val_1.fq.gz", "")
        out_bam = output_path / f"{sample}_sorted.bam"

        print(f"Aligning {sample}...")

        # Align reads, convert to BAM, and sort using a shell pipeline
        subprocess.run(
            f"bwa-mem2 mem -t {threads} {reference} {r1} {r2} | "
            f"samtools view -bS | "
            f"samtools sort -o {out_bam}",
            shell=True, check=True
        )

        # Index the sorted BAM file
        subprocess.run(["samtools", "index", str(out_bam)], check=True)

    print("All samples aligned and indexed.")

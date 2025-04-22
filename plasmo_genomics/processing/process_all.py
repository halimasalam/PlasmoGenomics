# Main Pipeline Script for Plasmodium Genomics Workflow
# This script runs the full genomic data processing pipeline: QC, mapping, deduplication,
# read group addition, and variant calling. It organizes outputs step-by-step into subfolders.

import argparse
from pathlib import Path

# Import each stage of the pipeline
from plasmo_genomics.processing.qc import run_trim_galore
from plasmo_genomics.processing.mapping import run_alignment
from plasmo_genomics.processing.dedup import run_dedup
from plasmo_genomics.processing.readgroups import add_read_groups
from plasmo_genomics.processing.variant_calling import run_variant_calling
from plasmo_genomics.processing.annotation import run_snpeff

def run_pipeline(
    input_dir: str,
    output_dir: str,
    reference: str,
    threads: int = 4
):
    # Convert inputs to full path objects
    input_path = Path(input_dir).expanduser()
    output_path = Path(output_dir).expanduser()
    reference_path = Path(reference).expanduser()

    # Check that the input directory exists
    if not input_path.exists():
        print(f"Input directory '{input_path}' does not exist.")
        return

    # Step 1: Quality control / trimming
    print("Step 1: Running QC...")
    qc_output = output_path / "01_qc"
    run_trim_galore(input_path, qc_output)

    # Step 2: Align trimmed reads to the reference genome
    print("Step 2: Running Mapping...")
    mapping_output = output_path / "02_mapping"
    run_alignment(qc_output, mapping_output, reference_path, threads)

    # Step 3: Remove PCR duplicates from the mapped BAM files
    print("Step 3: Removing Duplicates...")
    dedup_output = output_path / "03_dedup"
    run_dedup(mapping_output, dedup_output)

    # Step 4: Add read group information to BAM files
    print("Step 4: Adding Read Groups...")
    read_group_output = output_path / "04_readgroups"
    add_read_groups(dedup_output, read_group_output, reference_path)

    # Step 5: Perform variant calling
    print("Step 5: Variant Calling...")
    variant_output = output_path / "high_quality_variants.vcf"
    run_variant_calling(read_group_output, variant_output, reference_path)

    # Step 6: Annotate variants with SnpEff
    print("Step 6: Annotating Variants...")
    annotation_output = output_path / "high_quality_annotated_variants.vcf"
    run_snpeff(variant_output, annotation_output)


    print("Pipeline completed successfully!")

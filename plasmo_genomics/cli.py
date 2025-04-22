# plasmo_genomics/cli.py

import argparse
from .processing import qc, mapping, dedup, readgroups, variant_calling, annotation, process_all
import subprocess
import os
from pathlib import Path

def run_pca(args):
    import subprocess
    script_path = Path(__file__).parent / "scripts" / "run_pca.sh"
    
    if not os.access(script_path, os.X_OK):
        print("⚠️ Script not executable, fixing permissions...")
        os.chmod(script_path, 0o755)

    print(f"Running PCA shell pipeline...")
    subprocess.run(["bash", str(script_path), "-i", args.input, "-o", args.output], check=True)

    print("✅ PCA complete. All files and plots saved in:", args.output)

def run_admixture(input_prefix, chr_list_path, chr_map_path, output_dir=None):
    script_path = Path(__file__).parent /  "scripts" / "run_admixture_pipeline.sh"

    if not os.access(script_path, os.X_OK):
        print("⚠️ Script not executable, fixing permissions...")
        os.chmod(script_path, 0o755)

    print(f"Running admixture shell pipeline...")

    # Run the shell script with necessary environment vars
    cmd = [
        "bash",
        str(script_path),
        str(input_prefix),
        str(chr_list_path),
        str(chr_map_path)
    ]

    if output_dir:
        cmd.append(str(output_dir))

    subprocess.run(cmd, check=True)


def main():
    parser = argparse.ArgumentParser(description="PlasmoGenomics CLI Tool")
    subparsers = parser.add_subparsers(dest="command")

    # QC Subcommand
    parser_qc = subparsers.add_parser("qc", help="Run TrimGalore on paired-end reads")
    parser_qc.add_argument("-i", "--input", required=True, help="Input directory containing *_1.fastq.gz files")
    parser_qc.add_argument("--output", required=True, help="Output directory for trimmed files")

    # Alignment Subcommand
    parser_align = subparsers.add_parser("align", help="Run BWA-MEM2 on trimmed reads")
    parser_align.add_argument("-i", "--input", required=True, help="Input directory with trimmed reads")
    parser_align.add_argument("-o", "--output", required=True, help="Output directory for BAMs")
    parser_align.add_argument("-r", "--ref", required=True, help="Reference genome FASTA file")
    parser_align.add_argument("-t", "--threads", type=int, default=2, help="Number of threads for BWA")

    # Dedup
    parser_dedup = subparsers.add_parser("dedup", help="Remove duplicate reads using samtools")
    parser_dedup.add_argument("-i", "--input", required=True, help="Input dir with sorted BAMs")
    parser_dedup.add_argument("-o", "--output", required=True, help="Output dir for deduped BAMs")

    # Read Groups
    parser_rg = subparsers.add_parser("readgroups", help="Add read groups using GATK")
    parser_rg.add_argument("-i", "--input", required=True, help="Input dir with markdup BAMs")
    parser_rg.add_argument("-o", "--output", required=True, help="Output dir for BAMs with read groups")

    # Variant Calling
    parser_vc = subparsers.add_parser("variantcall", help="Run FreeBayes and merge VCFs")
    parser_vc.add_argument("-i", "--input", required=True, help="Input dir with *_rg.bam")
    parser_vc.add_argument("-o", "--output", required=True, help="Output dir for per-sample VCFs")
    parser_vc.add_argument("-r", "--ref", required=True, help="Reference genome FASTA")

    #Annotation
    parser_annotate = subparsers.add_parser("annotate", help="Annotate variants with SnpEff")
    parser_annotate.add_argument("-i", "--input", required=True, help="Input VCF file")
    parser_annotate.add_argument("-o", "--output", required=True, help="Output annotated VCF file")

    # Process all 
    parser_process_all = subparsers.add_parser("process_all", help="Run full PlasmoGenomics pipeline")
    parser_process_all.add_argument("-i", "--input", required=True, help="Input directory with raw FASTQ files")
    parser_process_all.add_argument("-o", "--output", required=True, help="Output directory")
    parser_process_all.add_argument("-r", "--reference", required=True, help="Reference genome FASTA file")
    parser_process_all.add_argument("-t", "--threads", type=int, default=4, help="Number of threads (default: 4)")

    # PCA command
    parser_pca = subparsers.add_parser("pca", help="Run PLINK-based PCA pipeline")
    parser_pca.add_argument("-i", "--input", required=True, help="VCF input file")
    parser_pca.add_argument("-o", "--output", required=True, help="Output directory")

    # Admixture command
    admix_parser = subparsers.add_parser("admixture", help="Run ADMIXTURE pipeline")
    admix_parser.add_argument("-i", "--input", required=True, help="Prefix for PLINK .bed/.bim/.fam (pruned files)")
    admix_parser.add_argument("-c", "--chr-list", required=True, help="File with allowed chromosomes")
    admix_parser.add_argument("-m", "--chr-map", required=True, help="Chromosome name → ID map")
    admix_parser.add_argument("-o", "--output", help="Optional output directory")

    

    args = parser.parse_args()

    if args.command == "qc":
        qc.run_trim_galore(args.input, args.output)
    elif args.command == "align":
        mapping.run_alignment(args.input, args.output, args.ref, args.threads)
    elif args.command == "dedup":
        dedup.run_dedup(args.input, args.output)
    elif args.command == "readgroups":
        readgroups.add_read_groups(args.input, args.output)
    elif args.command == "variantcall":
        variant_calling.run_variant_calling(args.input, args.output, args.ref)
    elif args.command == "annotate":
        annotation.run_snpeff(args.input, args.output)
    elif args.command == "process_all":
        process_all.run_pipeline(args.input, args.output, args.reference, args.threads)
    elif args.command == "pca":
        run_pca(args)
    elif args.command == "admixture":
        run_admixture(args.input, args.chr_list, args.chr_map, args.output)
    else:
        parser.print_help()

if __name__ == "__main__":
    main()

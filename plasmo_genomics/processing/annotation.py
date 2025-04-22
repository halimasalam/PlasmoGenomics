import subprocess
from pathlib import Path
from plasmo_genomics.utils.helper import ensure_dir, check_tool_installed
import sys

def run_snpeff(vcf_path, output_path, species="Plasmodium_malariae"):
    """
    Annotates a VCF file using SnpEff with a specified genome database.

    Parameters:
    - vcf_path (str or Path): Path to the input VCF file (unannotated).
    - output_path (str or Path): Path to save the annotated VCF output (file, not directory).
    - species (str): Name of the genome database defined in snpEff.config.
                     Default is 'Plasmodium_malariae'.
    """
    vcf_path = Path(vcf_path).expanduser()
    output_path = Path(output_path).expanduser()

    # Check if VCF file exists
    if not vcf_path.exists():
        print(f"❌ Input file '{vcf_path}' does not exist.")
        sys.exit(1)

    # Create parent directory for the output file
    ensure_dir(output_path.parent)

    # Confirm snpEff is available
    check_tool_installed("snpEff")

    print(f"🛠️ Annotating {vcf_path.name} with SnpEff database '{species}'...")

    command = [
        "snpEff", "ann", species, str(vcf_path)
    ]

    try:
        with open(output_path, "w") as out_vcf:
            subprocess.run(command, stdout=out_vcf, stderr=subprocess.PIPE, check=True)
        print(f"✅ Annotation complete. Output saved to: {output_path}")
    except subprocess.CalledProcessError as e:
        print(f"❌ SnpEff annotation failed:\n{e.stderr.decode()}")
        raise

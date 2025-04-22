#!/bin/bash

# === SETUP CUSTOM SNPEFF DATABASE FOR PLASMODIUM MALARIAE ===
# This script builds a Plasmodium malariae genome database for use with SnpEff

# ----------------------------------------
# ❓ What is $CONDA_PREFIX?
# It is an environment variable that points to your current conda environment path.
# For example: /Users/yourname/miniconda3/envs/bio_env
# It allows this script to work inside your current conda environment.
# ----------------------------------------

# USAGE:
# ./setup_snpeff.sh path/to/genome.gff path/to/genome.fasta

set -e  # Exit if any command fails

# Input files
GFF=$1
FA=$2
DB_NAME="Plasmodium_malariae"

if [[ -z "$GFF" || -z "$FA" ]]; then
    echo "Usage: $0 path/to/genome.gff path/to/genome.fasta"
    exit 1
fi

echo "[*] Using conda prefix at: $CONDA_PREFIX"

# === Add entry to snpEff.config if it doesn't exist ===
CONFIG_PATH=$(find $CONDA_PREFIX/share -name "snpEff.config" | head -n1)
if ! grep -q "$DB_NAME.genome" "$CONFIG_PATH"; then
    echo "[*] Adding '$DB_NAME' entry to snpEff.config..."
    echo "${DB_NAME}.genome : Plasmodium malariae" >> "$CONFIG_PATH"
fi

# === Create genome folder ===
GENOME_DIR=$(find $CONDA_PREFIX/share -type d -name "snpeff*" | head -n1)/data/$DB_NAME
mkdir -p "$GENOME_DIR"

# === Copy input files ===
echo "[*] Copying GFF and FASTA..."
cp "$GFF" "$GENOME_DIR/genes.gff"
cp "$FA" "$GENOME_DIR/sequences.fa"

# === Extract CDS and protein sequences using gffread ===
echo "[*] Extracting CDS and proteins with gffread..."
gffread -g "$GENOME_DIR/sequences.fa" \
        -x "$GENOME_DIR/cds.fa" \
        -y "$GENOME_DIR/protein.fa" \
        "$GENOME_DIR/genes.gff"

# === Build the custom database ===
SNPEFF_DIR=$(dirname "$CONFIG_PATH")
cd "$SNPEFF_DIR"
echo "[*] Building the database..."
java -jar snpEff.jar build -gff3 -v "$DB_NAME"

# === Confirm database ===
echo "[*] Confirming database is available..."
java -jar snpEff.jar databases | grep "$DB_NAME"

echo -e "\n✅ Done. You can now annotate using:"
echo "snpEff ann $DB_NAME your_variants.vcf > annotated_variants.vcf"

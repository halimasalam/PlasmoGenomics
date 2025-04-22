#!/bin/bash

set -e

# ====== INPUT FILES FROM CLI ======
VCF_BFILE="$1"                  # e.g., pruned_variants
CHR_LIST="$2"                   # e.g., data/chr_list.txt
CHR_MAP_FILE="$3"               # e.g., data/chr_name_to_int.txt
OUTPUT_DIR="${4:-.}"            # optional 4th arg, defaults to current directory

mkdir -p "$OUTPUT_DIR"

# ====== Filenames without paths ======
BASENAME=$(basename "$VCF_BFILE")
FIXED_BIM="${OUTPUT_DIR}/${BASENAME}_fixed.bim"
FINAL_BFILE="${OUTPUT_DIR}/${BASENAME}_final"
FILTERED_BFILE="${OUTPUT_DIR}/${BASENAME}_filtered"
INT_BFILE="${OUTPUT_DIR}/${BASENAME}_integer"

# ====== 1. Fix .bim file ======
echo "Fixing missing variant IDs in BIM file..."
awk '{if ($2 == ".") $2 = $1":"$4}1' OFS='\t' ${VCF_BFILE}.bim > "$FIXED_BIM"

# ====== 2. Recreate PLINK files ======
echo "Rebuilding PLINK files..."
plink --bfile "$VCF_BFILE" --bim "$FIXED_BIM" --make-bed --out "$FINAL_BFILE" --allow-extra-chr

# ====== 3. Exclude unwanted chromosomes ======
echo "Filtering chromosomes..."
CHRS=$(tr '\n' ' ' < "$CHR_LIST")
plink --bfile "$FINAL_BFILE" --chr $CHRS --make-bed --out "$FILTERED_BFILE" --allow-extra-chr

# ====== 4. Map chromosome names to integers ======
if [ ! -f "$CHR_MAP_FILE" ]; then
    echo "ERROR: Chromosome mapping file not found at $CHR_MAP_FILE"
    exit 1
fi

echo "Applying chromosome name to integer ID mapping..."
awk -v mapfile="$CHR_MAP_FILE" '
BEGIN {
    while ((getline < mapfile) > 0) {
        chr_map[$1] = $2
    }
}
{
    if ($1 in chr_map) $1 = chr_map[$1];
    print
}' "${FILTERED_BFILE}.bim" > "${INT_BFILE}.bim"

cp "${FILTERED_BFILE}.bed" "${INT_BFILE}.bed"
cp "${FILTERED_BFILE}.fam" "${INT_BFILE}.fam"

# ====== 5. Run ADMIXTURE ======
echo "Running ADMIXTURE for K=1 to 10..."
cd "$OUTPUT_DIR"
rm -f log*.out cv_errors.txt
for K in {1..10}; do
  echo "Running K=$K..."
  admixture --cv "${INT_BFILE}.bed" $K -j4 | tee "log${K}.out"
done

# ====== 6. Extract CV errors ======
echo "Extracting CV errors..."
grep -h "CV error" log*.out > cv_errors.txt

echo "All done ✅"

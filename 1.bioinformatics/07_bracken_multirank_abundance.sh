#!/bin/bash
#SBATCH --job-name=taxa_multirank_nonhuman
#SBATCH --account=uoo03549
#SBATCH --partition=milan
#SBATCH --qos=normal
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --array=0-95%24
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

# 07_bracken_multirank_abundance.sh
# Extracts kingdom/phylum/family summaries from Kraken2 reports and estimates genus abundance with Bracken.

set -euo pipefail

module purge
module load Bracken/2.7-GCC-11.3.0

BASE="/home/soths/00_nesi_projects/uoo03549_nobackup/Sereyboth"
KRAKEN_REPORT_DIR="${BASE}/Host_Removal/Kraken2_nonhuman/reports"
DB="/home/soths/00_nesi_projects/lincoln04415_nobackup/Metagenomics/Databases/Kraken2"
OUT_BASE="${BASE}/Host_Removal/Taxonomy_MultiRank"
OUT_K="${OUT_BASE}/kingdom_kraken"
OUT_P="${OUT_BASE}/phylum_kraken"
OUT_F="${OUT_BASE}/family_kraken"
OUT_G="${OUT_BASE}/genus_bracken"
READLEN=150
THRESH=10

mkdir -p "$OUT_K" "$OUT_P" "$OUT_F" "$OUT_G" logs
[[ -d "$DB" ]] || { echo "ERROR: Bracken/Kraken DB directory not found: ${DB}"; exit 1; }

mapfile -t REP_LIST < <(find "$KRAKEN_REPORT_DIR" -maxdepth 1 -type f -name "*.kreport.txt" | sort)
N="${#REP_LIST[@]}"
TASK_ID="${SLURM_ARRAY_TASK_ID:-0}"

[[ "$N" -gt 0 ]] || { echo "ERROR: No kreport files found in ${KRAKEN_REPORT_DIR}"; exit 1; }
[[ "$TASK_ID" -ge 0 && "$TASK_ID" -lt "$N" ]] || { echo "ERROR: Task ${TASK_ID} outside detected report range 0-$((N - 1))"; exit 1; }

REPORT="${REP_LIST[$TASK_ID]}"
SAMPLE="$(basename "$REPORT" .kreport.txt)"

echo "[$(date)] Processing taxonomic abundance for ${SAMPLE}"

awk -F'\t' '$4=="K" {gsub(/^ +/, "", $6); print $6 "\t" $1}' "$REPORT" > "${OUT_K}/${SAMPLE}.kingdom.tsv"
awk -F'\t' '$4=="P" {gsub(/^ +/, "", $6); print $6 "\t" $1}' "$REPORT" > "${OUT_P}/${SAMPLE}.phylum.tsv"
awk -F'\t' '$4=="F" {gsub(/^ +/, "", $6); print $6 "\t" $1}' "$REPORT" > "${OUT_F}/${SAMPLE}.family.tsv"

bracken \
  -d "$DB" \
  -i "$REPORT" \
  -o "${OUT_G}/${SAMPLE}.genus.bracken.tsv" \
  -r "$READLEN" \
  -l G \
  -t "$THRESH"

[[ -s "${OUT_G}/${SAMPLE}.genus.bracken.tsv" ]] || { echo "ERROR: Bracken genus output missing or empty for ${SAMPLE}"; exit 1; }

echo "[$(date)] Done ${SAMPLE}"

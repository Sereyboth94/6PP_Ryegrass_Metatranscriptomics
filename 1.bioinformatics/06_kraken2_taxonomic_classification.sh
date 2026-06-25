#!/bin/bash
#SBATCH --job-name=kraken2_nonhuman
#SBATCH --account=lincoln04415
#SBATCH --partition=milan
#SBATCH --qos=normal
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --array=0-95%10
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

# 06_kraken2_taxonomic_classification.sh
# Kraken2 taxonomic classification of plant/human-depleted paired-end reads.

set -euo pipefail

module purge
module load Kraken2/2.1.6-GCC-12.3.0

BASE="/home/soths/00_nesi_projects/uoo03549_nobackup/Sereyboth"
IN_DIR="${BASE}/Host_Removal/nonhuman"
KRAKEN_DB="/home/soths/00_nesi_projects/lincoln04415_nobackup/Metagenomics/Databases/Kraken2"
OUT_BASE="${BASE}/Host_Removal/Kraken2_nonhuman"
REP_DIR="${OUT_BASE}/reports"
KRAK_DIR="${OUT_BASE}/kraken"
LOG_DIR="${OUT_BASE}/logs"

mkdir -p "$OUT_BASE" "$REP_DIR" "$KRAK_DIR" "$LOG_DIR" logs
[[ -d "$KRAKEN_DB" ]] || { echo "ERROR: Kraken2 DB directory not found: ${KRAKEN_DB}"; exit 1; }

mapfile -t R1_LIST < <(find "$IN_DIR" -maxdepth 1 -type f -name "*.nonhuman.fastq.1.gz" | sort)
N="${#R1_LIST[@]}"
TASK_ID="${SLURM_ARRAY_TASK_ID:-0}"

[[ "$N" -gt 0 ]] || { echo "ERROR: No R1 files found in ${IN_DIR}"; exit 1; }
[[ "$TASK_ID" -ge 0 && "$TASK_ID" -lt "$N" ]] || { echo "ERROR: Task ${TASK_ID} outside detected sample range 0-$((N - 1))"; exit 1; }

R1="${R1_LIST[$TASK_ID]}"
R2="${R1/.nonhuman.fastq.1.gz/.nonhuman.fastq.2.gz}"
SAMPLE="$(basename "$R1" .nonhuman.fastq.1.gz)"

REPORT="${REP_DIR}/${SAMPLE}.kreport.txt"
KRAKEN_OUT="${KRAK_DIR}/${SAMPLE}.kraken.txt"
THREADS="${SLURM_CPUS_PER_TASK:-16}"

echo "[$(date)] Running Kraken2 for ${SAMPLE}"
[[ -f "$R1" && -f "$R2" ]] || { echo "ERROR: Missing input pair for ${SAMPLE}"; exit 1; }

if [[ -s "$REPORT" && -s "$KRAKEN_OUT" ]]; then
  echo "[$(date)] Kraken2 outputs already exist for ${SAMPLE}; skipping."
  exit 0
fi

kraken2 \
  --db "$KRAKEN_DB" \
  --paired \
  --gzip-compressed \
  --threads "$THREADS" \
  --report "$REPORT" \
  --output "$KRAKEN_OUT" \
  "$R1" "$R2"

[[ -s "$REPORT" ]] || { echo "ERROR: Kraken2 report missing or empty for ${SAMPLE}"; exit 1; }

echo "[$(date)] Finished ${SAMPLE}"

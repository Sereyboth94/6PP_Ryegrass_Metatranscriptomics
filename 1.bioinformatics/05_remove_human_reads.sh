#!/bin/bash
#SBATCH --job-name=host_rm_GRCh38
#SBATCH --account=lincoln04266
#SBATCH --partition=milan
#SBATCH --qos=normal
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --array=0-95%12
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

# 05_remove_human_reads.sh
# Additional depletion of possible human reads from plant-host-depleted FASTQs.

set -euo pipefail

module purge
module load Bowtie2/2.5.4-GCC-12.3.0

BASE="/home/soths/00_nesi_projects/uoo03549_nobackup/Sereyboth"
FASTQ_DIR="${BASE}/Host_Removal/unmapped"
OUT_DIR="${BASE}/Host_Removal/nonhuman"
HUMAN_INDEX="/home/soths/00_nesi_projects/lincoln04415_nobackup/Metagenomics/Databases/Human_GRCh38/bowtie2_index/grch38_primary"

mkdir -p "$OUT_DIR" logs

mapfile -t R1_LIST < <(find "$FASTQ_DIR" -maxdepth 1 -type f -name "*.unmapped.fastq.1.gz" | sort)
N="${#R1_LIST[@]}"
TASK_ID="${SLURM_ARRAY_TASK_ID:-0}"

[[ "$N" -gt 0 ]] || { echo "ERROR: No plant-host-depleted R1 files found in ${FASTQ_DIR}"; exit 1; }
[[ "$TASK_ID" -ge 0 && "$TASK_ID" -lt "$N" ]] || { echo "ERROR: Task ${TASK_ID} outside detected sample range 0-$((N - 1))"; exit 1; }
[[ -f "${HUMAN_INDEX}.1.bt2" || -f "${HUMAN_INDEX}.1.bt2l" ]] || { echo "ERROR: Bowtie2 index not found for prefix ${HUMAN_INDEX}"; exit 1; }

R1="${R1_LIST[$TASK_ID]}"
R2="${R1/.unmapped.fastq.1.gz/.unmapped.fastq.2.gz}"
SAMPLE="$(basename "$R1" .unmapped.fastq.1.gz)"

OUT_R1="${OUT_DIR}/${SAMPLE}.nonhuman.fastq.1.gz"
OUT_R2="${OUT_DIR}/${SAMPLE}.nonhuman.fastq.2.gz"
LOG="${OUT_DIR}/${SAMPLE}.human_bowtie2.log"
THREADS="${SLURM_CPUS_PER_TASK:-16}"

echo "[$(date)] Removing possible human reads for ${SAMPLE}"
[[ -f "$R1" && -f "$R2" ]] || { echo "ERROR: Missing input pair for ${SAMPLE}"; exit 1; }

if [[ -s "$OUT_R1" && -s "$OUT_R2" ]]; then
  echo "[$(date)] Nonhuman FASTQs already exist for ${SAMPLE}; skipping."
  exit 0
fi

bowtie2 \
  -x "$HUMAN_INDEX" \
  -1 "$R1" \
  -2 "$R2" \
  --very-sensitive \
  -p "$THREADS" \
  --un-conc-gz "${OUT_DIR}/${SAMPLE}.nonhuman.fastq.gz" \
  -S /dev/null \
  2> "$LOG"

mv "${OUT_DIR}/${SAMPLE}.nonhuman.fastq.1.gz" "$OUT_R1"
mv "${OUT_DIR}/${SAMPLE}.nonhuman.fastq.2.gz" "$OUT_R2"

[[ -s "$OUT_R1" && -s "$OUT_R2" ]] || { echo "ERROR: nonhuman outputs missing or empty for ${SAMPLE}"; exit 1; }

echo "[$(date)] Finished ${SAMPLE}"

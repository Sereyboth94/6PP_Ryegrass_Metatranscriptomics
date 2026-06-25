#!/bin/bash
#SBATCH --job-name=fastp_trim
#SBATCH --account=uoo03549
#SBATCH --partition=milan
#SBATCH --array=0-95%8
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

# 01_fastp_quality_trim.sh
# Quality trimming of paired-end metatranscriptomic FASTQ files with fastp.
# Processes one biological sample per SLURM array task across all listed flowcells and lanes.

set -euo pipefail

# -------------------------
# User-editable settings
# -------------------------
BASE="/home/soths/00_nesi_projects/uoo03549_nobackup/Sereyboth"
FLOWCELLS=("V350259628" "V350291364")
LANES=("L01" "L02" "L03" "L04")

# Input filename pattern expected:
#   ${BASE}/${FLOWCELL}/${LANE}/${FLOWCELL}_${LANE}_${SAMPLE}_1.fq.gz
#   ${BASE}/${FLOWCELL}/${LANE}/${FLOWCELL}_${LANE}_${SAMPLE}_2.fq.gz
# Output directory per lane:
#   ${BASE}/${FLOWCELL}/${LANE}/Quality_trimmed/

# Optional module load; edit for your HPC environment if needed.
# module purge
module load fastp/1.1.0-GCC-12.3.0

mkdir -p logs

# -------------------------
# Sample list A01..H12
# -------------------------
mapfile -t SAMPLES < <(
  for r in A B C D E F G H; do
    for c in $(seq -w 1 12); do
      echo "${r}${c}"
    done
  done
)

IDX="${SLURM_ARRAY_TASK_ID:-0}"
[[ "$IDX" -ge 0 && "$IDX" -lt "${#SAMPLES[@]}" ]] || {
  echo "ERROR: SLURM_ARRAY_TASK_ID=$IDX is outside sample range 0-$((${#SAMPLES[@]} - 1))"
  exit 1
}

SAMPLE="${SAMPLES[$IDX]}"
THREADS="${SLURM_CPUS_PER_TASK:-16}"

echo "[$(date)] Starting fastp for sample ${SAMPLE}"

processed=0
for FC in "${FLOWCELLS[@]}"; do
  for LANE in "${LANES[@]}"; do
    INDIR="${BASE}/${FC}/${LANE}"
    OUTDIR="${INDIR}/Quality_trimmed"

    R1="${INDIR}/${FC}_${LANE}_${SAMPLE}_1.fq.gz"
    R2="${INDIR}/${FC}_${LANE}_${SAMPLE}_2.fq.gz"

    OUT_R1="${OUTDIR}/${FC}_${LANE}_${SAMPLE}_quality_trimmed_1.fq.gz"
    OUT_R2="${OUTDIR}/${FC}_${LANE}_${SAMPLE}_quality_trimmed_2.fq.gz"
    JSON="${OUTDIR}/${FC}_${LANE}_${SAMPLE}.fastp.json"
    HTML="${OUTDIR}/${FC}_${LANE}_${SAMPLE}.fastp.html"

    if [[ ! -f "$R1" || ! -f "$R2" ]]; then
      echo "[$(date)] Missing pair for ${FC} ${LANE} ${SAMPLE}; skipping."
      continue
    fi

    mkdir -p "$OUTDIR"

    if [[ -s "$OUT_R1" && -s "$OUT_R2" ]]; then
      echo "[$(date)] Existing fastp outputs found for ${FC} ${LANE} ${SAMPLE}; skipping."
      processed=$((processed + 1))
      continue
    fi

    echo "[$(date)] Running fastp: ${FC} ${LANE} ${SAMPLE}"
    fastp \
      -i "$R1" \
      -I "$R2" \
      -o "$OUT_R1" \
      -O "$OUT_R2" \
      --length_required 30 \
      --qualified_quality_phred 20 \
      --thread "$THREADS" \
      --json "$JSON" \
      --html "$HTML"

    [[ -s "$OUT_R1" && -s "$OUT_R2" ]] || {
      echo "ERROR: fastp output missing or empty for ${FC} ${LANE} ${SAMPLE}"
      exit 1
    }

    processed=$((processed + 1))
  done
done

if [[ "$processed" -eq 0 ]]; then
  echo "ERROR: No input pairs were processed for sample ${SAMPLE}. Check BASE/FLOWCELLS/LANES and filename pattern."
  exit 1
fi

echo "[$(date)] Finished fastp for sample ${SAMPLE}; processed ${processed} file pair(s)."

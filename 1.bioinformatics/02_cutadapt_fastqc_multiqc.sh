#!/bin/bash
#SBATCH --job-name=cutadapt_qc
#SBATCH --account=uoo03549
#SBATCH --partition=milan
#SBATCH --array=0-95%8
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

# 02_cutadapt_fastqc_multiqc.sh
# Adapter/position trimming after fastp, followed by FastQC.
# MultiQC mode can be run once after all array tasks finish:
#   sbatch --export=MODE=multiqc 02_cutadapt_fastqc_multiqc.sh

set -euo pipefail

# -------------------------
# User-editable settings
# -------------------------
BASE="/home/soths/00_nesi_projects/uoo03549_nobackup/Sereyboth"
FLOWCELLS=("V350259628" "V350291364")
LANES=("L01" "L02" "L03" "L04")

FWD_ADAPTER="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
REV_ADAPTER="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
TRIM_5P_R1=6
TRIM_5P_R2=8
MODE="${MODE:-trim}"   # trim or multiqc

# Optional module loads; edit for your HPC environment if needed.
# module purge
module load cutadapt/5.2-foss-2023a-Python-3.11.6
module load FastQC/0.12.1
module load MultiQC/1.24.1-foss-2023a-Python-3.11.6

mkdir -p logs
THREADS="${SLURM_CPUS_PER_TASK:-16}"
CUTADAPT_CORES="$(( THREADS > 4 ? 4 : THREADS ))"

run_multiqc() {
  echo "[$(date)] Running MultiQC over all Final_trimmed QC directories."
  QC_DIRS=()
  for FC in "${FLOWCELLS[@]}"; do
    for LANE in "${LANES[@]}"; do
      d="${BASE}/${FC}/${LANE}/Quality_trimmed/Final_trimmed/QC_after_trimming"
      [[ -d "$d" ]] && QC_DIRS+=("$d")
    done
  done

  [[ "${#QC_DIRS[@]}" -gt 0 ]] || { echo "ERROR: No QC_after_trimming directories found."; exit 1; }
  OUT="${BASE}/MultiQC_after_trimming"
  mkdir -p "$OUT"
  multiqc "${QC_DIRS[@]}" -o "$OUT"
  echo "[$(date)] MultiQC report written to ${OUT}"
}

if [[ "$MODE" == "multiqc" ]]; then
  run_multiqc
  exit 0
fi

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
echo "[$(date)] Starting cutadapt/FastQC for sample ${SAMPLE}"

processed=0
for FC in "${FLOWCELLS[@]}"; do
  for LANE in "${LANES[@]}"; do
    QDIR="${BASE}/${FC}/${LANE}/Quality_trimmed"
    FINAL_DIR="${QDIR}/Final_trimmed"
    QC_DIR="${FINAL_DIR}/QC_after_trimming"

    R1="${QDIR}/${FC}_${LANE}_${SAMPLE}_quality_trimmed_1.fq.gz"
    R2="${QDIR}/${FC}_${LANE}_${SAMPLE}_quality_trimmed_2.fq.gz"

    # Important: these names match 03_merge_two_flowcells.sh
    OUT_R1="${FINAL_DIR}/${FC}_${LANE}_${SAMPLE}_trimmed_1.fq.gz"
    OUT_R2="${FINAL_DIR}/${FC}_${LANE}_${SAMPLE}_trimmed_2.fq.gz"
    LOG="${FINAL_DIR}/${FC}_${LANE}_${SAMPLE}.cutadapt.log"

    if [[ ! -f "$R1" || ! -f "$R2" ]]; then
      echo "[$(date)] Missing fastp pair for ${FC} ${LANE} ${SAMPLE}; skipping."
      continue
    fi

    mkdir -p "$FINAL_DIR" "$QC_DIR"

    if [[ -s "$OUT_R1" && -s "$OUT_R2" ]]; then
      echo "[$(date)] Existing cutadapt outputs found for ${FC} ${LANE} ${SAMPLE}; skipping trimming."
    else
      echo "[$(date)] Running cutadapt: ${FC} ${LANE} ${SAMPLE}"
      cutadapt \
        -u "$TRIM_5P_R1" \
        -U "$TRIM_5P_R2" \
        --cores "$CUTADAPT_CORES" \
        -a "$FWD_ADAPTER" \
        -A "$REV_ADAPTER" \
        -o "$OUT_R1" \
        -p "$OUT_R2" \
        "$R1" "$R2" \
        > "$LOG"
    fi

    [[ -s "$OUT_R1" && -s "$OUT_R2" ]] || {
      echo "ERROR: cutadapt output missing or empty for ${FC} ${LANE} ${SAMPLE}"
      exit 1
    }

    echo "[$(date)] Running FastQC: ${FC} ${LANE} ${SAMPLE}"
    fastqc -t "$THREADS" -o "$QC_DIR" "$OUT_R1" "$OUT_R2"
    processed=$((processed + 1))
  done
done

if [[ "$processed" -eq 0 ]]; then
  echo "ERROR: No fastp outputs were processed for sample ${SAMPLE}."
  exit 1
fi

echo "[$(date)] Finished cutadapt/FastQC for sample ${SAMPLE}; processed ${processed} file pair(s)."

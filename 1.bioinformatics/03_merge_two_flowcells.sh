#!/bin/bash
#SBATCH --job-name=merge_fc
#SBATCH --account=uoo03549
#SBATCH --partition=milan
#SBATCH --array=0-95%8
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:30:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

# 03_merge_two_flowcells.sh
# Concatenates final-trimmed paired-end FASTQ files from multiple flowcells/lanes into one R1/R2 pair per sample.

set -euo pipefail

BASE="/home/soths/00_nesi_projects/uoo03549_nobackup/Sereyboth"
FLOWCELLS=("V350259628" "V350291364")
LANES=("L01" "L02" "L03" "L04")
OUT="${BASE}/merged_flowcells"

mkdir -p "$OUT" logs

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
echo "[$(date)] Merging sample ${SAMPLE}"

R1_FILES=()
R2_FILES=()

for FC in "${FLOWCELLS[@]}"; do
  for LANE in "${LANES[@]}"; do
    f1="${BASE}/${FC}/${LANE}/Quality_trimmed/Final_trimmed/${FC}_${LANE}_${SAMPLE}_trimmed_1.fq.gz"
    f2="${BASE}/${FC}/${LANE}/Quality_trimmed/Final_trimmed/${FC}_${LANE}_${SAMPLE}_trimmed_2.fq.gz"

    if [[ -f "$f1" && -f "$f2" ]]; then
      R1_FILES+=("$f1")
      R2_FILES+=("$f2")
      echo "  found ${FC} ${LANE}"
    elif [[ -f "$f1" || -f "$f2" ]]; then
      echo "ERROR: Incomplete pair for ${FC} ${LANE} ${SAMPLE}"
      echo "R1=$f1"
      echo "R2=$f2"
      exit 1
    fi
  done
done

if [[ "${#R1_FILES[@]}" -eq 0 ]]; then
  echo "ERROR: No trimmed input files found for ${SAMPLE}."
  exit 1
fi

if [[ "${#R1_FILES[@]}" -ne "${#R2_FILES[@]}" ]]; then
  echo "ERROR: R1/R2 file count mismatch for ${SAMPLE}."
  exit 1
fi

TMP_R1="${OUT}/${SAMPLE}_R1.fastq.gz.tmp"
TMP_R2="${OUT}/${SAMPLE}_R2.fastq.gz.tmp"

cat "${R1_FILES[@]}" > "$TMP_R1"
cat "${R2_FILES[@]}" > "$TMP_R2"

mv -f "$TMP_R1" "${OUT}/${SAMPLE}_R1.fastq.gz"
mv -f "$TMP_R2" "${OUT}/${SAMPLE}_R2.fastq.gz"

[[ -s "${OUT}/${SAMPLE}_R1.fastq.gz" && -s "${OUT}/${SAMPLE}_R2.fastq.gz" ]] || {
  echo "ERROR: merged outputs are missing or empty for ${SAMPLE}."
  exit 1
}

echo "[$(date)] Done ${SAMPLE}: merged ${#R1_FILES[@]} input pair(s)."

#!/bin/bash
#SBATCH --job-name=hostrm_lp
#SBATCH --account=uoo03549
#SBATCH --partition=milan
#SBATCH --array=0-95%12
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=06:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

# 04_remove_plant_host.sh
# Removes perennial ryegrass host reads by Bowtie2 mapping to the Lolium perenne reference.

set -euo pipefail

module purge
module load Bowtie2/2.5.4-GCC-12.3.0
module load SAMtools/1.23.1-GCC-12.3.0

BASE="/home/soths/00_nesi_projects/uoo03549_nobackup/Sereyboth"
IN="${BASE}/merged_flowcells"
OUT="${BASE}/Host_Removal"
IDX="/home/soths/00_nesi_projects/lincoln04415_nobackup/Metagenomics/Fastq/Host_Reference/Lolium_perenne/Lolium_perenne_Kyuss_2.0"

mkdir -p "${OUT}/unmapped" "${OUT}/host_bam" "${OUT}/logs" logs

mapfile -t SAMPLES < <(
  for r in A B C D E F G H; do
    for c in $(seq -w 1 12); do
      echo "${r}${c}"
    done
  done
)

TASK_ID="${SLURM_ARRAY_TASK_ID:-0}"
[[ "$TASK_ID" -ge 0 && "$TASK_ID" -lt "${#SAMPLES[@]}" ]] || {
  echo "ERROR: SLURM_ARRAY_TASK_ID=$TASK_ID is outside sample range."
  exit 1
}

SAMPLE="${SAMPLES[$TASK_ID]}"
R1="${IN}/${SAMPLE}_R1.fastq.gz"
R2="${IN}/${SAMPLE}_R2.fastq.gz"

HOST_BAM="${OUT}/host_bam/${SAMPLE}.hostmapped.bam"
UNMAPPED_PREFIX="${OUT}/unmapped/${SAMPLE}.unmapped.fastq"
LOG="${OUT}/logs/${SAMPLE}.plant_host_bowtie2.log"
THREADS="${SLURM_CPUS_PER_TASK:-8}"

echo "[$(date)] Removing plant host reads for ${SAMPLE}"
[[ -f "$R1" && -f "$R2" ]] || { echo "ERROR: Missing merged input FASTQs for ${SAMPLE}"; exit 1; }
[[ -f "${IDX}.1.bt2" || -f "${IDX}.1.bt2l" ]] || { echo "ERROR: Bowtie2 index not found for prefix ${IDX}"; exit 1; }

if [[ -s "${UNMAPPED_PREFIX}.1.gz" && -s "${UNMAPPED_PREFIX}.2.gz" ]]; then
  echo "[$(date)] Plant-host-depleted FASTQs already exist for ${SAMPLE}; skipping."
  exit 0
fi

bowtie2 --very-sensitive \
  -x "$IDX" \
  -1 "$R1" \
  -2 "$R2" \
  -p "$THREADS" \
  --un-conc-gz "${UNMAPPED_PREFIX}.gz" \
  2> "$LOG" \
  | samtools view -@ "$THREADS" -bS - \
  | samtools sort -@ "$THREADS" -o "$HOST_BAM" -

samtools index "$HOST_BAM"

[[ -s "${UNMAPPED_PREFIX}.1.gz" && -s "${UNMAPPED_PREFIX}.2.gz" ]] || {
  echo "ERROR: Plant-host-depleted FASTQs are missing or empty for ${SAMPLE}."
  exit 1
}

echo "[$(date)] Done ${SAMPLE}"
echo "Unmapped: ${UNMAPPED_PREFIX}.1.gz ${UNMAPPED_PREFIX}.2.gz"
echo "Host BAM: ${HOST_BAM}"

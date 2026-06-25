#!/bin/bash
#SBATCH --job-name=humann4_nonhuman
#SBATCH --account=lincoln04266
#SBATCH --partition=milan
#SBATCH --qos=normal
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=180G
#SBATCH --array=0-95%4
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

# 08_humann4_functional_profiling.sh
# HUMAnN4 functional profiling of plant/human-depleted metatranscriptomic reads.
# This script keeps the project-specific MetaPhlAn compatibility wrapper used in the analysis.

set -euo pipefail

module purge
module load Miniconda3/23.10.0-1

# -------------------------
# User-editable settings
# -------------------------
BASE="/home/soths/00_nesi_projects/uoo03549_nobackup/Sereyboth"
ENV="/home/soths/00_nesi_projects/lincoln04415_nobackup/Metatrans_Test/envs/biobakery4_py312_clean"
INDIR="${BASE}/Host_Removal/nonhuman_merged_cat"
OUTBASE="${BASE}/Humann4_Output_nonhuman"

DBROOT="/home/soths/00_nesi_projects/lincoln04415_nobackup/Metatrans_Test/db"
HUMANN_NUC_DB="${DBROOT}/chocophlan"
HUMANN_PROT_DB="${DBROOT}/uniref"
HUMANN_UTIL_DB="${DBROOT}/utility_mapping"

METAPHLAN_DB="${DBROOT}/mpa_vJan25_CHOCOPhlAnSGB_202503"
METAPHLAN_INDEX="mpa_vJan25_CHOCOPhlAnSGB_202503"
HUMANN_EXPECTED_DB="vOct22_CHOCOPhlAnSGB_202403"

mkdir -p "$OUTBASE" logs

# -------------------------
# Activate environment
# -------------------------
source activate "$ENV"
export PATH="${ENV}/bin:${PATH}"
hash -r

# -------------------------
# MetaPhlAn compatibility wrapper
# -------------------------
WRAPDIR="${OUTBASE}/_wrappers"
mkdir -p "$WRAPDIR"

cat > "${WRAPDIR}/metaphlan" <<EOF_WRAP
#!/usr/bin/env bash
set -euo pipefail
REAL="${ENV}/bin/metaphlan"

if [[ "\${1:-}" == "--version" ]]; then
  "\$REAL" --version 2>/dev/null | head -n 1
  exit 0
fi

args=()
i=1
while [[ \$i -le \$# ]]; do
  a="\${!i}"
  case "\$a" in
    --bowtie2db)  i=\$((i+1)); args+=(--db_dir "\${!i}") ;;
    --bowtie2out) i=\$((i+1)); args+=(--mapout "\${!i}") ;;
    --add_viruses|--unclassified_estimation)
      ;;
    *) args+=("\$a") ;;
  esac
  i=\$((i+1))
done

exec "\$REAL" "\${args[@]}"
EOF_WRAP

chmod +x "${WRAPDIR}/metaphlan"
export PATH="${WRAPDIR}:${PATH}"
hash -r

# -------------------------
# Sanity checks
# -------------------------
echo "humann:        $(which humann)"
humann --version || true
echo "humann_config: $(which humann_config)"
echo "metaphlan:     $(which metaphlan)"
metaphlan --version || true

[[ -d "$HUMANN_NUC_DB" ]] || { echo "ERROR: missing ${HUMANN_NUC_DB}"; exit 1; }
[[ -d "$HUMANN_PROT_DB" ]] || { echo "ERROR: missing ${HUMANN_PROT_DB}"; exit 1; }
[[ -d "$HUMANN_UTIL_DB" ]] || { echo "ERROR: missing ${HUMANN_UTIL_DB}"; exit 1; }
[[ -f "$HUMANN_UTIL_DB/metacyc_reactions_level4ec_only.uniref.bz2" ]] || { echo "ERROR: missing MetaCyc reactions file in utility_mapping"; exit 1; }
[[ -e "$HUMANN_UTIL_DB/metacyc_pathways_structured_filtered_v24_subreactions" ]] || { echo "ERROR: missing MetaCyc pathways structure file in utility_mapping"; exit 1; }
[[ -d "$METAPHLAN_DB" ]] || { echo "ERROR: missing ${METAPHLAN_DB}"; exit 1; }
[[ -f "$METAPHLAN_DB/${METAPHLAN_INDEX}.1.bt2l" || -f "$METAPHLAN_DB/${METAPHLAN_INDEX}.1.bt2" ]] || { echo "ERROR: missing MetaPhlAn Bowtie2 index for ${METAPHLAN_INDEX}"; exit 1; }

# Provide a valid HUMAnN config to avoid accidental use of a stale/broken config file.
CFGDIR="${OUTBASE}/_humann_cfg"
mkdir -p "$CFGDIR"
CFGFILE="${CFGDIR}/humann.cfg"
cat > "$CFGFILE" <<EOF_CFG
[database_folders]
nucleotide = ${HUMANN_NUC_DB}
protein = ${HUMANN_PROT_DB}
utility_mapping = ${HUMANN_UTIL_DB}
EOF_CFG
export HUMANN_CONFIG_FILE="$CFGFILE"
echo "Using HUMANN_CONFIG_FILE=${HUMANN_CONFIG_FILE}"

# -------------------------
# Sample list from merged single-end/nonhuman cat files
# Expected input: ${INDIR}/${SAMPLE}_merged.fq.gz
# -------------------------
[[ -d "$INDIR" ]] || { echo "ERROR: input directory not found: ${INDIR}"; exit 1; }
mapfile -t SAMPLES < <(find "$INDIR" -maxdepth 1 -type f -name "*_merged.fq.gz" -printf "%f\n" | sed 's/_merged\.fq\.gz$//' | sort)
N="${#SAMPLES[@]}"
TASK_ID="${SLURM_ARRAY_TASK_ID:-0}"

[[ "$N" -gt 0 ]] || { echo "ERROR: No *_merged.fq.gz files found in ${INDIR}"; exit 1; }
[[ "$TASK_ID" -ge 0 && "$TASK_ID" -lt "$N" ]] || { echo "ERROR: Task ${TASK_ID} outside detected sample range 0-$((N - 1))"; exit 1; }

ID="${SAMPLES[$TASK_ID]}"
IN="${INDIR}/${ID}_merged.fq.gz"
OUTDIR="${OUTBASE}/${ID}"
THREADS="${SLURM_CPUS_PER_TASK:-32}"
MP_THREADS="$(( THREADS > 16 ? 16 : THREADS ))"

mkdir -p "$OUTDIR"

echo "[$(date)] Sample:  ${ID}"
echo "[$(date)] Input:   ${IN}"
echo "[$(date)] Outdir:  ${OUTDIR}"
echo "[$(date)] Threads: ${THREADS} (MetaPhlAn: ${MP_THREADS})"
[[ -s "$IN" ]] || { echo "ERROR: Missing or empty input ${IN}"; exit 1; }

FINAL_GF="${OUTDIR}/${ID}_merged_2_genefamilies.tsv"
FINAL_RXN="${OUTDIR}/${ID}_merged_3_reactions.tsv"
FINAL_PWY="${OUTDIR}/${ID}_merged_4_pathabundance.tsv"

if [[ -s "$FINAL_GF" && -s "$FINAL_RXN" && -s "$FINAL_PWY" ]]; then
  if [[ "$(wc -l < "$FINAL_PWY")" -gt 5 ]]; then
    echo "[$(date)] HUMAnN outputs already exist and look OK. Skipping ${ID}."
    exit 0
  fi
  echo "[$(date)] Existing pathabundance file appears too small; re-running ${ID}."
fi

# -------------------------
# MetaPhlAn profile
# -------------------------
MP_RAW="${OUTDIR}/${ID}_metaphlan_rel_ab_w_read_stats.tsv"
MP_OK="${OUTDIR}/${ID}_metaphlan_for_humann.tsv"
MP_MAPOUT="${OUTDIR}/${ID}_metaphlan_mapout.bz2"

echo "[$(date)] Running MetaPhlAn..."
metaphlan "$IN" \
  --input_type fastq \
  --bowtie2db "$METAPHLAN_DB" \
  -x "$METAPHLAN_INDEX" \
  --nproc "$MP_THREADS" \
  --offline \
  -t rel_ab_w_read_stats \
  --bowtie2out "$MP_MAPOUT" \
  -o "$MP_RAW"

[[ -s "$MP_RAW" ]] || { echo "ERROR: MetaPhlAn output empty: ${MP_RAW}"; exit 1; }

# Project-specific compatibility patch: HUMAnN alpha expected a specific DB label.
sed "1s/${METAPHLAN_INDEX}/${HUMANN_EXPECTED_DB}/" "$MP_RAW" > "$MP_OK"
[[ -s "$MP_OK" ]] || { echo "ERROR: patched MetaPhlAn output empty: ${MP_OK}"; exit 1; }

# -------------------------
# HUMAnN
# -------------------------
echo "[$(date)] Running HUMAnN..."
humann \
  -i "$IN" \
  -o "$OUTDIR" \
  --threads "$THREADS" \
  --nucleotide-database "$HUMANN_NUC_DB" \
  --protein-database "$HUMANN_PROT_DB" \
  --utility-database "$HUMANN_UTIL_DB" \
  --taxonomic-profile "$MP_OK"

PWY_LINES="$(wc -l < "$FINAL_PWY" || echo 0)"
echo "[$(date)] pathabundance lines: ${PWY_LINES}"
if [[ "$PWY_LINES" -le 5 ]]; then
  echo "ERROR: pathabundance looks wrong (<=5 lines). Check ${OUTDIR}/${ID}_merged_0.log"
  exit 2
fi

# -------------------------
# Cleanup intermediates; keep key HUMAnN outputs and logs.
# -------------------------
echo "[$(date)] Cleaning intermediates..."
rm -f "$MP_RAW" "$MP_OK" "$MP_MAPOUT" 2>/dev/null || true
rm -f \
  "$OUTDIR"/*_bowtie2out*.bz2 \
  "$OUTDIR"/*_bowtie2out*.txt \
  "$OUTDIR"/*diamond*.daa \
  "$OUTDIR"/*tmp* \
  "$OUTDIR"/*temp* 2>/dev/null || true

echo "[$(date)] Done ${ID}"

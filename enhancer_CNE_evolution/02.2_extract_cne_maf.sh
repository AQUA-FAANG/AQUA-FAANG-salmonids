#!/bin/bash
# =============================================================================
# extract_cne_maf.sh
# SLURM array job: extract MAF alignment blocks for CNE regions using pike
# genome coordinates as the reference species for mafExtractor.
#
# Each array task processes one line from the relevant TSV list.
# The list type (exclusive | shared | alignable) is passed via the environment
# variable LIST_TYPE, set by submit_pipeline.sh at submission time.
#
# Input TSV columns (tab-separated, no header):
#   1  salmon_chr   e.g. ssa10
#   2  pike_chr     e.g. LG19
#   3  pike_start   e.g. 10716213
#   4  pike_end     e.g. 10716295
#   5  maf_filename e.g. ssa10_113462995_113463677.maf
#
# Output filename format (inside the per-list output folder):
#   {maf_filename}_{pike_chr}_{pike_start}_{pike_end}.maf
#   e.g.  exclusive/ssa10_113462995_113463677.maf_LG19_10716213_10716295.maf
# =============================================================================
#SBATCH --output=logs/extract_%A_%a.out
#SBATCH --error=logs/extract_%A_%a.err
#SBATCH --time=04:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
# NOTE: --array and --job-name are set dynamically by submit_pipeline.sh

set -euo pipefail

# ── Configuration ────────────────────────────────────────────────────────────
# Path to the combined MAF file produced by download_cat_maf.sh
COMBINED_MAF="all_alignments.maf"

# Pike species name prefix as it appears in the MAF species strings.
# In salmon-trout-pike MAF files the pike (Esox lucius) sequence lines are
# typically labelled  "Eluc.LGxx" — verify against your MAF headers and
# adjust if needed (e.g. "EsoxLuc", "pike", etc.)
PIKE_SPECIES_PREFIX="Eluc"

# Path to mafExtractor binary
MAFEXTRACTOR="softwares/mafTools/bin/mafExtractor"
# ────────────────────────────────────────────────────────────────────────────

# ── Validate environment ─────────────────────────────────────────────────────
echo "[$(date '+%F %T')] === Array task ${SLURM_ARRAY_TASK_ID} started ==="
echo "LIST_TYPE : ${LIST_TYPE:?'ERROR: LIST_TYPE env variable not set. Use submit_pipeline.sh.'}"
echo "Node      : $(hostname)"
echo ""

# ── Resolve input list and output directory ──────────────────────────────────
case "${LIST_TYPE}" in
    exclusive) INPUT_LIST="exclusive_CNE_extraction_list.tsv" ; OUT_DIR="exclusive"  ;;
    shared)    INPUT_LIST="shared_CNE_extraction_list.tsv"    ; OUT_DIR="shared"     ;;
    alignable) INPUT_LIST="alignable_CNE_extraction_list.tsv" ; OUT_DIR="alignable"  ;;
    *)
        echo "ERROR: Unknown LIST_TYPE '${LIST_TYPE}'. Must be exclusive | shared | alignable."
        exit 1
        ;;
esac

# ── Validate required files ──────────────────────────────────────────────────
[[ -f "${INPUT_LIST}" ]]   || { echo "ERROR: input list not found: ${INPUT_LIST}";   exit 1; }
[[ -f "${COMBINED_MAF}" ]] || { echo "ERROR: combined MAF not found: ${COMBINED_MAF}"; exit 1; }
[[ -x "${MAFEXTRACTOR}" ]] || { echo "ERROR: mafExtractor not found/executable: ${MAFEXTRACTOR}"; exit 1; }

mkdir -p "${OUT_DIR}"

# ── Retrieve this task's input line ─────────────────────────────────────────
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${INPUT_LIST}")

if [[ -z "${LINE}" ]]; then
    echo "ERROR: No line ${SLURM_ARRAY_TASK_ID} in ${INPUT_LIST}."
    exit 1
fi

# Parse the five TSV columns
read -r sal_chr pike_chr pike_start pike_end maf_file <<< "${LINE}"

# ── Build mafExtractor species string and output path ────────────────────────
# mafExtractor expects the species string exactly as written in the MAF 's' lines,
# e.g.  s Eluc.LG19  10716213  82  +  ...
SPECIES_STRING="${PIKE_SPECIES_PREFIX}.${pike_chr}"

# Output filename: {original_maf_file}_{pike_chr}_{pike_start}_{pike_end}.maf
OUTPUT_FILE="${OUT_DIR}/${maf_file}_${pike_chr}_${pike_start}_${pike_end}.maf"

echo "Input line    : ${LINE}"
echo "Salmon chr    : ${sal_chr}"
echo "Pike region   : ${pike_chr}:${pike_start}-${pike_end}"
echo "Species string: ${SPECIES_STRING}"
echo "Source MAF    : ${COMBINED_MAF}"
echo "Output file   : ${OUTPUT_FILE}"
echo ""

# ── Run mafExtractor ─────────────────────────────────────────────────────────
"${MAFEXTRACTOR}" \
    --maf "${COMBINED_MAF}" \
    --seq "${SPECIES_STRING}" \
    --start "${pike_start}" \
    --stop  "${pike_end}" \
    > "${OUTPUT_FILE}"

# Sanity check: warn if output is empty (no overlapping alignment block found)
if [[ ! -s "${OUTPUT_FILE}" ]]; then
    echo "WARNING: Output file is empty — no alignment block found for"
    echo "         ${SPECIES_STRING} ${pike_start}-${pike_end}"
    echo "         Check PIKE_SPECIES_PREFIX and coordinate system."
fi

echo ""
echo "[$(date '+%F %T')] === Task ${SLURM_ARRAY_TASK_ID} finished: ${OUTPUT_FILE} ==="

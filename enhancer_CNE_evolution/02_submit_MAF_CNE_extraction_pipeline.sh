#!/bin/bash
# =============================================================================
# submit_pipeline.sh
# Master submission script for CNE MAF extraction pipeline.
#
# Usage:
#   bash submit_pipeline.sh
#
# Steps:
#   1. Submits download + concatenation job (download_cat_maf.sh)
#   2. Submits three array extraction jobs (extract_cne_maf.sh),
#      one per CNE list type, dependent on step 1 completing successfully.
#
# Outputs (created by extraction jobs):
#   exclusive/   shared/   alignable/
# =============================================================================

set -euo pipefail

# ── User-configurable paths ──────────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DOWNLOAD_SCRIPT="${SCRIPT_DIR}/download_cat_maf.sh"
EXTRACT_SCRIPT="${SCRIPT_DIR}/extract_cne_maf.sh"

EXCL_LIST="exclusive_CNE_extraction_list.tsv"
SHAR_LIST="shared_CNE_extraction_list.tsv"
ALIG_LIST="alignable_CNE_extraction_list.tsv"
# ────────────────────────────────────────────────────────────────────────────

mkdir -p logs

# Sanity-check that all three input lists are present before submitting
for f in "${EXCL_LIST}" "${SHAR_LIST}" "${ALIG_LIST}"; do
    [[ -f "${f}" ]] || { echo "ERROR: input list not found: ${f}"; exit 1; }
done

# Count lines (skip blank lines just in case)
N_EXCL=$(grep -c '[^[:space:]]' "${EXCL_LIST}" || true)
N_SHAR=$(grep -c '[^[:space:]]' "${SHAR_LIST}" || true)
N_ALIG=$(grep -c '[^[:space:]]' "${ALIG_LIST}" || true)

echo "Input sizes:"
echo "  exclusive  : ${N_EXCL} entries"
echo "  shared     : ${N_SHAR} entries"
echo "  alignable  : ${N_ALIG} entries"

# ── Step 1: download + concatenation ────────────────────────────────────────
DL_JOB=$(sbatch --parsable "${DOWNLOAD_SCRIPT}")
echo "Submitted download job: ${DL_JOB}"

# ── Step 2: three array extraction jobs (run only after step 1 succeeds) ────
EXCL_JOB=$(sbatch --parsable \
    --dependency=afterok:"${DL_JOB}" \
    --array="1-${N_EXCL}" \
    --export=ALL,LIST_TYPE=exclusive \
    --job-name=extract_exclusive \
    "${EXTRACT_SCRIPT}")
echo "Submitted exclusive extraction job: ${EXCL_JOB}  (array 1-${N_EXCL})"

SHAR_JOB=$(sbatch --parsable \
    --dependency=afterok:"${DL_JOB}" \
    --array="1-${N_SHAR}" \
    --export=ALL,LIST_TYPE=shared \
    --job-name=extract_shared \
    "${EXTRACT_SCRIPT}")
echo "Submitted shared extraction job: ${SHAR_JOB}  (array 1-${N_SHAR})"

ALIG_JOB=$(sbatch --parsable \
    --dependency=afterok:"${DL_JOB}" \
    --array="1-${N_ALIG}" \
    --export=ALL,LIST_TYPE=alignable \
    --job-name=extract_alignable \
    "${EXTRACT_SCRIPT}")
echo "Submitted alignable extraction job: ${ALIG_JOB}  (array 1-${N_ALIG})"

echo ""
echo "Pipeline submitted. Dependency chain:"
echo "  [${DL_JOB}] download_cat_maf"
echo "    └─ [${EXCL_JOB}] extract_exclusive   (${N_EXCL} tasks)"
echo "    └─ [${SHAR_JOB}] extract_shared      (${N_SHAR} tasks)"
echo "    └─ [${ALIG_JOB}] extract_alignable   (${N_ALIG} tasks)"
echo ""
echo "Monitor with:  squeue -u \$USER"
echo "Logs in:       logs/"

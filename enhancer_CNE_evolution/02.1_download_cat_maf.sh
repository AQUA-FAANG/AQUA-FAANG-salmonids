#!/bin/bash
# =============================================================================
# download_cat_maf.sh
# SLURM job: download all MAF files from SalmoBase and concatenate them into
# a single all_alignments.maf file for downstream mafExtractor use.
# =============================================================================
#SBATCH --job-name=download_cat_maf
#SBATCH --output=logs/download_%j.out
#SBATCH --error=logs/download_%j.err
#SBATCH --time=06:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4

set -euo pipefail

# ── Configuration ────────────────────────────────────────────────────────────
MAF_URL="https://test.salmobase.org/datafiles/datasets/Aqua-Faang/alignments/salmon-trout-pike/maf_ref_Ssal_A/"
MAF_DIR="maf_files"
COMBINED_MAF="all_alignments.maf"
# ────────────────────────────────────────────────────────────────────────────

echo "[$(date '+%F %T')] === Download + Concatenation job started ==="
echo "Destination dir : ${MAF_DIR}"
echo "Combined output : ${COMBINED_MAF}"
echo ""

mkdir -p "${MAF_DIR}"

# ── Download ─────────────────────────────────────────────────────────────────
echo "[$(date '+%F %T')] Downloading MAF files from:"
echo "  ${MAF_URL}"

# -r        recursive
# -l1       one level deep (the directory listing only)
# -nd       no subdirectories — save all files flat in MAF_DIR
# -np       no-parent (don't crawl up)
# -A "*.maf" only download .maf files
# -P        save to MAF_DIR
wget \
    --recursive \
    --level=1 \
    --no-directories \
    --no-parent \
    --accept "*.maf" \
    --directory-prefix="${MAF_DIR}" \
    --progress=dot:mega \
    "${MAF_URL}"

N_FILES=$(ls "${MAF_DIR}"/*.maf 2>/dev/null | wc -l)
echo ""
echo "[$(date '+%F %T')] Downloaded ${N_FILES} MAF files."

if [[ ${N_FILES} -eq 0 ]]; then
    echo "ERROR: No .maf files found in ${MAF_DIR}. Aborting."
    exit 1
fi

# ── Concatenation ────────────────────────────────────────────────────────────
# MAF format:
#   Line 1 of every file: "##maf version=1 ..."   <- keep only from first file
#   Remaining lines: alignment blocks (a/s/i/e/q lines + blank separators)
#
# Strategy:
#   1. Write the full first file (including its header) to COMBINED_MAF
#   2. For every subsequent file, skip the ##maf header line and append the rest
# ────────────────────────────────────────────────────────────────────────────
echo "[$(date '+%F %T')] Concatenating MAF files into ${COMBINED_MAF} ..."

# Sort for reproducibility
mapfile -t MAF_FILES < <(ls "${MAF_DIR}"/*.maf | sort)

# Write first file in full
cp "${MAF_FILES[0]}" "${COMBINED_MAF}"
echo "  [1/${N_FILES}] ${MAF_FILES[0]} (full, with header)"

# Append the rest, stripping the ##maf header line from each
for (( i=1; i<${#MAF_FILES[@]}; i++ )); do
    maf="${MAF_FILES[$i]}"
    # grep -v exits 1 if no lines match (i.e. file has no ##maf line); allow that
    grep -v '^##maf' "${maf}" >> "${COMBINED_MAF}" || true
    echo "  [$((i+1))/${N_FILES}] ${maf}"
done

COMBINED_SIZE=$(du -sh "${COMBINED_MAF}" | cut -f1)
echo ""
echo "[$(date '+%F %T')] Combined MAF written: ${COMBINED_MAF} (${COMBINED_SIZE})"
echo "[$(date '+%F %T')] === Download + Concatenation job finished ==="

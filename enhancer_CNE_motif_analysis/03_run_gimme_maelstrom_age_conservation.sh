#!/bin/bash
#SBATCH --job-name=gimme_maelstrom_age
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH --time=24:00:00
#SBATCH --output=logs/gimme_maelstrom_age_%j.out
#SBATCH --error=logs/gimme_maelstrom_age_%j.err

###############################################################################
# GimmeMotifs Maelstrom Analysis (Age + Conservation category)
#
# Purpose:
#   Run Maelstrom motif enrichment on:
#     1. Age × conservation category interaction clusters
#
# Input:
#   Maelstrom formatted region file (age or interaction-based)
#
# Output:
#   Motif enrichment results per cluster type
#
###############################################################################

set -euo pipefail

# =============================================================================
# INPUT ARGUMENTS
# =============================================================================

INPUT_FILE=${1:Provide Maelstrom input file} # maelstrom_age_category_interaction_200bp_input.txt
OUTPUT_DIR=${2: Provide output directory} # ssal_CNE_age_cons_int_output

# =============================================================================
# CONFIGURATION
# =============================================================================

GENOME="./genomes/Ssal_v3.1_lex/Ssal_v3.1_lex.fa"
MOTIF_DB="gimme.vertebrate.v5.0"

# =============================================================================
# LOGGING
# =============================================================================

echo "================================================="
echo "GimmeMotifs Maelstrom (Age-based workflow)"
echo "================================================="
echo "Host        : $(hostname)"
echo "Job ID      : ${SLURM_JOB_ID:-NA}"
echo "Date        : $(date)"
echo "Input       : $INPUT_FILE"
echo "Output      : $OUTPUT_DIR"
echo "Threads     : ${SLURM_CPUS_PER_TASK}"
echo "================================================="

# =============================================================================
# ENVIRONMENT
# =============================================================================

source ~/.bashrc
conda activate gimme

# =============================================================================
# INPUT CHECKS
# =============================================================================

if [[ ! -f "$INPUT_FILE" ]]; then
    echo "ERROR: Input file not found: $INPUT_FILE"
    exit 1
fi

if [[ ! -f "$GENOME" ]]; then
    echo "ERROR: Genome FASTA not found: $GENOME"
    exit 1
fi

# =============================================================================
# OUTPUT SETUP
# =============================================================================

mkdir -p logs
mkdir -p "$OUTPUT_DIR"

# =============================================================================
# GIMMEMOTIFS CACHE HANDLING
# =============================================================================

CACHE_BASE="${TMPDIR:-/tmp}"
LOCAL_CACHE="${CACHE_BASE}/${USER}_gimme_cache"

mkdir -p "$LOCAL_CACHE"

if [[ -z "${XDG_CACHE_HOME:-}" ]]; then
    XDG_CACHE_HOME="$HOME/.cache"
fi

if [[ -d "$XDG_CACHE_HOME/gimmemotifs" ]]; then
    echo "Copying GimmeMotifs cache to local node storage..."
    cp -r "$XDG_CACHE_HOME/gimmemotifs" "$LOCAL_CACHE/"
fi

export XDG_CACHE_HOME="$LOCAL_CACHE"

echo "Using cache: $XDG_CACHE_HOME"

# =============================================================================
# RUN MAELSTROM
# =============================================================================

echo
echo "Running Maelstrom..."
echo

gimme maelstrom \
    "$INPUT_FILE" \
    "$GENOME" \
    "$OUTPUT_DIR" \
    --no-filter \
    -p "$MOTIF_DB" \
    -N "${SLURM_CPUS_PER_TASK}"

# =============================================================================
# CLEANUP
# =============================================================================

echo
echo "Completed successfully at $(date)"

rm -rf "$LOCAL_CACHE"

echo "================================================="
#!/bin/bash
#SBATCH --job-name=gimme_maelstrom
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH --time=24:00:00
#SBATCH --output=logs/gimme_maelstrom_%j.out
#SBATCH --error=logs/gimme_maelstrom_%j.err

###############################################################################
# GimmeMotifs Maelstrom Analysis
#
# Purpose:
#   Run GimmeMotifs Maelstrom motif enrichment analysis on processed
#   CNE regions generated from the preprocessing workflow.
#
# Inputs:
#   - Maelstrom region file
#   - Reference genome FASTA
#
# Outputs:
#   - Maelstrom motif enrichment results
#
# Usage:
#   sbatch run_gimme_maelstrom.sh \
#       maelstrom_input_age_200bp.txt \
#       ssal_CNE_age_output
#
###############################################################################

set -euo pipefail

# =============================================================================
# INPUT ARGUMENTS
# =============================================================================

INPUT_FILE=${1:Maelstrom input file} #maelstrom_input_200bp.txt
OUTPUT_DIR=${2:output directory} #ssal_CNE_age_output

# =============================================================================
# CONFIGURATION
# =============================================================================

GENOME_FASTA="./genomes/Ssal_v3.1_lex/Ssal_v3.1_lex.fa"
MOTIF_DATABASE="gimme.vertebrate.v5.0"

# =============================================================================
# START LOGGING
# =============================================================================

echo "================================================="
echo "GimmeMotifs Maelstrom Analysis"
echo "================================================="
echo "Host          : $(hostname)"
echo "Job ID        : ${SLURM_JOB_ID:-NA}"
echo "Date          : $(date)"
echo "Input file    : $INPUT_FILE"
echo "Output folder : $OUTPUT_DIR"
echo "Genome        : $GENOME_FASTA"
echo "Threads       : ${SLURM_CPUS_PER_TASK}"
echo "================================================="

# =============================================================================
# LOAD ENVIRONMENT
# =============================================================================

source ~/.bashrc

conda activate gimme

echo "Active environment:"
conda info --envs | grep "*"

# =============================================================================
# CHECK INPUTS
# =============================================================================

if [[ ! -f "$INPUT_FILE" ]]; then
    echo "ERROR: Input file not found:"
    echo "$INPUT_FILE"
    exit 1
fi

if [[ ! -f "$GENOME_FASTA" ]]; then
    echo "ERROR: Genome FASTA not found:"
    echo "$GENOME_FASTA"
    exit 1
fi

# =============================================================================
# CREATE OUTPUT DIRECTORIES
# =============================================================================

mkdir -p logs
mkdir -p "$OUTPUT_DIR"

# =============================================================================
# GIMMEMOTIFS CACHE WORKAROUND
# =============================================================================
#
# Some HPC systems experience filesystem contention when multiple
# jobs access the shared GimmeMotifs cache simultaneously.
#
# Copy cache to node-local storage and redirect access there.
#
# =============================================================================

CACHE_ROOT="${TMPDIR:-/tmp}"
LOCAL_CACHE="${CACHE_ROOT}/${USER}_gimme_cache"

mkdir -p "$LOCAL_CACHE"

if [[ -z "${XDG_CACHE_HOME:-}" ]]; then
    XDG_CACHE_HOME="$HOME/.cache"
fi

if [[ -d "$XDG_CACHE_HOME/gimmemotifs" ]]; then

    echo "Copying GimmeMotifs cache..."

    rsync -a \
        "$XDG_CACHE_HOME/gimmemotifs/" \
        "$LOCAL_CACHE/gimmemotifs/"
fi

export XDG_CACHE_HOME="$LOCAL_CACHE"

echo "Using cache:"
echo "$XDG_CACHE_HOME"

# =============================================================================
# RUN MAELSTROM
# =============================================================================

echo
echo "Starting Maelstrom..."
echo

gimme maelstrom \
    "$INPUT_FILE" \
    "$GENOME_FASTA" \
    "$OUTPUT_DIR" \
    --no-filter \
    -p "$MOTIF_DATABASE" \
    -N "${SLURM_CPUS_PER_TASK}"

# =============================================================================
# CLEANUP
# =============================================================================

echo
echo "Analysis completed successfully."

if [[ -d "$LOCAL_CACHE" ]]; then
    rm -rf "$LOCAL_CACHE"
fi

echo "Finished: $(date)"
echo "================================================="
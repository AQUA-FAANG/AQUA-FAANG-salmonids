#!/bin/bash

#$ -V -cwd
#$ -l h_rt=50:00:00
#$ -l h_vmem=5G
#$ -t 1-29
#$ -tc 29

###############################################################################
# Script: 01_prepare_peak_subsets.sh
#
# Purpose:
#   1. Extract ATAC consensus peaks for a single chromosome.
#   2. Create a chromosome-specific output directory.
#   3. Split the chromosome peak file into five smaller chunks to
#      accelerate downstream MAF extraction.
#
# Input:
#   salmon_3field_consensus_peaks.bed
#
# Output:
#   ssa<CHR>/
#   ├── <CHR>_extracted_peaks.bed
#   ├── <CHR>_subset_extracted_peaksaa
#   ├── <CHR>_subset_extracted_peaksab
#   ├── ...
#
# Example:
#   Task 1 processes chromosome 1
#   Task 2 processes chromosome 2
#   ...
###############################################################################

set -euo pipefail

CHR=${SGE_TASK_ID}

INPUT_BED="salmon_3field_consensus_peaks.bed"
OUTPUT_DIR="ssa${CHR}"

echo "Processing chromosome ${CHR}"

# Create chromosome-specific directory
mkdir -p "${OUTPUT_DIR}"

# Extract peaks belonging to this chromosome
awk -v chr="${CHR}" '$1 == chr' "${INPUT_BED}" \
    > "${CHR}_extracted_peaks.bed"

# Split into five approximately equal chunks
split \
    -n l/5 \
    -a 2 \
    "${CHR}_extracted_peaks.bed" \
    "${CHR}_subset_extracted_peaks"

echo "Finished chromosome ${CHR}"
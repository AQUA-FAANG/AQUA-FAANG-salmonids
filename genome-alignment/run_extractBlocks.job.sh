#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem=20G
#SBATCH --constraint=avx2
#SBATCH --job-name=extractBlock

module load SAMtools

set -x

for script in $(dir data/blocks/*/extractSequences.sh); do
    source "$script"
done



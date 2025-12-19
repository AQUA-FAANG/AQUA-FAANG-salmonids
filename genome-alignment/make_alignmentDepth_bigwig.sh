#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem=4G
#SBATCH --constraint=avx2
#SBATCH --job-name=halAlignmentDepth
#SBATCH --array 1-29

set -e

# Usage: sbatch make_alignmentDepth_bigwig <OUTDIR> <chrom.sizes file> <Extra parameters for halAlignmentDepth>
OUTDIR=$1
CHROMSIZES=$2
# Shift the positional parameters to the left, so that remaining parameters can be passed as $@
shift 2

# Example parameters:
# * Ssal to all encestors: Ssal_A --step 1 --targetGenomes Anc0,Anc1,Anc2
# * Ssal to duplicate: Ssal_A --step 1 --noAncestors --targetGenomes Ssal_B
# * Omyk to duplicate: Omyk_A --step 1 --noAncestors --targetGenomes Omyk_B

CACTUS="singularity exec --contain --bind data:/data --pwd / cactus_v2.7.1.sif"

mkdir -p $OUTDIR


# Loop through chromosomes 1 to 29
#for CHROM_NUM in {1..29}
#do

#  CHROM="$(printf "ssa%02d" $CHROM_NUM)"
  CHROM="$(printf "ssa%02d" $SLURM_ARRAY_TASK_ID)"

  echo "Processing $CHROM..."
  
  $CACTUS halAlignmentDepth data/blocks/$CHROM/output.hal $@ \
    | python fix_wiggle.py stdin > $OUTDIR/$CHROM.wig

  wigToBigWig $OUTDIR/$CHROM.wig $CHROMSIZES $OUTDIR/$CHROM.bw

  rm $OUTDIR/$CHROM.wig
#done

echo "Done."

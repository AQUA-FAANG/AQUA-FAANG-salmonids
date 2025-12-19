#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem=5G
#SBATCH --constraint=avx2
#SBATCH --job-name=maf_fix
#SBATCH --array=1-29%10

# if $1 is empty, use Ssal_A as the reference genome
if [ -z "$1" ]; then
  REFGENOME="Ssal_A"
else
  REFGENOME=$1
fi

BLOCK_DIR="data/blocks/ssa$(printf "%02d" $SLURM_ARRAY_TASK_ID)"
MAF_FILE=$BLOCK_DIR/$REFGENOME.maf
OUTDIR=data/maf_ref_$REFGENOME
mkdir -p $OUTDIR
MAF_FIXED="$OUTDIR/ssa$(printf "%02d" $SLURM_ARRAY_TASK_ID)_AB.maf"

echo Converting $MAF_FILE to $MAF_FIXED...

python fix_maf.py --keep_subgenome $MAF_FILE chromSizes.tsv > $MAF_FIXED

echo Compressing...

gzip $MAF_FIXED

echo Done.
#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem=10G
#SBATCH --constraint=avx2
#SBATCH --job-name=mergepafs
#SBATCH --array=1-4

ALL_SUFFIX=("Omyk_vs_Omyk" "Omyk_vs_Ssal" "Ssal_vs_Omyk" "Ssal_vs_Ssal")
ALL_PREFIX=("ref_Omyk_A" "ref_Omyk_A" "ref_Ssal_A" "ref_Ssal_A")
SUFFIX=${ALL_SUFFIX[$SLURM_ARRAY_TASK_ID-1]}
PREFIX=${ALL_PREFIX[$SLURM_ARRAY_TASK_ID-1]}
echo Task: $SLURM_ARRAY_TASK_ID
echo Suffix: $SUFFIX
INDIR="data/paf/$PREFIX"
echo InDir: $INDIR

OUTDIR="data/paf/merged"
mkdir -p $OUTDIR
cat $INDIR/*_$SUFFIX.paf | sort -k 6,8n -t$'\t' > $OUTDIR/$SUFFIX.paf

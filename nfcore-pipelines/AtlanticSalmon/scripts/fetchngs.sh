#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem=4G
#SBATCH --job-name=fetchngs
#SBATCH --output=log-fetchngs-%j.job

module load Nextflow/21.03
export NXF_SINGULARITY_CACHEDIR="/cluster/software/nf-core/CACHE"

SCRIPTS=/mnt/project/Aqua-Faang/nfcore/AtlanticSalmon/scripts

$SCRIPTS/nextflow run nf-core/fetchngs \
  -r 1.10.0 -profile singularity \
  -config $SCRIPTS/orion.config \
  --input seq_ids.txt \
  --outdir seq_data \
  --force_sratools_download \
  -resume
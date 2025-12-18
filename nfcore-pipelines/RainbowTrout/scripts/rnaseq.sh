#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem=4G
#SBATCH --job-name=rnaseq
#SBATCH --output=log-rnaseq-%j.job

module load Nextflow/21.03
export NXF_SINGULARITY_CACHEDIR="/cluster/software/nf-core/CACHE"
export NXF_OPTS='-Xms1g -Xmx4g'

SAMPLES=design.csv
SCRIPTS=/mnt/project/Aqua-Faang/nfcore/RainbowTrout/scripts
ORIONCONFIG=$SCRIPTS/orion.config
VERSION=3.8.1
GENOME=/mnt/project/Aqua-Faang/nfcore/RainbowTrout/genome
FASTA=$GENOME/Oncorhynchus_mykiss-GCA_013265735.3-softmasked.fa
GTF=$GENOME/Oncorhynchus_mykiss-GCA_013265735.3-2020_12-genes.gtf
NAME=$(pwd . | sed s/.*nfcore\\/// | sed s/\\//_/g)
WORK=/mnt/ScratchProjects/Aqua-Faang/nfcore/$NAME

$SCRIPTS/nextflow run nf-core/rnaseq \
  -r $VERSION -profile singularity \
  --input $SAMPLES \
  --outdir results \
  --multiqc_title $NAME \
  -c $ORIONCONFIG \
  --fasta $FASTA \
  --gtf $GTF \
  -resume \
  --validate_params 0 \
  -w $WORK

#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem=4G
#SBATCH --job-name=chipseq
#SBATCH --output=log-chipseq-%j.job

module load Nextflow/21.03
module list

export NXF_SINGULARITY_CACHEDIR="/cluster/software/nf-core/CACHE"
export NXF_OPTS='-Xms1g -Xmx4g'

SAMPLES=design.csv
SCRIPTS=/mnt/project/Aqua-Faang/seq_results/RainbowTrout/scripts
ORIONCONFIG=$SCRIPTS/orion.config
VERSION=1.2.2
GENOME=/mnt/project/Aqua-Faang/seq_results/RainbowTrout/genome
FASTA=$GENOME/Oncorhynchus_mykiss-GCA_013265735.3-softmasked.fa
GTF=$GENOME/Oncorhynchus_mykiss-GCA_013265735.3-2020_12-genes.gtf
GSIZE=2070645207
NAME=$(pwd . | sed s/.*seq_results\\/// | sed s/\\//_/g)
WORK=/mnt/ScratchProjects/Aqua-Faang/nfcore/$NAME

# Broad peaks
$SCRIPTS/nextflow run nf-core/chipseq \
  --input $SAMPLES \
  -profile singularity \
  -r $VERSION \
  -c $ORIONCONFIG \
  --fasta $FASTA \
  --gtf $GTF \
  --macs_gsize $GSIZE \
  --name $NAME \
  -resume \
  -w $WORK

# Narrow peaks
$SCRIPTS/nextflow run nf-core/chipseq \
  --input $SAMPLES \
  -profile singularity \
  -r $VERSION \
  -c $ORIONCONFIG \
  --fasta $FASTA \
  --gtf $GTF \
  --macs_gsize $GSIZE \
  --name $NAME \
  -resume \
  --narrow_peak \
  -w $WORK

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
SCRIPTS=/mnt/project/Aqua-Faang/nfcore/AtlanticSalmon/scripts
ORIONCONFIG=$SCRIPTS/orion.config
VERSION=1.2.2
GENOME=/mnt/project/Aqua-Faang/nfcore/AtlanticSalmon/genome
FASTA=$GENOME/Salmo_salar-GCA_905237065.2-softmasked.fa
GTF=$GENOME/Salmo_salar-GCA_905237065.2-softmasked_genes.gtf
BWA_INDEX=$GENOME/BWAIndex/Salmo_salar.Ssal_v3.1.dna_sm.toplevel.fa
GSIZE=2756563003
NAME=$(pwd . | sed s/.*nfcore\\/// | sed s/\\//_/g)
WORK=/mnt/ScratchProjects/Aqua-Faang/nfcore/$NAME

# Broad peaks
$SCRIPTS/nextflow run nf-core/chipseq \
  --input $SAMPLES \
  -profile singularity \
  -r $VERSION \
  -c $ORIONCONFIG \
  --fasta $FASTA \
  --gtf $GTF \
  --bwa_index $BWA_INDEX \
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
  --bwa_index $BWA_INDEX \
  --macs_gsize $GSIZE \
  --name $NAME \
  -resume \
  -w $WORK \
  --narrow_peak

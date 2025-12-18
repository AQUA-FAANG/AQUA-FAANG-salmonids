#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem=4G
#SBATCH --job-name=atacseq
#SBATCH --output=log-atacseq-%j.out

module load Nextflow/21.03

export NXF_SINGULARITY_CACHEDIR="/cluster/software/nf-core/CACHE"
export NXF_OPTS='-Xms1g -Xmx4g'

SAMPLES=design.csv
SCRIPTS=/path/to/project/Aqua-Faang/nfcore/AtlanticSalmon/scripts
VERSION=1.2.1
GENOME=/path/to/project/Aqua-Faang/nfcore/AtlanticSalmon/genome
FASTA=$GENOME/Salmo_salar-GCA_905237065.2-softmasked.fa
GTF=$GENOME/Salmo_salar-GCA_905237065.2-softmasked_genes.gtf
GSIZE=2756563003
NAME=$(pwd . | sed s/.*nfcore\\/// | sed s/\\//_/g)
WORK=/path/to/ScratchProjects/Aqua-Faang/nfcore/$NAME

$SCRIPTS/nextflow run nf-core/atacseq \
   -r $VERSION -profile singularity \
   --input $SAMPLES \
   --outdir results \
   --fasta $FASTA \
   --gtf $GTF \
   --macs_gsize $GSIZE \
   -c $SCRIPTS/orion.config \
   --name $NAME \
   -resume \
   -w $WORK \
   --narrow_peak

#!/bin/bash
#SBATCH --ntasks=12
#SBATCH --nodes=1
#SBATCH --mem=100G
#SBATCH --job-name=maelstrom
#SBATCH --output=maelstrom-%j.log

## Create and activate gimme environment
# module load Miniconda3
# conda activate micromamba
## If not created: 
#     mamba create -c conda-forge -c bioconda -c defaults -n gimme gimmemotifs
# mamba activate gimme

echo "Running gimme maelstom"

## AtlanticSalmon
  
echo "Atlantic salmon DevMap promoters"
gimme maelstrom \
  AtlanticSalmon_DevMap_promoters.tsv \
  Salmo_salar-GCA_905237065.2-softmasked.fa \
  AtlanticSalmon_DevMap_promoters \
  -p JASPAR2022_vertebrates \
  -N 12

echo "Atlantic salmon DevMap enhancers"
gimme maelstrom \
  AtlanticSalmon_DevMap_enhancers.tsv \
  Salmo_salar-GCA_905237065.2-softmasked.fa \
  AtlanticSalmon_DevMap_enhancers \
  -p JASPAR2022_vertebrates \
  -N 12
  
echo "Atlantic salmon BodyMap promoters"
gimme maelstrom \
  AtlanticSalmon_BodyMap_promoters.tsv \
  Salmo_salar-GCA_905237065.2-softmasked.fa \
  AtlanticSalmon_BodyMap_promoters \
  -p JASPAR2022_vertebrates \
  -N 12
  
echo "Atlantic salmon BodyMap enhancers"
gimme maelstrom \
  AtlanticSalmon_BodyMap_enhancers.tsv \
  Salmo_salar-GCA_905237065.2-softmasked.fa \
  AtlanticSalmon_BodyMap_enhancers \
  -p JASPAR2022_vertebrates \
  -N 12
  
## RainbowTrout
  
echo "Rainbow trout DevMap promoters"
gimme maelstrom \
  RainbowTrout_DevMap_promoters.tsv \
  Oncorhynchus_mykiss-GCA_013265735.3-softmasked.fa \
  RainbowTrout_DevMap_promoters \
  -p JASPAR2022_vertebrates \
  -N 12

echo "Rainbow trout DevMap enhancers"
gimme maelstrom \
  RainbowTrout_DevMap_enhancers.tsv \
  Oncorhynchus_mykiss-GCA_013265735.3-softmasked.fa \
  RainbowTrout_DevMap_enhancers \
  -p JASPAR2022_vertebrates \
  -N 12
  
echo "Rainbow trout BodyMap promoters"
gimme maelstrom \
  RainbowTrout_BodyMap_promoters.tsv \
  Oncorhynchus_mykiss-GCA_013265735.3-softmasked.fa \
  RainbowTrout_BodyMap_promoters \
  -p JASPAR2022_vertebrates \
  -N 12
  
echo "Rainbow trout BodyMap enhancers"
gimme maelstrom \
  RainbowTrout_BodyMap_enhancers.tsv \
  Oncorhynchus_mykiss-GCA_013265735.3-softmasked.fa \
  RainbowTrout_BodyMap_enhancers \
  -p JASPAR2022_vertebrates \
  -N 12
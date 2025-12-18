#!/bin/bash
#SBATCH --array=1-208%30
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --constraint=avx2
#SBATCH --mem=3G
#SBATCH --job-name=bamtobed
#SBATCH --output=logs/bamtobed-%j_%a.log

# For an input folder of bams going into chromHMM:
# Convert bam files from nfcore results to bed files, 
# while filtering out reads mapping to backlisted regions of the genome

FILE_LIST=( $( ls bams/*bam | sed 's/bams\///' | sed 's/.bam//' ) )
FILE=${FILE_LIST[$SLURM_ARRAY_TASK_ID - 1]}
GROUP=$( echo $FILE | sed 's/_R..mLb.clN.sorted//' | sed 's/[^_]*_[^_]*_//' )
BLACKLIST=$(ls *blacklist.bed )

module load BEDTools

echo "Convert bams/$FILE.bam to bed/$GROUP/$FILE.bed, using $BLACKLIST"

mkdir -p bed/$GROUP/

bedtools bamtobed -i bams/$FILE.bam | bedtools subtract -a stdin -b $BLACKLIST > bed/$GROUP/$FILE.bed

grep $GROUP design.tsv  | sed 's/bam/bed/g' > bed/$GROUP/design.tsv

echo "Done."

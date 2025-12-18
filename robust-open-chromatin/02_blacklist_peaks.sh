#!/bin/bash
#SBATCH --array=1-125%30
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --constraint=avx2
#SBATCH --job-name=blacklist_peaks
#SBATCH --mem=1G
#SBATCH --partition=smallmem
#SBATCH --out=logs/blacklist_peaks_%A_%a.log

FILE_LIST=( $( find ATAC_narrowPeak -name "*Peak" -type f ) )
FILE=${FILE_LIST[$SLURM_ARRAY_TASK_ID - 1]}
OUT=$( echo $FILE | sed s/ATAC_narrowPeak// )

SPECIES="AtlanticSalmon"
if [[ "$FILE" == *"RainbowTrout"* ]]; then
      SPECIES="RainbowTrout"
fi

MAP="DevMap"
if [[ "$FILE" == *"BodyMap"* ]]; then
      MAP="BodyMap"
fi

DIR=blacklisted_peaks/$SPECIES/$MAP
if [[ ! -d $DIR ]]; then
  mkdir -p $DIR
fi

OUT=$DIR/$( echo $FILE | sed 's/.*\///g' ).blacklisted

BLACKLIST=$(ls ../blacklist/$SPECIES/${SPECIES}_blacklist.bed)

module load BEDTools

cat $FILE | \
  # Remove Atlantic salmon contigs
  grep -v "CAJNNT" | \
  # Remove Rainbow trout contigs
  grep -v "JAAXML" | \
  # Subtract regions of peaks overlapping the blacklist regions
  bedtools intersect -v -f 1 -F 1 -e -a stdin -b $BLACKLIST > $OUT

#!/bin/bash
#SBATCH --array=1-42
#SBATCH --ntasks=1               
#SBATCH --nodes=1                
#SBATCH --job-name=merge_peaks      
#SBATCH --mem=1G                    
#SBATCH --partition=smallmem          
#SBATCH --out=logs/merge_peaks_%A_%a.log

FILE_LIST=( $( find IDR_peaks -type f -name "*narrowPeak" | sed 's/_R.*//' | sort -u ) )
FILE=${FILE_LIST[$SLURM_ARRAY_TASK_ID - 1]}
NAME=$( echo $FILE | sed 's/.*\///g' )

SPECIES="AtlanticSalmon"
if [[ "$FILE" == *"RainbowTrout"* ]]; then
      SPECIES="RainbowTrout"
fi

MAP="DevMap"
if [[ "$FILE" == *"BodyMap"* ]]; then
      MAP="BodyMap"
fi

OUT=merged_peaks/$SPECIES/$MAP/
if [[ ! -d "$OUT" ]]; then
  mkdir -p $OUT
fi

module load BEDTools

cat $FILE*narrowPeak | sort -k1,1 -k2,2n | \
  bedtools merge -d -1 -c 4,5,6 -o distinct,max,distinct -i stdin > $OUT/$NAME.bed


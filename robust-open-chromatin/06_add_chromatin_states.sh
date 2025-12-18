#!/bin/bash
#SBATCH --array=1-38
#SBATCH --ntasks=1               
#SBATCH --nodes=1                
#SBATCH --job-name=add_chromatin_states   
#SBATCH --mem=1G             
#SBATCH --partition=smallmem
#SBATCH --out=logs/add_chromatin_states_%J_%a.log

FILE_LIST=( $( find summits -type f -name "*summits" | grep -v ".*RainbowTrout.*Gonad.*" | sort ) )
FILE=${FILE_LIST[$SLURM_ARRAY_TASK_ID - 1]}
NAME=$( echo $FILE | sed 's/.*\///g' | sed 's/.summits//' )

SPECIES="AtlanticSalmon"
if [[ "$FILE" == *"RainbowTrout"* ]]; then
      SPECIES="RainbowTrout"
fi

MAP="DevMap"
if [[ "$FILE" == *"BodyMap"* ]]; then
      MAP="BodyMap"
fi

OUT=add_chromatin_states/$SPECIES/$MAP
if [[ ! -d "$OUT" ]]; then
  mkdir -p $OUT
fi

STATES=../chromatin_states/$SPECIES/$MAP/annotated_bed/$( echo $NAME".bed" | sed 's/ATAC_//' )
PEAKS=merged_peaks/$SPECIES/$MAP/$NAME.bed

module load BEDTools

bedtools intersect -a $PEAKS -b $FILE -wa -wb | bedtools intersect -a stdin -b $STATES -wa -wb -loj > $OUT/$NAME.bed




#!/bin/bash
#SBATCH --array=1-190%20
#SBATCH --ntasks=1               
#SBATCH --nodes=1
#SBATCH --job-name=expression
#SBATCH --partition=orion
#SBATCH --mem=3G
#SBATCH --out=logs/expression_%A_%a.log

FILE_LIST=( $( find ../chromatin_states/{AtlanticSalmon,RainbowTrout}/{BodyMap,DevMap}/bams -name "*.bam" | grep -vi "input" | sed 's/_R.*bam//' | sort -u ) )
FILE=${FILE_LIST[$SLURM_ARRAY_TASK_ID - 1]}
NAME=$( echo $FILE | sed 's/.*\///g' )

SPECIES="AtlanticSalmon"
if [[ "$FILE" == *"RainbowTrout"* ]]; then
      SPECIES="RainbowTrout"
fi

PEAKS=unified_annotated_peaks/$SPECIES"_unified_peaks.saf"

MAP="DevMap"
if [[ "$FILE" == *"BodyMap"* ]]; then
      MAP="BodyMap"
fi

OUT=unified_annotated_peaks/assay_expression/counts/$SPECIES/$MAP
if [[ ! -d "$OUT" ]]; then
  mkdir -p $OUT
fi

scripts/subread-2.0.6/bin/featureCounts -O -p --countReadPairs -a $PEAKS -F SAF -o $OUT/$NAME'_counts'.txt $FILE*
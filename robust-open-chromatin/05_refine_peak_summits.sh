#!/bin/bash
#SBATCH --array=1-42
#SBATCH --ntasks=1               
#SBATCH --nodes=1                
#SBATCH --job-name=refine_peaks    
#SBATCH --mem=3G                
#SBATCH --constraint=avx2
#SBATCH --partition=smallmem          
#SBATCH --out=logs/refine_peaks_%A_%a.log

FILE_LIST=( $( find merged_peaks -type f -name "*bed" ) )
FILE=${FILE_LIST[$SLURM_ARRAY_TASK_ID - 1]}
NAME=$( echo $FILE | sed 's/.*\///g' | sed 's/.bed//' )

SPECIES="AtlanticSalmon"
if [[ "$FILE" == *"RainbowTrout"* ]]; then
      SPECIES="RainbowTrout"
fi

MAP="DevMap"
if [[ "$FILE" == *"BodyMap"* ]]; then
      MAP="BodyMap"
fi

OUT=summits/$SPECIES/$MAP
if [[ ! -d "$OUT" ]]; then
  mkdir -p $OUT
fi

BAM=$( find ../../seq_results/$SPECIES/$MAP/ -type f -name "$NAME.mRp.clN.sorted.bam" )

# conda activate micromamba
# micromamba activate macs2

macs2 refinepeak -b $FILE -i $BAM -o $OUT/$NAME.summits



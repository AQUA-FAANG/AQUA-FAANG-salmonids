#!/bin/bash
#SBATCH --array=1-42
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --job-name=idr
#SBATCH --mem=1G
#SBATCH --partition=smallmem
#SBATCH --out=logs/idr_%A_%a.log

# On IDR thresholds:
# The following IDR thresholds are recommended for the different types of IDR analyses.
# 
# For self-consistency and comparison of true replicates
# 
#     If you started with ~150 to 300K relaxed pre-IDR peaks for large genomes (human/mouse), 
#     then threshold of 0.01 or 0.02 generally works well. 
#     
#     If you started with < 100K pre-IDR peaks for large genomes (human/mouse), then threshold 
#     of 0.05 is more appropriate. This is because the IDR sees a smaller noise component and 
#     the IDR scores get weaker. This is typically for use with peak callers that are unable to 
#     be adjusted to call large number of peaks (eg. PeakSeq or QuEST)
# 
#     For smaller genomes such as worm, if you start with ~15K to 40K peaks then once again IDR 
#     thresholds of 0.01 or 0.02 work well.
#     
#     For self-consistency analysis of datasets with shallow sequencing depths, you can use an 
#     IDR threshold as relaxed as 0.1 if you start with < 100K pre-IDR peaks.
# 
THRESHOLD=0.1

# Get unique group names from peak files
FILE_LIST=( $( find blacklisted_peaks -type f -name "*narrowPeak.blacklisted" | sed 's/_R.*//' | sort -u ) )
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

OUT=IDR_peaks/$SPECIES/$MAP/
if [[ ! -d "$OUT" ]]; then
  mkdir -p $OUT
fi

FILE1=$( ls ${FILE}_R1* )
if [ -f "$FILE1" ]; then

  FILE2=$( ls ${FILE}_R2* )
  if [ -f "$FILE2" ]; then

    echo "Compare replicate 1 vs 2"
    idr --samples $FILE1 $FILE2 \
      --input-file-type narrowPeak \
      --rank p.value \
      --soft-idr-threshold $THRESHOLD \
      --output-file $OUT/$NAME"_R1vR2_idr.narrowPeak" \
      --output-file-type narrowPeak \
      --plot \
      --log-output-file $OUT/$NAME"_R1vR2_idr.log"
      
    sed -i -E 's/^(([^\t]+\t){3})\./\1R1_R2/' $OUT/$NAME"_R1vR2_idr.narrowPeak"
    
    FILE3=$( ls ${FILE}_R3* )
    if [ -f "$FILE3" ]; then
      
      echo "Compare replicate 1 vs 3"
      idr --samples $FILE1 $FILE3 \
        --input-file-type narrowPeak \
        --rank p.value \
        --soft-idr-threshold $THRESHOLD \
        --output-file $OUT/$NAME"_R1vR3_idr.narrowPeak" \
        --output-file-type narrowPeak \
        --plot \
        --log-output-file $OUT/$NAME"_R1vR3_idr.log"
      
      sed -i -E 's/^(([^\t]+\t){3})\./\1R1_R3/' $OUT/$NAME"_R1vR3_idr.narrowPeak"
        
      echo "Compare replicate 2 vs 3"
      idr --samples $FILE2 $FILE3 \
        --input-file-type narrowPeak \
        --rank p.value \
        --soft-idr-threshold $THRESHOLD \
        --output-file $OUT/$NAME"_R2vR3_idr.narrowPeak" \
        --output-file-type narrowPeak \
        --plot \
        --log-output-file $OUT/$NAME"_R2vR3_idr.log"
        
      sed -i -E 's/^(([^\t]+\t){3})\./\1R2_R3/' $OUT/$NAME"_R2vR3_idr.narrowPeak"
      
      echo "NICE: 3 replicates for $FILE"
    else 
      echo "OK: 2 replicates, FILE3 does not exist"
    fi
  else 
    echo "WARNING: 1 replicate, FILE2 does not exist"
  fi
else 
  echo "ERROR: no replicates, FILE1 does not exist"
fi

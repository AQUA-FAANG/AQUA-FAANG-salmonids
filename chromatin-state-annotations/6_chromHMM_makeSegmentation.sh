#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem=10G
#SBATCH --constraint=avx2
#SBATCH --partition=orion
#SBATCH --job-name=chromHMM_makeSegmentation
#SBATCH --output=logs/chromHMM_makeSegmentation-%J.log

# Run ChromHMM makeSegmentation to segment the genome into 200bp bins,
# assinging each bin to one of the learned chromatin states.
# Select learn model to use e.g. 12 state model.

# For rainbow trout BodyMap, separate models were used for some tissues,
# based on differences in the strength of signal coming from some ChIP assays.
# Described in the README.


unset DISPLAY

OUTDIR=$(pwd)
ASSEMBLY=$(cat assembly.txt)
CHROMHMM=/mnt/project/Aqua-Faang/chromatin_states/scripts/ChromHMM

# Select model
MODEL=$OUTDIR/models/model_12/model_12.txt

echo "Make segmentations"
CMD="MakeSegmentation $MODEL $(pwd)/binary $OUTDIR/segmentations"
echo $CMD
cd $CHROMHMM
java -mx10000M -jar ChromHMM.jar $CMD
echo "Done."

# Run only for specific tissues for rainbow trout BodyMap
# TISSUE=""
# MODEL=$OUTDIR/$TISSUE/models/model_10/model_10.txt
# echo "Make segmentations"
# CMD="MakeSegmentation $MODEL $(pwd)/$TISSUE/binary $OUTDIR/segmentations"
# echo $CMD
# cd $CHROMHMM
# java -mx10000M -jar ChromHMM.jar $CMD
# echo "Done."


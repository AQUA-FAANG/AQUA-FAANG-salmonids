#!/bin/bash
#SBATCH --array=1-4
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem=100G
#SBATCH --constraint=avx2
#SBATCH --partition=hugemem
#SBATCH --job-name=binarizeBed
#SBATCH --output=logs/binarizeBed-%A_%a.log

# Binarize the blacklisted bed files of assay data going into chromHMM

unset DISPLAY

ASSEMBLY=$(cat assembly.txt)
CHROMHMM=/mnt/project/Aqua-Faang/chromatin_states/scripts/ChromHMM
SIZES=$CHROMHMM/CHROMSIZES/$ASSEMBLY.txt

GROUP_LIST=( $( ls bed ) )
GROUP=${GROUP_LIST[$SLURM_ARRAY_TASK_ID - 1]}

echo "Binarize bed files"
CMD="BinarizeBed $SIZES bed/$GROUP bed/$GROUP/design.tsv binary"
echo $CMD
java -mx100000M -jar $CHROMHMM/ChromHMM.jar $CMD
echo "Done."


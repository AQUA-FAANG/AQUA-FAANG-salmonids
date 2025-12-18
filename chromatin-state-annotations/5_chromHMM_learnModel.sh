#!/bin/bash
#SBATCH --array=10-13
#SBATCH --ntasks=20
#SBATCH --nodes=1
#SBATCH --mem=100G
#SBATCH --constraint=avx2
#SBATCH --job-name=chromHMM_learnModel
#SBATCH --output=logs/chromHMM_learnModel-%j_%a.log

# Run ChromHMM learnModel on binarized input assay data.
# Number of states to model can be set as a range with the SLURM array option.
# e.g. 10-13 states

unset DISPLAY

OUTDIR=$(pwd)/models
ASSEMBLY=$(cat assembly.txt)
CHROMHMM=/mnt/project/Aqua-Faang/chromatin_states/scripts/ChromHMM
NSTATES=$SLURM_ARRAY_TASK_ID
CPUS=20

echo "Learn model"
CMD="LearnModel -p $CPUS $(pwd)/binary $OUTDIR/model_$NSTATES $NSTATES $ASSEMBLY"
echo $CMD
cd $CHROMHMM
java -mx100000M -jar ChromHMM.jar $CMD
echo "Done."

#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem=5G
#SBATCH --constraint=avx2
#SBATCH --job-name=hal2maf
#SBATCH --array=1-29%10

CACTUS_SIF="cactus_v2.4.4.sif"

# if $1 is empty, use Ssal_A as the reference genome
if [ -z "$1" ]; then
  REFGENOME="Ssal_A"
else
  REFGENOME=$1
fi


BLOCK_DIR="data/blocks/ssa$(printf "%02d" $SLURM_ARRAY_TASK_ID)"
HAL_FILE=$BLOCK_DIR/output.hal
MAF_FILE=$BLOCK_DIR/$REFGENOME.maf
echo Generating $MAF_FILE

mkdir -p data/salmon-trout-pike
MAF_FIXED="data/salmon-trout-pike/ssa$(printf "%02d" $SLURM_ARRAY_TASK_ID).maf"


singularity exec --contain --bind .:/data -H /data  $CACTUS_SIF hal2maf \
  --maxBlockLen 1000000 --noAncestors --noDupes \
  --refGenome $REFGENOME \
  $HAL_FILE $MAF_FILE

# Note: the binaries in the image are in the /home directory which is bound by default on Orion.
# the --contain option overrides default bindings.


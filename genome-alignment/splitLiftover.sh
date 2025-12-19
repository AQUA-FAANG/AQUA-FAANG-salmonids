#!/bin/bash

set -e

CACTUS_SIF="cactus_v2.7.1.sif"
SYNTENY_TBL="Salmonid_Synteny_for_alignments_2023.04.12.xlsx"

# Usage: splitLiftover.sh <in_bed> <out_bed> <src_genome> <target_genome> 
# 
# PSL support: If <out_bed> ends with ".psl" it will produce psl (with names) files instead of bed
# 
# Note: Assumes the in/out bed files are under the data/ directory.
# Also that the hal files are data/blocks/ssa[1-29]/output.hal
#
#
in_bed=$1; out_bed=$2; src_genome=$3; target_genome=$4


# Get extension of output file:
out_ext=${out_bed##*.}



echo "#### split the input bed file (${in_bed}) ####"

SPLIT_DIR=$in_bed.split.$src_genome
mkdir $SPLIT_DIR

module load R/4.3.3

Rscript split_bed_to_blocks.R $in_bed $SYNTENY_TBL $src_genome $SPLIT_DIR/

module unload R



echo "#### liftover ####"


SPLITLIFT_DIR=$out_bed.splitlift.$target_genome
mkdir $SPLITLIFT_DIR

# check if output is PSL
if [[ "$out_ext" == "psl" ]];then
  outPSLWithName="--outPSLWithName"
else
  outPSLWithName=""
fi

# iterate over split input files
for INPUT_FILE in $SPLIT_DIR/*.bed; do
  filename=$(basename "$INPUT_FILE")
  SSABLOCK="${filename%.*}"
  OUTPUT_FILE=$SPLITLIFT_DIR/$SSABLOCK.$out_ext
  HAL_FILE=data/blocks/$SSABLOCK/output.hal
  # halLiftover [Options] <halFile> <srcGenome> <srcBed> <tgtGenome> <tgtBed>
  singularity exec --contain --bind data:/data --pwd /  $CACTUS_SIF halLiftover \
    $outPSLWithName \
    $HAL_FILE $src_genome $INPUT_FILE $target_genome $OUTPUT_FILE
done



echo "#### Fix coordinates and merge ####"


# check if output is PSL
if [[ "$out_ext" == "psl" ]];then
  for pslfile in $SPLITLIFT_DIR/*.psl; do
    python fix_psl.py $pslfile >> $out_bed
  done
else
  for bedfile in $SPLITLIFT_DIR/*.bed; do
    python fix_bed.py $bedfile >> $out_bed
  done
fi

echo "#### Cleanup.. remove temporary files ####"

rm -r $SPLIT_DIR
rm -r $SPLITLIFT_DIR

echo "DONE."



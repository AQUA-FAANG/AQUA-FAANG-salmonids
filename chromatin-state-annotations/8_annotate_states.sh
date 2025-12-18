#!/bin/bash
#SBATCH --array=1-38
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem=3G
#SBATCH --job-name=annotate_states
#SBATCH --output=logs/annotate_states-%A-%a.log

# Annotate ChromHMM segmentation bed files with curated state functional labels and colours.

cd /mnt/project/Aqua-Faang/chromatin_states

FILE_LIST=( $( find AtlanticSalmon/DevMap/segmentations/ \
                    AtlanticSalmon/BodyMap/segmentations/ \
                    RainbowTrout/DevMap/segmentations/ \
                    RainbowTrout/BodyMap/segmentations/ \
                    -type f -name "*segments.bed" | sort ) )
FILE=${FILE_LIST[$SLURM_ARRAY_TASK_ID - 1]}
NAME=$( echo $FILE | sed 's/.*\///g' | sed -E 's/_[^_]+_segments.bed//' )

SPECIES="AtlanticSalmon"
if [[ "$FILE" == *"RainbowTrout"* ]]; then
      SPECIES="RainbowTrout"
fi

MAP="DevMap"
if [[ "$FILE" == *"BodyMap"* ]]; then
      MAP="BodyMap"
fi

OUT=$SPECIES/$MAP/annotated_bed
if [[ ! -d "$OUT" ]]; then
  mkdir -p $OUT
fi

if [[ "$FILE" == *"RainbowTrout/BodyMap"* ]]; then
      MAP=$( echo $FILE | sed 's/.*BodyMap\/segmentations/BodyMap/g' | sed 's/_.*//g' )
fi

ANNOTATION=$SPECIES/$MAP/state_annotation.tsv
COLOURS=$SPECIES/$MAP/annotation_colours.tsv

java -mx3000M -jar scripts/ChromHMM/ChromHMM.jar MakeBrowserFiles -c $COLOURS -m $ANNOTATION $FILE $NAME $OUT/$SPECIES"_"$NAME

sed -E 's/^(([^\t]+\t){3})[0-9]+_/\1/' $OUT/$SPECIES"_"$NAME*dense.bed > $OUT/$SPECIES"_"$NAME".bed"

rm $OUT/$SPECIES"_"$NAME*expanded.bed
rm $OUT/$SPECIES"_"$NAME*dense.bed

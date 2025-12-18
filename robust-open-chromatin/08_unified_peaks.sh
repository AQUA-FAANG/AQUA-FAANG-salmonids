#!/bin/bash
#SBATCH --ntasks=1               
#SBATCH --nodes=1                
#SBATCH --job-name=unified_peaks      
#SBATCH --mem=3G                    
#SBATCH --out=logs/unified_peaks_%J.log

OUT=unified_annotated_peaks

if [[ ! -d "$OUT" ]]; then
  mkdir -p $OUT
fi

module load BEDTools

for species in AtlanticSalmon RainbowTrout; do
  unified_bed="$OUT/${species}_unified_peaks.bed"
  unified_saf="$OUT/${species}_unified_peaks.saf"

  # Build unified BED
  find "annotated_peaks/${species}" -type f -name "*bed" -print0 | \
    xargs -0 awk -F $'\t' 'BEGIN{OFS="\t"} {print $1,$2,$3,FILENAME":"$4,$5,$6}' | \
    sed 's/ /_/g' | \
    sed 's@.*/@@g' | \
    sed "s/${species}_ATAC_//" | \
    sed 's/\.bed$//' | \
    sort -k1,1 -k2,2n | \
    bedtools merge -d -1 -c 4,5,6 -o distinct,max,distinct -i stdin \
    > "$unified_bed"

  # Convert BED → SAF
  awk -F'\t' 'BEGIN{OFS="\t"; print "GeneID","Chr","Start","End","Strand"}
              {print NR, $1, $2, $3, $6}' \
    "$unified_bed" > "$unified_saf"
done
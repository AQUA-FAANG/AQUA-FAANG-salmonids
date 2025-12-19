#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem=1G
#SBATCH --constraint=avx2
#SBATCH --job-name=maf2pafs
#SBATCH --array=1-29%5

# check that there are atleast 2 arguments
if [ $# -lt 2 ]; then
    echo "Usage: $0 <ref_genome> <genome2> [genome3] [genome4] ..."
    exit 1
fi

REFGENOME=$1

CHROM="$(printf "ssa%02d" $SLURM_ARRAY_TASK_ID)"
MAF_FILE="data/maf_ref_$REFGENOME/${CHROM}_AB.maf.gz"
OUTDIR="data/paf/ref_$REFGENOME"
OUTPREFIX=$OUTDIR/${CHROM}_

echo MAF file: $MAF_FILE   Output prefix: $OUTPREFIX    Genomes included: $@

mkdir -p $OUTDIR
zcat $MAF_FILE | python maf2pafs.py --use_reference --remove_subgenome --out_prefix $OUTPREFIX - $@
OUTDIR=data/unified_peaks_liftover

mkdir -p $OUTDIR

cd $OUTDIR
# download unified peaks
wget https://salmobase.org/datafiles/datasets/Aqua-Faang/robust_ATAC_peaks/unified_annotated_peaks/AtlanticSalmon_unified_peaks.bed
wget https://salmobase.org/datafiles/datasets/Aqua-Faang/robust_ATAC_peaks/unified_annotated_peaks/RainbowTrout_unified_peaks.bed
# Make bed files that only contains a peak ID in forth column (so we know what peak that was lifted over)
awk '{print $1, $2, $3, "OmykUniPeak_" $1 ":" $2 "-" $3}' OFS='\t' RainbowTrout_unified_peaks.bed > OmykUniPeak.bed
awk '{print $1, $2, $3, "SsalUniPeak_" $1 ":" $2 "-" $3}' OFS='\t' AtlanticSalmon_unified_peaks.bed > SsalUniPeak.bed
cd -

# splitLiftover.sh <in_bed> <out_bed> <src_genome> <target_genome> 
bash splitLiftover.sh $OUTDIR/SsalUniPeak.bed $OUTDIR/SsalUniPeak_liftover_to_Omyk_A.bed Ssal_A Omyk_A

Rscript get_UniPeak_Ssal_overlapping_Omyk.R
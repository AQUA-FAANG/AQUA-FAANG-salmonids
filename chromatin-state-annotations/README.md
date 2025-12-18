# Chromatin state annotations

*Summary*

Chromatin states were annotated across the Atlantic salmon and rainbow trout genomes using ChromHMM Input signal data were BAM alignments for the ATAC-seq and ChIP-seq (H3K27ac, H3K4me1, H3K4me3, and H3K27me3 marks) data generated as described above. Sampled matched input controls were given to estimate background signals for the ChIP-seq data. Prior to ChromHMM, BAM alignments were converted to BED format and filtered to keep reads only in chromosomes and exclude blacklist regions. The alignments were binarized using the ChromHMM BinarizeBed function with default bin size (200 bp), before running the ChromHMM LearnModel function.

After testing models with various state numbers, with the aim to capture the most biologically relevant chromatin states we proceeded with: two separate 12-state models for Atlantic salmon DevMaps (i.e. embryo stages) and BodyMaps (i.e. tissue samples; brain, gonad, liver, muscle); and for rainbow trout, a 12-state model for DevMaps, and three separate 10-state models for each included tissue (brain, liver, muscle).

Note that each rainbow trout tissue was treated separately, as data quality was more variable between tissues compared to Atlantic salmon. After generating the ChromHMM models, the non-blacklisted BAM files representing ATAC-seq and multiple histone mark ChIP-seq data were binarized using the ChromHMM BinarizeBam function with default bin size (200 bp). The combination of binarized signals from ATAC-seq and ChIP-seq marks were used together with the ChromHMM model for the respective species and condition to annotate regions across the genome representing different chromatin states at 200 bp resolution.

*Scripts*

1.  make_genome_data.R

2.  make_design.R

3.  bamtobed_blacklist_filtered.sh

4.  chromHMM_binarizeBed.sh

5.  chromHMM_learnModel.sh

6.  chromHMM_makeSegmentation.sh

7.  state_annotations_to_tsv.R

8.  annotate_states.sh

*State annotation reports*

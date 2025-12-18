# Robust open chromatin definition and functional annotation

*Summary*

A set of robust open chromatin regions, defined as those that are highly reproducible across biological replicates per condition, was generated using the ATAC-seq narrow peak files from the nfcore-pipelines, post-filtering to remove peaks either not found on chromosome sequences, or completely overlapping blacklist regions. The reproducibility of peaks between all possible pairs of biological replicates per condition was determined using the Irreproducible Discovery Rate (IDR) approach as done by the ENCODE project. IDR measures the reproducibility of the scoring of peaks between replicates to determine an optimal cut-off for peak significance. IDR was run using options --rank='p.value' and --soft-idr-threshold=0.1. This gave a set of significantly reproducible peaks between pairs of replicates. To deal with having three biological replicates per condition, the resulting reproducible peaks from each pair were merged together with bedtools merge (=\>1bp overlap), taking the maximum IDR score for the resulting peak. With one set of reproducible peaks per condition, the summits of those peaks were re-calculated using the MACS2 refinepeak function and all replicate ATAC-seq BAM alignments were merged together.

Finally, the reproducible ATAC-seq peaks with summits were assigned with chromatin states per condition by intersecting ChromHMM annotations using bedtools intersect. When a peak overlapped more than one ChromHMM annotation in a given condition, the assigned chromatin state was chosen in order of importance, based on how confidently we could call the state given the combination of assay signals required. A single unified set of robust open chromatin regions (reproducible ATAC-seq peaks) was generated per species by merging all overlapping peak regions (=\>1bp overlap) across all conditions. Read counts from each ATAC-seq and ChIP-seq assays were counted with featureCounts using the BAM alignments for each biological replicate against the unified peak set of each species. This table of assay counts within robust open chromatin peaks for all samples was used for downstream analyses.

*Scripts*

1.  Collect peaks - Link to all ATAC narrow peak files from nfcore results

-   script: 01_collect_all_peaks.R
-   results: ./ATAC_narrowPeak/

2.  Blacklist peaks - Filter out ATAC peaks not in chromosome sequences and that are completely overlapped by a species specific blacklist region

-   script: 02_blacklist_peaks.sh
-   results: ./blacklisted_peaks/

3.  IDR peaks - IDR testing of ATAC narrow peaks between pairs of biological replicates, using options --rank='p.value', and --soft-idr-threshold=0.1

-   script: 03_run_idr.sh
-   results: ./IDR_peaks/

4.  Merge IDR peaks - Merging of overlapping IDR peaks (\>=1bp overlap) per biological condition, taking max IDR score for the resulting peak

-   script: 04_merge_peaks.sh
-   results: ./merged_peaks/

5.  Refine peak summits - Recalculate summits of merged IDR peaks using together the bam alignments per biological condition

-   script: 05_refine_peak_summits.sh
-   results: ./summits/

6.  Add chromatin states - Annotate merged IDR peaks with the intersecting refined summits, and chromatin states of the same biological condition

-   script: 06_add_chromatin_states.sh
-   results: ./add_chromatin_states/

7.  Finalising annotated peaks - Take the annotated peaks and produce a final set of bed and narrowPeak files with proper column orders, and assign a single best represented chromatin state for the name

-   script: 07_annotate_peaks.R
-   results: ./annotated_peaks/

8.  Unified peaks per species - Merging of overlapping final peak regions into a per species unified file of peaks, taking for merged peaks the best score and the names are listing the biological conditions from which peaks were merged, together with the condition specific chromatin states

-   script: 08_unified_peaks.sh
-   results: ./unified_annotated_peaks/

9.  Count assay expression of unified peaks - Count read coverage from assay BAM files over unified peak regions.

-   script: 09_unified_expression.sh
-   results: ./unified_annotated_peaks/assay_expression/

10. Join and normalise assay expression of unified peaks - Combine peak assay counts into a single dataset and normalise counts to CPM.

-   script: 10_unified_expression_combined.R
-   results: ./unified_annotated_peaks/assay_expression/

# Transcriptome analysis scripts

*Summary*

Raw count files generated with the nf-core RNAseq pipeline were processed, replicates merged, filtered, TPM values scaled and prepared for subsequent analysis. Scaled expression files were clustered using SOM, and visualized through heatmaps, umaps and line plots. Further, clusters were compared for enriched GO terms and analysed for ortholog enrichment across SOM clusters between atlantic salmon and rainbow trout (Fig 2).
Duplicated gene expression divergence was assessed across embryonic stages and tissues for conserved AORE and LORE ohnologs as well as duplicated and singleton orthologs.
*Scripts*

1 - Load count data, merge replicates and normalise expression, perform SOM clustering

-   script: Data_Loading_and_Normalization.Rmd
-   functions: Data_Loading_and_Normalization_Functions.R

2 - Visualise normalized expression data using Heatmaps and UMAPS

-   script: Clustering_and_Visualization_Fig2.Rmd
-   functions: clustering_and_visualization_functions.R
-   results: results/

3 - Calculate and visualise shared orthologs across Atlantic salmon and rainbow trout SOM clusters

-   script: observed_expected.Rmd

4 - Create custom Atlantic salmon annotation package

-   script: makeOrgDb.Rmd

5 - Create custom rainbow trout annotation package

-   script: makeOrbDb_trout.Rmd

6 - Perform GO enrichment test across Atlantic salmon SOM clusters

-   script: GO.Rmd

7 - Perform GO enrichment test across rainbow trout SOM clusters

-   script: GO_trout.Rmd


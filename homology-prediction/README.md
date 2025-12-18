# Homology prediction

Identification of high confidence ohnologs and orthologs in salmon and trout

## Method

We used Ensembl compara gene trees (Compara 106 downloaded from <http://ftp.ensembl.org/pub/misc/aqua-faang/ensembl_trees/gene_trees/Compara.106.protein_default.nhx.emf.gz>) to identify ohnolog (i.e. duplicates from the Ss4R WGD) and ortholog relationship between Atlantic salmon and Rainbow trout genes. However, these relationships can be obscured by incorrect gene-tree estimation, gene-loss or delayed re-diploidisation. To verify that gene duplicates derive from the WGD, i.e. they are ohnologs, we take advantage of the retained syntenic blocks that the vast majority of duplicates reside in. The process starts by extracting the salmonid+pike (Protacanthopterygii) clades from the gene trees. If both duplicates have been retained each clade should contain two Atlantic salmon and two Rainbow trout genes, i.e. 2:2 orthologs. Using the 2:2 ortholog clades we identify the syntenic blocks by selecting the commonly observed combinations of the four chromosomes (\>20 clades with same combination of chromosomes). We assign the ortholog relationship of the four chromosome using the most commonly observed phylogeny of the clades that are in agreement with Ancestral Ortholog Resolution (AORe), i.e. there are two clades, each with one Atlantic salmon and one Rainbow trout ortholog. Since there are regions that have undergone delayed re-diploidization there are clades with Lineage specific Ortholog Resolution (LORe) phylogeny, i.e. both copies in each species are in separate clades. The syntenic relationship of the 2:2 genes is then used to verify the synteny of incomplete clades (i.e. 2:1, 2:0, 1:1).

### Files and scripts

-   `download.sh` - downloads compara gene trees and species tree from Ensembl
-   `selected-speces.txt` - Defines the subset of species to extract from the full compare gene trees
-   `read-compara-trees.R` - Extracts trees with a subset of species and generates ortholog group table (OGtbl)
-   `R/*` - Some helper functions used by `read-compara-trees.R`
-   `Salmon-Trout_quartets_Ens106.Rmd` - Uses the OGtbl and trees to generate a table of salmon and trout orthologs
-   `suppFigure_orthologs.svg` - Manually edited composite figure for the supplementary

### Running the scripts

Make sure you have all the R packages installed and run in this order:

1.  download.sh
2.  read-compara-trees.R
3.  Salmon-Trout_quartets_Ens106.Rmd

### Final output

-   `SalmonTroutOrthologs.tsv` - Table of orthologs
-   `Salmon-Trout_quartets_Ens106.html` - report with some explanations of the columns in the table.

> Note: These files are available at https://salmobase.org/datafiles/datasets/Aqua-Faang/salmon-trout-orthologs/ 

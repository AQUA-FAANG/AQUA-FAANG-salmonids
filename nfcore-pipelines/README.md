## Nf-core pipelines: rnaseq, chipseq, atacseq

*RNA-seq alignment and quantification*

RNA-seq data was processed using a standardised approach through nf-core pipeline rnaseq v3.8.1 with standard parameters. The pipeline involves adapter and quality trimming of reads (Trim Galore!), alignment of reads to the reference genome using STAR before projection of alignments to the transcriptome for BAM-level quantification of gene and transcript level read counts using RSEM. MultiQC reports were produced and assessed for all samples.

*ATAC-seq alignment and peak calling*

ATAC-seq data was processed using a standardised approach through nf-core pipeline atacseq v1.2.1 with standard parameters. The pipeline involves adapter and quality trimming of reads (Trim Galore!), alignment of reads to the reference genome using the Burrows-Wheeler Alignment (BWA) tool and (narrow) peak calling using MACS2. The estimated genome size used was: 2,756,563,003 for Atlantic salmon, and 2,070,645,207 for rainbow trout. MultiQC reports were produced and assessed for all samples.

*ChIP-seq alignment and peak calling*

ChIP-seq data was processed using a standardised approach through the nf-core pipeline chipseq v1.2.2 with standard parameters. The pipeline involves adapter and quality trimming of reads (Trim Galore!), alignment of reads to the reference genome using BWA, and peak calling using MACS2 with representative ChIP input control samples used to establish background signal. The estimated genome size used was as detailed above for ATAC-seq. Broad peaks were called for H3K4me1, H3K27me3 and narrow peaks for H3K4me3, H3K27ac. MultiQC reports were produced and assessed for all samples.

### Reference genomes

The reference genomes and annotations used were downloaded from the Ensembl rapid release database; for Atlantic salmon: assembly Ssal_v3.1, softmasked fasta file (Salmo_salar-GCA_905237065.2-softmasked.fa) and gene annotation GTF file (Salmo_salar-GCA_905237065.2-2021_07-genes.gtf.gz). For rainbow trout: assembly USDA_OmykA_1.1, softmasked fasta file (Oncorhynchus_mykiss-GCA_013265735.3-softmasked.fa) and gene annotation GTF file (Oncorhynchus_mykiss-GCA_013265735.3-2020_12-genes.gtf).

Genomes and annotations were saved in this relative folder structure:

-   `AtlanticSalmon/genome`
    -   Genome: `Salmo_salar-GCA_905237065.2-softmasked.fa`
    -   Genes: `Salmo_salar-GCA_905237065.2-softmasked_genes.gtf`\
-   `RainbowTrout/genome`
    -   Genome: `Oncorhynchus_mykiss-GCA_013265735.3-softmasked.fa`
    -   Genes: `Oncorhynchus_mykiss-GCA_013265735.3-2020_12-genes.gtf`

### Scripts

All script files were saved under two relative folder structures for Atlantic salmon and rainbowtrout

-   `AtlanticSalmon/scripts`
-   `RainbowTrout/scripts`

*Scripts to create sample designs*

These R scripts, under both Atlantic salmon and rainbow trout directories, created all nf-core sample design tables for all pipelines, separately for BodyMap and DevMap samples.

-   BodyMap: `make_designs_bodymap.R`
-   DevMap: `make_designs_devmap.R`

*Nfcore scripts*

These bash scripts, under both Atlantic salmon and rainbow trout directories, ran all nextflow nf-core pipelines (nextflow verison 22.04.5.5708), including fetching raw fastq data from ENA, then processing RNA-seq, ATAC-seq, and ChIP-seq data with the respective pipeline. The config file was specific to the local HPC (Orion).

-   Download sequence data from ENA: `fetchngs.sh`
-   RNA-seq pipeline: `rnaseq.sh`
-   ATAC-seq pipeline: `atacseq.sh`
-   ChIP-seq pipeline: `chipseq.sh`
-   Config file: `orion.config`

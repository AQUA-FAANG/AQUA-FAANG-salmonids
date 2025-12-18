# AQUA-FAANG-salmonids

This repository contains the analysis scripts used for the manuscript titled "Salmonids reveal principles of regulatory evolution following autotetraploidization"

The study presents large-scale regulatory omics datasets (RNA-seq, ATAC-seq, ChIP-seq) in Atlantic salmon and rainbow trout, spanning embryonic development (DevMap) and adult tissues (BodyMap), to investigate the evolution of gene expression and regulatory elements following the salmonid-specific Whole Genome Duplication (4R autotetraploidization).

## Repository Structure

The analysis scripts are organized into sub-directories, with each major section of the manuscript having its own folder. For clarity and reproducibility, scripts related to the main figures (Fig. 1 to Fig. 6) are placed within their corresponding analysis directories.

```
.
├── README.md
├── blacklist-regions/         # Generation of genome-specific blacklist tracks to mask artefacts
├── chromatin-state-annotations/ # ChromHMM scripts and reports for genome-wide chromatin state calls
├── homology-prediction/       # Scripts for identifying ohnologs and orthologs
└── nfcore-pipelines/          # Pipeline configurations and resources for RNA-seq, ATAC-seq, and ChIP-seq processing
```

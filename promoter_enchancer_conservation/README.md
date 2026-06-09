# Promoter and Enhancer Conservation Across the Atlantic Salmon Genome

*This repository contains the code used to identify and analyse the conservation of regulatory elements across the Atlantic salmon (Salmo salar) genome. The workflow combines whole-genome alignments, ATAC-seq peak annotations, active promoter and enhancer datasets, and rediploidization histories to investigate how regulatory element conservation varies across ontogeny and between early- and late-rediploidized genomic regions.*


## Workflow
```

ATAC peaks
    ↓
Split peaks by chromosome
     ↓
Extract MAF alignments
     ↓
Overlay peaks onto alignments
     ↓
Merge chromosome-level results
     ↓
Curate alignments
     ↓
Generate conservation datasets
     ↓
Overlay promoter and active enhancer annotations
     ↓
Analyse conservation across:
    • ontogeny
    • rediploidization history
    • conservation categories
     ↓
Generate figures and summary statistics

Generate circos visualization of rediploidization regions
```
## Repository Structure

### 01_prepare_peak_subset.sh  
*Prepares chromosome-specific subsets of ATAC-seq peaks.*

### 02_extract_peak_mafs.sh  
*Extracts MAF alignment blocks overlapping each chromosome-specific peak set.*

### 03_overlay_maf_with_atac.R  
*Projects ATAC peaks onto MAF alignments and identifies conserved peak regions.*

### 04_global_ATAC_conservation_analysis.R  
*Performs genome-wide conservation analyses and classifies regulatory elements into conservation categories.*

### 05_promoter_active_enhancer_dynamics_across_ontogeny.R  
*Integrates conservation datasets with active promoter and enhancer annotations.*

### 06_circos_figure_4B.R  
*Generates the circos plot shown in Figure 4B.*  
*This script is independent from the main conservation pipeline and can be run separately.*




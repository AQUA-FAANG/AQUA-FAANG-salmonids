# Gimme Maelstrom Enhancer CNE Motif Enrichment Pipeline
*This repository contains a complete workflow for Enchancer CNE (Conserved Non-coding Element) motif enrichment analysis using GimmeMotifs Maelstrom, including preprocessing, motif enrichment, visualization, and annotation of transcription factor binding motifs.*

## Overview
The workflow consists of five modular steps:
1. CNE preprocessing and region construction
2. Maelstrom motif enrichment analysis (Age and Age × Conservation)
3. Post-processing and motif ranking
4. Visualization of enriched motif families
5. Automated motif annotation (TF mapping + logos)

```
01_prepare_CNE_regions_for_maelstrom.R
02_run_gimme_maelstrom_age.sh
03_run_gimme_maelstrom_age_conservation.sh
04_plot_maelstrom_results.R
05_add_motif_logos.py
```
## Workflow description

### 1. CNE preprocessing (01_prepare_CNE_regions_for_maelstrom.R)
**Purpose:**  
Prepares conserved non-coding elements (CNEs) for motif analysis by:  
*Merging nearby CNEs within a defined distance threshold*  
*Ensuring minimum region size*  
*Preventing biologically invalid merges across annotation classes*  
*Generating Maelstrom-compatible region files*  

### 2. Maelstrom motif enrichment – Age analysis Script (02_run_gimme_maelstrom_age.sh)  
**Purpose:**  
*Runs GimmeMotifs Maelstrom to identify motif enrichment across evolutionary age categories.*  

### 3. Maelstrom motif enrichment – Age × Conservation Script (03_run_gimme_maelstrom_age_conservation.sh)  
**Purpose:**  
*Identifies motif enrichment patterns across combined evolutionary age and conservation context, enabling detection of interaction effects.*  

### 4. Motif visualization and enrichment analysis Script(04_plot_maelstrom_results.R)  
**Purpose:**  
Processes Maelstrom outputs to generate:  
*Top enriched motifs per category*  
*Non-redundant motif families*  
*Heatmaps of motif activity (z-scores)*  
*Comparative plots across:*  
*Evolutionary age*  
*Age × conservation interactions*  

### 5. Motif annotation and logo generation Script (05_add_motif_logos.py)  
**Purpose:**  
Automatically annotates motif enrichment results with:  
*Transcription factor (TF) names*  
*Motif family classification*  
*DNA binding logos from GimmeMotifs database*  

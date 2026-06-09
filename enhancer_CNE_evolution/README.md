# Evolution of Enhancer-CNEs Following the Salmonid Whole Genome Duplication
*Workflow used to generate the conserved non-coding element (CNE) and enhancer integration analyses presented in Figure 6*

## Overview
*The objective of this workflow is to investigate how enhancer CNE activity evolves following the salmonid-specific whole genome duplication (Ss4R). We quantify how phylogenetic age, enhancer conservation state, developmental stage, and rediploidization history interact with enchancer CNE elements.*

## Workflow
### Script : Prepare CNEs for MAF extraction
```
01_get_CNE_maf_extraction_list.R
```
**Purpose:**  
*Starting from the full set of pike-derived CNE coordinates, applies multiple quality-control filters. Removes problematic alignments. Generates a curated set of high-confidence CNEs suitable for downstream comparative analyses.*  
  
### Scripts: Extract MAF alignments for curated CNEs  
```
02_submit_MAF_CNE_extraction_pipeline.sh  
02.1_download_cat_maf.sh  
02.2_extract_cne_maf.sh  
```
**Purpose:**   
*These scripts automate extraction of whole-genome multiple alignment (MAF) blocks corresponding to curated CNE coordinates.*

**Workflow:**  

02.1_download_cat_maf.sh  
*Downloads or prepares genome-wide MAF alignments.*  
*Concatenates chromosome-specific alignment files when required.*  

02.2_extract_cne_maf.sh  
*Extracts alignment blocks overlapping curated CNE coordinates.*  
*Generates individual CNE alignment files.*  

02_submit_MAF_CNE_extraction_pipeline.sh  
*HPC submission wrapper.*  
*Launches large-scale parallel extraction jobs.*  
  
### Step 3 – Compute CNE conservation statistics  
### Script  
```
03_compute_CNE_stats_summary.R
```
**Purpose:**  
Processes extracted MAF alignments and generates summary stats etc  

### Step 4 – Integrate active enhancer data  
### Script  
```
04_enhancer_CNE_overlap_ontogeny.R
```
**Purpose**  
Generate figures in Fig6 from the manuscript.  
*Distribution of enhancer-associated CNEs across:*  
*Phylogenetic age categories, Shared, Alignable, Exclusive conservation classes*  
*including Fisher's Exact Test enrichment analyses.*  
  
*Normalised proportions of enhancer-associated CNEs across:*
*Embryonic developmental stages, Adult tissues, AORe regions, LORe regions, Conservation categories, Phylogenetic age classes*  

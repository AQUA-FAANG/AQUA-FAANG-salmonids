# Genome alignment with cactus by syntenic blocks

This directory contains the scripts used to generate salmonid whole genome alignments. We found that feeding the entire genome directly into the cactus genome aligner performed poorly for the duplicated genomes. We therefore apply a strategy for improving the alignments by splitting genomes into __blocks__ of duplicate syntenic regions. The genomes of salmon, trout, and pike (Ssal, Omyk, Eluc) are aligned, with the duplicated salmonid genomes split into two __sub-genomes__ (suffixed _A and _B).

The __blocks__ are manually curated and can be found in the table `Salmonid_Synteny_for_alignments_2023.04.12.xlsx`. The Atlantic salmon chromosomes are used as reference, i.e. there is one block for each salmon chromosome (this is Ssal_A in the table and is used as block ID). For each block the coordinates for the corresponding syntenic regions in the other (sub-)genomes are specified (typically multiple regions from multiple chromosomes). Note that every part of each genome is represented twice in the table since a duplicated genome is used as reference.

The coordinates from the alignments output from cactus are relative to the extracted syntenic regions, a.k.a. __block coordinates__, and must be converted back to original coordinates.

## Overview

This repository contains a set of scripts to:

1. Main workflow for generating the alignments. (Extract blocks, run cactus, HAL to MAF, coordinate conversion)
2. Tool to run halLiftover on the 
3. Tool to convert MAF to PAF files.
4. Tool for making alignmentDepth bigwig
5. Scripts for mapping unified peaks between salmon and trout using liftover


## Prerequisites

- R
- samtools
- cactus (includes halLiftover and hal2maf)

TODO: Add exact versions: R, samtools, cactus (and Singularity image tag), bedGraphToBigWig/bigWigMerge, Python deps (BioPython, pandas, etc.).

## Setup

For now the scripts are hardcoded to use my local setup, including the use of slurm jobs and environment modules, so some modification is needed to make them run elsewhere. Cactus was run using singularity with image pulled from the docker image at https://github.com/ComparativeGenomicsToolkit/cactus/releases

TODO: Make the code more reproducible

Also the input genome fasta are hardcoded to local paths in `extractBlocks.R`. 
TODO: make a download script that downloads genomes from Ensembl and modify extractBlocks.R to point to these.

## Workflow

### Sequence Extraction

* `extractBlocks.R`: Reads the synteny table `Salmonid_Synteny_for_alignments_2023.04.12.xlsx` and generates a directory structure (one for each salmon chromosome) and the commands (samtools faidx) needed to extract the corresponding blocks for each (sub-)genome.
* `run_extractBlocks.job.sh`: Executes the commands generated in the previous step. (This was necessary because samtools wouldn't run directly from R).

### Cactus Alignment
* `run_cactus.job.sh`: Slurm array that runs Cactus for each of the directories (one for each salmon chromosome).

###  Coordinate Conversion

* `run_hal2maf.job.sh`: Slurm array that runs hal2maf ($1 is reference genome "Ssal_A" by default).
* `run_fix_maf.job.sh`: Slurm array that runs `fix_maf.py`, only keeping the A/B suffix.
* `fix_maf.py`: Converts coordinates from block coordinates to full genome coordinates (optionally keeps the _A/_B suffix).
* `makeChromSizes.R`: Generates a table of chromosome sizes needed for `fix_maf.py`.

## Tools

### splitLiftover.sh - halLiftover with coordinate conversion

Usage: splitLiftover.sh <in_bed> <out_bed> <src_genome> <target_genome> 

This script lets you to perform liftover with synteny block alignments by converting coordinates before and after running halLiftover. It is performing the following steps:

1. Convert coordinates in BED file to block coordinates and split with `split_bed_to_blocks.R`.
2. Run halLiftover for each block.
3. Convert coordinates in BED/PSL files from block coordinates back to sequence coordinates using `fix_bed.py` or `fix_psl.py`.

### MAF to PAF

* `maf2pafs.py`: Converts MAF files to pairwise .PAF files that can be viewed in JBrowse.
* `run_maf2pafs.py.job.sh`: Slurm array that runs `maf2pafs.py`
* `run_merge_pafs.py.job.sh`: Merge and sort pafs to get whole genome in one file

### Alignment depth to bigwig

Several scripts that allows to create bigwig files with alignment depth from the split hal files.

* `make_alignmentDepth_bigwig.sh` - runs halAlignmentDepth on each block, converts block coordinates with fix_wiggle.py and converts to bigwig
* `merge_bigwigs.sh` - uses bigWigMerge and bedGraphToBigWig to merge and calls `fixMergedBigWig.R` to fill the gaps
* `fixMergedBigWig.R` - Fills all gaps with 0's (Since bigWigMerge removes all region with 0's)
* `run_all_alignmentdepths.sh` - runs make_alignmentDepth_bigwig.sh for multiple combinations of genomes

## Additional scripts

### Mapping unified peaks between salmon and pike

* `unified_peaks_liftover.sh` - Script that downloads peaks and runs liftover
* `get_UniPeak_Ssal_overlapping_Omyk.R` - Script that generates a mapping table from liftover results

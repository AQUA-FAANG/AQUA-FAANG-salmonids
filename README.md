# AQUA-FAANG salmonids — regulatory evolution after autotetraploidization

Analysis code for the manuscript **"Salmonids reveal principles of regulatory
evolution following autotetraploidization"** (Baudement, Perojil Morata, Gillard
*et al.*; corresponding authors Daniel J. Macqueen and Sigbjørn Lien).

## What this study is

Salmonids underwent a salmonid-specific whole-genome duplication (the **4R**
autotetraploidization, ~80–100 Mya). Because that duplication is comparatively
recent and rediploidized *asynchronously*, salmonid genomes preserve a graded
record of how duplicated genes (**ohnologs**) and their regulatory elements
diverge over time. This project uses large-scale multiomics to read that record.

We generated species-matched **RNA-seq, ATAC-seq and ChIP-seq** (histone marks
H3K27ac, H3K4me1, H3K4me3, H3K27me3) for two salmonids — **Atlantic salmon**
(*Salmo salar*) and **rainbow trout** (*Oncorhynchus mykiss*) — across two
sample panels:

- **DevMap** — a developmental time course spanning key stages of embryogenesis.
- **BodyMap** — an adult tissue panel (brain, gonad, liver, muscle, and more),
  for both sexes and at two maturity stages.

From these data the study builds a functional genome annotation (ChromHMM
chromatin states, a robust open-chromatin atlas, active promoters and
enhancers), a salmon–trout whole-genome alignment, and a high-confidence
ohnolog/ortholog map, then asks how expression and regulatory-element activity
evolve following the 4R WGD as a function of **ontogeny**, **tissue context**,
and **rediploidization history** (Early vs. Late). The headline result is a
period of *maximal regulatory constraint at advanced stages of embryogenesis*,
mirroring constraints seen across vertebrate species.

Key resource numbers (Atlantic salmon / rainbow trout):

| Resource | Salmon | Trout |
|---|---|---|
| Robust open-chromatin regions (IDR) | 319,718 | 194,002 |
| Active promoters (ChromHMM) | 45,802 | 28,138 |
| Active enhancers (ChromHMM) | 78,424 | 42,164 |

High-confidence salmon–trout homology: **10,521 ohnolog pairs** retained from
the 4R WGD plus **5,832 singletons**; the comparative atlas covers 26,157 genes
and 194,528 unique promoters/enhancers active in at least one sample across both
species.

## Data availability

The repository contains **code and small derived/figure-ready assets only** —
not the raw sequencing data or the large intermediate files. Inputs come from:

- **Raw sequencing data** — the FAANG data portal / ENA (771 sequencing
  datasets; see manuscript and `nfcore-pipelines/`).
- **Processed datasets and derived tables** — [salmobase.org](https://salmobase.org)
  under `datafiles/datasets/Aqua-Faang/` (nf-core results per species/map/assay,
  the salmon–trout ortholog table, the whole-genome Cactus alignments, etc.).
  Several scripts download directly from these endpoints.
- **Reference genomes** — Ensembl rapid release:
  Atlantic salmon `Ssal_v3.1` (GCA_905237065.2), rainbow trout `USDA_OmykA_1.1`
  (GCA_013265735.3), northern pike `fEsoLuc1` (GCA_011004845.1) as the
  non-duplicated outgroup.

## Repository structure

Each subdirectory corresponds to one stage of the analysis and has its **own
README** with the detailed method and script-by-script description. They are
listed below roughly in pipeline order; the "Figure" column shows the main
manuscript figure each one feeds.

| Directory | What it does | Related figure(s) |
|---|---|---|
| [`nfcore-pipelines/`](nfcore-pipelines/) | Raw RNA-seq / ATAC-seq / ChIP-seq processing via nf-core (alignment, quantification, peak calling) and sample-design generation. | Upstream of all; design schematic Fig. 1a,b |
| [`blacklist-regions/`](blacklist-regions/) | Build genome-specific blacklist tracks (Umap mappability + Blacklist) to mask signal-artefact regions. | Upstream masking (not a plotted panel) |
| [`chromatin-state-annotations/`](chromatin-state-annotations/) | ChromHMM genome-wide chromatin-state models and per-condition segmentations. | Fig. 1c,d (states/colours also in 1e,1g); trout = Ext. Data Fig. 1 |
| [`robust-open-chromatin/`](robust-open-chromatin/) | Define reproducible ATAC peaks (IDR), refine summits, assign chromatin states, build the unified robust open-chromatin atlas and assay-count tables. | Fig. 1e,f (feeds 1g); peak/promoter/enhancer sets used throughout |
| [`homology-prediction/`](homology-prediction/) | Identify high-confidence 4R ohnologs and salmon–trout orthologs from Ensembl Compara trees + synteny. | Foundational — example in Fig. 1g; underpins Fig. 2c, 3, 5, 6 |
| [`genome-alignment/`](genome-alignment/) | Whole-genome alignment of salmon, trout and pike with Cactus, split by syntenic blocks; coordinate conversion, liftover and PAF/bigWig tooling. | Foundational — Fig. 1g (track), 2d, 4, 5, 6 |
| [`TFBS/`](TFBS/) | Validate active promoters/enhancers by TF-binding-site enrichment (GimmeMotifs maelstrom, JASPAR). | Supplementary (TFBS validation; Supp. Figs. 10–12) |
| [`transcriptome-analysis/`](transcriptome-analysis/) | RNA-seq clustering (SOMs/UMAP), GO enrichment, ortholog/ohnolog SOM co-clustering. | Fig. 2a,c; Ext. Data Figs. 4a,c, 5, 6; GO = Suppl. |
| [`chromatin-accesibility-clustering-and-visualization/`](chromatin-accesibility-clustering-and-visualization/) | ATAC-seq clustering (SOMs/UMAP) and the chromatin panels of Figure 2. | Fig. 2b,d; Ext. Data Fig. 4b,d |
| [`promoter_enchancer_conservation/`](promoter_enchancer_conservation/) | Promoter/enhancer conservation categories (Shared / Alignable / Exclusive) across ontogeny and rediploidization; circos and TE-overlap panels; conservation-vs-JSD correlations. | Fig. 4 (4B circos, 4C/4D proportions + TE/JSD-correlation panels) |
| [`enhancer_CNE_evolution/`](enhancer_CNE_evolution/) | Evolutionary fate of enhancer-associated conserved non-coding elements (CNEs) by phylogenetic age and conservation class. | Fig. 6a,b |
| [`enhancer_CNE_motif_analysis/`](enhancer_CNE_motif_analysis/) | TF-motif enrichment in enhancer-CNEs across age and conservation (GimmeMotifs maelstrom) with motif logos. | Ext. Data Fig. 8 (supports Fig. 6) |

> Note on figure attribution: Figures are assembled in Illustrator from several
> sources, so no single script produces a whole figure.

## How to run

The code falls into two tiers with very different reproducibility profiles.

### Tier 1 — turnkey figure reproduction (runs on a laptop)

Two directories ship the deposited, manuscript-exact derived assets and
re-render their figures with **no HPC and almost no download**:

- [`transcriptome-analysis/`](transcriptome-analysis/) → Figure 2A
- [`chromatin-accesibility-clustering-and-visualization/`](chromatin-accesibility-clustering-and-visualization/) → Figure 2B

These are R-markdown notebooks. They auto-install their R package stack on first
run and default to a `reproduce` mode that reads the small gzipped TSV / RDS
assets tracked in this repo, so the published panels reproduce deterministically.
For example:

```r
# transcriptome-analysis/ — open in RStudio and run 00..06 in order, same session
rmarkdown::render("00_setup.Rmd")   # selects reproduce mode, builds paths
rmarkdown::render("01_data_loading.Rmd")
# ... 02..06

# chromatin-accesibility-clustering-and-visualization/ — knit in order
rmarkdown::render("00_setup.Rmd"); rmarkdown::render("01_atac_loading.Rmd")
rmarkdown::render("02_atac_soms.Rmd"); rmarkdown::render("03_fig2B_atac_soms.Rmd")
```

Each of those directories' READMEs documents the `reproduce` vs `recompute`
modes, a `TEST_MODE` smoke-test switch, expected run times, and the
SOM/BLAS reproducibility caveats. Start there.

### Tier 2 — upstream and comparative workflow scripts (HPC)

Everything else (`nfcore-pipelines`, `blacklist-regions`,
`chromatin-state-annotations`, `robust-open-chromatin`, `genome-alignment`,
`TFBS`, `homology-prediction`, `promoter_enchancer_conservation`,
`enhancer_CNE_evolution`, `enhancer_CNE_motif_analysis`) is the production
pipeline as it was actually run on an HPC cluster. These scripts:

- are **numbered** (`01_…`, `02_…`, …) and meant to run **in order** within each
  directory — follow the per-directory README;
- contain **hardcoded local paths** (e.g. `/mnt/project/Aqua-Faang/…`), **SLURM**
  directives (`#SBATCH`), and environment-module / Singularity calls specific to
  the original cluster;
- expect the large raw and intermediate files that are **not** shipped here.

To re-run them elsewhere you will need to **edit the hardcoded paths**, adapt the
job-submission and module-loading lines to your environment, and stage the
public input data (raw reads from FAANG/ENA, processed inputs from Salmobase,
reference genomes from Ensembl). Treat them as a precise, runnable record of the
methods rather than a one-command pipeline. The
[`genome-alignment/`](genome-alignment/) README in particular flags remaining
hardcoded inputs and reproducibility TODOs.

## Software requirements

R-notebook directories install their own R/Bioconductor packages on first run
(see each `00_setup`/`R/setup.R`). The upstream/comparative scripts assume the
following external tools were available on the cluster (versions as reported in
the manuscript Methods):

- **Pipelines:** Nextflow 22.04.5, nf-core `rnaseq` v3.8.1, `atacseq` v1.2.1,
  `chipseq` v1.2.2 (bundling Trim Galore!, STAR, RSEM, BWA, MACS2, MultiQC).
- **Peak / signal processing:** ChromHMM, IDR, bedtools, MACS2 `refinepeak`,
  featureCounts (Subread), Umap, Blacklist.
- **Alignment / comparative genomics:** Cactus (progressive, via Singularity),
  HAL tools (`hal2maf`, `halLiftover`, `halAlignmentDepth`), mafTools
  (`mafStats`, `mafExtractor`), minimap2 v2.18, samtools,
  bedGraphToBigWig / bigWigMerge, BioPython, pandas.
- **Motifs / TFBS:** GimmeMotifs v0.17.1 (`gimme maelstrom`) with the
  JASPAR2022_vertebrates and `gimme.vertebrate.v5.0` motif databases.
- **Homology:** Ensembl Compara 106 gene trees; R tree-handling helpers
  (`homology-prediction/R/`).
- **Validation / CNEs:** Bowtie2 v2.5.0, SAMtools v1.5, CAGEr v2.11.4 (CAGE
  validation); GERP++, Multiz, RepeatMasker (CNE / TE analyses).
- **Key R packages** (across the analysis directories): tidyverse, data.table,
  kohonen, uwot, ComplexHeatmap, preprocessCore, phylentropy, clusterProfiler,
  AnnotationForge, biomaRt, circlize, ggalluvial, UpSetR, mclust, ggpubr.


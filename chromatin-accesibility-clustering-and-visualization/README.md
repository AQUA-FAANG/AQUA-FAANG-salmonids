# Chromatin accessibility clustering and visualization

R notebooks for the ATAC-seq clustering analyses in Atlantic salmon (*Salmo
salar*) and rainbow trout (*Oncorhynchus mykiss*), across the DevMap (embryonic
development) and BodyMap (adult tissue) contexts. Starting from the deposited
per-species consensus-expression tables and per-context SOM assets, the pipeline
builds scaled fold-enrichment matrices, assembles the 5×5 hexagonal
self-organising maps, and reproduces the ATAC panels of Figure 2 (Figure 2B) of
the manuscript.

## Notebooks

| File | Purpose |
|---|---|
| `00_setup.Rmd` | Packages, on-disk paths, presence checks, helper sourcing. |
| `01_atac_loading.Rmd` | Read the consensus-expression tables and the deposited per-context peak subsets, write scaled fold-enrichment matrices to `data/derived/`. |
| `02_atac_soms.Rmd` | Load the deposited SOMs, super-clusters and UMAP embeddings; build the ribbon-plot summaries. |
| `03_fig2B_atac_soms.Rmd` | Heatmaps, UMAPs and the Figure 2B ATAC composite. |
| `99_helpers.R` | Deposited-asset reader, scaled-FE matrix builder, SOM summary, palettes and plot builders. |

## How to run

Knit the notebooks in order, in the same R session:

```r
rmarkdown::render("00_setup.Rmd")
rmarkdown::render("01_atac_loading.Rmd")
rmarkdown::render("02_atac_soms.Rmd")
rmarkdown::render("03_fig2B_atac_soms.Rmd")
```

The analysis runs **fully offline** from the deposited assets; nothing is
downloaded. Every expensive chunk is cached under `cache/`.

## Data

All inputs are the deposited manuscript assets shipped in this repository (the
ATAC counterpart of the transcriptome `data/norm/` + `data/soms/`), so knitting
reproduces the published Figure 2B exactly — the SOMs are **loaded**, not re-fit:

- `data/atac/{AtlanticSalmon,RainbowTrout}_consensus_expression.rds` — the exact
  manuscript ATAC values (the unified robust-peak set × per-condition ATAC
  fold-enrichment, the SOM input).
- `data/atac_soms/{devmap,bodymap}_{salmon,trout}_{idx,unit_classif,superclasses,umap}.tsv.gz`
  — per context: the manuscript peak subset (`idx`), the SOM unit assignment per
  peak (`unit_classif`), the PAM(`k = 4`) super-clusters, and the UMAP
  coordinates.

`01` loads the deposited `idx` to subset each context; `02` loads the deposited
SOMs, super-clusters and UMAPs. `data/derived/`, `results/` and `cache/` are
git-ignored; the deposited `data/atac/*_consensus_expression.rds` tables and the
`data/atac_soms/` assets are tracked.

## Reproducibility

The published SOMs, super-clusters and UMAP embeddings are deposited directly, so
the figures reproduce bit-for-bit without re-fitting (SOM/UMAP fitting is not
bit-stable across BLAS/OS builds). Every expensive chunk uses `cache = TRUE`;
delete the `cache/` directory to force a clean re-run.

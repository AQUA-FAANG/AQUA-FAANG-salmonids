# Chromatin accessibility clustering and visualization

R notebooks for the ATAC-seq clustering analyses in Atlantic salmon (*Salmo
salar*) and rainbow trout (*Oncorhynchus mykiss*), across the DevMap
(embryonic development) and BodyMap (adult tissue) contexts. Starting from
per-species consensus expression tables and per-context IDR consensus BED
files, the pipeline builds scaled fold-enrichment matrices, fits 5×5
hexagonal self-organising maps, derives PAM super-clusters and UMAP
embeddings, and assembles the ATAC composite that feeds Figure 2 of the
manuscript.

## Notebooks

| File | Purpose |
|---|---|
| `00_setup.Rmd` | Packages, on-disk paths, helper sourcing. |
| `01_atac_loading.Rmd` | Read consensus expression RDS and IDR BEDs, build context unions, write scaled fold-enrichment matrices to `data/derived/`. |
| `02_atac_soms.Rmd` | Fit four SOMs (DevMap and BodyMap × salmon and trout), derive PAM super-clusters and UMAP embeddings. |
| `03_figure2_atac.Rmd` | Heatmaps, UMAPs and the Figure 2 ATAC composite. |
| `99_helpers.R` | I/O wrappers, SOM wrappers, palettes and plot builders. |

## How to run

Knit the notebooks in order:

```r
rmarkdown::render("00_setup.Rmd")
rmarkdown::render("01_atac_loading.Rmd")
rmarkdown::render("02_atac_soms.Rmd")
rmarkdown::render("03_figure2_atac.Rmd")
```

The IDR BED imports in `01_atac_loading.Rmd` and the four SOM fits in
`02_atac_soms.Rmd` are the dominant chunks; every expensive chunk is cached.

## Data

Two modes are supported, serving two different reproducibility goals.

### `paths$mode = "legacy"` — identical figure reproduction

Reads the exact derived assets used to produce the manuscript figures:

- `data/atac/AtlanticSalmon_consensus_expression.RDS`
- `data/atac/RainbowTrout_consensus_expression.RDS`
- `data/atac/{bodymap,devmap}_idr/{salmon,trout}/`

Knitting in this mode reproduces every ATAC panel from the fixed analysis
state used in the manuscript.

### `paths$mode = "salmobase"` — independent biological replication

The Salmobase ATAC trees publish only the upstream nf-core alignment
outputs:

- <https://salmobase.org/datafiles/datasets/Aqua-Faang/nfcore/AtlanticSalmon/DevMap/ATAC/results/>
- <https://salmobase.org/datafiles/datasets/Aqua-Faang/nfcore/AtlanticSalmon/BodyMap/ATAC/>
- <https://salmobase.org/datafiles/datasets/Aqua-Faang/nfcore/RainbowTrout/DevMap/ATAC/results/>
- <https://salmobase.org/datafiles/datasets/Aqua-Faang/nfcore/RainbowTrout/BodyMap/ATAC/>

The IDR consensus BEDs and the consensus-expression RDS files are derived
in the upstream `robust-open-chromatin/` pipeline (see the sibling
directory in this repository). A user redoing the analysis from raw
alignments runs that pipeline first to regenerate the derived assets, then
points `paths$legacy` at them. In this path the SOM fits are redone and
unit numbering will differ from the published figures, but the qualitative
results should reproduce.

## Reproducibility

Master seed `12345`, applied inside `fit_atac_som` for `somInit` and
`kohonen::som` and inside `02_atac_soms.Rmd` for every `uwot::umap` call.
Every expensive chunk uses `cache = TRUE`; delete the per-notebook
`*_cache/` directory to force a clean re-run.

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

The per-context IDR-consensus peak sets are fetched on the fly from
Salmobase as bigBed (`.bb`) files; URLs are baked into
`00_setup.Rmd::paths$idr_*_urls`. The downloads are cached under
`cache/atac_peaks/` (git-ignored) so re-knits do not re-hit the
network. Canonical file lists:

- <https://salmobase.org/datafiles/datasets/Aqua-Faang/robust_ATAC_peaks/salmon_files.txt>
- <https://salmobase.org/datafiles/datasets/Aqua-Faang/robust_ATAC_peaks/trout_files.txt>

The per-species consensus-expression RDS files are **not** redistributed
in this repository — they are several hundred MB each and the lab is
still discussing where to host them (Zenodo DOI vs regeneration from
Salmobase BAMs via `robust-open-chromatin/09_unified_expression.sh` +
`10_unified_expression_combined.R`). In the meantime, the manuscript
state of these tables lives on OneDrive as `Ssal_Figure_2.RData` and
`Omyk_Figure_2.RData`. To populate `data/atac/` for a local run:

```r
options(
  aquafaang.ssal_rdata = "<path>/Ssal_Figure_2.RData",
  aquafaang.omyk_rdata = "<path>/Omyk_Figure_2.RData",
  aquafaang.atac_out   = "data/atac"
)
source("scripts/extract_legacy_consensus_expression_rdata.R")
```

`data/atac/`, `data/derived/`, `results/`, and `cache/` are git-ignored.

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

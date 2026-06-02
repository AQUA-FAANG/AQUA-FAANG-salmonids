# Transcriptome analysis

R notebooks for the bulk RNA-seq analyses in Atlantic salmon (*Salmo salar*)
and rainbow trout (*Oncorhynchus mykiss*), covering a 14-stage developmental
time course (DevMap) and an adult tissue panel (BodyMap, eight tissues across
sex and sexual maturity). The pipeline merges technical replicates,
quantile-normalises and z-scales TPM matrices, fits self-organising maps
(4×4 hexagonal for DevMap, 6×6 for BodyMap, both species), and analyses the
resulting clusters for GO term enrichment, salmon–trout ortholog
co-clustering against a permutation null, and salmon ohnolog co-clustering
stratified by rediploidization timing.

## Notebooks

| Notebook | Purpose |
|---|---|
| `00_setup.Rmd` | Entry point: choose the data source, install/load packages, build `paths`, print the configuration. Run first. |
| `01_data_loading.Rmd` | Load (reproduce mode) or re-derive (recompute mode) the four normalised expression matrices. |
| `02_soms.Rmd` | Load the deposited SOMs (reproduce) or re-fit them (recompute); build per-cluster annotation tables. |
| `03_fig2A_expression_soms.Rmd` | The **transcriptome** panels of Figure 2 (= Fig 2A: expression-SOM heatmaps + UMAPs) plus the supplementary full SOMs. Each panel is saved to `results/` as a PDF. |
| `04_ohnolog_coclustering.Rmd` | Per-cluster composition and per-pair co-clustering of AORe and LORe salmon ohnologs across the DevMap and BodyMap SOMs (supplementary figures). |
| `05_observed_expected.Rmd` | Observed vs expected ortholog co-clustering against a 10,000-permutation null — this produces the **Fig 2A observed/expected ratio bands**. RNA and ATAC Sankeys. |
| `06_go.Rmd` | GO biological-process enrichment via `compareCluster` for the four SOMs. |
| `99_helpers.R` | I/O, deposited-asset loaders, normalisation, SOM fit/annotation, palettes, heatmap builders. |
| `99a_make_orgdb_salmon.Rmd` | Build and install `org.Ssalar.eg.db` from Ensembl BioMart. |
| `99b_make_orgdb_trout.Rmd` | Build and install `org.Omykiss.eg.db` from Ensembl BioMart. |
| `R/setup.R` | Sourced by every notebook: packages, constants, the `paths` list, mode handling. |
| `R/url_or_path.R` | Loader dispatch that accepts a local path or a Salmobase URL, caching downloads under `cache/`. |
| `scripts/` | R scripts sourced by `04_ohnolog_coclustering.Rmd` for the ohnolog analysis and its figure panels. |

### Where Figure 2 comes from

The manuscript's Figure 2 is assembled in Illustrator from several sources —
no single notebook produces the whole figure:

| Figure 2 piece | Produced by |
|---|---|
| 2A heatmaps + UMAPs (RNA expression SOMs) | `03_fig2A_expression_soms.Rmd` |
| 2A observed/expected ratio bands (RNA) | `05_observed_expected.Rmd` |
| 2B heatmaps + UMAPs (ATAC peak intensity) | `../chromatin-accesibility-clustering-and-visualization/03_fig2B_atac_soms.Rmd` |
| 2B synteny bands (ATAC) | chromatin notebooks |

## How to run

Open `00_setup.Rmd` in RStudio and run its chunks top to bottom. The first
chunk already selects the default `reproduce` mode, so you can simply run it as
is — no edit required. Sourcing `R/setup.R` then loads packages, declares
`paths`, creates `./data`, `./results` and `./cache`, and prints the resolved
configuration. Then run `01..06` **in order, in the same R session** (later
notebooks reuse objects created by earlier ones).

To use the alternative mode or the fast smoke-test shortcut, change the options
in `00_setup.Rmd` (or set them before sourcing `R/setup.R`):

```r
options(aquafaang.mode = "recompute")   # re-derive from raw RNA-seq; default is "reproduce"
options(aquafaang.test = TRUE)          # short permutations / bootstraps (smoke test)
```

The OrgDB notebooks (`99a`, `99b`) build and install `org.Ssalar.eg.db` /
`org.Omykiss.eg.db` from Ensembl BioMart and only need to be run once per
environment.

Rough wall times in `reproduce` mode:

- `01_data_loading.Rmd`: seconds (loads the deposited normalised matrices).
- `02_soms.Rmd`: seconds (loads the deposited SOM assignments).
- `03_fig2A_expression_soms.Rmd`: ~1 minute.
- `04_ohnolog_coclustering.Rmd`: ~1–2 minutes (downloads salmon raw TPM once for
  the fold-change panel, then cached).
- `05_observed_expected.Rmd`: ~3 hours per permutation block in production,
  about a minute with `TEST_MODE = TRUE` (100 permutations instead of 10,000).
- `06_go.Rmd`: ~15 minutes once the OrgDB packages are installed.

## Data

Two modes, serving two reproducibility goals. The mode is set with
`options(aquafaang.mode = ...)` (see `00_setup.Rmd`).

### `reproduce` (default) — identical figure reproduction

Reads the small deposited derived assets shipped in this repo:

- normalised expression matrices: `data/norm/{sd,td,sb,tb}_norm.tsv.gz`
- SOM unit assignments + summaries: `data/soms/*.tsv.gz` + `*_stage_levels.txt`

This skips RNA-seq loading and normalisation entirely and feeds the saved
matrices / unit assignments straight into the pipeline, reproducing every panel
bit-for-bit (within RNG-free chunks). The only thing downloaded is the
orthogroup table (from Salmobase, cached). The large raw TPM matrices are **not**
shipped here — they are unchanged on Salmobase; use `recompute` mode to derive
the matrices from them.

### `recompute` — independent biological replication

Reads the raw RNA-seq from the public Salmobase endpoints. DevMap is one merged
TSV per species:

- <https://salmobase.org/datafiles/datasets/Aqua-Faang/nfcore/AtlanticSalmon/DevMap/RNA/results/star_rsem/rsem.merged.gene_tpm.tsv>
- <https://salmobase.org/datafiles/datasets/Aqua-Faang/nfcore/RainbowTrout/DevMap/RNA/results/star_rsem/rsem.merged.gene_tpm.tsv>

BodyMap is published per tissue under
`{Species}/BodyMap/RNA/{Tissue}/results/star_rsem/rsem.merged.gene_tpm.tsv`,
with `Tissue` in `{Brain, DistalIntestine, Gill, Gonad, HeadKidney, Liver,
Muscle}`. `process_tissue_data()` pulls the seven per-tissue TSVs and
concatenates them; downloads are cached under `cache/`.

In this mode the SOMs are re-fit from scratch, so the unit numbering — and
therefore the clustering and figures — can differ slightly from the manuscript
(kohonen/BLAS ordering is not bit-stable across platforms; see Reproducibility).
The re-derived matrices and SOMs are written to sibling `*_regen.tsv.gz` files
so the deposited assets are never overwritten.

## Reproducibility

The SOMs whose unit assignments are deposited under `data/soms/` were fit on
R 4.1.3 (Windows) with `set.seed(387334)`. The `kohonen` RNG path is identical
on modern R, but BLAS-level floating-point ordering means a fresh fit on macOS
or Linux produces different unit-to-gene assignments. The notebooks therefore
default (reproduce mode) to reconstructing the deposited SOMs from the gzipped
TSVs; the re-fit chunks in `02_soms.Rmd` are gated behind `recompute` mode and
write a sibling `*_regen` set so the deposited files stay intact.

The observed/expected permutation in `05_observed_expected.Rmd` does not set a
seed, so the empirical p-values drift slightly at the 1e-6 floor across re-runs.

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
| `00_setup.Rmd` | Package stack, `paths` configuration (local or Salmobase), shared constants, helper sourcing. |
| `01_data_loading.Rmd` | Load DevMap and BodyMap TPM matrices, merge replicates, quantile-normalise. |
| `02_soms.Rmd` | Fit (or load) the four SOMs; build per-cluster annotation tables. |
| `03_figure2.Rmd` | Heatmaps, UMAP panels and the Figure 2 composite. |
| `04_ohnolog_coclustering.Rmd` | Per-cluster composition and per-pair co-clustering of AORe and LORe salmon ohnologs across the DevMap and BodyMap SOMs. |
| `05_observed_expected.Rmd` | Observed vs expected ortholog co-clustering against a 10,000-permutation null. RNA and ATAC Sankeys. |
| `06_go.Rmd` | GO biological-process enrichment via `compareCluster` for the four SOMs. |
| `99_helpers.R` | I/O, normalisation, SOM fit and annotation, palettes, heatmap builders. |
| `99a_make_orgdb_salmon.Rmd` | Build and install `org.Ssalar.eg.db` from Ensembl BioMart. |
| `99b_make_orgdb_trout.Rmd` | Build and install `org.Omykiss.eg.db` from Ensembl BioMart. |
| `R/url_or_path.R` | Dispatch for the data loaders to accept either a local path or a Salmobase URL, with caching under `cache/`. |
| `scripts/` | R scripts sourced by `04_ohnolog_coclustering.Rmd` for the ohnolog co-clustering analysis and its figure panels. |

## How to run

Open `00_setup.Rmd` in RStudio and run its chunks first; that sources
`R/setup.R`, which loads the packages, declares `paths`, creates the
expected directories under `./data`, `./results` and `./cache`, and prints
the configuration in scope. Then run the chunks of `01..06` in order.

To switch between data sources or to enable the test-mode shortcuts, set
the relevant option before sourcing `R/setup.R` (or in `00_setup.Rmd`
itself):

```r
options(aquafaang.mode = "salmobase")   # pull from Salmobase; default is "legacy"
options(aquafaang.test = TRUE)          # short permutations / bootstraps
```

The OrgDB notebooks (`99a`, `99b`) build and install
`org.Ssalar.eg.db` / `org.Omykiss.eg.db` from Ensembl BioMart and only
need to be run once per environment.

Rough wall times in `legacy` mode:

- `01_data_loading.Rmd`: seconds with the deposited normalised matrices;
  a few minutes when regenerating.
- `02_soms.Rmd`: seconds when loading; ~10 minutes when re-fitting in
  `salmobase` mode.
- `03_figure2.Rmd`: ~1 minute.
- `04_ohnolog_coclustering.Rmd`: ~1 minute.
- `05_observed_expected.Rmd`: ~3 hours per permutation block in production,
  about a minute with `TEST_MODE = TRUE` (100 permutations instead of
  10,000).
- `06_go.Rmd`: ~15 minutes once the OrgDB packages are installed.

## Data

Two modes are supported, serving two different reproducibility goals.

### `paths$mode = "legacy"` — identical figure reproduction

Reads the deposited normalised expression matrices
(`data/norm/{sd,td,sb,tb}_norm.tsv.gz`) and the deposited SOM unit
assignments + per-unit summary tables (`data/soms/*.tsv.gz` plus the
per-SOM `*_stage_levels.txt`). Knitting in this mode skips the RNA-seq
loading + normalisation step entirely and feeds the saved matrices /
unit assignments straight into the downstream pipeline, reproducing
every panel bit-for-bit (within RNG-free chunks). The orthogroup table
is still fetched from Salmobase on first use; everything else is local.

The raw TPM TSVs that produced these matrices are not redistributed
here — they are large and unchanged on Salmobase. Run `salmobase` mode
if you want to re-derive the matrices from the upstream TPMs.

### `paths$mode = "salmobase"` — independent biological replication

Reads from the public Salmobase endpoints. DevMap is a single merged TSV
per species:

- <https://salmobase.org/datafiles/datasets/Aqua-Faang/nfcore/AtlanticSalmon/DevMap/RNA/results/star_rsem/rsem.merged.gene_tpm.tsv>
- <https://salmobase.org/datafiles/datasets/Aqua-Faang/nfcore/RainbowTrout/DevMap/RNA/results/star_rsem/rsem.merged.gene_tpm.tsv>

BodyMap is published per tissue, one merged TSV per tissue under
`{Species}/BodyMap/RNA/{Tissue}/results/star_rsem/rsem.merged.gene_tpm.tsv`,
with `Tissue` in `{Brain, DistalIntestine, Gill, Gonad, HeadKidney, Liver,
Muscle}`. `process_tissue_data()` pulls the seven per-tissue TSVs and
concatenates them; downloads are cached under `cache/` after the first
fetch.

In this mode the SOM fits are redone from scratch and the resulting unit
numbering will differ from the published figures (kohonen RNG and BLAS
ordering are not bit-stable across platforms; see Reproducibility). Use
this if you want a biological replication: cluster assignments, panel
appearance and exact significance values will differ slightly.

## Reproducibility

The SOMs whose unit assignments are deposited under `data/soms/` were
fit on R 4.1.3 (Windows) with `set.seed(387334)`. The `kohonen` RNG path
is identical on modern R but BLAS-level floating-point ordering means
that a fresh fit on macOS or Linux will produce different unit-to-gene
assignments. The notebooks therefore default to reconstructing the
deposited SOMs from the gzipped TSVs; the regeneration chunks in
`02_soms.Rmd` are gated behind `salmobase` mode and write a sibling
`*_regen` set so the deposited files stay intact.

The observed/expected permutation in `05_observed_expected.Rmd` does not
set a seed, so the empirical p-values drift slightly at the 1e-6 floor
across re-runs.

# Shared setup. Sourced from every notebook in this directory.
#
# -----------------------------------------------------------------------------
# USER TOGGLES — edit either inline here, or override before sourcing this
# file (e.g. options(aquafaang.mode = "salmobase") in an interactive session).
# -----------------------------------------------------------------------------
#
# DATA_MODE controls where the input data is read from:
#   "legacy"    — read the deposited normalised matrices (data/norm/) and
#                 the deposited SOM unit assignments (data/soms/) that
#                 produced the manuscript figures. Skips RNA-seq loading
#                 entirely; everything downstream of normalisation is
#                 reproduced bit-for-bit. The orthogroup table is fetched
#                 from Salmobase on first use (cached locally).
#   "salmobase" — read the per-tissue / per-stage RNA-seq TPM matrices
#                 directly from Salmobase URLs and re-fit the SOMs from
#                 scratch. SOM unit numbering will differ slightly from
#                 the manuscript because kohonen::som is not bit-stable
#                 across BLAS builds.
#
# TEST_MODE = TRUE shortens the slow chunks (permutations, bootstraps) so
# the pipeline can be exercised end-to-end in minutes rather than hours.
# Switch to FALSE for the production runs.

DATA_MODE <- getOption("aquafaang.mode", "legacy")
TEST_MODE <- getOption("aquafaang.test", FALSE)

stopifnot(DATA_MODE %in% c("legacy", "salmobase"))

# -----------------------------------------------------------------------------
# Packages — auto-install on a fresh R session.
# -----------------------------------------------------------------------------
# Union of every CRAN/Bioc dependency used anywhere in this notebook tree
# (00..06, 99a/99b, 99_helpers, scripts/). Listed once at the entry point
# so a user with a clean R install can knit the pipeline without first
# hand-installing dependencies.
.cran_pkgs <- c(
  "tidyverse", "data.table", "foreach", "kohonen", "uwot", "circlize",
  "viridis", "viridisLite", "gridExtra", "ggalluvial", "scales", "broom",
  "cowplot", "ggridges", "ggtext", "magick", "patchwork", "UpSetR"
)
.bioc_pkgs <- c(
  "ComplexHeatmap", "preprocessCore", "AnnotationForge", "biomaRt",
  "clusterProfiler"
)
.ensure_pkgs <- function(cran = character(), bioc = character()) {
  miss_cran <- cran[!vapply(cran, requireNamespace,
                            logical(1), quietly = TRUE)]
  miss_bioc <- bioc[!vapply(bioc, requireNamespace,
                            logical(1), quietly = TRUE)]
  if (length(miss_cran)) {
    message("[setup] installing CRAN packages: ",
            paste(miss_cran, collapse = ", "))
    install.packages(miss_cran, repos = "https://cloud.r-project.org")
  }
  if (length(miss_bioc)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager", repos = "https://cloud.r-project.org")
    }
    message("[setup] installing Bioconductor packages: ",
            paste(miss_bioc, collapse = ", "))
    BiocManager::install(miss_bioc, ask = FALSE, update = FALSE)
  }
}
.ensure_pkgs(.cran_pkgs, .bioc_pkgs)

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(foreach)
  library(kohonen)
  library(uwot)
  library(ComplexHeatmap)
  library(circlize)
  library(viridis)
  library(gridExtra)
  library(ggalluvial)
  library(scales)
  library(preprocessCore)
})

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
source("99_helpers.R")
source("R/url_or_path.R")

# -----------------------------------------------------------------------------
# Constants
# -----------------------------------------------------------------------------
stage_names <- c("LC",
                 "EB", "MB", "LB",
                 "EG", "MG", "LG",
                 "ES", "EMS", "MS", "MLS", "LS",
                 "EP", "LP")

indices_order_ssal <- c(1, 2, 3, 4, 10,
                        12, 9, 11, 21, 22, 23, 24, 17, 18,
                        19, 20, 13, 14, 15, 16, 5, 6,
                        7, 8, 25, 26, 27, 28)

indices_order_omyk <- c(1, 2, 3, 4, 10, 9, 11,
                        20, 21, 22, 23, 16,
                        17, 18, 19, 12, 13,
                        14, 15, 5, 6, 7, 8, 24, 25, 26, 27)

bm_tissues <- c("Brain", "DistalIntestine", "Gill", "Gonad",
                "HeadKidney", "Liver", "Muscle")

# Test-mode iteration counts. The notebooks read these directly.
N_PERMUTATIONS  <- if (TEST_MODE) 100L else 10000L
SOM_LOAD_LEGACY <- TRUE   # always load deposited SOMs in this notebook tree;
                          # set FALSE in salmobase mode to force a fresh fit.
if (DATA_MODE == "salmobase") SOM_LOAD_LEGACY <- FALSE

# -----------------------------------------------------------------------------
# Paths — both modes available; .active is whichever DATA_MODE selects.
# -----------------------------------------------------------------------------
paths <- list(
  mode = DATA_MODE,

  # Local deposited derived assets. Both modes use these for reproducible
  # figure reproduction; salmobase mode also overwrites them with re-fit
  # versions when re-running from raw RNA-seq.
  legacy = list(
    norm_dir = "data/norm",
    soms_dir = "data/soms"
  ),

  # Remote RNA-seq inputs + orthogroups. Always fetched from Salmobase —
  # we do not ship the raw TPM matrices or the orthogroup table in the
  # repository.
  salmobase = list(
    devmap_salmon = "https://salmobase.org/datafiles/datasets/Aqua-Faang/nfcore/AtlanticSalmon/DevMap/RNA/results/star_rsem/rsem.merged.gene_tpm.tsv",
    devmap_trout  = "https://salmobase.org/datafiles/datasets/Aqua-Faang/nfcore/RainbowTrout/DevMap/RNA/results/star_rsem/rsem.merged.gene_tpm.tsv",
    bodymap_salmon_urls = setNames(
      sprintf(paste0("https://salmobase.org/datafiles/datasets/Aqua-Faang/",
                     "nfcore/AtlanticSalmon/BodyMap/RNA/%s/results/star_rsem/",
                     "rsem.merged.gene_tpm.tsv"), bm_tissues),
      bm_tissues),
    bodymap_trout_urls = setNames(
      sprintf(paste0("https://salmobase.org/datafiles/datasets/Aqua-Faang/",
                     "nfcore/RainbowTrout/BodyMap/RNA/%s/results/star_rsem/",
                     "rsem.merged.gene_tpm.tsv"), bm_tissues),
      bm_tissues),
    orthogroups = "https://salmobase.org/datafiles/datasets/Aqua-Faang/salmon-trout-orthologs/SalmonTroutOrthologs.tsv"
  )
)
# Both modes need the salmobase URLs (orthogroups is always fetched from
# there); legacy mode additionally exposes the local deposited assets.
paths$active <- if (DATA_MODE == "legacy") {
  c(paths$legacy, paths$salmobase["orthogroups"])
} else {
  paths$salmobase
}

# -----------------------------------------------------------------------------
# Workspace bootstrap — create only the dirs notebooks WRITE into.
# -----------------------------------------------------------------------------
for (d in c("data/derived", "results", "cache")) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

message(sprintf("[setup] DATA_MODE = %s   TEST_MODE = %s   N_PERMUTATIONS = %d",
                DATA_MODE, TEST_MODE, N_PERMUTATIONS))

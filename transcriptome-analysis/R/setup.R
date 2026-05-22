# Shared setup. Sourced from every notebook in this directory.
#
# -----------------------------------------------------------------------------
# USER TOGGLES — edit either inline here, or override before sourcing this
# file (e.g. options(aquafaang.mode = "salmobase") in an interactive session).
# -----------------------------------------------------------------------------
#
# DATA_MODE controls where the input data is read from:
#   "legacy"    — read the local RData and per-tissue TSVs that produced
#                 the manuscript figures (deterministic, fast, requires the
#                 derived assets to be in data/ already).
#   "salmobase" — read directly from Salmobase. DevMap RNA is a single TSV
#                 per species; BodyMap RNA is per-tissue and the loader
#                 pulls the seven tissue TSVs and concatenates them.
#                 SOM fits are then redone from scratch and the resulting
#                 unit numbering will differ slightly from the manuscript.
#
# TEST_MODE = TRUE shortens the slow chunks (permutations, bootstraps) so
# the pipeline can be exercised end-to-end in minutes rather than hours.
# Switch to FALSE for the production runs.

DATA_MODE <- getOption("aquafaang.mode", "legacy")
TEST_MODE <- getOption("aquafaang.test", FALSE)

stopifnot(DATA_MODE %in% c("legacy", "salmobase"))

# -----------------------------------------------------------------------------
# Packages
# -----------------------------------------------------------------------------
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

  legacy = list(
    soms               = "data/original_soms/soms.RData",
    norm               = "data/original_norm_files/norm.RData",
    devmap_salmon      = "data/devmap/salmon/rsem.merged.gene_tpm.tsv",
    devmap_trout       = "data/devmap/trout/rsem.merged.gene_tpm.tsv",
    bodymap_salmon_dir = "data/bodymap/salmon",
    bodymap_trout_dir  = "data/bodymap/trout",
    orthogroups        = "data/orthogroups/SalmonTroutOrthologs.tsv"
  ),

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
paths$active <- paths[[DATA_MODE]]

# -----------------------------------------------------------------------------
# Workspace bootstrap — create the directories the notebooks write into.
# -----------------------------------------------------------------------------
for (d in c("data", "data/orthogroups",
            "data/original_soms", "data/original_norm_files",
            "data/devmap/salmon", "data/devmap/trout",
            "data/bodymap/salmon", "data/bodymap/trout",
            "data/derived",
            "results", "cache")) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# Seed the orthogroups TSV from Salmobase if it is not already in place.
# 12 MB; downloads once, written to data/orthogroups/ as-is so legacy mode
# does not need to know about the URL after the first run.
if (!file.exists(paths$legacy$orthogroups)) {
  message("[setup] orthogroups TSV missing locally; fetching from Salmobase")
  utils::download.file(paths$salmobase$orthogroups,
                       paths$legacy$orthogroups,
                       mode = "wb", quiet = FALSE)
}

message(sprintf("[setup] DATA_MODE = %s   TEST_MODE = %s   N_PERMUTATIONS = %d",
                DATA_MODE, TEST_MODE, N_PERMUTATIONS))

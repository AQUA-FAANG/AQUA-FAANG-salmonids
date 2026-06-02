# Shared setup for the transcriptome notebooks. Every notebook (00..06,
# 99a/99b) sources this file as its first step. It (1) auto-installs the
# package stack, (2) defines the project constants, and (3) builds the `paths`
# list that tells the loaders where to read data from. Sourcing it twice in a
# session is harmless.

# -----------------------------------------------------------------------------
# USER TOGGLES
# -----------------------------------------------------------------------------
# There are exactly two switches, both set with options() BEFORE this file is
# sourced (00_setup.Rmd shows how). Defaults are chosen so a first-time user
# who changes nothing gets the published figures.
#
#   aquafaang.mode  — where the input data comes from:
#     "reproduce"  (default) Use the small deposited derived assets shipped in
#                  this repo (data/norm/ = normalised expression matrices,
#                  data/soms/ = SOM unit assignments). These are the exact
#                  objects behind the manuscript figures, so the downstream
#                  notebooks reproduce those figures. Fast, deterministic, and
#                  needs almost no download (only the orthogroup table, which
#                  is pulled from Salmobase once and cached).
#     "recompute"  Re-derive everything from the raw RNA-seq TPM matrices,
#                  which are downloaded from Salmobase. The SOMs are then re-fit
#                  from scratch; because kohonen::som is not bit-stable across
#                  BLAS/OS builds, the SOM unit numbering — and therefore the
#                  exact clustering and figures — can differ slightly from the
#                  manuscript depending on the system used.
#
#   aquafaang.test  — TRUE shortens the slow stochastic chunks (permutations,
#                  bootstraps) so the whole pipeline runs in minutes for a
#                  smoke test. Leave FALSE (default) for production runs.

# Read the data-source switch; fall back to "reproduce" when unset.
DATA_MODE <- getOption("aquafaang.mode", "reproduce")
# Read the test-mode switch; fall back to FALSE (full-length run) when unset.
TEST_MODE <- getOption("aquafaang.test", FALSE)

# Guard against typos in the mode string — fail loudly rather than silently
# behaving like one of the two valid modes.
stopifnot(DATA_MODE %in% c("reproduce", "recompute"))

# -----------------------------------------------------------------------------
# Packages — auto-install on a fresh R session
# -----------------------------------------------------------------------------
# Union of every CRAN and Bioconductor package used ANYWHERE in this notebook
# tree (00..06, 99a/99b, 99_helpers.R, scripts/). Declared once here so a user
# with a clean R install can knit the pipeline without hand-installing anything.
.cran_pkgs <- c(
  "tidyverse", "data.table", "foreach", "kohonen", "uwot", "circlize",
  "viridis", "viridisLite", "gridExtra", "ggalluvial", "scales", "broom",
  "cowplot", "ggridges", "ggtext", "magick", "patchwork", "UpSetR"
)
.bioc_pkgs <- c(
  "ComplexHeatmap", "preprocessCore", "AnnotationForge", "biomaRt",
  "clusterProfiler"
)

# .ensure_pkgs — install any of the listed packages that are not already
# available. CRAN packages go through install.packages(); Bioconductor packages
# go through BiocManager::install() (installing BiocManager itself first if
# needed). Packages already present are left untouched, so this is a no-op on a
# configured machine.
.ensure_pkgs <- function(cran = character(), bioc = character()) {
  # Keep only the requested CRAN packages that cannot currently be loaded.
  miss_cran <- cran[!vapply(cran, requireNamespace, logical(1), quietly = TRUE)]
  # Same for the requested Bioconductor packages.
  miss_bioc <- bioc[!vapply(bioc, requireNamespace, logical(1), quietly = TRUE)]
  # Install the missing CRAN packages from a fixed mirror (reproducible URL).
  if (length(miss_cran)) {
    message("[setup] installing CRAN packages: ", paste(miss_cran, collapse = ", "))
    install.packages(miss_cran, repos = "https://cloud.r-project.org")
  }
  # Install the missing Bioconductor packages via BiocManager.
  if (length(miss_bioc)) {
    # BiocManager is itself a CRAN package; install it if it is not present.
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager", repos = "https://cloud.r-project.org")
    }
    message("[setup] installing Bioconductor packages: ", paste(miss_bioc, collapse = ", "))
    # ask=FALSE / update=FALSE: do not prompt and do not touch other packages.
    BiocManager::install(miss_bioc, ask = FALSE, update = FALSE)
  }
}
# Run the installer for this notebook tree's full dependency set.
.ensure_pkgs(.cran_pkgs, .bioc_pkgs)

# Attach the packages the helpers and notebooks call unqualified. (The OrgDB /
# annotation packages used only by 06_go.Rmd and 99a/99b are attached in those
# notebooks, not here.) suppressPackageStartupMessages keeps the knit log clean.
suppressPackageStartupMessages({
  library(tidyverse)       # dplyr/tidyr/ggplot2/readr/stringr/tibble/purrr
  library(data.table)      # fast fread + reshape used by the TPM loaders
  library(foreach)         # %do% loop in process_tissue_data()
  library(kohonen)         # self-organising map fitting (recompute mode)
  library(uwot)            # UMAP embeddings for the Figure 2 panels
  library(ComplexHeatmap)  # SOM expression heatmaps
  library(circlize)        # colorRamp2() colour scales for the heatmaps
  library(viridis)         # perceptually-uniform palettes
  library(gridExtra)       # arrange multiple grobs into a composite
  library(ggalluvial)      # alluvial/sankey geoms
  library(scales)          # axis formatting helpers
  library(preprocessCore)  # quantile normalisation in normalize()
})

# -----------------------------------------------------------------------------
# Helpers — source the shared function library and the URL/path dispatcher.
# -----------------------------------------------------------------------------
source("99_helpers.R")     # process_tissue_data, normalize, get_som, loaders, ...
source("R/url_or_path.R")  # read_table_anywhere(), fetch_to_cache()

# -----------------------------------------------------------------------------
# Constants
# -----------------------------------------------------------------------------
# The 14 DevMap developmental stages, in chronological order. Used to rename
# and reorder the DevMap expression columns after replicate-merging.
stage_names <- c("LC",
                 "EB", "MB", "LB",
                 "EG", "MG", "LG",
                 "ES", "EMS", "MS", "MLS", "LS",
                 "EP", "LP")

# Column permutation that reorders the salmon BodyMap matrix into the canonical
# tissue x maturity x sex order used by every BodyMap figure in the manuscript.
indices_order_ssal <- c(1, 2, 3, 4, 10,
                        12, 9, 11, 21, 22, 23, 24, 17, 18,
                        19, 20, 13, 14, 15, 16, 5, 6,
                        7, 8, 25, 26, 27, 28)

# The equivalent column permutation for the trout BodyMap matrix.
indices_order_omyk <- c(1, 2, 3, 4, 10, 9, 11,
                        20, 21, 22, 23, 16,
                        17, 18, 19, 12, 13,
                        14, 15, 5, 6, 7, 8, 24, 25, 26, 27)

# The seven BodyMap tissues. Drives both the Salmobase per-tissue URL list and
# the column layout that process_tissue_data() assembles.
bm_tissues <- c("Brain", "DistalIntestine", "Gill", "Gonad",
                "HeadKidney", "Liver", "Muscle")

# Permutation count for the observed/expected analysis (05). 100 is enough to
# exercise the code in test mode; 10,000 is the manuscript value.
N_PERMUTATIONS <- if (TEST_MODE) 100L else 10000L

# Whether 02_soms.Rmd should LOAD the deposited SOM assignments (reproduce mode)
# rather than RE-FIT the SOMs from scratch (recompute mode). The 02 chunks are
# gated on this flag so only one of the two paths runs.
SOM_LOAD_DEPOSITED <- (DATA_MODE == "reproduce")

# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------
# `paths` carries BOTH sub-lists at all times so any notebook can reach either
# source explicitly (e.g. 04's fold-change script always needs the raw TPM URLs
# regardless of mode). `paths$active` is the one matching DATA_MODE and is what
# the loaders read by default.
paths <- list(
  mode = DATA_MODE,

  # reproduce: the small derived assets shipped in this repo. Nothing here is
  # downloaded except the orthogroup table (spliced into $active below).
  reproduce = list(
    norm_dir = "data/norm",   # {sd,td,sb,tb}_norm.tsv.gz   (normalised matrices)
    soms_dir = "data/soms"    # {tag}_unit_classif/_unit_stats/_stage_levels
  ),

  # recompute: the raw RNA-seq TPM matrices + the orthogroup table, all served
  # from Salmobase. We never ship these large files in the repository.
  recompute = list(
    # DevMap RNA: one merged TPM matrix per species.
    devmap_salmon = "https://salmobase.org/datafiles/datasets/Aqua-Faang/nfcore/AtlanticSalmon/DevMap/RNA/results/star_rsem/rsem.merged.gene_tpm.tsv",
    devmap_trout  = "https://salmobase.org/datafiles/datasets/Aqua-Faang/nfcore/RainbowTrout/DevMap/RNA/results/star_rsem/rsem.merged.gene_tpm.tsv",
    # BodyMap RNA: one merged TPM matrix per tissue; named vector tissue -> URL.
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
    # Salmon x trout sub-orthogroup table — the cross-species join key. Always
    # from Salmobase (used by BOTH modes), never shipped in the repo.
    orthogroups = "https://salmobase.org/datafiles/datasets/Aqua-Faang/salmon-trout-orthologs/SalmonTroutOrthologs.tsv"
  )
)

# Select the active source set. reproduce mode also needs the orthogroup URL
# (the only thing it downloads), so splice it onto the deposited-asset list.
paths$active <- if (DATA_MODE == "reproduce") {
  c(paths$reproduce, paths$recompute["orthogroups"])
} else {
  paths$recompute
}

# -----------------------------------------------------------------------------
# Workspace bootstrap — create only the directories the notebooks WRITE into
# (read locations are version-controlled and already present).
# -----------------------------------------------------------------------------
for (d in c("data/derived", "results", "cache")) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# Echo the resolved configuration so it is visible at the top of every knit log.
message(sprintf("[setup] DATA_MODE = %s   TEST_MODE = %s   N_PERMUTATIONS = %d",
                DATA_MODE, TEST_MODE, N_PERMUTATIONS))

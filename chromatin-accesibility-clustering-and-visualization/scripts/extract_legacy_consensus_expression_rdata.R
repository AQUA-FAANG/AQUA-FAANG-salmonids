# Materialise the two per-species consensus-expression RDS files that
# 01_atac_loading.Rmd expects under data/atac/, from the manuscript-state
# Ssal_Figure_2.RData / Omyk_Figure_2.RData bundles.
#
# The IDR consensus peak BEDs are no longer extracted here — those come
# from Salmobase as bigBed files on the fly (see 00_setup.Rmd /
# 99_helpers.R::load_idr_urls). This script only writes the
# consensus-expression tables, which are several hundred MB each and have
# not (yet) been redistributed by Salmobase or Zenodo.
#
# Source data location (OneDrive):
#   OneDrive-ImperialCollegeLondon/salmonid_paper_analysis/
#     AQUA-FAANG_Salmon_RNA-seq/Transcriptome_analysis/
#       Ssal_Figure_2.RData
#       Omyk_Figure_2.RData
#
# Override paths via options() before sourcing, e.g.
#   options(aquafaang.ssal_rdata = "/path/Ssal_Figure_2.RData")
#   options(aquafaang.omyk_rdata = "/path/Omyk_Figure_2.RData")
#   options(aquafaang.atac_out   = "data/atac")

ssal_rdata <- getOption(
  "aquafaang.ssal_rdata",
  "~/Library/CloudStorage/OneDrive-ImperialCollegeLondon/salmonid_paper_analysis/AQUA-FAANG_Salmon_RNA-seq/Transcriptome_analysis/Ssal_Figure_2.RData"
)
omyk_rdata <- getOption(
  "aquafaang.omyk_rdata",
  "~/Library/CloudStorage/OneDrive-ImperialCollegeLondon/salmonid_paper_analysis/AQUA-FAANG_Salmon_RNA-seq/Transcriptome_analysis/Omyk_Figure_2.RData"
)
out_root <- getOption("aquafaang.atac_out", "data/atac")

ssal_rdata <- path.expand(ssal_rdata)
omyk_rdata <- path.expand(omyk_rdata)
out_root   <- path.expand(out_root)

stopifnot(file.exists(ssal_rdata), file.exists(omyk_rdata))
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

extract_one <- function(rdata_path, out_rds, candidate_names) {
  cat("Loading", rdata_path, "\n")
  e <- new.env()
  load(rdata_path, envir = e)
  nm <- candidate_names[candidate_names %in% ls(e)][1]
  if (is.na(nm)) {
    stop("None of ", paste(candidate_names, collapse = ", "),
         " found in ", rdata_path)
  }
  expr <- get(nm, envir = e)
  cat("  Writing", out_rds, "(", nrow(expr), "x", ncol(expr), ")\n")
  saveRDS(expr, out_rds)
}

extract_one(ssal_rdata,
            file.path(out_root, "AtlanticSalmon_consensus_expression.RDS"),
            c("ssal_consensus_expression"))
extract_one(omyk_rdata,
            file.path(out_root, "RainbowTrout_consensus_expression.RDS"),
            c("omyk_consensus_expression", "ssal_consensus_expression"))

cat("\nDone. consensus_expression RDS files written to", out_root, "\n")

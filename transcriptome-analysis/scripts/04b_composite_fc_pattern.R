#!/usr/bin/env Rscript
# =============================================================================
# 04b_composite_fc_pattern.R
# -----------------------------------------------------------------------------
# Purpose
#   Assemble the two supplementary composite figures (DevMap + BodyMap) for
#   the FC × pattern mechanism analysis (03f). Layout mirrors 04:
#
#   Panel A : fc_pattern_hex        (2D hex of mag_gap × pattern_r, faceted
#                                    by same_cluster × proxiPhylogeny)
#   Panel B : fc_pattern_ridges     (3 metrics, split by class × same_cluster)
#   Panel C : mechanism_stacked  |  logit_forest        (side by side)
#   Panel D : example_trajectories  (representative pair profiles)
#
#   Inputs : results/ohnolog/panels/C1_*.png (falls back to .pdf if PNG missing)
#   Outputs: results/ohnolog/composite/Figure_Cx_{devmap,bodymap}_fc_pattern
#            .{pdf,png}
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(cowplot)
  library(magick)
})

font_family <- "Helvetica"
theme_set(theme_bw(base_size = 11, base_family = font_family))

panels_dir    <- "results/ohnolog/panels"
composite_dir <- "results/ohnolog/composite"
dir.create(composite_dir, showWarnings = FALSE, recursive = TRUE)

# Read a panel by filename and wrap it as ggdraw. Routes by filename into
# panels/{pdf|png}/{devmap|bodymap}/{const|noconst}/ (same pattern as 04).
read_panel <- function(fname) {
  base_noext <- tools::file_path_sans_ext(fname)
  map        <- if (grepl("_bodymap_", base_noext)) "bodymap" else "devmap"
  variant    <- if (grepl("_noconst$", base_noext)) "noconst" else "const"
  png_path   <- file.path(panels_dir, "png", map, variant, paste0(base_noext, ".png"))
  pdf_path   <- file.path(panels_dir, "pdf", map, variant, fname)
  path <- if (file.exists(png_path)) png_path else pdf_path
  ggdraw() + draw_image(image_read(path))
}

# Build one composite (DevMap or BodyMap).
build_composite <- function(prefix, out_base) {
  A   <- read_panel(sprintf("C1_%s_fc_pattern_hex.pdf",         prefix))
  B   <- read_panel(sprintf("C1_%s_fc_pattern_ridges.pdf",      prefix))
  C_l <- read_panel(sprintf("C1_%s_mechanism_stacked.pdf",      prefix))
  C_r <- read_panel(sprintf("C1_%s_logit_forest.pdf",           prefix))
  D   <- read_panel(sprintf("C1_%s_example_trajectories.pdf",   prefix))

  panelC <- plot_grid(C_l, C_r, ncol = 2, rel_widths = c(1, 1),
                      labels = NULL)

  top_row <- plot_grid(A, B, ncol = 2, rel_widths = c(1, 1),
                       labels = c("A", "B"), label_size = 16)
  mid_row <- plot_grid(panelC, ncol = 1, labels = c("C"), label_size = 16)
  bot_row <- plot_grid(D,      ncol = 1, labels = c("D"), label_size = 16)

  composite <- plot_grid(top_row, mid_row, bot_row, ncol = 1,
                         rel_heights = c(1, 0.9, 1.3)) +
    theme(plot.background = element_rect(fill = "white", colour = NA))

  w <- 14; h <- 18
  ggsave(file.path(composite_dir, paste0(out_base, ".pdf")),
         composite, width = w, height = h, device = cairo_pdf, bg = "white")
  ggsave(file.path(composite_dir, paste0(out_base, ".png")),
         composite, width = w, height = h, dpi = 300, bg = "white")
  cat("  Saved", out_base, "(.pdf and .png)\n")
}

cat("=== DevMap FC-pattern composite ===\n")
build_composite("devmap",  "Figure_Cx_devmap_fc_pattern")
cat("=== BodyMap FC-pattern composite ===\n")
build_composite("bodymap", "Figure_Cx_bodymap_fc_pattern")

cat("\nDone.\n")

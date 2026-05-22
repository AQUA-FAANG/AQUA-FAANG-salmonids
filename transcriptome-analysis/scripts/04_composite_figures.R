#!/usr/bin/env Rscript
# =============================================================================
# 04_composite_figures.R
# -----------------------------------------------------------------------------
# Purpose
#   Assemble the two supplementary composite figures (DevMap + BodyMap) from
#   the individual panel PNG/PDFs produced by 03c / 03d / 03e.
#
# Layout (identical for both maps)
#   Panel A (top-left)  : (i) coclustering_by_{stage,tissue}  (stacked)
#                         (ii) lore_proportion_noconst
#   Panel B (top-right) : sankey_early_noconst | sankey_late_noconst
#   Panel C (bottom)    : upset_early_noconst  /  upset_late_noconst (stacked)
#
# Inputs  : results/ohnolog/panels/C1_*.png (falls back to .pdf if PNG missing)
# Outputs : results/ohnolog/composite/Figure_Cx_{devmap,bodymap}_ohnolog_coclustering
#           .{pdf,png}
#
# Why PNG > PDF as the source
#   UpSetR draws with grid/base graphics. pdftools::pdf_render_page (used by
#   magick::image_read_pdf) rasterizes these PDFs inconsistently — under some
#   paths the whole UpSet panel comes out blank. The PNGs that 03e writes
#   alongside the PDFs always render correctly, so we prefer PNGs here.
#
# Why explicit white background
#   cowplot / draw_image produce transparent regions where the child image's
#   aspect ratio doesn't fill the allocated cell; many viewers render
#   transparent PNGs as black. Forcing plot.background = white on the final
#   canvas avoids that.
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

# Read a panel by filename and wrap it as a ggdraw object that can be
# composed with plot_grid. Accepts a .pdf filename for API symmetry but
# prefers the sibling .png when present (see header note on UpSetR).
# Args:  fname — panel filename (with .pdf extension) in `panels_dir`.
# Returns: a ggdraw/ggplot.
read_panel <- function(fname) {
  # Route by filename to panels/{pdf|png}/{devmap|bodymap}/{const|noconst}/.
  # Prefer the PNG because UpSetR PDFs rasterise unreliably via pdftools.
  base_noext <- tools::file_path_sans_ext(fname)
  map        <- if (grepl("_bodymap_", base_noext)) "bodymap" else "devmap"
  variant    <- if (grepl("_noconst$", base_noext)) "noconst" else "const"
  png_path   <- file.path(panels_dir, "png", map, variant, paste0(base_noext, ".png"))
  pdf_path   <- file.path(panels_dir, "pdf", map, variant, fname)
  path <- if (file.exists(png_path)) png_path else pdf_path
  img <- image_read(path)
  ggdraw() + draw_image(img)
}

# Build one composite (DevMap or BodyMap) and write the paired PDF + PNG.
# Args:
#   prefix   — "devmap" or "bodymap"; drives panel filename substitution.
#   out_base — composite filename without extension.
# Returns: invisible NULL (saves files as side effect).
build_composite <- function(prefix, out_base) {
  by_group_file <- if (prefix == "devmap")
    "C1_devmap_coclustering_by_stage.pdf"
  else
    "C1_bodymap_coclustering_by_tissue.pdf"

  A_top    <- read_panel(by_group_file)
  A_bot    <- read_panel(sprintf("C1_%s_composition_proportion_noconst.pdf", prefix))
  B_left   <- read_panel(sprintf("C1_%s_sankey_early_noconst.pdf", prefix))
  B_right  <- read_panel(sprintf("C1_%s_sankey_late_noconst.pdf",  prefix))
  C_top    <- read_panel(sprintf("C1_%s_upset_early_noconst.pdf",  prefix))
  C_bot    <- read_panel(sprintf("C1_%s_upset_late_noconst.pdf",   prefix))

  # Sub-grids
  panelA <- plot_grid(A_top, A_bot, ncol = 1, rel_heights = c(1, 1),
                      labels = c("i", "ii"), label_size = 11,
                      label_fontface = "italic")
  panelB <- plot_grid(B_left, B_right, ncol = 2, rel_widths = c(1, 1),
                      labels = c("Early", "Late"), label_size = 11,
                      label_fontface = "plain", label_x = 0.5, hjust = 0.5)
  panelC <- plot_grid(C_top, C_bot, ncol = 1, rel_heights = c(1, 1),
                      labels = c("Early", "Late"), label_size = 11,
                      label_fontface = "plain", label_x = 0.02)

  top_row    <- plot_grid(panelA, panelB, ncol = 2, rel_widths = c(1, 1.2),
                          labels = c("A", "B"), label_size = 16)
  composite  <- plot_grid(top_row, panelC, ncol = 1, rel_heights = c(1, 1.1),
                          labels = c("", "C"), label_size = 16)

  composite <- composite +
    theme(plot.background = element_rect(fill = "white", colour = NA))

  w <- 14; h <- 16
  ggsave(file.path(composite_dir, paste0(out_base, ".pdf")),
         composite, width = w, height = h, device = cairo_pdf, bg = "white")
  ggsave(file.path(composite_dir, paste0(out_base, ".png")),
         composite, width = w, height = h, dpi = 300, bg = "white")
  cat("  Saved", out_base, "(.pdf and .png)\n")
}

cat("=== Building DevMap composite ===\n")
build_composite("devmap", "Figure_Cx_devmap_ohnolog_coclustering")

cat("=== Building BodyMap composite ===\n")
build_composite("bodymap", "Figure_Cx_bodymap_ohnolog_coclustering")

cat("\nDone.\n")

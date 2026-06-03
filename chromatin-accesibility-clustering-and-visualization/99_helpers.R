# Shared helpers for the ATAC clustering and visualisation notebooks.

# === I/O ===

#' Read a deposited ATAC SOM asset, `{tag}_{what}.tsv.gz` under
#' `dir` (default data/atac_soms/). These are the exact manuscript SOM products
#' loaded by 01/02 (the ATAC analogue of the transcriptome data/soms/).
#' `tag` is one of devmap_salmon / bodymap_salmon / devmap_trout /
#' bodymap_trout; `what` is one of:
#'   "idx"          peak subset (row indices into consensus_expression)
#'   "unit_classif" SOM unit per subset peak (row-aligned to the FE matrix)
#'   "superclasses" PAM(k=4) super-cluster per SOM unit (length 25)
#'   "umap"         2-D UMAP coordinates per subset peak (matrix, cols V1,V2)
read_atac_som_asset <- function(tag, what, dir = "data/atac_soms") {
  tbl <- readr::read_tsv(file.path(dir, paste0(tag, "_", what, ".tsv.gz")),
                         show_col_types = FALSE)
  switch(what,
         idx          = as.integer(tbl$idx),
         unit_classif = as.integer(tbl$unit),
         superclasses = as.integer(tbl$superclass),
         umap         = as.matrix(tbl[, c("V1", "V2")]),
         stop("unknown asset '", what, "'"))
}

#' Build the scaled fold-enrichment matrix used as SOM input. ATAC columns are
#' divided by their column mean (fold enrichment over the per-sample average),
#' then each region is z-scaled across samples (center = FALSE, scale = SD).
build_scaled_fe_matrix <- function(expression_df, atac_col_idx) {
  fe <- expression_df %>%
    dplyr::select(dplyr::contains("ATAC")) %>%            # keep only ATAC columns
    dplyr::select(dplyr::all_of(atac_col_idx)) %>%        # subset to this context
    dplyr::mutate(dplyr::across(dplyr::everything(), ~ .x / mean(.x)))  # fold enrichment vs column mean
  # Transpose so scale() works per region, z-scale (SD only, no centering), transpose back.
  fe %>% as.matrix() %>% t() %>% scale(center = FALSE) %>% t()
}

# === SOM summary ===

#' Ported from SOMbodymapATACfunc.R. Builds the per-unit label vector, long
#' tidy frame and per-unit mean/sd table used by the ribbon plots, then
#' returns them together with the input SOM and a default ribbon ggplot.
get_som <- function(in_mat, xdim, ydim, somvar) {
  som_obj <- somvar

  # Per-unit facet labels of the form "Class <id> (<n regions>)".
  label_unit <- table(som_obj$unit.classif)
  tmp_name <- names(label_unit)
  label_unit <- stringr::str_c("Class ", tmp_name, " (", label_unit, ")")
  names(label_unit) <- tmp_name

  # Long tidy frame: one row per (region, sample), tagged with its SOM unit;
  # stage/class as ordered factors so facets/axes read in the right order.
  tmp_df <- in_mat %>% as.data.frame() %>%
    {cbind(class = som_obj$unit.classif, .)} %>%
    tidyr::pivot_longer(cols = !class,
                        names_to = "stage",
                        values_to = "value") %>%
    dplyr::mutate(stage = factor(stage, levels = colnames(in_mat)),
                  class = factor(class,
                                 levels = stringr::str_sort(unique(class),
                                                            numeric = TRUE)))

  # Per-(unit, sample) mean and SD, the ribbon-plot summary table.
  sum_stat <- tmp_df %>% dplyr::group_by(class, stage) %>%
    dplyr::summarize(mean_value = mean(value),
                     sd = stats::sd(value),
                     .groups = "drop")

  out_plot <- ggplot2::ggplot(sum_stat,
                              ggplot2::aes(x = stage, y = mean_value,
                                           group = class)) +
    ggplot2::geom_line(linewidth = 1.5) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = mean_value - 2 * sd,
                                      ymax = mean_value + 2 * sd),
                         alpha = 0.2) +
    ggplot2::facet_wrap(~ class, ncol = xdim,
                        labeller = ggplot2::labeller(class = label_unit),
                        scale = "free_y", as.table = FALSE) +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   legend.position = "none")

  list(som = som_obj,
       label_unit = label_unit,
       tmp_df = tmp_df,
       stat = sum_stat,
       plot = out_plot)
}

# === palettes ===

#' Build the inferno-plus-grey palette used for the DevMap UMAPs. The first
#' `n_dynamic` SOM units (i.e. the dynamic / non-constitutive units listed in
#' `splits_levels`) are coloured with viridis::inferno restricted to the
#' 0.285-0.627 range; the remaining `n_constitutive` units share the
#' "#becacaff" grey. The returned vector is reordered so that
#' `palette[i]` is the colour for SOM unit `i`.
#'
#' @param n_dynamic number of dynamic SOM units (e.g. 17 for salmon DevMap).
#' @param n_constitutive number of constitutive SOM units (e.g. 8).
#' @param order_perm character vector of length `n_dynamic + n_constitutive`
#'   giving the SOM unit labels in the order in which the palette colours
#'   should be applied (the union of `splits_levels` and the constitutive
#'   units, in the order they appear in the figure).
inferno_plus_grey_palette <- function(n_dynamic, n_constitutive, order_perm) {
  dyn  <- viridis::inferno(n_dynamic, begin = 0.285, end = 0.627)  # dynamic-unit colours
  cons <- rep("#becacaff", n_constitutive)                          # constitutive grey
  final_colors <- c(dyn, cons)                                      # in figure order
  total <- n_dynamic + n_constitutive
  # Invert order_perm so the returned vector is keyed by raw SOM unit id
  # (palette[i] = colour for unit i), as the ggplot scales expect.
  color_order_for_objects <- match(seq_len(total), order_perm)
  final_colors[color_order_for_objects]
}

#' Hand-curated tissue palette for the salmon BodyMap UMAP (25 SOM units).
salmon_bodymap_palette <- function() {
  c("#8a0000ff", "skyblue",   "forestgreen", "skyblue",    "#ffa300ff",
    "#8a0000ff", "#8a0000ff", "forestgreen", "skyblue",    "skyblue",
    "#becacaff", "#becacaff", "#dce648ff",   "#becacaff",  "#becacaff",
    "skyblue",   "#ffa300ff", "skyblue",     "#dce648ff",  "#becacaff",
    "salmon",    "salmon",    "salmon",      "#dce648ff",  "#dce648ff")
}

#' Hand-curated tissue palette for the trout BodyMap UMAP (25 SOM units).
trout_bodymap_palette <- function() {
  c("skyblue",   "salmon",    "salmon",    "skyblue",   "skyblue",
    "skyblue",   "skyblue",   "skyblue",   "skyblue",   "skyblue",
    "skyblue",   "skyblue",   "#becacaff", "#becacaff", "#becacaff",
    "#8a0000ff", "#8a0000ff", "skyblue",   "skyblue",   "#dce648ff",
    "#8a0000ff", "#becacaff", "#becacaff", "#becacaff", "#dce648ff")
}

# === plots ===

#' Render the ATAC SOM heatmap used in Figure 2. Rows are restricted to the
#' SOM units listed in `splits_levels` (the dynamic units) and split by that
#' factor in the given order. Columns are kept in input order so the per-
#' species sample axis is preserved.
make_atac_heatmap <- function(in_mat, som, splits_levels, palette,
                              show_column_names = TRUE,
                              show_heatmap_legend = FALSE,
                              use_raster = TRUE) {
  # Restrict rows to regions assigned to the listed (dynamic) SOM units.
  keep <- which(som$unit.classif %in% splits_levels)
  ComplexHeatmap::Heatmap(
    in_mat[keep, ],
    cluster_rows = TRUE,
    row_split = factor(som$unit.classif[keep], levels = splits_levels),
    show_row_dend = FALSE,
    cluster_columns = FALSE,
    cluster_row_slices = FALSE,
    show_column_names = show_column_names,
    col = palette,
    name = "Normalized peak intensity",
    raster_by_magick = TRUE,
    show_heatmap_legend = show_heatmap_legend,
    use_raster = use_raster
  )
}

#' Render the SOM-coloured ATAC UMAP. `flip_y` mirrors the y-axis (used for
#' the DevMap panels so the colour gradient matches the heatmap row order).
make_atac_umap <- function(umap_df, som, palette, flip_y = FALSE) {
  df <- as.data.frame(umap_df)
  y_expr <- if (flip_y) ggplot2::aes(V1, V2 * (-1)) else ggplot2::aes(V1, V2)
  ggplot2::ggplot(df, y_expr) +
    ggplot2::geom_point(ggplot2::aes(colour = factor(som$unit.classif)),
                        alpha = 0.4, size = 0.1) +
    ggplot2::scale_color_manual(values = palette) +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "none",
                   axis.line = ggplot2::element_blank(),
                   axis.text = ggplot2::element_blank(),
                   axis.title = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank())
}

#' Ribbon panel summarising mean +/- 2 SD per SOM unit across samples. Used
#' for the four DevMap / BodyMap ribbon strips in Figure 2.
make_ribbon_panel <- function(tmp_df, label_unit, ncol = 5) {
  ggplot2::ggplot(tmp_df, ggplot2::aes(stage, value, group = "Class")) +
    ggplot2::stat_summary(fun.data = "mean_sdl",
                          fun.args = list(mult = 2), geom = "ribbon") +
    ggplot2::stat_summary(fun = mean, geom = "line") +
    ggplot2::facet_wrap(~ class, ncol = ncol,
                        labeller = ggplot2::labeller(class = label_unit),
                        as.table = FALSE) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

# Shared helpers for the ATAC clustering and visualisation notebooks.

# === I/O ===

#' Load IDR consensus BEDs from a directory into a named list of GRanges.
#'
#' Wrapper around rtracklayer::import that strips the
#' "<Species>_ATAC_" prefix and the "_consensus_peaks.bed" suffix from the
#' file names so the resulting list is keyed by tissue / stage labels.
#'
#' @param dir directory containing the IDR consensus BED files.
#' @param species_prefix prefix to remove from each file basename (e.g.
#'   "AtanticSalmon_ATAC_" or "RainbowTrout_ATAC_"). The misspelling of
#'   "AtanticSalmon" reproduces the legacy file naming.
#' @param reader function used to read each BED. Defaults to
#'   rtracklayer::import; pass a tsv reader for trout BodyMap BEDs that have
#'   extra annotation columns.
load_idr_dir <- function(dir, species_prefix = "", reader = rtracklayer::import) {
  files <- list.files(dir, pattern = "\\.bed$", full.names = TRUE)
  if (!length(files)) {
    stop("No .bed files found in ", dir)
  }
  out <- lapply(files, reader)
  names(out) <- sub("_consensus_peaks\\.bed$", "",
                    sub(species_prefix, "", basename(files), fixed = TRUE))
  out
}

#' Reader for the trout BodyMap IDR BEDs, which carry extra ChromHMM columns
#' beyond the standard BED6.
read_trout_bodymap_bed <- function(x) {
  tbl <- readr::read_tsv(
    x,
    col_names = c("seqnames", "start", "end", "name",
                  "score", "strand", "tickStart", "tickEnd", "itemRgb",
                  "no_idea_1", "no_idea_2", "no_idea_3", "Replicates",
                  "ChromHMM_annotation"),
    show_col_types = FALSE
  )
  GenomicRanges::makeGRangesFromDataFrame(
    tbl, keep.extra.columns = TRUE, starts.in.df.are.0based = TRUE
  )
}

#' Reduce a list of GRanges into a single union (non-overlapping intervals).
union_granges <- function(gr_list) {
  as(gr_list, "GRangesList") %>% unlist() %>% GenomicRanges::reduce()
}

#' Build the scaled fold-enrichment matrix used as SOM input. ATAC columns are
#' divided by their column mean (fold enrichment over the per-sample average),
#' then each region is z-scaled across samples (center = FALSE, scale = SD).
build_scaled_fe_matrix <- function(expression_df, atac_col_idx) {
  fe <- expression_df %>%
    dplyr::select(dplyr::contains("ATAC")) %>%
    dplyr::select(dplyr::all_of(atac_col_idx)) %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), ~ .x / mean(.x)))
  fe %>% as.matrix() %>% t() %>% scale(center = FALSE) %>% t()
}

# === SOM wrappers ===

#' Fit a 5 x 5 hexagonal SOM with the legacy somInit initialisation. The
#' wrapper sets the seed for reproducibility and returns the kohonen SOM
#' object.
#'
#' @param in_mat scaled fold-enrichment matrix (regions x samples).
#' @param xdim,ydim grid dimensions; default 5 x 5.
#' @param seed integer seed; default 12345 to match the manuscript.
fit_atac_som <- function(in_mat, xdim = 5, ydim = 5, seed = 12345) {
  set.seed(seed)
  init <- aweSOM::somInit(traindat = in_mat, xdim, ydim)
  set.seed(seed)
  kohonen::som(in_mat,
               grid = kohonen::somgrid(xdim, ydim, "hexagonal"),
               init = init)
}

#' Ported from SOMbodymapATACfunc.R. Builds the per-unit label vector, long
#' tidy frame and per-unit mean/sd table used by the ribbon plots, then
#' returns them together with the input SOM and a default ribbon ggplot.
get_som <- function(in_mat, xdim, ydim, somvar) {
  som_obj <- somvar

  label_unit <- table(som_obj$unit.classif)
  tmp_name <- names(label_unit)
  label_unit <- stringr::str_c("Class ", tmp_name, " (", label_unit, ")")
  names(label_unit) <- tmp_name

  tmp_df <- in_mat %>% as.data.frame() %>%
    {cbind(class = som_obj$unit.classif, .)} %>%
    tidyr::pivot_longer(cols = !class,
                        names_to = "stage",
                        values_to = "value") %>%
    dplyr::mutate(stage = factor(stage, levels = colnames(in_mat)),
                  class = factor(class,
                                 levels = stringr::str_sort(unique(class),
                                                            numeric = TRUE)))

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

#' Convenience: PAM super-clustering on a SOM's prototype codes.
pam_superclusters <- function(som_obj, k = 4) {
  cluster::pam(som_obj$codes[[1]], k)$clustering
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
  dyn  <- viridis::inferno(n_dynamic, begin = 0.285, end = 0.627)
  cons <- rep("#becacaff", n_constitutive)
  final_colors <- c(dyn, cons)
  total <- n_dynamic + n_constitutive
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

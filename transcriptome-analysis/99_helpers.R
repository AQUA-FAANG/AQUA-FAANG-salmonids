# Shared helpers for the transcriptome notebooks.

# === I/O =====================================================================

# process_tissue_data — concatenate the seven per-tissue BodyMap TSVs into
# one 84-column raw matrix in the order merge_reps_bodymap() expects. Accepts
# either `file_path` (directory of local TSVs) or `urls` (named vector of
# tissue → URL, e.g. Salmobase per-tissue endpoints). URLs are downloaded
# once via fetch_to_cache(); the cached files are then read by data.table.
process_tissue_data <- function(fish, file_path = NULL, urls = NULL) {
  if (!is.null(urls)) {
    raw_files <- vapply(sort(urls), fetch_to_cache, character(1),
                        USE.NAMES = FALSE)
  } else {
    if (is.null(file_path)) file_path <- paste0("../data/bodymap/", fish)
    raw_files <- list.files(path = file_path,
                            pattern = "gene_tpm.tsv",
                            full.names = TRUE)
  }

  temp <- lapply(raw_files, fread, sep = "\t")
  nams <- vector(mode = "character", length = 84)
  nseq <- seq(1, 73, 12)

  raw_data <- data.frame(matrix(nrow = nrow(temp[[1]]), ncol = 84))
  row.names(raw_data) <- temp[[1]]$gene_id

  foreach(i = seq_along(temp), l = nseq) %do% {
    raw_data[, l:(l + 11)] <- temp[[i]][, 3:14]
    nams[l:(l + 11)] <- colnames(temp[[i]][, 3:14])
  }

  colnames(raw_data) <- nams
  raw_data
}

# === normalization ===========================================================

merge_reps_devmap <- function(unmerged) {
  setDT(unmerged)

  unmerged <- unmerged[, `transcript_id(s)` := NULL]
  unmerged <- melt(unmerged,
                   id.vars = "gene_id",
                   variable.name = "sample",
                   value.name = "tpm")
  unmerged[, c("stage", "replicate") := tstrsplit(sample, "_R", fixed = TRUE)]

  merged <- unmerged[, .(tpm = mean(tpm, na.rm = TRUE)),
                     by = .(gene_id, stage)]

  merged <- as_tibble(merged)
  merged <- merged %>%
    pivot_wider(names_from = "stage", values_from = "tpm") %>%
    column_to_rownames(var = "gene_id")

  # Reorder columns into chronological stage order, then rename to `stage_names`.
  merged <- merged[, c(12, 13, 14, 6, 1, 7, 2, 8, 3, 9, 4, 10, 5, 11)]
  colnames(merged) <- stage_names
  merged
}

merge_reps_bodymap <- function(unmerged) {
  merged <- unmerged %>% rownames_to_column("gene_id")
  setDT(merged)

  merged <- melt(merged,
                 id.vars = "gene_id",
                 variable.name = "sample",
                 value.name = "tpm")

  merged[, c("Species", "Seq", "Tissue", "Maturity", "Sex", "Replicate") :=
           tstrsplit(sample, "_", fixed = TRUE)]
  merged[, sample := NULL]
  merged[, Name := paste(Tissue, Maturity, Sex, sep = "_")]

  merged <- merged[, .(tpm = mean(tpm, na.rm = TRUE)), by = .(gene_id, Name)]
  merged <- as_tibble(merged)

  merged %>%
    pivot_wider(names_from = "Name", values_from = "tpm") %>%
    column_to_rownames(var = "gene_id") %>%
    rename_with(~ gsub("Gonad", "Ova", .x),
                starts_with("Gonad_Immature_Female")) %>%
    rename_with(~ gsub("Gonad", "Testes", .x),
                starts_with("Gonad_Immature_Male")) %>%
    rename_with(~ gsub("Gonad", "Ova", .x),
                starts_with("Gonad_Mature_Female")) %>%
    rename_with(~ gsub("Gonad", "Testes", .x),
                starts_with("Gonad_Mature_Male"))
}

normalize <- function(merged) {
  row_names <- merged %>% filter(if_any(everything(), ~ .x > 1))
  row_names <- rownames(row_names)

  norm <- merged %>%
    filter(if_any(everything(), ~ .x > 1)) %>%
    data.matrix() %>%
    preprocessCore::normalize.quantiles() %>%
    t() %>% scale(center = FALSE, scale = TRUE) %>% t()

  rownames(norm) <- row_names
  colnames(norm) <- colnames(merged)
  norm
}

# === SOM fit =================================================================

get_som <- function(in_mat, xdim, ydim) {
  som_obj <- kohonen::som(in_mat,
                          grid = somgrid(xdim, ydim, topo = "hexagonal"))

  label_unit <- table(som_obj$unit.classif)
  tmp_name <- names(label_unit)
  label_unit <- str_c("Class ", tmp_name, " (", label_unit, ")")
  names(label_unit) <- tmp_name

  tmp_df <- in_mat %>% as.data.frame() %>%
    {cbind(class = som_obj$unit.classif, .)} %>%
    pivot_longer(cols = !class,
                 names_to = "stage",
                 values_to = "value") %>%
    mutate(stage = factor(stage, levels = colnames(in_mat)),
           class = factor(class, levels = str_sort(unique(class), numeric = TRUE)))

  sum_stat <- tmp_df %>% group_by(class, stage) %>%
    summarize(mean_value = mean(value),
              sd = sd(value),
              .groups = "drop") %>% ungroup()

  out_plot <- ggplot(sum_stat, aes(x = stage, y = mean_value, group = class)) +
    geom_line(size = 1.5) +
    geom_ribbon(aes(ymin = mean_value - 2 * sd,
                    ymax = mean_value + 2 * sd),
                alpha = 0.2) +
    facet_wrap(~ class, ncol = xdim, labeller = labeller(class = label_unit),
               scale = "free_y", as.table = FALSE) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")

  list(som = som_obj, stat = sum_stat, plot = out_plot)
}

get_som_bodymap <- function(in_mat, xdim, ydim, order_indices = NULL) {
  if (is.null(order_indices)) order_indices <- seq_len(ncol(in_mat))
  ordered_stages <- colnames(in_mat)[order_indices]

  som_obj <- kohonen::som(in_mat,
                          grid = somgrid(xdim, ydim, topo = "hexagonal"))

  label_unit <- table(som_obj$unit.classif)
  tmp_name <- names(label_unit)
  label_unit <- str_c("Class ", tmp_name, " (", label_unit, ")")
  names(label_unit) <- tmp_name

  tmp_df <- in_mat %>% as.data.frame() %>%
    {cbind(class = som_obj$unit.classif, .)} %>%
    pivot_longer(cols = !class,
                 names_to = "stage",
                 values_to = "value") %>%
    mutate(stage = factor(stage, levels = ordered_stages),
           class = factor(class, levels = str_sort(unique(class), numeric = TRUE)))

  sum_stat <- tmp_df %>% group_by(class, stage) %>%
    summarize(mean_value = mean(value), sd = sd(value), .groups = "drop") %>%
    ungroup()

  out_plot <- ggplot(sum_stat, aes(x = stage, y = mean_value, group = class)) +
    geom_line(size = 1.5) +
    geom_ribbon(aes(ymin = mean_value - 2 * sd,
                    ymax = mean_value + 2 * sd), alpha = 0.2) +
    facet_wrap(~ class, ncol = xdim, labeller = labeller(class = label_unit),
               scale = "free_y", as.table = FALSE) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")

  list(som = som_obj, stat = sum_stat, plot = out_plot)
}

# === SOM annotation ==========================================================

som_ids <- function(norm_file, som_file) {
  listid <- list()
  for (i in seq_along(unique(som_file$som$unit.classif))) {
    listid[[i]] <- rownames(norm_file[which(som_file$som$unit.classif == i), ])
  }
  listid
}

som_file_info_d <- function(som_file) {
  class_max <- vector(); constitutive_class <- vector()
  class_max_mean <- vector(); class_number <- c(1:16)

  for (i in seq_along(unique(som_file$som$unit.classif))) {
    class_max[i] <-
      stage_names[which.max(som_file$stat$mean_value[which(som_file$stat$class == i)])]
    class_max_mean[i] <-
      som_file$stat$mean_value[
        which.max(som_file$stat$mean_value[which(som_file$stat$class == i)]) +
          (length(stage_names) * (i - 1))]

    if (sd(som_file$stat$mean_value[which(som_file$stat$class == i)]) < 0.4) {
      constitutive_class[i] <- "TRUE"
    } else {
      constitutive_class[i] <- "FALSE"
    }
  }

  class_info <- as.data.frame(cbind(class_number, class_max,
                                    class_max_mean, constitutive_class))
  class_info$class_max <- factor(class_info$class_max, levels = stage_names)

  arrange(class_info, constitutive_class, class_max, desc(class_max_mean))
}

som_file_info_b <- function(som_file) {
  class_max <- vector(); constitutive_class <- vector()
  stage_specific <- vector(); class_max_mean <- vector()
  class_number <- c(1:36)
  tissue_names <- as.vector(unique(som_file$stat$stage))

  for (i in seq_along(unique(som_file$som$unit.classif))) {
    class_max[i] <-
      tissue_names[which.max(som_file$stat$mean_value[which(som_file$stat$class == i)])]
    class_max_mean[i] <-
      som_file$stat$mean_value[
        which.max(som_file$stat$mean_value[which(som_file$stat$class == i)]) +
          (length(unique(som_file$stat$stage)) * (i - 1))]

    if (sd(som_file$stat$mean_value[which(som_file$stat$class == i)]) < 0.4) {
      constitutive_class[i] <- "TRUE"
    } else {
      constitutive_class[i] <- "FALSE"
    }

    # 2.165 is the per-cluster peak-mean threshold above which the SOM unit
    # is treated as tissue-specific (matched manually against the BodyMap
    # heatmaps; see Clustering_and_Visualization_Fig2.Rmd for the manual
    # adjustments applied on top).
    if (class_max_mean[i] > 2.165) {
      stage_specific[i] <- "TRUE"
    } else {
      stage_specific[i] <- "FALSE"
    }
  }

  class_info <- as.data.frame(cbind(class_number, class_max, class_max_mean,
                                    constitutive_class, stage_specific))
  class_info$class_max <- factor(class_info$class_max, levels = tissue_names)
  class_info$tissue <- sub("\\_.*", "", class_info$class_max)

  arrange(class_info, constitutive_class, desc(stage_specific),
          class_max, desc(class_max_mean))
}

# === palettes ================================================================

palette_devmap <- function(n_stage_specific, n_constitutive, som_order) {
  stage_specific_cols <- inferno(n_stage_specific, begin = 0.15, end = 0.9)
  constitutive <- c("#becacaff")
  final_colors <- c(stage_specific_cols, rep(constitutive, n_constitutive))

  color_order <- match(1:16, som_order)
  final_colors[color_order]
}

palette_bodymap <- function(som_info) {
  brain_col       <- "#dce648ff"
  testes_col      <- "forestgreen"
  ova_col         <- "#ffa300ff"
  muscle_col      <- "salmon"
  liver_col       <- "#8a0000ff"
  head_kidney_col <- "black"
  gill_col        <- "red2"
  intestine_col   <- "hotpink"
  multi_col       <- "skyblue"
  constitutive_col <- "#becacaff"

  pal <- vector()
  for (i in seq_len(nrow(som_info))) {
    if (som_info$tissue[i] == "Brain" & som_info$stage_specific[i] == TRUE) {
      pal[i] <- brain_col
    } else if (som_info$tissue[i] == "Testes" & som_info$stage_specific[i] == TRUE) {
      pal[i] <- testes_col
    } else if (som_info$tissue[i] == "Ova" & som_info$stage_specific[i] == TRUE) {
      pal[i] <- ova_col
    } else if (som_info$tissue[i] == "Muscle" & som_info$stage_specific[i] == TRUE) {
      pal[i] <- muscle_col
    } else if (som_info$tissue[i] == "Liver" & som_info$stage_specific[i] == TRUE) {
      pal[i] <- liver_col
    } else if (som_info$tissue[i] == "HeadKidney" & som_info$stage_specific[i] == TRUE) {
      pal[i] <- head_kidney_col
    } else if (som_info$tissue[i] == "Gill" & som_info$stage_specific[i] == TRUE) {
      pal[i] <- gill_col
    } else if (som_info$tissue[i] == "SuppDistalIntestine" & som_info$stage_specific[i] == TRUE) {
      pal[i] <- intestine_col
    } else if (som_info$stage_specific[i] == FALSE & som_info$constitutive_class[i] == FALSE) {
      pal[i] <- multi_col
    } else if (som_info$constitutive_class[i] == TRUE) {
      pal[i] <- constitutive_col
    }
  }

  pal
}

# === heatmaps ================================================================

make_heatmap <- function(norm_file, som_file, order) {
  splits <- factor(som_file$som$unit.classif, levels = rev(order)) %>% na.omit()
  Heatmap(norm_file,
          name = "Scaled TPM",
          cluster_column_slices = FALSE,
          cluster_row_slices = FALSE,
          split = splits,
          show_row_dend = FALSE,
          show_row_names = FALSE,
          col = colorRamp2(0:3, cividis(5)[2:5]),
          show_column_names = FALSE,
          show_column_dend = FALSE,
          use_raster = TRUE,
          raster_by_magick = TRUE,
          cluster_columns = FALSE,
          show_heatmap_legend = FALSE)
}

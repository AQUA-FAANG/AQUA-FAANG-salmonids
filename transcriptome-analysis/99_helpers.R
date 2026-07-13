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

# === deposited-asset loaders (reproduce mode) ================================
# These read the small derived assets shipped in the repo and inject named
# objects into the CALLER's environment, exactly like base load() does for an
# .RData file. They are the single source of truth for "load the deposited
# SOMs / normalised matrices", used by BOTH 02_soms.Rmd and the ohnolog scripts
# under scripts/, so the two can never drift apart.

# load_deposited_norm — read the four deposited normalised expression matrices
# (data/norm/{sd,td,sb,tb}_norm.tsv.gz) and inject them into `envir` as the
# matrices sd_norm, td_norm, sb_norm, tb_norm. Each file is a gzipped TSV whose
# first column is the gene_id (-> matrix rownames) and whose remaining columns
# are the per-sample normalised values.
load_deposited_norm <- function(norm_dir = "data/norm", envir = parent.frame()) {
  # Read one tag's matrix: gene_id column -> rownames, the rest -> numeric body.
  read_one <- function(tag) {
    df <- readr::read_tsv(file.path(norm_dir, paste0(tag, "_norm.tsv.gz")),
                          show_col_types = FALSE)
    m <- as.matrix(df[, -1L])     # drop the gene_id column to get the numeric body
    rownames(m) <- df[[1L]]       # restore gene IDs as matrix row names
    m
  }
  # Inject sd_norm/td_norm/sb_norm/tb_norm into the caller's environment.
  for (tag in c("sd", "td", "sb", "tb")) {
    assign(paste0(tag, "_norm"), read_one(tag), envir = envir)
  }
  invisible(NULL)
}

# load_deposited_soms — reconstruct the four SOM objects from the deposited
# per-gene unit assignments (data/soms/{tag}_unit_classif.tsv.gz) and per-unit
# summary tables (data/soms/{tag}_unit_stats.tsv.gz + {tag}_stage_levels.txt),
# injecting sd_som, td_som, sb_som, tb_som into `envir`. Each reconstructed
# object is a list carrying exactly the two fields the rest of the pipeline
# reads: $som$unit.classif (integer SOM unit per gene) and $stat (the per-unit
# mean/sd summary tibble). The unit vector is aligned to the matching deposited
# norm matrix's row order so the downstream positional zip
# tibble(geneID = rownames(<norm>), cluster = <som>$som$unit.classif) is correct
# regardless of file ordering.
load_deposited_soms <- function(soms_dir = "data/soms",
                                norm_dir = "data/norm",
                                envir = parent.frame()) {
  for (tag in c("sd", "td", "sb", "tb")) {
    # Per-gene unit assignments: columns gene_id, unit.
    uc <- readr::read_tsv(file.path(soms_dir, paste0(tag, "_unit_classif.tsv.gz")),
                          show_col_types = FALSE)
    # Per-unit summary table: columns class, stage, mean_value, sd.
    st <- readr::read_tsv(file.path(soms_dir, paste0(tag, "_unit_stats.tsv.gz")),
                          show_col_types = FALSE)
    # The canonical stage/tissue factor order for this SOM's summary table.
    stage_levels <- readLines(file.path(soms_dir, paste0(tag, "_stage_levels.txt")))
    # The gene order of the matching norm matrix (its first/gene_id column only).
    norm_ids <- readr::read_tsv(file.path(norm_dir, paste0(tag, "_norm.tsv.gz")),
                                show_col_types = FALSE, col_select = 1L)[[1L]]
    # Reorder the unit assignments into norm's gene order via an explicit join.
    unit <- uc$unit[match(norm_ids, uc$gene_id)]
    # Rebuild the class factor in natural numeric order (1,2,...,K not 1,10,11).
    class_levels <- stringr::str_sort(unique(as.character(st$class)), numeric = TRUE)
    st$class <- factor(st$class, levels = class_levels)
    # Rebuild the stage factor in the deposited canonical order.
    st$stage <- factor(st$stage, levels = stage_levels)
    # Assemble the minimal SOM object and inject it as <tag>_som.
    assign(paste0(tag, "_som"),
           list(som  = list(unit.classif = as.integer(unit)),
                stat = tibble::as_tibble(st)),
           envir = envir)
  }
  invisible(NULL)
}

# === normalization ===========================================================

# merge_reps_devmap — collapse the per-replicate DevMap TPM columns to one
# mean-TPM column per stage, then reorder + rename the 14 columns into
# chronological `stage_names` order. Returns a gene x stage matrix.
merge_reps_devmap <- function(unmerged) {
  setDT(unmerged)

  unmerged <- unmerged[, `transcript_id(s)` := NULL]
  unmerged <- melt(unmerged,
                   id.vars = "gene_id",
                   variable.name = "sample",
                   value.name = "tpm")
  # Strip the trailing replicate tag (_R1/_R2/...) to recover the stage. Robust to
  # a sample prefix such as "AtlanticSalmon_RNA_" that itself contains "_R".
  unmerged[, stage := sub("_R[0-9]+$", "", as.character(sample))]

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

# merge_reps_bodymap — collapse the per-replicate BodyMap columns to one mean-TPM
# column per (Tissue x Maturity x Sex) sample, and rename the gonad samples to
# Ova / Testes by sex. Returns a gene x sample matrix.
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

# normalize — turn a merged (gene x sample) mean-TPM matrix into the SOM input:
# drop genes never expressed above 1 TPM in any sample, quantile-normalise
# across samples, then z-scale each gene (SD only, no centering). Returns a
# matrix with the surviving gene IDs as rownames.
normalize <- function(merged) {
  # Gene IDs of rows kept (any sample > 1 TPM).
  row_names <- merged %>% filter(if_any(everything(), ~ .x > 1))
  row_names <- rownames(row_names)

  norm <- merged %>%
    filter(if_any(everything(), ~ .x > 1)) %>%        # drop never-expressed genes
    data.matrix() %>%
    preprocessCore::normalize.quantiles() %>%          # quantile-normalise across samples
    t() %>% scale(center = FALSE, scale = TRUE) %>% t() # per-gene z-scale (SD only)

  rownames(norm) <- row_names                           # quantile.normalise drops names; restore
  colnames(norm) <- colnames(merged)
  norm
}

# === SOM fit =================================================================

# get_som — fit a DevMap SOM on `in_mat` (xdim x ydim hexagonal grid) and return
# a list(som, stat, plot): the kohonen object, a per-(unit, stage) mean/SD
# summary table, and a ribbon ggplot. Used in recompute mode (02_soms.Rmd).
get_som <- function(in_mat, xdim, ydim) {
  # Fit the SOM.
  som_obj <- kohonen::som(in_mat,
                          grid = somgrid(xdim, ydim, topo = "hexagonal"))

  # Per-unit facet labels "Class <id> (<n genes>)".
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

# get_som_bodymap — like get_som, but the summary's stage factor is ordered by
# `order_indices` (the canonical tissue x maturity x sex column order) so the
# BodyMap ribbon/heatmap columns read correctly.
get_som_bodymap <- function(in_mat, xdim, ydim, order_indices = NULL) {
  if (is.null(order_indices)) order_indices <- seq_len(ncol(in_mat))  # default: input order
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

# som_ids — return a list where element i is the vector of gene IDs assigned to
# SOM unit i (gene IDs taken from the norm matrix rownames in unit.classif order).
som_ids <- function(norm_file, som_file) {
  listid <- list()
  for (i in seq_along(unique(som_file$som$unit.classif))) {
    listid[[i]] <- rownames(norm_file[which(som_file$som$unit.classif == i), ])
  }
  listid
}

# som_file_info_d — per-unit annotation table for a DevMap SOM. For each of the
# 16 units it records: class_max (the peak stage), class_max_mean (that peak's
# mean value), and constitutive_class (TRUE if the unit's across-stage SD < 0.4,
# i.e. a flat/housekeeping unit). Rows are sorted constitutive -> peak stage ->
# peak mean, the manuscript display order.
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

# som_file_info_b — per-unit annotation table for a BodyMap SOM (36 units). Like
# som_file_info_d but also adds: a `stage_specific` flag (peak mean > 2.165; see
# inline note) and a `tissue` column (the peak sample's tissue prefix). Sorted
# constitutive -> stage_specific -> peak tissue -> peak mean.
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

# palette_devmap — inferno colours for the stage-specific DevMap units + grey for
# the constitutive ones, then reordered so palette[i] is the colour for raw SOM
# unit i (the ggplot scales key on raw unit id).
palette_devmap <- function(n_stage_specific, n_constitutive, som_order) {
  stage_specific_cols <- inferno(n_stage_specific, begin = 0.15, end = 0.9)
  constitutive <- c("#becacaff")
  final_colors <- c(stage_specific_cols, rep(constitutive, n_constitutive))

  color_order <- match(1:16, som_order)   # invert som_order -> key by raw unit id
  final_colors[color_order]
}

# palette_bodymap — map each BodyMap SOM unit to a colour by its annotation:
# a distinct hue per peak tissue for stage-specific units, one shared "multi"
# blue for non-specific non-constitutive units, and grey for constitutive units.
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

# make_heatmap — ComplexHeatmap of a normalised expression matrix, rows split by
# SOM unit in the given `order` (reversed so the first unit sits at the top), on
# a cividis scaled-TPM colour ramp.
make_heatmap <- function(norm_file, som_file, order) {
  # Row split factor: SOM unit per gene, in reversed display order.
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

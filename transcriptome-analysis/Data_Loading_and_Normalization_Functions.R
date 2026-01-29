process_tissue_data <- function(fish) {
  file_path <- paste0("../data/bodymap/", fish)
  raw_files <- list.files(path = file_path, 
                          pattern = "gene_tpm.tsv", 
                          full.names = TRUE)
  
  temp <- lapply(raw_files, fread, sep ="\t")
  nams <- vector(mode = "character", length = 84)
  nseq <- seq(1, 73, 12)
  
  raw_data <- data.frame(matrix(
    nrow = nrow(temp[[1]]),
    ncol = 84
  ))
  
  row.names(raw_data) <- temp[[1]]$gene_id
  
  foreach(i = 1:length(temp), l = nseq) %do% {
    raw_data[,l:(l+11)] <- temp[[i]][,3:14]
    nams[l:(l+11)] <- colnames(temp[[i]][,3:14])
  }
  
  colnames(raw_data) <- nams
  
  return(raw_data)
}

merge_reps_devmap <- function(unmerged) {
  # Convert to data.table
  setDT(unmerged)
  
  # Perform the required transformations before summarise
  unmerged <- unmerged[, `transcript_id(s)` := NULL]
  unmerged <- melt(unmerged, 
                   id.vars = "gene_id", 
                   variable.name = "sample", 
                   value.name = "tpm")
  unmerged[, c("stage", "replicate") := tstrsplit(sample, "_R", fixed = TRUE)]
  
  # Efficient summarisation with data.table
  merged <- unmerged[, .(tpm = mean(tpm, na.rm = TRUE)), 
                     by = .(gene_id, stage)]
  
  # Convert back to tibble
  merged <- as_tibble(merged)
  
  merged <- merged %>%
    pivot_wider(names_from = "stage", # put back into wide format
                values_from = "tpm") %>%
    column_to_rownames(var = "gene_id") # gene_id col to rownames
  
  merged <- merged[,c(12,13,14,6,1,7,2,8,3,9,4,10,5,11)]
  colnames(merged) <- stage_names
  unmerged <- merged
}


merge_reps_bodymap <- function(unmerged) {
  
  # Rownames to column
  merged <- unmerged %>% rownames_to_column("gene_id")
  
  # Convert to data.table
  setDT(merged)
  
  # Perform the initial transformations
  merged <- melt(merged, 
                 id.vars = "gene_id", 
                 variable.name = "sample", 
                 value.name = "tpm")
  
  merged[, c("Species", "Seq", "Tissue", "Maturity", "Sex", "Replicate") := 
           tstrsplit(sample, "_", fixed = TRUE)]
  
  merged[, sample := NULL]
  merged[, Name := paste(Tissue, Maturity, Sex, sep = "_")]
  
  # Summarize using data.table
  merged <- merged[, .(tpm = mean(tpm, na.rm = TRUE)), by = .(gene_id, Name)]
  
  # Convert back to tibble for further manipulation
  merged <- as_tibble(merged)
  
  # Pivot wider and rename columns using dplyr
  unmerged <- merged %>%
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
    data.matrix() %>% normalize.quantiles() %>%
    t() %>% scale(center = F, scale = T) %>% t()
  
  rownames(norm) <- row_names
  colnames(norm) <- colnames(merged)
  
  merged <- norm
}

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
    mutate(stage = factor(stage,
                          levels = colnames(in_mat)),
           class = factor(class,
                          levels = str_sort(unique(class), numeric = T)))
  
  sum_stat <- tmp_df %>% group_by(class, stage) %>%
    summarize(mean_value = mean(value), 
              sd = sd(value), 
              .groups = 'drop') %>% ungroup()
  
  out_plot <- ggplot(sum_stat, aes(x = stage, y = mean_value, group = class)) +
    geom_line(size = 1.5) +   # added color aesthetic and size
    geom_ribbon(aes(ymin = mean_value - 2*sd, 
                    ymax = mean_value + 2*sd), 
                alpha = 0.2) +  # plot points with the same color as the line
    facet_wrap(~ class, ncol = xdim, labeller = labeller(class = label_unit),
               scale = "free_y", as.table = F) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), 
          legend.position = "none")
  
  list(som = som_obj,
       stat = sum_stat,
       plot = out_plot)
}

get_som_bodymap <- function(in_mat, xdim, ydim, order_indices = NULL) {
  
  if (is.null(order_indices)) {
    order_indices <- 
      1:ncol(in_mat) # Default to the original order if none is provided.
  }
  
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
           class = factor(class, levels = str_sort(unique(class), numeric = T)))
  
  sum_stat <- tmp_df %>% group_by(class, stage) %>%
    summarize(mean_value = mean(value), sd = sd(value), .groups = 'drop') %>% 
    ungroup()
  
  out_plot <- ggplot(sum_stat, aes(x = stage, y = mean_value, group = class)) +
    geom_line(size = 1.5) + 
    geom_ribbon(aes(ymin = mean_value - 2*sd, 
                    ymax = mean_value + 2*sd), alpha = 0.2) +
    facet_wrap(~ class, ncol = xdim, labeller = labeller(class = label_unit),
               scale = "free_y", as.table = F) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), 
          legend.position = "none")
  
  list(som = som_obj, stat = sum_stat, plot = out_plot)
}


som_ids <- function(norm_file, som_file) {
  listid <- list()
  for (i in 1:length(unique(som_file$som$unit.classif))) {
    listid[[i]] <- rownames(norm_file[which(som_file$som$unit.classif == i),])
  }
  norm_file <- listid
}

som_file_info_d <- function(som_file) {
  class_max <- vector()
  constitutive_class <- vector()
  class_max_mean <- vector()
  class_number <- c(1:16)
  
  for (i in 1:length(unique(som_file$som$unit.classif))) {
    class_max[i] <- 
      stage_names[which.max(som_file$stat$mean_value[which(
        som_file$stat$class == i)])]
    class_max_mean[i] <- 
      som_file$stat$mean_value[
        which.max(
          som_file$stat$mean_value[
            which(som_file$stat$class == i)
          ]
        ) + (length(stage_names) * (i - 1))
      ]
    
    if (sd(som_file$stat$mean_value[which(som_file$stat$class == i)]) < 0.4) {
      constitutive_class[i] <- "TRUE"
    }
    else {constitutive_class[i] <- "FALSE"}
  }
  
  class_info <- 
    as.data.frame(cbind(class_number,
                        class_max,
                        class_max_mean,
                        constitutive_class))
  
  class_info$class_max <- factor(class_info$class_max,levels=stage_names)
  
  som_file <- arrange(class_info,constitutive_class,class_max,
                      desc(class_max_mean))
}

som_file_info_b <- function(som_file) {
  class_max <- vector()
  constitutive_class <- vector()
  stage_specific <- vector()
  class_max_mean <- vector()
  class_number <- c(1:36)
  tissue_names <- as.vector(unique(som_file$stat$stage))
  
  for (i in 1:length(unique(som_file$som$unit.classif))) {
    class_max[i] <- 
      tissue_names[which.max(som_file$stat$mean_value[which(
        som_file$stat$class == i)])]
    class_max_mean[i] <- 
      som_file$stat$mean_value[
        which.max(
          som_file$stat$mean_value[
            which(som_file$stat$class == i)
          ]
        ) + (length(unique(som_file$stat$stage)) * (i - 1))
      ]
    
    if (sd(som_file$stat$mean_value[which(som_file$stat$class == i)]) < 0.4) {
      constitutive_class[i] <- "TRUE"
    }
    else {constitutive_class[i] <- "FALSE"}
    
    if (class_max_mean[i] > 2.165) {
      stage_specific[i] <- "TRUE"
    }
    else {stage_specific[i] <- "FALSE"}
  }
  
  class_info <- 
    as.data.frame(cbind(class_number,
                        class_max,
                        class_max_mean,
                        constitutive_class,
                        stage_specific))
  
  class_info$class_max <- factor(class_info$class_max,levels=tissue_names)
  class_info$tissue <- sub("\\_.*", "", class_info$class_max)
  
  som_file <- arrange(class_info,constitutive_class,desc(stage_specific),class_max,
                      desc(class_max_mean))
}

palette_devmap <- function(n_stage_specific,n_constitutive,som_order) {
  stage_specific_cols <- inferno(n_stage_specific,begin = 0.15, end = 0.9)
  constitutive <- c("#becacaff")
  final_colors <- c(stage_specific_cols, rep(constitutive,n_constitutive))
  
  order <- som_order
  color_order <- match(1:16,order)
  ordered_palette <- final_colors[color_order]
  return(ordered_palette)
}

palette_bodymap <- function(som_info) {
  brain_col <- "#dce648ff"
  testes_col <- "forestgreen"
  ova_col <- "#ffa300ff"
  muscle_col <- "salmon"
  liver_col <- "#8a0000ff"
  head_kidney_col <- "black"
  gill_col <- "red2"
  intestine_col <-   "hotpink"
  multi_col <- "skyblue"
  constitutive_col <- "#becacaff"
  
  pal <- vector()
  class_num <- som_info$class_number
  
  for (i in 1:nrow(som_info)) {
    if (som_info$tissue[i] == "Brain" & som_info$stage_specific[i] == TRUE) {
      pal[i] <- brain_col
    }
    else if (som_info$tissue[i] == "Testes" & som_info$stage_specific[i] == TRUE) {
      pal[i] <- testes_col
    }
    else if (som_info$tissue[i] == "Ova" & som_info$stage_specific[i] == TRUE) {
      pal[i] <- ova_col
    }
    else if (som_info$tissue[i] == "Muscle" & som_info$stage_specific[i] == TRUE) {
      pal[i] <- muscle_col
    }
    else if (som_info$tissue[i] == "Liver" & som_info$stage_specific[i] == TRUE) {
      pal[i] <- liver_col
    }
    else if (som_info$tissue[i] == "HeadKidney" & som_info$stage_specific[i] == TRUE) {
      pal[i] <- head_kidney_col
    }
    else if (som_info$tissue[i] == "Gill" & som_info$stage_specific[i] == TRUE) {
      pal[i] <- gill_col
    }
    else if (som_info$tissue[i] == "SuppDistalIntestine" & som_info$stage_specific[i] == TRUE) {
      pal[i] <- intestine_col
    }
    else if (som_info$stage_specific[i] == FALSE & som_info$constitutive_class[i] == FALSE) {
      pal[i] <- multi_col
    }
    else if (som_info$constitutive_class[i] == TRUE) {
      pal[i] <- constitutive_col
    }
  }
  
  return(pal)
}

  
  
  
make_heatmap <- function(norm_file,som_file,order) {
  splits <- factor(som_file$som$unit.classif,levels = rev(order)) %>% na.omit()
  ht <- Heatmap(norm_file,
                name = "Scaled TPM",
                cluster_column_slices = F,
                cluster_row_slices = F,
                split = splits,
                show_row_dend = F,
                show_row_names = F,
                col = colorRamp2(0:3,cividis(5)[2:5]),
                show_column_names = F,
                show_column_dend = F,
                use_raster = T,
                raster_by_magick = T,
                cluster_columns = F,
                show_heatmap_legend = F)
  norm_file <- ht
}



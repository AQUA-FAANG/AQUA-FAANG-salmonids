get_som <- function(in_mat, xdim, ydim, somvar) {
  
  som_obj <- somvar
  
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

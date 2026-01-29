library(circlize)
library(tidyverse)
library(viridis)
library(ComplexHeatmap)
library(data.table)
library(gtable)

sda_mat <- devMap_in_mat
sda_som <- devMap_atac.som
sda_umap <- devMap_umap_atac

pal <- colorRamp2(c(0:2),cividis(3)[0:3])
splits_levels <- rev(c("25","24","14","20",
                   "10","5","4","8",
                   "2","1","6","7",
                   "11","17","16","22","21"))

sda_heatmap <- Heatmap(sda_mat[which(sda_som$unit.classif %in% splits_levels),],
        cluster_rows = T,
        row_split = factor(sda_som$unit.classif[which(sda_som$unit.classif %in% splits_levels)],
                           levels = splits_levels),
        show_row_dend = F,
        cluster_columns = F,
        cluster_row_slices = F,
        show_column_names = T,
        col = pal,
        name = "Normalized peak intensity",
        raster_by_magick = T,
        show_heatmap_legend = F,
        use_raster = T
      )

 



splits_levels_c <- c(splits_levels,c("3","9","13","12","15","18","19","23"))

umap_palette_quantnorm <- inferno(17,begin = 0.285, end = 0.627)
constitutive <- rep(c("#becacaff"),8)
final_colors <- c(umap_palette_quantnorm, constitutive)
palette_matrix <- matrix(1:25, nrow = 1) # Creating a matrix of the 25 numbers
image(palette_matrix, col = final_colors, axes = FALSE, 
      main = "Custom Color Palette")
color_order_for_objects <- match(1:25, splits_levels_c)
ordered_palette_salmon_devmap <- final_colors[color_order_for_objects]

sda_umap <- devMap_umap_atac
sda_som <- devMap_atac.som

sda_umap <- sda_umap %>% as.data.frame() %>% 
  ggplot(aes(V1, V2*(-1))) + #geom_density_2d() + 
  geom_point(aes(colour = factor(sda_som$unit.classif)), 
             alpha = .4, size = .1) + 
  scale_color_manual(values = ordered_palette_salmon_devmap) +
  theme_classic() +# guides(colour = guide_legend(override.aes = list(size=10)))
  theme(legend.position = "none",
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank()) 
#################################################################################
colnames_bodymap <- colnames(bodyMap_in_mat)
colnames_bodymap[5] <- "Ovary_Inmature_Female_ATAC"
colnames_bodymap[6] <- "Testes_Inmature_Male_ATAC"
colnames_bodymap[7] <- "Ovary_Mature_Female_ATAC"
colnames_bodymap[8] <- "Testes_Mature_Male_ATAC"
colnames(bodyMap_in_mat) <- colnames_bodymap
bodymap_in_mat_reorder <- bodyMap_in_mat[,c(1:4,6,8,5,7,9:16)]


splits_levels_b <- rev(c("24","25","19","13","8","3","17","5","1","6","7","23","22"))
splits_levels_b_c <- c(splits_levels_b,c("20","2","14","15","12","11","18","16","21","9","10","12"))

sba_heatmap <- Heatmap(bodymap_in_mat_reorder[which(bodyMap_atac.som$unit.classif %in% splits_levels_b),],
        cluster_rows = T,
        row_split = factor(bodyMap_atac.som$unit.classif[which(bodyMap_atac.som$unit.classif %in% splits_levels_b)],
                           levels = splits_levels_b),
        show_row_dend = F,
        cluster_columns = F,
        cluster_row_slices = F,
        show_column_names = T,
        col = pal,
        name = "Normalized peak intensity"
)

ht <- Heatmap(bodymap_in_mat_reorder,
                       cluster_rows = T,
                       row_split = factor(bodyMap_atac.som$unit.classif,levels = c(1:25)),
                       show_row_dend = F,
                       cluster_columns = F,
                       cluster_row_slices = F,
                       show_column_names = F,
                       col = pal,
                       name = "Normalized peak intensity"
)



bodymap_palette_salmon <- c("#8a0000ff","skyblue","forestgreen","skyblue","#ffa300ff",
                     "#8a0000ff","#8a0000ff","forestgreen","skyblue","skyblue",
                     "#becacaff","#becacaff","#dce648ff","#becacaff","#becacaff",
                     "skyblue","#ffa300ff","skyblue","#dce648ff","#becacaff",
                     "salmon","salmon","salmon","#dce648ff","#dce648ff")

color_order_for_objects_b <- match(1:25, splits_levels_b_c)
ordered_palette_salmon_bodymap <- bodymap_palette[color_order_for_objects_b]

atac_salmon <- bodyMap_umap_atac

sba_umap <- bodyMap_umap_atac %>% as.data.frame() %>% 
  ggplot(aes(V1, V2)) + #geom_density_2d() + 
  geom_point(aes(colour = factor(bodyMap_atac.som$unit.classif)), 
             alpha = .4, size = .1) + 
  scale_color_manual(values = bodymap_palette_salmon) +
  theme_classic() +# guides(colour = guide_legend(override.aes = list(size=10)))
  theme(legend.position = "none",
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank()) 
  
g1 <- sba_umap
  stages_cols <- viridis(14)
  
  plot(NULL, xlim=c(0,length(stages_cols)), ylim=c(0,1), 
       +      xlab="", ylab="", xaxt="n", yaxt="n") + rect(0:(length(stages_cols)-1), 0, 1:length(stages_cols), 1, col=stages_cols)

  

  
   salmonatacdev <- devMap_tmp_df %>%
    ggplot(aes(stage, value,
               group = "Class")) +
    stat_summary(fun.data  = "mean_sdl",
                 fun.args = list(mult = 2), geom = "ribbon") + 
    stat_summary(fun = mean, geom = "line") +
    facet_wrap(~ class, ncol = 5, 
               labeller = labeller(class = devMap_label_unit),
               as.table = F) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

    
   troutatacdev <- devMap_tmp_df %>%
     ggplot(aes(stage, value,
                group = "Class")) +
     stat_summary(fun.data  = "mean_sdl",
                  fun.args = list(mult = 2), geom = "ribbon") + 
     stat_summary(fun = mean, geom = "line") +
     facet_wrap(~ class, ncol = 5, 
                labeller = labeller(class = devMap_label_unit),
                as.table = F) +
     theme_bw() +
     theme(axis.text.x = element_text(angle = 45, hjust = 1))
   
   salmonatacbody <- bodyMap_tmp_df %>%
     ggplot(aes(stage, value,
                group = "Class")) +
     stat_summary(fun.data  = "mean_sdl",
                  fun.args = list(mult = 2), geom = "ribbon") + 
     stat_summary(fun = mean, geom = "line") +
     facet_wrap(~ class, ncol = 5, 
                labeller = labeller(class = devMap_label_unit),
                as.table = F) +
     theme_bw() +
     theme(axis.text.x = element_text(angle = 45, hjust = 1))
   
   troutatacbody <- bodyMap_tmp_df %>%
     ggplot(aes(stage, value,
                group = "Class")) +
     stat_summary(fun.data  = "mean_sdl",
                  fun.args = list(mult = 2), geom = "ribbon") + 
     stat_summary(fun = mean, geom = "line") +
     facet_wrap(~ class, ncol = 5, 
                labeller = labeller(class = devMap_label_unit),
                as.table = F) +
     theme_bw() +
     theme(axis.text.x = element_text(angle = 45, hjust = 1))
   
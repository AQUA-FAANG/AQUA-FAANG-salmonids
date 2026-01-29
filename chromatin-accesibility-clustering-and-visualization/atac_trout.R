library(circlize)
library(tidyverse)
library(viridis)
library(ComplexHeatmap)
library(gridExtra)



pal_tr <- colorRamp2(c(0:2),cividis(3)[0:3])

splits_levels <- rev(c("25","15","14","20","24",
                       "23","18","5",
                       "3","2",
                       "19","1",
                       "7","6","21","17","16"))

Heatmap(devMap_in_mat,
        cluster_rows = T,
        row_split = devMap_atac.som$unit.classif,
        show_row_dend = F,
        cluster_columns = F,
        cluster_row_slices = F,
        show_column_names = F,
        col = pal,
        name = "Normalized peak intensity",
        raster_by_magick = T
)


tda_heatmap <- Heatmap(devMap_in_mat[which(devMap_atac.som$unit.classif %in% splits_levels),],
        cluster_rows = T,
        row_split = factor(devMap_atac.som$unit.classif[which(devMap_atac.som$unit.classif %in% splits_levels)],
                           levels = splits_levels),
        show_row_dend = F,
        cluster_columns = F,
        cluster_row_slices = F,
        show_column_names = F,
        col = pal,
        name = "Normalized peak intensity",
        raster_by_magick = T,
        show_heatmap_legend = T,
        use_raster = T
)



splits_levels_c <- c(splits_levels,c("8","9","10","11","12","13","22","4"))

umap_palette_quantnorm <- inferno(17,begin = 0.285, end = 0.627)
constitutive <- rep(c("#becacaff"),8)
final_colors <- c(umap_palette_quantnorm, constitutive)
palette_matrix <- matrix(1:25, nrow = 1) # Creating a matrix of the 25 numbers
image(palette_matrix, col = final_colors, axes = FALSE, 
      main = "Custom Color Palette")
color_order_for_objects <- match(1:25, splits_levels_c)
ordered_palette_salmon_devmap <- final_colors[color_order_for_objects]


tda_umap <- devMap_umap_atac %>% as.data.frame() %>% 
  ggplot(aes(V1, V2*(-1))) + #geom_density_2d() + 
  geom_point(aes(colour = factor(devMap_atac.som$unit.classif)), 
             alpha = .4, size = .1) + 
  scale_color_manual(values = ordered_palette_salmon_devmap) +
  theme_classic() +# guides(colour = guide_legend(override.aes = list(size=10)))
  theme(legend.position = "none",
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank()) 


sda_heatmap_specific_grob <- grid.grabExpr(draw(sda_heatmap))
tda_heatmap_specific_grob <- grid.grabExpr(draw(tda_heatmap))

atac_heatmaps <- grid.arrange(sda_heatmap_specific_grob,tda_heatmap_specific_grob,nrow = 1)
atac_umaps <- grid.arrange(sda_umap,tda_umap,nrow = 1)
#################################################################################
colnames_bodymap <- colnames(bodyMap_in_mat)
colnames(bodyMap_in_mat) <- colnames_bodymap
bodymap_in_mat_reorder <- bodyMap_in_mat[,c(1:4,5:8,9:12)]


splits_levels_b <- rev(c("25","20","21","16","17","2","3"))
splits_levels_b_c <- c(splits_levels_b,c("4","9","11","12","13","14","15","16","20"))


tba_heatmap <- Heatmap(bodymap_in_mat_reorder[which(bodyMap_atac.som$unit.classif %in% splits_levels_b),],
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
                       row_split = bodyMap_atac.som$unit.classif,
                       show_row_dend = F,
                       cluster_columns = T,
                       cluster_row_slices = F,
                       show_column_names = T,
                       col = pal,
                       name = "Normalized peak intensity"
)

bodymap_palette_trout <- c("skyblue","salmon","salmon","skyblue","skyblue",
                     "skyblue","skyblue","skyblue","skyblue","skyblue",
                     "skyblue","skyblue","#becacaff","#becacaff","#becacaff",
                     "#8a0000ff","#8a0000ff","skyblue","skyblue","#dce648ff",
                     "#8a0000ff","#becacaff","#becacaff","#becacaff","#dce648ff")

color_order_for_objects_b <- match(1:25, splits_levels_b_c)
ordered_palette_salmon_bodymap <- bodymap_palette[color_order_for_objects_b]



tba_umap <- bodyMap_umap_atac %>% as.data.frame() %>% 
  ggplot(aes(V1, V2)) + #geom_density_2d() + 
  geom_point(aes(colour = factor(bodyMap_atac.som$unit.classif)), 
             alpha = 0.4, size = 0.1) + 
  scale_color_manual(values = bodymap_palette_trout) +
  theme_classic() +# guides(colour = guide_legend(override.aes = list(size=10)))
  theme(legend.position = "none",
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank()) 

g2 <- tba_umap

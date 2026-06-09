#######################################
###TE overlap across Active enhancers
#######################################

library(tidyverse)
library(RColorBrewer)

#load data

salmon_TE_overlap<-readRDS("enhancer_TE_overlap_by_group_redip_salmon.RDS")

# Define colors
blue_shades <- brewer.pal(3, "Blues")[3:1]
red_shades <- brewer.pal(3, "Reds")[3:1]

# Map conservation categories to colors
fill_colors_aore<- c(
  "exclusive" = blue_shades[3],
  "alignable" = blue_shades[2],
  "shared" = blue_shades[1]
)

gap_df <- expand.grid(
  group = "gap",
  redip = unique(salmon_TE_overlap$redip),
  enhancer_type=unique(salmon_TE_overlap$enhancer_type)) %>%
  mutate(prop = 0) 

salmon_TE_overlap_aug <- bind_rows(salmon_TE_overlap, gap_df)

group_levels <- c("LateBlastulation", "MidGastrulation", "EarlySomitogenesis", "MidSomitogenesis", "LateSomitogenesis","gap",  "Brain_Immature", "Brain_Mature",
                  "Testis_Immature", "Testis_Mature", "Ovary_Immature", "Ovary_Mature", "Liver_Immature", "Liver_Mature", "Muscle_Immature", "Muscle_Mature")    

enhancer_type_levels <- c("exclusive","alignable", "shared")



##### merge them and plot

fill_colors_all<- c(
  "exclusive_LORe" = red_shades[3],
  "alignable_LORe" = red_shades[2],
  "shared_LORe" = red_shades[1],
  "exclusive_AORe" = blue_shades[3],
  "alignable_AORe" = blue_shades[2],
  "shared_AORe" = blue_shades[1]
)


# Custom labels for facet strip
facet_labels <- c(
  "AORe" = "Early rediploidization",
  "LORe" = "Late rediploidization"
)

type_redip_levels<-c("exclusive_AORe","exclusive_LORe","alignable_AORe",
              "alignable_LORe","shared_AORe","shared_LORe")

salmon_TE_overlap_aug %>%
  unite(c("enhancer_type","redip"),col="type_redip",sep="_",remove = F) %>%
  mutate(group = factor(group, levels = group_levels)) %>%
  mutate(type_redip = factor(type_redip, levels = type_redip_levels)) %>%
  ggplot(aes(x = group, y = TE_proportion, fill =type_redip )) +
  geom_bar(position = "stack", stat = "identity", color = "black", linewidth = 0.07, width = 0.85) +
  facet_wrap(~ redip, nrow = 1, strip.position = "top",labeller = as_labeller(facet_labels)) +
  theme_classic(base_size = 9)+
  theme(
    strip.background = element_rect(fill = "white", color = "black"),
    panel.border = element_rect(color = "black", fill = NA),
    axis.text.x = element_text(vjust = 0.5, size = 6),
    axis.text.y = element_text(size = 6),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 8),
    legend.position = "none"
  )+
  scale_y_continuous(limits = c(0, 0.5), breaks = seq(0, 0.5,0.1))+
  scale_fill_manual(values = fill_colors_all) +
  scale_color_manual(values = fill_colors_all) +
  theme(legend.position = "none") +
  scale_x_discrete(labels=c("LB","MG",'ES','MS', 'LS','', 'Br\n(I)','Br\n(M)','Mu\n(I)','Mu\n(M)','Li\n(I)','Li\n(M)', 
                            'Ov\n(I)','Ov\n(M)','Te\n(I)','Te\n(M)'))+
  ylab(label = "Active enhancer overlap with repeats")
  

ggsave("F4_TE_overlap_active_enhancers.tiff", dpi = 600, width = 11.9, height = 6, units = "cm")

ggsave(filename = "F4_TE_overlap_active_enhancers.pdf",  plot = last_plot(), device = cairo_pdf, width = 11.9, height = 6, dpi = 600,units = "cm")

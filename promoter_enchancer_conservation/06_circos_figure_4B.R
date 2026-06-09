library(circlize)
library(tidyverse)
library(openxlsx)
library(readxl)
library(colorRamp2)
library(ComplexHeatmap)

circos.clear()

cytoband_df <- read.xlsx(
  "salmon_data.xlsx",
  sheet = "Sheet1",
  colNames = FALSE
) %>%
  mutate(
    X1 = factor(as.character(X1), levels = unique(X1))
  )

link_df_1 <- read.xlsx("salmon_data.xlsx", sheet = "Sheet5")
link_df_2 <- read.xlsx("salmon_data.xlsx", sheet = "Sheet3")

identity_df <- read_excel(
  "one_mb_identity.xlsx",
  col_types = c("numeric", "numeric", "numeric", "text")
) %>%
  mutate(
    Interval_end = as.numeric(gsub(",", "", Interval_end)),
    chromosome   = as.numeric(str_extract(Chr, "\\d+"))
  ) %>%
  arrange(chromosome, Interval_end) %>%
  group_by(chromosome) %>%
  mutate(
    Interval_start = lag(Interval_end, default = 0)
  ) %>%
  ungroup() %>%
  select(
    chromosome,
    Interval_start,
    Interval_end,
    Ident
  )

identity_df <- as.data.frame(identity_df)

# initialise circos plot

circos.clear()

circos.par(
  start.degree = 90,
  gap.degree = 1.5,
  cell.padding = c(0.001, 0, 0.001, 0),
  canvas.xlim = c(-1.2, 1.2),
  canvas.ylim = c(-1.2, 1.2)
)

circos.initializeWithIdeogram(
  cytoband_df,
  plotType = NULL
)

# identity track

identity_col_fun <- colorRamp2(
  c(80, 85, 90, 95, 100),
  c(
    "green",
    "darkgreen",
    "yellow3",
    "tomato2",
    "darkred"
  )
)

circos.genomicTrack(
  identity_df,
  ylim = range(identity_df$Ident),
  panel.fun = function(region, value, ...) {
    
    circos.genomicLines(
      region,
      value,
      type = "h",
      col = identity_col_fun(value[[1]]),
      border = NA,
      ...
    )
  },
  track.height = 0.09,
  track.margin = c(0, 0),
  bg.border = "lightgrey"
)

# chromosome labels

circos.genomicTrack(
  cytoband_df,
  ylim = c(0, 1),
  panel.fun = function(region, value, ...) {
    
    circos.genomicText(
      region,
      value,
      y = 0.5,
      labels.column = 1,
      cex = 0.4,
      niceFacing = TRUE,
      ...
    )
  },
  track.height = 0.02,
  track.margin = c(0.005, 0.005),
  bg.border = NA
)

# ideogram

circos.genomicIdeogram(
  cytoband = cytoband_df,
  track.height = 0.04
)

# homologous links
circos.genomicLink(
  region1 = link_df_1,
  region2 = link_df_2,
  col = adjustcolor(link_df_1$colour, alpha.f = 0.6),
  border = adjustcolor(link_df_1$colour, alpha.f = 0.7),
  lwd = 0.1
)

#legend
identity_legend <- Legend(
  title = "Identity (%)",
  col_fun = identity_col_fun,
  at = seq(80, 100, by = 5),
  direction = "vertical",
  legend_width = unit(4, "cm"),
  title_gp = gpar(fontsize = 10, fontface = "bold"),
  labels_gp = gpar(fontsize = 9)
)

pdf(
  "identity_legend.pdf",
  width = 6,
  height = 2
)

draw(
  identity_legend,
  x = unit(0.5, "npc"),
  y = unit(0.5, "npc"),
  just = c("center", "center")
)

dev.off()

#Identity statistics for AORe and LORe windows

link_regions <- link_df_1 %>%
  rename(
    chromosome = Chr,
    start_link = Start,
    end_link = End
  )

overlap_identity <- identity_df %>%
  rename(
    start_ident = Interval_start,
    end_ident = Interval_end
  ) %>%
  inner_join(link_regions, by = "chromosome") %>%
  filter(
    start_ident < end_link,
    end_ident > start_link
  )

identity_summary_sig <- overlap_identity %>%
  group_by(colour) %>%
  summarise(
    mean_ident = mean(Ident),
    sd_ident = sd(Ident),
    n = n(),
    .groups = "drop"
  )

identity_summary_sig

#Genome proportion occupied by LORe

GENOME_SIZE_BP <- 2673458927

salmon_LORe <- read.table(
  "../Ssal_late_rediploidized_regions.tsv",
  header = TRUE
)

LORe_summary <- salmon_LORe %>%
  mutate(
    window_size = end - start
  ) %>%
  summarise(
    total_bp = sum(window_size)
  ) %>%
  mutate(
    genome_bp = GENOME_SIZE_BP,
    proportion = total_bp / genome_bp
  )

LORe_summary

#Identity in LORe vs AORe

LORe_regions <- salmon_LORe %>%
  rename(
    chromosome = chrom,
    start_LORe = start,
    end_LORe = end
  )

identity_windows <- identity_df %>%
  rename(
    start_ident = Interval_start,
    end_ident = Interval_end
  )

#Lore
LORe_overlap <- identity_windows %>%
  inner_join(
    LORe_regions,
    by = "chromosome"
  ) %>%
  filter(
    start_ident < end_LORe,
    end_ident > start_LORe
  ) %>%
  mutate(
    overlap_status = "LORe"
  )

#AORe

AORe_windows <- identity_windows %>%
  anti_join(
    LORe_overlap,
    by = c(
      "chromosome",
      "start_ident",
      "end_ident",
      "Ident"
    )
  ) %>%
  mutate(
    overlap_status = "AORe"
  )

#Summary
identity_redip_summary <- bind_rows(
  LORe_overlap,
  AORe_windows
) %>%
  group_by(overlap_status) %>%
  summarise(
    mean_ident = mean(Ident),
    median_ident = median(Ident),
    sd_ident = sd(Ident),
    n = n(),
    .groups = "drop"
  )

identity_redip_summary





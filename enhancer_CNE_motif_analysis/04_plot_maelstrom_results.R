############################################################
# GIMME MAELSTROM – MOTIF ENRICHMENT VISUALISATION
#
# Kept analyses:
#   1. Age
#   2. Age × Conservation interaction
#
#
############################################################

# ============================================================
# LIBRARIES
# ============================================================

library(tidyverse)
library(readxl)
library(writexl)
library(stringr)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

# ============================================================
# HELPER FUNCTIONS
# ============================================================

extract_family <- function(motif_id) {
  str_extract(motif_id, "(?<=GM\\.5\\.0\\.)([^.]+)")
}

get_top_motifs <- function(df, cols, top_n = 60) {
  
  motif_col <- colnames(df)[1]
  colnames(df)[1] <- "motif"
  
  lapply(cols, function(col) {
    
    df %>%
      arrange(desc(.data[[col]])) %>%
      slice_head(n = top_n) %>%
      mutate(top_category = col,
             rank = row_number())
  }) %>% bind_rows()
}

get_top10_nonredundant <- function(df) {
  
  df <- df %>%
    mutate(family = extract_family(motif))
  
  df %>%
    group_by(top_category) %>%
    group_modify(~ {
      
      .x %>%
        group_by(family) %>%
        slice_min(rank, n = 1, with_ties = FALSE) %>%
        ungroup() %>%
        slice_min(rank, n = 10, with_ties = FALSE) %>%
        mutate(all_motifs = paste(unique(motif), collapse = " | "),
               n_candidates = n_distinct(motif))
    }) %>%
    ungroup() %>%
    mutate(row_label = paste0(family, " | ", top_category))
}

# ============================================================
# AGE ANALYSIS
# ============================================================

df_age <- read_excel("age.xlsx")

age_order <- c(
  "Euteleostomi", "Actinopteri", "Neopterygii",
  "Teleostei", "Clupeocephala",
  "Euteleosteomorpha", "Protacanthopterygii"
)

age_labels <- c("Et", "Ac", "Ne", "Te", "Cl", "Eu", "Pr")

df_age <- df_age %>%
  select(motif = 1, all_of(age_order))

colnames(df_age) <- c("motif", age_labels)

top_age <- get_top_motifs(df_age, age_labels)
top10_age <- get_top10_nonredundant(top_age)

write_xlsx(top10_age,
           "maelstrom_age_top10_families.xlsx")

# ============================================================
# AGE × CONSERVATION ANALYSIS
# ============================================================

df_int <- read_excel("age_cons_int.xlsx")

cons <- c("shared", "alignable", "exclusive")

age_full <- c(
  "Euteleostomi","Actinopteri","Neopterygii",
  "Teleostei","Clupeocephala",
  "Euteleosteomorpha","Protacanthopterygii"
)

age_short <- c("Et","Ac","Ne","Te","Cl","Eu","Pr")

col_order <- outer(cons, age_full, paste, sep = "_")
col_labels <- outer(cons, age_short, paste, sep = "_")

existing <- col_order[col_order %in% colnames(df_int)]
labels   <- col_labels[col_order %in% colnames(df_int)]

df_int <- df_int %>%
  select(motif = 1, all_of(existing))

colnames(df_int) <- c("motif", labels)

top_int <- get_top_motifs(df_int, labels)
top10_int <- get_top10_nonredundant(top_int)

write_xlsx(top10_int,
           "maelstrom_age_conservation_top10.xlsx")

# ============================================================
# AGE HEATMAP (ComplexHeatmap)
# ============================================================

plot_heatmap <- function(df, title, file) {
  
  mat <- df %>%
    select(-c(top_category, rank, family, motif)) %>%
    column_to_rownames("row_label") %>%
    as.matrix()
  
  max_abs <- max(abs(mat), na.rm = TRUE)
  
  col_fun <- colorRamp2(
    c(-max_abs, 0, max_abs),
    c("#2166AC", "white", "#B2182B")
  )
  
  ht <- Heatmap(
    mat,
    col = col_fun,
    name = "Z-score",
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_names_side = "right",
    show_row_names = TRUE,
    column_title = title
  )
  
  png(file, width = 5000, height = 3500, res = 200)
  draw(ht)
  dev.off()
}

plot_heatmap(
  top10_age,
  "Top 10 Motifs – Age",
  "maelstrom_age_heatmap.png"
)

plot_heatmap(
  top10_int,
  "Top 10 Motifs – Age × Conservation",
  "maelstrom_age_conservation_heatmap.png"
)

# ============================================================
# EXPORT FINAL TABLES
# ============================================================

write_xlsx(top10_age, "FINAL_age_top10.xlsx")
write_xlsx(top10_int, "FINAL_age_conservation_top10.xlsx")

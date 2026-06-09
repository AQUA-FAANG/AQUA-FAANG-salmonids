############################################################
## CNE MAF processing + conservation + percent identity
## Input: MAF alignments per CNE category
## Output: curated alignment metrics + identity estimates + summary stats + supplementary figures
############################################################

library(tidyverse)
library(fuzzyjoin)
library(openxlsx)
library(tidygenomics)
library(janitor)
library(pbapply)
library(stringr)
library(broom)
library(purrr)

############################################################
## 1. MAF PARSER (basic block extractor)
############################################################

MAF2tbl <- function(mafFile) {
  
  x <- readLines(mafFile)
  
  # MAF blocks start with "a", alignment lines start with "s"
  aLines <- grepl("^a", x)
  sLines <- grepl("^s", x)
  
  tibble(idx = cumsum(aLines), txt = x) %>%
    filter(sLines) %>%
    separate(
      txt,
      into = c("s", "src", "start", "size", "strand", "srcSize", "text"),
      sep = "[ \t]+",
      fill = "right"
    ) %>%
    select(-s) %>%
    mutate(across(c(start, size, srcSize), as.integer))
}

############################################################
## 2. Load MAF files for each CNE category
############################################################

shared_files    <- list.files("./shared_CNEs/", pattern = "\\.maf$", full.names = TRUE)
exclusive_files <- list.files("./exclusive_CNEs/", pattern = "\\.maf$", full.names = TRUE)
alignable_files <- list.files("./alignable_CNEs/", pattern = "\\.maf$", full.names = TRUE)

shared_maf_list    <- setNames(pblapply(shared_files, MAF2tbl), shared_files)
exclusive_maf_list <- setNames(pblapply(exclusive_files, MAF2tbl), exclusive_files)
alignable_maf_list <- setNames(pblapply(alignable_files, MAF2tbl), alignable_files)

############################################################
## 3. CURATION FUNCTION (Salmon + Eluc alignment stats)
############################################################

curate_maf <- function(maf_table) {
  
  maf_table %>%
    mutate(
      n_nucleotides = str_count(text, "[ACGTacgtN]"),
      n_gaps        = str_count(text, "-")
    ) %>%
    separate(src, into = c("species", "chromosome"), sep = "\\.", remove = FALSE) %>%
    filter(species %in% c("Ssal_A", "Ssal_B", "Eluc")) %>%
    mutate(
      end = ifelse(strand == "+", start + size, start + size),
      
      # genome-relative coordinates (strand-aware approximation)
      pos_start = ifelse(strand == "+", start, srcSize - (start + size)),
      pos_end   = ifelse(strand == "+", start + size, srcSize - start),
      
      total_sites = n_nucleotides + n_gaps
    ) %>%
    group_by(species, chromosome) %>%
    summarise(
      n_nucleotides = sum(n_nucleotides),
      n_gaps        = sum(n_gaps),
      total_sites   = sum(total_sites),
      
      pos_start = min(pos_start),
      pos_end   = max(pos_end),
      
      .groups = "drop"
    ) %>%
    mutate(
      conservation = n_nucleotides / total_sites,
      span_length  = pos_end - pos_start
    )
}

############################################################
## 4. Run curation across all MAF files
############################################################

curate_all_maf <- function(maf_list) {
  out <- pblapply(maf_list, curate_maf)
  names(out) <- names(maf_list)
  out
}

shared_curated    <- curate_all_maf(shared_maf_list)
exclusive_curated <- curate_all_maf(exclusive_maf_list)
alignable_curated <- curate_all_maf(alignable_maf_list)

shared_CNEs_ssal_coords <- bind_rows(shared_curated, .id = "Origin") %>%
  mutate(Origin = basename(Origin))

exclusive_CNEs_ssal_coords <- bind_rows(exclusive_curated, .id = "Origin") %>%
  mutate(Origin = basename(Origin))

alignable_CNEs_ssal_coords <- bind_rows(alignable_curated, .id = "Origin") %>%
  mutate(Origin = basename(Origin))

############################################################
## 5. QUALITY CHECK: alignment completeness vs expected size
############################################################

check_alignment_quality <- function(df, threshold = 1) {
  
  df %>%
    mutate(cne_info = str_extract(Origin, "(?<=\\.maf_).*")) %>%
    separate(cne_info, into = c("chr", "start", "end"), sep = "_", convert = TRUE) %>%
    mutate(
      end = as.integer(str_remove(end, "\\.maf$")),
      true_length = end - start,
      coverage_ratio = total_sites / true_length
    ) %>%
    filter(coverage_ratio < threshold | conservation < 1) %>%
    distinct(Origin)
}

shared_CNEs_ssal_coords %>% check_alignment_quality(0.95)
exclusive_CNEs_ssal_coords %>% check_alignment_quality(0.95)
alignable_CNEs_ssal_coords %>% check_alignment_quality(0.95)

############################################################
## 6. Percent identity calculation (Eluc vs Ssal)
############################################################

compute_percent_identity <- function(df,
                                     eluc = "Eluc",
                                     ssal = c("Ssal_A", "Ssal_B")) {
  
  stitched <- df %>%
    separate(src, into = c("species", "chromosome"), sep = "\\.", remove = FALSE) %>%
    filter(species %in% c(eluc, ssal)) %>%
    group_by(species) %>%
    summarise(full_seq = paste0(toupper(text), collapse = ""), .groups = "drop")
  
  if (!eluc %in% stitched$species) {
    return(tibble())
  }
  
  eluc_seq <- strsplit(stitched %>% filter(species == eluc) %>% pull(full_seq), "")[[1]]
  
  # Eluc vs Ssal comparison
  out <- stitched %>%
    filter(species %in% ssal) %>%
    rowwise() %>%
    mutate(
      other_seq = list(strsplit(full_seq, "")[[1]]),
      aligned_len = min(length(other_seq), length(eluc_seq)),
      
      matches = sum(other_seq[1:aligned_len] == eluc_seq[1:aligned_len]),
      mismatches = aligned_len - matches,
      
      percent_identity = 100 * matches / aligned_len,
      seq_length = aligned_len
    ) %>%
    ungroup() %>%
    select(
      Eluc = eluc,
      Ssal_ID = species,
      percent_identity,
      seq_length,
      mismatches
    )
  
  # Ssal_A vs Ssal_B identity (optional summary metric)
  if (all(ssal %in% stitched$species)) {
    a <- strsplit(stitched %>% filter(species == ssal[1]) %>% pull(full_seq), "")[[1]]
    b <- strsplit(stitched %>% filter(species == ssal[2]) %>% pull(full_seq), "")[[1]]
    
    L <- min(length(a), length(b))
    out$Ssal_identity <- 100 * sum(a[1:L] == b[1:L]) / L
  }
  
  out
}

############################################################
## 7. Run identity analysis
############################################################

shared_pid    <- pblapply(shared_maf_list, compute_percent_identity)
exclusive_pid <- pblapply(exclusive_maf_list, compute_percent_identity)
alignable_pid <- pblapply(alignable_maf_list, compute_percent_identity)

shared_CNEs_percent_identity <- bind_rows(shared_pid, .id = "Origin") %>%
  mutate(Origin = basename(Origin))

exclusive_CNEs_percent_identity <- bind_rows(exclusive_pid, .id = "Origin") %>%
  mutate(Origin = basename(Origin))

alignable_CNEs_percent_identity <- bind_rows(alignable_pid, .id = "Origin") %>%
  mutate(Origin = basename(Origin))

############################################################
## 8. Merge identity back to curated CNE tables
############################################################

attach_identity <- function(cne_df, pid_df) {
  
  cne_df %>%
    left_join(pid_df, by = c("Origin", "species" = "Ssal_ID")) %>%
    mutate(
      cne_info = str_extract(Origin, "(?<=\\.maf_).*")
    ) %>%
    separate(cne_info, into = c("chr", "start", "end"), sep = "_", convert = TRUE) %>%
    mutate(
      end = as.integer(str_remove(end, "\\.maf$")),
      true_length = end - start,
      alignment_capture = total_sites / true_length
    )
}

shared_final    <- attach_identity(shared_CNEs_ssal_coords, shared_CNEs_percent_identity)
exclusive_final <- attach_identity(exclusive_CNEs_ssal_coords, exclusive_CNEs_percent_identity)
alignable_final <- attach_identity(alignable_CNEs_ssal_coords, alignable_CNEs_percent_identity)

############################################################
## 9. Final QC filters
############################################################

shared_final %>%
  filter(alignment_capture < 0.90 | conservation < 1 | percent_identity < 60) %>%
  distinct(Origin)

exclusive_final %>%
  filter(alignment_capture < 0.90 | conservation < 1 | percent_identity < 60) %>%
  distinct(Origin)

alignable_final %>%
  filter(alignment_capture < 0.90 | conservation < 1 | percent_identity < 60) %>%
  distinct(Origin)

############################################################
## 10. Load original Pike CNE annotation table
############################################################

pike_cne_file <- "../Pike_CNE_details.txt"

pike_cne_coords <- read.table(pike_cne_file, header = TRUE)

############################################################
## 11. Attach Pike annotation to curated CNE tables
############################################################

attach_pike_annotation <- function(df) {
  df %>%
    left_join(
      pike_cne_coords,
      by = c(
        "cne_chr" = "chromosome",
        "cne_start" = "start",
        "cne_end" = "end"
      )
    )
}

shared_CNEs_annotated    <- attach_pike_annotation(shared_CNEs_2)
alignable_CNEs_annotated <- attach_pike_annotation(alignable_CNEs_2)
exclusive_CNEs_annotated <- attach_pike_annotation(exclusive_CNEs_2)

############################################################
## 12. SUMMARY STATISTICS (by age + species)
############################################################

summarise_cne_set <- function(df, category_label, species_filter = NULL) {
  
  df %>%
    group_by(age, species) %>%
    summarise(
      mean_seq_length  = mean(seq_length, na.rm = TRUE),
      mean_identity    = mean(percent_identity, na.rm = TRUE),
      sd_identity      = sd(percent_identity, na.rm = TRUE),
      mean_conservation = mean(conservation, na.rm = TRUE),
      counts = n(),
      .groups = "drop"
    ) %>%
    { if (!is.null(species_filter)) filter(., species %in% species_filter) else . } %>%
    mutate(category = category_label)
}

shared_summary <- summarise_cne_set(shared_CNEs_annotated, "Shared")
alignable_summary <- summarise_cne_set(alignable_CNEs_annotated, "Alignable")
exclusive_summary <- summarise_cne_set(
  exclusive_CNEs_annotated,
  "Exclusive",
  species_filter = "Ssal_A"
)

cne_age_summary <- bind_rows(shared_summary, alignable_summary, exclusive_summary)

write.xlsx(cne_age_summary, "CNE_age_summary.xlsx")

############################################################
## 13. FILTERING: high-confidence conserved CNEs
############################################################

filter_cnes <- function(df, category_label, min_identity = 80) {
  
  df %>%
    filter(species != "Eluc") %>%
    filter(cne_alignment_capture > 0.90, conservation == 1) %>%
    group_by(Origin) %>%
    filter(
      all(c("Ssal_A", "Ssal_B") %in% unique(species)),
      any(percent_identity > min_identity)
    ) %>%
    ungroup() %>%
    mutate(category = category_label)
}

shared_CNEs_filtered <- filter_cnes(shared_CNEs_annotated, "shared")
alignable_CNEs_filtered <- filter_cnes(alignable_CNEs_annotated, "alignable")

exclusive_CNEs_filtered <- exclusive_CNEs_annotated %>%
  filter(species == "Ssal_A") %>%
  filter(cne_alignment_capture > 0.90, conservation == 1) %>%
  group_by(Origin) %>%
  filter(all("Ssal_A" %in% unique(species))) %>%
  ungroup() %>%
  mutate(category = "exclusive")

############################################################
## 14. EXPORT FILTERED SETS
############################################################

write.table(shared_CNEs_filtered, "shared_CNEs_post_filtering.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(alignable_CNEs_filtered, "alignable_CNEs_post_filtering.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(exclusive_CNEs_filtered, "exclusive_CNEs_post_filtering.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

############################################################
## 15. MERGE ALL FILTERED CNEs
############################################################

total_cne_post_filtering <- bind_rows(
  shared_CNEs_filtered,
  alignable_CNEs_filtered,
  exclusive_CNEs_filtered
)

write.table(total_cne_post_filtering,
            "total_CNEs_post_filtering.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

############################################################
## 16. Compute Ssal_A vs Ssal_B identity ratio
############################################################

compute_identity_ratio <- function(df) {
  
  ratio_tbl <- df %>%
    select(Origin, species, percent_identity) %>%
    pivot_wider(names_from = species, values_from = percent_identity) %>%
    mutate(identity_ratio = pmin(Ssal_A, Ssal_B) / pmax(Ssal_A, Ssal_B)) %>%
    select(Origin, identity_ratio)
  
  df %>%
    left_join(ratio_tbl, by = "Origin") %>%
    mutate(identity_ratio = if_else(species == "Ssal_A", identity_ratio, NA_real_))
}

shared_CNEs_filtered <- compute_identity_ratio(shared_CNEs_filtered)
alignable_CNEs_filtered <- compute_identity_ratio(alignable_CNEs_filtered)

############################################################
## 17. ADD SALMON COORDINATES + SIZE CONSISTENCY CHECK
############################################################

add_salmon_overlap_check <- function(df, salmon_coords) {
  
  df %>%
    left_join(
      salmon_coords %>%
        select(
          Origin, species, chromosome,
          n_nucleotides, n_gaps,
          nucleotides_gaps_total,
          conservation,
          pos_start, pos_end, pos_length
        ),
      by = c(
        "Origin", "species", "chromosome",
        "n_nucleotides", "n_gaps",
        "nucleotides_gaps_total",
        "conservation"
      )
    ) %>%
    mutate(true_vs_salmon_ratio = cne_true_length / pos_length)
}

shared_CNE_final    <- add_salmon_overlap_check(shared_CNEs_filtered, shared_CNEs_ssal_coords)
alignable_CNE_final <- add_salmon_overlap_check(alignable_CNEs_filtered, alignable_CNEs_ssal_coords)
exclusive_CNE_final <- add_salmon_overlap_check(exclusive_CNEs_filtered, exclusive_CNEs_ssal_coords)

############################################################
## 18. FILTER BY COORDINATE CONSISTENCY (important QC step)
############################################################

filter_coordinate_consistency <- function(df) {
  df %>%
    group_by(Origin) %>%
    filter(all(true_vs_salmon_ratio >= 0.8 & true_vs_salmon_ratio <= 1.2)) %>%
    ungroup()
}

shared_CNE_final <- filter_coordinate_consistency(shared_CNE_final)
alignable_CNE_final <- filter_coordinate_consistency(alignable_CNE_final)
exclusive_CNE_final <- filter_coordinate_consistency(exclusive_CNE_final)

############################################################
## 19. EXPORT BED-LIKE FILES FOR GIMME ANALYSIS
############################################################

prepare_bed <- function(df) {
  df %>%
    select(chromosome, pos_start, pos_end, age, category)
}

all_cne_salmon_coords <- bind_rows(
  prepare_bed(shared_CNE_final),
  prepare_bed(alignable_CNE_final),
  prepare_bed(exclusive_CNE_final)
)

write.table(all_cne_salmon_coords,
            "salmon_CNE_list_raw_input_gimme.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

############################################################
## 20. GIMME MAELSTROM INPUT FORMAT
############################################################

gimme_maelstrom_input <- all_cne_salmon_coords %>%
  mutate(
    loc = paste0(chromosome, ":", pos_start, "-", pos_end),
    cluster = age
  ) %>%
  select(loc, cluster)

write.table(gimme_maelstrom_input,
            "gimme_maelstrom_input.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

############################################################
## 21. Age ordering (used across all plots)
############################################################

age_order <- c(
  "Euteleostomi", "Actinopteri", "Neopterygii", "Teleostei",
  "Clupeocephala", "Euteleosteomorpha", "Protacanthopterygii"
)

############################################################
## 22. GLOBAL SUMMARY: counts + size per category & age
############################################################

cne_age_summary <- bind_rows(shared_CNEs_4, alignable_CNEs_4, exclusive_CNEs_4) %>%
  distinct(Origin, .keep_all = TRUE) %>%
  group_by(age, category) %>%
  summarise(
    n = n(),
    total_length = sum(cne_true_length, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = category,
    values_from = c(n, total_length),
    names_sep = "_"
  )

write.xlsx(
  cne_age_summary,
  "cne_age_summary_post_filtering_counts_length.xlsx"
)

############################################################
## 23. HELPER: base dataset for identity comparisons
############################################################

identity_df <- bind_rows(shared_CNEs_4, alignable_CNEs_4)

############################################################
## 24. OVERALL Ssal identity comparison (shared vs alignable)
############################################################

identity_df %>%
  filter(!is.na(Ssal_identity)) %>%
  mutate(category = factor(category, levels = c("shared", "alignable"))) %>%
  ggplot(aes(x = category, y = Ssal_identity, fill = category)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.8) +
  scale_fill_viridis_d(option = "E", end = 0.85) +
  scale_y_continuous(limits = c(70, 103), breaks = seq(70, 100, 10)) +
  theme_classic(base_size = 12) +
  labs(
    x = "Conservation category",
    y = "Ssal duplicate identity (%)"
  ) +
  theme(
    panel.border = element_rect(color = "black", fill = NA),
    axis.text = element_text(size = 7),
    axis.title = element_text(size = 8),
    legend.position = "none"
  )

ggsave("F6_shared_alignable_Ssal_identity_boxplot.jpeg",
       dpi = 600, width = 8, height = 8, units = "cm")

############################################################
## 25. SUMMARY STATISTICS (overall)
############################################################

identity_df %>%
  filter(!is.na(Ssal_identity)) %>%
  group_by(category) %>%
  summarise(
    median = median(Ssal_identity),
    mean   = mean(Ssal_identity),
    .groups = "drop"
  )

############################################################
## 26. STATISTICAL TEST: shared vs alignable
############################################################

wilcox.test(Ssal_identity ~ category, data = identity_df %>% filter(!is.na(Ssal_identity)))

############################################################
## 27. AGE-SPLIT BOXPLOT (Ssal identity)
############################################################

identity_df %>%
  filter(!is.na(Ssal_identity)) %>%
  mutate(
    category = factor(category, levels = c("shared", "alignable")),
    age = factor(age, levels = age_order)
  ) %>%
  ggplot(aes(x = age, y = Ssal_identity, fill = category)) +
  geom_boxplot(
    position = position_dodge(0.7),
    outlier.shape = NA,
    width = 0.6
  ) +
  scale_fill_viridis_d(option = "E", end = 0.85) +
  scale_y_continuous(limits = c(70, 101), breaks = seq(70, 100, 10)) +
  theme_classic(base_size = 12) +
  labs(
    x = "Age group",
    y = "Ssal duplicate identity (%)"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    legend.position = "right"
  )

ggsave("F6_identity_by_age_boxplot.jpeg",
       dpi = 600, width = 14, height = 8, units = "cm")

############################################################
## 28. MEDIAN + MEAN BY AGE
############################################################

identity_df %>%
  filter(!is.na(Ssal_identity)) %>%
  group_by(age, category) %>%
  summarise(
    median = median(Ssal_identity),
    mean   = mean(Ssal_identity),
    .groups = "drop"
  )

############################################################
## 29. WILCOX TEST BY AGE (shared vs alignable)
############################################################

wilcox_by_age_identity <- identity_df %>%
  filter(!is.na(Ssal_identity)) %>%
  group_by(age) %>%
  filter(n_distinct(category) == 2) %>%
  summarise(
    test = list(wilcox.test(Ssal_identity ~ category)),
    .groups = "drop"
  ) %>%
  mutate(tidy = map(test, broom::tidy)) %>%
  unnest(tidy) %>%
  transmute(
    age,
    p.value,
    statistic,
    p.adj = p.adjust(p.value, method = "fdr"),
    significant = p.adj < 0.05
  )

print(wilcox_by_age_identity)

############################################################
## 30. Pike vs Ssal_A/Ssal_B identity (alignable only)
############################################################

alignable_CNEs_4 %>%
  filter(!is.na(percent_identity)) %>%
  mutate(species = factor(species, levels = c("Ssal_A", "Ssal_B"))) %>%
  ggplot(aes(x = species, y = percent_identity, fill = species)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.8) +
  scale_fill_viridis_d(option = "E", end = 0.85) +
  scale_y_continuous(limits = c(60, 101), breaks = seq(60, 100, 10)) +
  theme_classic(base_size = 12) +
  labs(
    x = "Ssal reference",
    y = "Percent identity vs Pike"
  ) +
  theme(legend.position = "none")

ggsave("F6_alignable_pike_identity.jpeg",
       dpi = 600, width = 8, height = 8, units = "cm")

wilcox.test(percent_identity ~ species,
            data = alignable_CNEs_4 %>% filter(!is.na(percent_identity)))

############################################################
## 31. IDENTITY RATIO (Ssal_A vs Ssal_B)
############################################################

identity_ratio_df <- identity_df %>%
  filter(!is.na(identity_ratio))

identity_ratio_df %>%
  ggplot(aes(x = category, y = identity_ratio, fill = category)) +
  geom_boxplot(outlier.shape = 1) +
  scale_fill_viridis_d(option = "E", end = 0.85) +
  scale_y_continuous(limits = c(0.8, 1.01)) +
  theme_classic(base_size = 12) +
  labs(
    x = "Category",
    y = "Ssal identity ratio"
  )

ggsave("F6_identity_ratio_boxplot.jpeg",
       dpi = 600, width = 8, height = 8, units = "cm")

############################################################
## 32. FINAL SUMMARY TABLE (all metrics)
############################################################

pid_median <- identity_df %>%
  filter(!is.na(percent_identity)) %>%
  group_by(age, category) %>%
  summarise(
    median = median(percent_identity),
    mean   = mean(percent_identity),
    .groups = "drop"
  )


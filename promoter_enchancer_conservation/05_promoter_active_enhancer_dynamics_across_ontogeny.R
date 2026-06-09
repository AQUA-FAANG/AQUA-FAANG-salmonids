############################################################
# TITLE: ATAC promoter conservation across ontogeny
# PURPOSE: Clean integration of ATAC peak conservation,
#          promoter activity, and rediploidization context
############################################################

############################
# 1. LIBRARIES
############################

library(tidyverse)
library(fuzzyjoin)
library(openxlsx)
library(tidygenomics)
library(RColorBrewer)
library(rstatix)
library(stringr)

############################
# 2. DATA INPUT : all this data is output of the previous script
############################
# NOTE FOR MANUSCRIPT:
# Describe each dataset briefly in README or methods

two_two_atac_peaks <- read.table(
  "salmon_atac_conservation_2to2_summit_cutoff_2024.txt",
  header = TRUE
)

two_two_peaks_summary <- read.table(
  "salmon_two_way_conservation_data_summit_overlap_2024.txt",
  header = TRUE
) %>% mutate(across(everything(), as.character))

one_one_atac_peaks <- read.table(
  "salmon_atac_conservation_1to1_summit_cutoff_2024.txt",
  header = TRUE
)

one_one_peaks_summary <- read.table(
  "salmon_one_way_conservation_data_2024.txt",
  header = TRUE
) %>% mutate(across(everything(), as.character))

two_one_atac_peaks <- read.table(
  "salmon_atac_conservation_2to1_2024.txt",
  header = TRUE
)

two_one_peaks_summary <- read.table(
  "salmon_twoseq_onepeak_conservation_data_2024.txt",
  header = TRUE
) %>% mutate(across(everything(), as.character))


############################
# 3. STANDARDISE SAMPLE METADATA
############################

atac_numbers <- read.table("ATAC_numbers_across_tissues_salmon.txt", header = TRUE) %>%
  mutate(stage = case_when(
    stage == "gonad_immature_male" ~ "testis_immature_male",
    stage == "gonad_immature_female" ~ "ovary_immature_female",
    stage == "gonad_mature_female" ~ "ovary_mature_female",
    stage == "gonad_mature_male" ~ "testis_mature_male",
    TRUE ~ stage
  )) %>%
  mutate(stage = as.character(stage))

############################
# 4. PROMOTER ANNOTATION
############################
# NOTE:
# This block defines promoter universe used across all analyses

active_prom <- read.table(
  "../../atac_merged_gareth/redo_salmon_new_annotation_2024/Active_TSS_1.bed",
  col.names = c("chromosome", "start", "end", "details", "peaksize", "strand")
) %>%
  unite("Origin", chromosome, start, end, sep = "_", remove = FALSE) %>%
  select(-peaksize, -strand) %>%
  mutate(
    chromosome = sprintf("%02d", chromosome),
    chromosome = paste0("Ssal", chromosome),
    Origin = paste0("ssa", Origin, ".maf"),
    details = as.character(details)
  )

# expand tissue annotations
active_prom_2 <- active_prom %>%
  mutate(details = strsplit(details, ",")) %>%
  unnest(details) %>%
  separate(details, into = c("tissue", "chrom"), sep = ":") %>%
  filter(chrom == "Active_TSS")

# harmonise naming
rep_str <- c(
  "Gonad_Immature_Male" = "testis_immature_male",
  "Gonad_Immature_Female" = "ovary_immature_female",
  "Gonad_Mature_Female" = "ovary_mature_female",
  "Gonad_Mature_Male" = "testis_mature_male"
)

active_prom_2 <- active_prom_2 %>%
  mutate(
    tissue = str_replace_all(tissue, rep_str),
    tissue = tolower(tissue)
  )

write.table(
  active_prom_2,
  "active_promoters_across_ontogeny_2024.bed",
  quote = FALSE, row.names = FALSE
)

############################
# 5. HELPER FUNCTION 
############################
# NOTE: to process peaks and conservation

process_peak_set <- function(atac_df, summary_df, type_label) {
  
  peak_merge <- atac_df %>%
    select(1:12) %>%
    left_join(
      active_prom_2,
      by = c(
        "chromosome.y" = "chromosome",
        "start.y" = "start",
        "end.y" = "end"
      )
    ) %>%
    filter(!is.na(Origin.y)) %>%
    mutate(tissue = tolower(tissue))
  
  filtered <- summary_df %>%
    left_join(
      peak_merge,
      by = c("Origin" = "Origin.x", "stage" = "tissue")
    ) %>%
    filter(!is.na(chromosome.x)) %>%
    select(-peaks_open)
  
  counts <- filtered %>%
    count(Origin, stage, name = "n")
  
  list(
    filtered = filtered,
    counts = counts,
    type = type_label
  )
}

############################
# 6. APPLY PIPELINE
############################

two_two <- process_peak_set(two_two_atac_peaks, two_two_peaks_summary, "shared")
one_one <- process_peak_set(one_one_atac_peaks, one_one_peaks_summary, "exclusive")
two_one <- process_peak_set(two_one_atac_peaks, two_one_peaks_summary, "alignable")

############################
# 7. CLASSIFY PROMOTERS
############################

dup_prom  <- filter(two_two$counts, n == 2)
sing_prom <- filter(one_one$counts, n == 1)
twoone_prom <- filter(two_one$counts, n == 1)

write.table(dup_prom, "duplicated_promoters_across_ontogeny_2024.bed", row.names = FALSE)
write.table(sing_prom, "singleton_promoters_across_ontogeny_2024.bed", row.names = FALSE)
write.table(twoone_prom, "two_one_promoters_across_ontogeny_2024.bed", row.names = FALSE)

############################
# 8. DETAILED TABLES
############################

dup_prom_detailed <- dup_prom %>%
  distinct(Origin) %>%
  left_join(two_two_atac_peaks, by = "Origin")

sing_prom_detailed <- sing_prom %>%
  distinct(Origin) %>%
  left_join(one_one_atac_peaks, by = "Origin")

twoone_prom_detailed <- twoone_prom %>%
  distinct(Origin) %>%
  left_join(two_one_atac_peaks, by = "Origin")

############################
# 9. LORe / AORe ANNOTATION - rediploidization
############################

salmon_LORe <- read.table("./Ssal_late_rediploidized_regions.tsv", header = TRUE) %>%
  rename(chromosome = chrom) %>%
  mutate(
    chromosome = sprintf("%02d", chromosome),
    chromosome = paste0("ssa", chromosome),
    start = as.numeric(start),
    end = as.numeric(end)
  )

# Combine promoter sets
promoters_list <- bind_rows(
  dup_prom,
  sing_prom,
  twoone_prom
) %>%
  mutate(
    Origin = str_remove(Origin, "\\.maf$")
  ) %>%
  separate(Origin, into = c("chromosome", "start", "end"), sep = "_", convert = TRUE)

promoters_list <- fuzzy_left_join(
  promoters_list,
  salmon_LORe,
  by = c("chromosome", "start", "end"),
  match_fun = list(`==`, `>=`, `<=`)
) %>%
  mutate(redip = if_else(!is.na(start.y), "LORe", "AORe")) %>%
  select(-ends_with(".y")) %>%
  rename_with(~ str_remove(., "\\.x$"))

############################
# 10. SUMMARY OUTPUT
############################

prom_activity_redip <- promoters_list %>%
  count(stage, redip)

write.xlsx(prom_activity_redip, "promoter_activity_rediploidization_summary.xlsx")

############################################################
# 11. LORe / AORe ANNOTATION PIPELINE (REFRACTORED)
############################################################
# This section assigns promoter sets to rediploidization regions
# using genomic interval overlap (fuzzy matching)

############################
# helper function
############################

annotate_lore_aore <- function(df, label) {
  
  df %>%
    group_by(Origin, stage) %>%
    summarise(n = n(), .groups = "drop") %>%
    mutate(Origin = str_remove(Origin, "\\.maf$")) %>%
    separate(Origin, into = c("chromosome", "start", "end"),
             sep = "_", convert = TRUE) %>%
    mutate(chromosome = sprintf("ssa%02d", as.integer(str_remove(chromosome, "ssa")))) %>%
    
    fuzzy_left_join(
      salmon_LORe,
      by = c("chromosome", "start", "end"),
      match_fun = list(`==`, `>=`, `<=`)
    ) %>%
    
    mutate(
      redip = if_else(!is.na(start.y), "LORe", "AORe"),
      set = label
    ) %>%
    select(-ends_with(".y")) %>%
    rename_with(~ str_remove(., "\\.x$"))
}

############################
# apply to all promoter classes
############################

prom_filter_counts_lore_aore <- annotate_lore_aore(prom_filter, "shared")
prom_filter_counts_sing_lore_aore <- annotate_lore_aore(prom_filter_sing, "exclusive")
prom_filter_counts_twoone_lore_aore <- annotate_lore_aore(prom_filter_twoone, "alignable")

############################################################
# 12. PROMOTER BREADTH OF ACTIVITY
############################################################

del4 <- prom_filter_counts_lore_aore %>%
  filter(n == 2) %>%
  mutate(category = "shared") %>%
  slice(rep(1:n(), each = 2)) %>%
  group_by(chromosome, start, end, stage) %>%
  mutate(copy_id = row_number()) %>%
  ungroup()

del5 <- prom_filter_counts_sing_lore_aore %>%
  filter(n == 1) %>%
  mutate(category = "exclusive", copy_id = 1)

del6 <- prom_filter_counts_twoone_lore_aore %>%
  filter(n == 1) %>%
  mutate(category = "alignable", copy_id = 1)

del7 <- bind_rows(del4, del5, del6) %>%
  separate(stage, into = c("tissue", "maturity", "sex"),
           sep = "_", remove = FALSE) %>%
  mutate(across(everything(), ~ replace_na(.x, "embryo"))) %>%
  unite(col = "tissue", tissue, maturity, sep = "_")

############################
# breadth metric
############################

breadth <- del7 %>%
  group_by(redip, chromosome, start, end, copy_id, category) %>%
  summarise(n_tissues = n_distinct(tissue), .groups = "drop")

############################
# summary stats
############################

breadth_summary <- breadth %>%
  group_by(redip, category) %>%
  summarise(
    mean_tissues = mean(n_tissues),
    median_tissues = median(n_tissues),
    sd_tissues = sd(n_tissues),
    n_peaks = n(),
    .groups = "drop"
  )

############################################################
# 13. STATISTICAL TESTING (Wilcoxon)
############################################################

breadth_wilcox <- breadth %>%
  group_by(redip) %>%
  wilcox_test(n_tissues ~ category, p.adjust.method = "BH") %>%
  add_significance() %>%
  select(redip, group1, group2, statistic, p, p.adj, p.adj.signif)

############################################################
# 14. PLOTS (BREADTH DISTRIBUTION)
############################################################

cat_cols <- c(
  "shared" = "#1b9e77",
  "alignable" = "#d95f02",
  "exclusive" = "#7570b3"
)

p_hist <- ggplot(breadth, aes(x = n_tissues, fill = category)) +
  geom_histogram(binwidth = 1, position = "dodge", color = "black", alpha = 0.85) +
  scale_fill_manual(values = cat_cols) +
  facet_wrap(~ redip) +
  labs(x = "Number of tissues active",
       y = "Number of promoters",
       fill = "Category") +
  theme_bw(base_size = 12) +
  theme(legend.position = "top")

ggsave("promoter_breadth_histogram.pdf", p_hist, width = 18, height = 10)

p_box <- ggplot(breadth, aes(x = category, y = n_tissues, fill = category)) +
  geom_boxplot(outlier.colour = "black") +
  scale_fill_manual(values = cat_cols) +
  facet_wrap(~ redip) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

ggsave("promoter_breadth_boxplot.pdf", p_box, width = 18, height = 10)

############################################################
# 15. TISSUE-SPECIFIC COMBINATIONS
############################################################

tissue_combos <- del7 %>%
  group_by(redip, chromosome, start, end, copy_id, category) %>%
  summarise(
    tissue_combo = str_c(sort(unique(tissue)), collapse = ","),
    .groups = "drop"
  )

combo_summary <- tissue_combos %>%
  group_by(redip, category, tissue_combo) %>%
  summarise(n_peaks = n(), .groups = "drop")

single_tissue_peaks <- combo_summary %>%
  filter(!str_detect(tissue_combo, ","))

############################################################
# 16. PLOT: TISSUE-SPECIFIC PROMOTERS
############################################################

stage_order <- c(
  "lateblastulation_embryo","midgastrulation_embryo","earlysomitogenesis_embryo",
  "midsomitogenesis_embryo","latesomitogenesis_embryo","gap",
  "brain_immature","brain_mature","muscle_immature","muscle_mature",
  "liver_immature","liver_mature","ovary_immature","ovary_mature",
  "testis_immature","testis_mature"
)

plot_data <- combo_summary %>%
  filter(!str_detect(tissue_combo, ",")) %>%
  mutate(
    redip = recode(redip,
                   AORe = "Early rediploidization",
                   LORe = "Late rediploidization"),
    tissue_combo = factor(tissue_combo, levels = stage_order),
    category = factor(category, levels = c("shared", "alignable", "exclusive"))
  )

p_tissue <- ggplot(plot_data,
                   aes(x = tissue_combo, y = n_peaks, fill = category)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = cat_cols) +
  facet_wrap(~ redip, scales = "free_x") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 70, hjust = 1))

ggsave("tissue_specific_promoters.pdf", p_tissue, width = 25, height = 12)


################
dups_prom_activity<-prom_filter_counts_lore_aore %>%
  filter(n == 2) %>%
  select(c("stage","redip")) %>%
  count(stage,redip) %>%
  mutate(nos = n*2 ) %>% ## account for duplication here
  select(-c(n))

sing_prom_activity<- prom_filter_counts_sing_lore_aore %>%
  filter(n == 1) %>%
  select(c("stage","redip")) %>%
  count(stage,redip) %>%
  mutate(nos = n)%>%
  select(-c(n))

twoone_prom_activity<-prom_filter_counts_twoone_lore_aore %>%
  filter(n == 1) %>%
  select(c("stage","redip")) %>%
  count(stage,redip) %>%
  mutate(nos = n)%>%
  select(-c(n))

prom_redip_activity<- bind_rows(dups_prom_activity,sing_prom_activity,twoone_prom_activity) %>% group_by(stage,redip) %>%
  summarise(n=sum(nos))

##### get proportions now

dups_conserved_props_prom_lore_aore<-prom_filter_counts_lore_aore %>%
  filter(n == 2) %>%
  select(c("stage","redip")) %>%
  count(stage,redip) %>%
  mutate(nos = n*2 ) %>%
  left_join(prom_redip_activity, by = c("stage","redip")) %>%
  mutate(perc_conserved = round((nos / n.y) * 100, 2))

sing_conserved_props_prom_lore_aore <- prom_filter_counts_sing_lore_aore %>%
  filter(n == 1) %>%
  select(c("stage","redip")) %>%
  count(stage,redip) %>%
  mutate(nos = n) %>%
  left_join(prom_redip_activity, by = c("stage","redip")) %>%
  mutate(perc_conserved = round((nos / n.y) * 100, 2))

twoone_conserved_props_prom_lore_aore <- prom_filter_counts_twoone_lore_aore %>%
  filter(n == 1) %>%
  select(c("stage","redip")) %>%
  count(stage,redip) %>%
  mutate(nos = n) %>%
  left_join(prom_redip_activity, by = c("stage","redip")) %>%
  mutate(perc_conserved = round((nos / n.y) * 100, 2))

#### get a detailed report here

dups_conserved_props_active_pro_lore_aore_list<-prom_filter_counts_lore_aore %>%
  filter(n == 2) %>% mutate(chromosome = str_replace(chromosome, "^ssa0*(\\d+)$", "ssa\\1")) %>%
  unite(c(chromosome,start,end),col="Origin", sep="_") %>% select(-c(stage)) %>%
  distinct(Origin,redip) %>%  mutate(Origin = paste0(Origin, ".maf")) %>% left_join(dup_prom_detailed,by="Origin")

write.table(dups_conserved_props_active_pro_lore_aore_list,"shared_active_promoters_aore_lore_detailed.txt",
            quote = F, na="NA", row.names = F, col.names = T)

dups_conserved_props_active_pro_lore_aore_list %>% distinct(Origin,redip) %>% count(redip)

sing_conserved_props_active_pro_lore_aore_list <- prom_filter_counts_sing_lore_aore %>%
  filter(n == 1) %>% mutate(chromosome = str_replace(chromosome, "^ssa0*(\\d+)$", "ssa\\1")) %>%
  unite(c(chromosome,start,end),col="Origin", sep="_") %>% select(-c(stage)) %>%
  distinct(Origin,redip) %>%  mutate(Origin = paste0(Origin, ".maf")) %>% left_join(sing_prom_detailed,by="Origin")

write.table(sing_conserved_props_active_pro_lore_aore_list,"exclusive_active_promoters_aore_lore_detailed.txt",
            quote = F, na="NA", row.names = F, col.names = T)

twoone_conserved_props_active_pro_lore_aore_list <- prom_filter_counts_twoone_lore_aore %>%
  filter(n == 1) %>% mutate(chromosome = str_replace(chromosome, "^ssa0*(\\d+)$", "ssa\\1")) %>%
  unite(c(chromosome,start,end),col="Origin", sep="_") %>% select(-c(stage)) %>%
  distinct(Origin,redip) %>%  mutate(Origin = paste0(Origin, ".maf")) %>% left_join(twoone_prom_detailed,by="Origin")

write.table(twoone_conserved_props_active_pro_lore_aore_list,"aligned_active_promoters_aore_lore_detailed.txt",
            quote = F, na="NA", row.names = F, col.names = T)


## plot this

summ2_p_lore_aore <- dups_conserved_props_prom_lore_aore %>% mutate(cons = "twoway")
summ1_p_lore_aore <- sing_conserved_props_prom_lore_aore %>% mutate(cons = "a_unique")
summ21_p_lore_aore <- twoone_conserved_props_prom_lore_aore %>% mutate(cons = "two_one")
summ12_p_lore_aore <- bind_rows(summ2_p_lore_aore, summ1_p_lore_aore, summ21_p_lore_aore)

summ12_p_lore_aore_actual_counts <- summ12_p_lore_aore


##run fishers test across all combinations, remember to filter for AORE and LORe regions separately.

summ12_p_lore_aore_ft<-summ12_p_lore_aore_actual_counts %>% separate(stage,into=c("tissue","maturity","sex"),sep="_",remove=F) %>% replace(is.na(.), "embryo") %>%
  group_by(redip,cons,tissue,maturity) %>% summarise(cons_peaks=round(mean(nos),0), total_peaks = round(mean(n.y),0)) %>% ungroup() %>%
  unite(col="tissue", tissue, maturity, sep ="_")

write.xlsx(summ12_p_lore_aore_ft,"promoter_conservation_actual_numbers_lore_aore_ontogeny_2025.xlsx") #to be used as supplementary data

cons_vec_active_pro_lore_aore<- unique(summ12_p_lore_aore_ft$cons)

result_fishers_active_pro_aore<-lapply(cons_vec_active_pro_lore_aore, function(a1) { ## for AORe regions
  
  t_group<- summ12_p_lore_aore_ft %>% filter(redip=="AORe") %>% filter(cons==a1)
  
  idx = t(combn(seq_along(t_group$tissue),2))
  
  lapply(1:nrow(idx),function(i){
    
    tissue1_idx <- idx[i, 1]
    tissue2_idx <- idx[i, 2]
    
    test = fisher.test(t_group[idx[i,],c("cons_peaks","total_peaks")])
    
    data.frame(
      tissue1 = t_group$tissue[tissue1_idx],
      conserved_peaks_1 = t_group$cons_peaks[tissue1_idx],
      total_number_peaks_1 = t_group$total_peaks[tissue1_idx],
      tissue2 = t_group$tissue[tissue2_idx],
      conserved_peaks_2 = t_group$cons_peaks[tissue2_idx],
      total_number_peaks_2 = t_group$total_peaks[tissue2_idx],
      odds_ratio = as.numeric(test$estimate),
      p = as.numeric(test$p.value)
      
      
    )
  })
  
})

names(result_fishers_active_pro_aore)<-cons_vec_active_pro_lore_aore

##merge all datasets 

# Create an empty list to store merged data frames
merged_dfs_active_pro_AORe <- list()

# Iterate over each element in my_list
for (i in seq_along(result_fishers_active_pro_aore)) {
  # Extract the sublists for each letter
  sublists <- result_fishers_active_pro_aore[[i]]
  
  # Merge the data frames within each sublist
  merged_df <- do.call(bind_rows, sublists)
  
  # Add an additional column with the name of i
  merged_df <- mutate(merged_df, tissue = names(result_fishers_active_pro_aore)[i])
  
  # Store the merged data frame in merged_dfs list
  merged_dfs_active_pro_AORe[[i]] <- merged_df
  
}

# Combine all merged data frames into a single data frame

final_merged_results_fishers_active_promoters_AORe_conservation <- do.call(bind_rows, merged_dfs_active_pro_AORe)

final_merged_results_fishers_active_promoters_AORe_conservation <- final_merged_results_fishers_active_promoters_AORe_conservation %>%
  mutate(
    element = "Active Promoter",
    redip = "Early rediploidization",
    p_bonferroni = p * n(),  # Bonferroni correction
    p_bonferroni = ifelse(p_bonferroni > 1, 1, p_bonferroni),  # cap at 1
    significant_pre_correction = ifelse(p < 0.05, "yes", "no"),
    significant_post_correction = ifelse(p_bonferroni < 0.05, "yes", "no")
  )

write.xlsx(final_merged_results_fishers_active_promoters_AORe_conservation,"active_promoters_AORE_conservation_across_stages_fishers_test_2025.xlsx")


##now same for LORE regions

result_fishers_active_pro_lore<-lapply(cons_vec_active_pro_lore_aore, function(a1) { ## for LORe regions
  
  t_group<- summ12_p_lore_aore_ft %>% filter(redip=="LORe") %>% filter(cons==a1)
  
  idx = t(combn(seq_along(t_group$tissue),2))
  
  lapply(1:nrow(idx),function(i){
    
    tissue1_idx <- idx[i, 1]
    tissue2_idx <- idx[i, 2]
    
    test = fisher.test(t_group[idx[i,],c("cons_peaks","total_peaks")])
    
    data.frame(
      tissue1 = t_group$tissue[tissue1_idx],
      conserved_peaks_1 = t_group$cons_peaks[tissue1_idx],
      total_number_peaks_1 = t_group$total_peaks[tissue1_idx],
      tissue2 = t_group$tissue[tissue2_idx],
      conserved_peaks_2 = t_group$cons_peaks[tissue2_idx],
      total_number_peaks_2 = t_group$total_peaks[tissue2_idx],
      odds_ratio = as.numeric(test$estimate),
      p = as.numeric(test$p.value)
      
      
    )
  })
  
})

names(result_fishers_active_pro_lore)<-cons_vec_active_pro_lore_aore

##merge all datasets 

# Create an empty list to store merged data frames
merged_dfs_active_pro_LORe <- list()

# Iterate over each element in my_list
for (i in seq_along(result_fishers_active_pro_lore)) {
  # Extract the sublists for each letter
  sublists <- result_fishers_active_pro_lore[[i]]
  
  # Merge the data frames within each sublist
  merged_df <- do.call(bind_rows, sublists)
  
  # Add an additional column with the name of i
  merged_df <- mutate(merged_df, tissue = names(result_fishers_active_pro_lore)[i])
  
  # Store the merged data frame in merged_dfs list
  merged_dfs_active_pro_LORe[[i]] <- merged_df
  
}

# Combine all merged data frames into a single data frame

final_merged_results_fishers_active_promoters_LORe_conservation <- do.call(bind_rows, merged_dfs_active_pro_LORe)

final_merged_results_fishers_active_promoters_LORe_conservation <- final_merged_results_fishers_active_promoters_LORe_conservation %>%
  mutate(
    element = "Active Promoter",
    redip = "Late rediploidization",
    p_bonferroni = p * n(),  # Bonferroni correction
    p_bonferroni = ifelse(p_bonferroni > 1, 1, p_bonferroni),  # cap at 1
    significant_pre_correction = ifelse(p < 0.05, "yes", "no"),
    significant_post_correction = ifelse(p_bonferroni < 0.05, "yes", "no")
  )

write.xlsx(final_merged_results_fishers_active_promoters_LORe_conservation,"active_promoters_LORE_conservation_across_stages_fishers_test_2025.xlsx")

############################################################
# 17. PROMOTER CONSERVATION ACROSS DEVELOPMENT (FIGURE)
############################################################
# This figure summarizes promoter conservation proportions
# across developmental stages and rediploidization regimes

############################
# 17.1 CLEAN + STRUCTURE INPUT DATA
############################

prom_conservation_df <- summ12_p_lore_aore_actual_counts %>%
  
  # split stage into interpretable components
  separate(stage,
           into = c("tissue", "maturity", "sex"),
           sep = "_",
           remove = FALSE,
           fill = "right") %>%
  
  mutate(
    sex = replace_na(sex, "embryo"),
    maturity = replace_na(maturity, "embryo")
  ) %>%
  
  # average conservation per group
  group_by(cons, tissue, maturity, redip) %>%
  summarise(avg_perc = mean(perc_conserved), .groups = "drop") %>%
  
  # annotate developmental vs adult mapping
  mutate(map = if_else(maturity == "embryo", "devmap", "bodymap")) %>%
  
  unite("tissue_maturity", tissue, maturity, sep = "_", remove = FALSE)

############################
# 17.2 NORMALISE WITHIN STAGE
############################

plot_df <- prom_conservation_df %>%
  mutate(
    tissue_maturity = factor(tissue_maturity, levels = stage_order),
    redip = factor(redip, levels = c("AORe", "LORe"))
  ) %>%
  
  group_by(redip, tissue_maturity) %>%
  mutate(prop = avg_perc / sum(avg_perc)) %>%
  ungroup()

############################
# 17.3 CREATE COLOR MAPPING
############################

blue_shades <- brewer.pal(3, "Blues")[3:1]
red_shades  <- brewer.pal(3, "Reds")[3:1]

fill_colors <- c(
  "a_unique"   = red_shades[3],
  "two_one"    = red_shades[2],
  "twoway"     = red_shades[1],
  "AORe.a_unique" = blue_shades[3],
  "AORe.two_one"  = blue_shades[2],
  "AORe.twoway"   = blue_shades[1]
)

facet_labels <- c(
  AORe = "Early rediploidization",
  LORe = "Late rediploidization"
)

############################
# 17.4 CREATE COMPOSITE CATEGORY ID
############################

plot_df <- plot_df %>%
  mutate(
    cons_redip = if_else(
      redip == "AORe",
      paste0("AORe.", cons),
      cons
    )
  )

############################
# 17.5 ADD VISUAL GAP (DEV vs BODY SEPARATION)
############################

gap_df <- expand.grid(
  tissue_maturity = "gap",
  redip = unique(plot_df$redip),
  cons_redip = unique(plot_df$cons_redip)
) %>%
  mutate(prop = 0)

plot_df_aug <- bind_rows(plot_df, gap_df) %>%
  mutate(
    tissue_maturity = factor(tissue_maturity, levels = stage_order)
  )

############################
# 17.6 FINAL PLOT
############################

p_promoter_conservation <- ggplot(
  plot_df_aug,
  aes(x = tissue_maturity, y = prop, fill = cons_redip)
) +
  
  geom_col(
    color = "black",
    width = 0.85,
    linewidth = 0.07
  ) +
  
  scale_fill_manual(values = fill_colors) +
  
  scale_x_discrete(labels = c(
    "LB", "MG", "ES", "MS", "LS", "",
    "Br\n(I)", "Br\n(M)",
    "Mu\n(I)", "Mu\n(M)",
    "Li\n(I)", "Li\n(M)",
    "Ov\n(I)", "Ov\n(M)",
    "Te\n(I)", "Te\n(M)"
  )) +
  
  facet_wrap(
    ~ redip,
    nrow = 1,
    strip.position = "top",
    labeller = as_labeller(facet_labels)
  ) +
  
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  
  labs(y = "Proportion of promoters") +
  
  theme_classic(base_size = 9) +
  theme(
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(face = "bold"),
    panel.border = element_rect(color = "black", fill = NA),
    axis.text.x = element_text(size = 6, vjust = 0.5),
    axis.text.y = element_text(size = 6),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 7),
    legend.position = "none",
    panel.spacing = unit(0.2, "lines")
  )

############################
# 17.7 SAVE OUTPUT ### Figure 4 C
############################

ggsave(
  filename = "F4_promoter_conservation_dev_body_redip_2025.tiff",
  plot = p_promoter_conservation,
  dpi = 600,
  width = 11.9,
  height = 6,
  units = "cm"
)

ggsave(
  filename = "F4_promoter_conservation_dev_body_redip_2025.pdf",
  plot = p_promoter_conservation,
  device = cairo_pdf,
  width = 11.9,
  height = 6,
  dpi = 600,
  units = "cm"
)


############################################################
# ACTIVE ENHANCERS — CONSERVATION PIPELINE
############################################################
# This section quantifies conservation of active enhancers
# across ontogeny and genomic duplication classes

############################################################
# 1. LOAD + STANDARDISE ACTIVE ENHANCERS
############################################################

active_enh <- read.table(
  "../../atac_merged_gareth/redo_salmon_new_annotation_2024/Active_enhancer_1.bed",
  col.names = c("chromosome", "start", "end", "details", "peaksize", "strand")
) %>%
  unite("Origin", chromosome, start, end, sep = "_", remove = FALSE) %>%
  select(-peaksize, -strand) %>%
  mutate(
    chromosome = sprintf("Ssal%02d", as.integer(chromosome)),
    Origin = paste0("ssa", Origin, ".maf"),
    details = as.character(details)
  )

active_enh_2 <- active_enh %>%
  mutate(details = strsplit(details, ",")) %>%
  unnest(details) %>%
  separate(details, into = c("tissue", "chrom"), sep = ":") %>%
  filter(chrom == "Active_enhancer")

rep_str <- c(
  "Gonad_Immature_Male" = "testis_immature_male",
  "Gonad_Immature_Female" = "ovary_immature_female",
  "Gonad_Mature_Female" = "ovary_mature_female",
  "Gonad_Mature_Male" = "testis_mature_male"
)

active_enh_2 <- active_enh_2 %>%
  mutate(tissue = str_replace_all(tissue, rep_str)) %>%
  mutate(tissue = tolower(tissue))

############################################################
# 2. GENERAL FUNCTION
############################################################

process_active_enhancers <- function(atac_df, summary_df, label) {
  
  peak_merge <- atac_df %>%
    select(1:12) %>%
    left_join(
      active_enh_2,
      by = c("chromosome.y" = "chromosome",
             "start.y" = "start",
             "end.y" = "end")
    ) %>%
    filter(!is.na(Origin.y)) %>%
    mutate(tissue = tolower(tissue))
  
  filtered <- summary_df %>%
    left_join(
      peak_merge,
      by = c("Origin" = "Origin.x", "stage" = "tissue")
    ) %>%
    filter(!is.na(chromosome.x)) %>%
    select(-peaks_open)
  
  counts <- filtered %>%
    count(Origin, stage, name = "n")
  
  list(
    filtered = filtered,
    counts = counts,
    label = label
  )
}

############################################################
# 3. APPLY TO ALL THREE CATEGORIES
############################################################

dup_enh  <- process_active_enhancers(two_two_atac_peaks, two_two_peaks_summary, "shared")
sing_enh <- process_active_enhancers(one_one_atac_peaks, one_one_peaks_summary, "exclusive")
twoone_enh <- process_active_enhancers(two_one_atac_peaks, two_one_peaks_summary, "alignable")

############################################################
# 4. CLASSIFY ENHANCERS
############################################################

dup_active_enh  <- filter(dup_enh$counts, n == 2)
sing_active_enh <- filter(sing_enh$counts, n == 1)
twoone_active_enh <- filter(twoone_enh$counts, n == 1)

############################################################
# 5. WRITE BASIC OUTPUT TABLES
############################################################

write.table(dup_active_enh, "duplicated_active_enhancers_across_ontogeny.bed", row.names = FALSE)
write.table(sing_active_enh, "singleton_active_enhancers_across_ontogeny.bed", row.names = FALSE)
write.table(twoone_active_enh, "two_one_active_enhancers_across_ontogeny.bed", row.names = FALSE)

############################################################
# 6. DETAILED OUTPUT FUNCTION
############################################################

make_detailed_table <- function(df, atac_original) {
  
  df %>%
    distinct(Origin) %>%
    left_join(atac_original, by = "Origin") %>%
    rename(
      Peak_ID = Origin,
      Alignment_chromosome = chromosome.x,
      Alignment_start = start.x,
      Alignment_end = end.x,
      ATAC_ID = atac_loc,
      ATAC_chromosome = chromosome.y,
      ATAC_start = start.y,
      ATAC_end = end.y
    ) %>%
    select(-any_of("species_b"))
}

dup_active_enh_detailed <- make_detailed_table(dup_active_enh, two_two_atac_peaks)
sing_active_enh_detailed <- make_detailed_table(sing_active_enh, one_one_atac_peaks)
twoone_active_enh_detailed <- make_detailed_table(twoone_active_enh, two_one_atac_peaks)

############################################################
# 7. WRITE DETAILED OUTPUT
############################################################

write.table(dup_active_enh_detailed,
            "duplicated_active_enhancers_across_ontogeny_detailed_2025.txt",
            row.names = FALSE)

write.table(sing_active_enh_detailed,
            "singleton_active_enhancers_across_ontogeny_detailed_2025.txt",
            row.names = FALSE)

write.table(twoone_active_enh_detailed,
            "aligned_active_enhancers_across_ontogeny_detailed_2025.txt",
            row.names = FALSE)

############################################################
# 8. SUMMARY ACTIVE ENHANCER ACTIVITY
############################################################

active_enh_activity <- bind_rows(
  dup_active_enh,
  sing_active_enh,
  twoone_active_enh
) %>%
  count(stage, name = "n_enhancers")


############################################
## LOAD LORe REGIONS
############################################

salmon_LORe <- read.table("./Ssal_late_rediploidized_regions.tsv", header = TRUE) %>%
  rename(chromosome = chrom) %>%
  mutate(
    chromosome = sprintf("ssa%02d", as.integer(chromosome)),
    start = as.numeric(start),
    end = as.numeric(end)
  )

############################################
##  LORe ANNOTATION PIPELINE
############################################

add_lore_annotation <- function(df, lore_df) {
  
  df %>%
    mutate(Origin1 = str_remove(Origin, "\\.maf$")) %>%
    separate(Origin1,
             into = c("chromosome", "start", "end"),
             sep = "_",
             convert = TRUE) %>%
    mutate(chromosome = sprintf("ssa%02d",
                                as.integer(str_remove(chromosome, "ssa")))) %>%
    
    fuzzy_left_join(
      lore_df,
      by = c("chromosome" = "chromosome",
             "start" = "start",
             "end" = "end"),
      match_fun = list(`==`, `>=`, `<=`)
    ) %>%
    
    mutate(
      redip = if_else(!is.na(start.y), "LORe", "AORe")
    ) %>%
    
    select(-ends_with(".y")) %>%
    rename_with(~ str_remove(., "\\.x$")) %>%
    select(-any_of(c("label", "assembly_name")))
}


############################################
## APPLY LORe/AORe CLASSIFICATION
## (DUPLICATE / SINGLE / TWO-ONE)
############################################

active_enh_filter_counts_lore_aore <- active_enh_filter %>%
  group_by(Origin, stage) %>%
  summarise(n = n(), .groups = "drop") %>%
  add_lore_annotation(salmon_LORe)

active_enh_filter_counts_sing_lore_aore <- active_enh_filter_sing %>%
  group_by(Origin, stage) %>%
  summarise(n = n(), .groups = "drop") %>%
  add_lore_annotation(salmon_LORe)

active_enh_filter_counts_twoone_lore_aore <- active_enh_filter_twoone %>%
  group_by(Origin, stage) %>%
  summarise(n = n(), .groups = "drop") %>%
  add_lore_annotation(salmon_LORe)

###########################
## test for breadth of activity in active enhancers
###########################

shared <- active_enh_filter_counts_lore_aore %>%
  filter(n == 2) %>%
  mutate(category = "shared") %>%
  slice(rep(1:n(), each = 2)) %>%
  group_by(chromosome, start, end, stage) %>%
  mutate(copy_id = row_number()) %>%
  ungroup()

exclusive <- active_enh_filter_counts_sing_lore_aore %>%
  filter(n == 1) %>%
  mutate(category = "exclusive", copy_id = 1)

alignable <- active_enh_filter_counts_twoone_lore_aore %>%
  filter(n == 1) %>%
  mutate(category = "alignable", copy_id = 1)

breadth_df <- bind_rows(shared, exclusive, alignable) %>%
  separate(stage, into = c("tissue", "maturity", "sex"),
           sep = "_", remove = FALSE) %>%
  mutate(
    tissue = replace_na(tissue, "embryo"),
    tissue = paste(tissue, maturity, sep = "_")
  )

#------------------------------------------------------------
# 7. Compute breadth
#------------------------------------------------------------

breadth <- breadth_df %>%
  group_by(redip, chromosome, start, end, copy_id, category) %>%
  summarise(n_tissues = n_distinct(tissue), .groups = "drop")

breadth_summary <- breadth %>%
  group_by(redip, category) %>%
  summarise(
    mean_tissues   = mean(n_tissues),
    median_tissues = median(n_tissues),
    sd_tissues     = sd(n_tissues),
    n_peaks        = n(),
    .groups = "drop"
  )

breadth_wilcox <- breadth %>%
  group_by(redip) %>%
  wilcox_test(n_tissues ~ category) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()

#------------------------------------------------------------
# 8. Plotting setup
#------------------------------------------------------------

cat_cols <- c(
  "shared"    = "#1b9e77",
  "alignable" = "#d95f02",
  "exclusive" = "#7570b3"
)

#------------------------------------------------------------
# 9. Histogram
#------------------------------------------------------------

ggplot(breadth, aes(n_tissues, fill = category)) +
  geom_histogram(binwidth = 1, color = "black", alpha = 0.85,
                 position = "dodge") +
  scale_fill_manual(values = cat_cols) +
  facet_wrap(~ redip) +
  theme_bw(base_size = 12) +
  labs(x = "Number of tissues active",
       y = "Number of active enhancers",
       fill = "Category")

ggsave("enhancer_breadth_histogram.pdf", width = 25, height = 12)


#------------------------------------------------------------
# 10. Boxplot
#------------------------------------------------------------

ggplot(breadth, aes(category, n_tissues, fill = category)) +
  geom_boxplot(outlier.shape = 21) +
  scale_fill_manual(values = cat_cols) +
  facet_wrap(~ redip) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none") +
  labs(x = "Category",
       y = "Breadth of enhancer activity")

ggsave("enhancer_breadth_boxplot.pdf", width = 25, height = 12)

#------------------------------------------------------------
# 11. Tissue-specific enhancer counts
#------------------------------------------------------------

stage_order <- c(
  "lateblastulation_embryo","midgastrulation_embryo",
  "earlysomitogenesis_embryo","midsomitogenesis_embryo",
  "latesomitogenesis_embryo","gap",
  "brain_immature","brain_mature",
  "muscle_immature","muscle_mature",
  "liver_immature","liver_mature",
  "ovary_immature","ovary_mature",
  "testis_immature","testis_mature"
)

combo_summary <- breadth_df %>%
  group_by(redip, chromosome, start, end, copy_id, category) %>%
  summarise(
    tissue_combo = str_c(sort(unique(tissue)), collapse = ","),
    .groups = "drop"
  ) %>%
  count(redip, category, tissue_combo, name = "n_peaks")

plot_data <- combo_summary %>%
  filter(!str_detect(tissue_combo, ",")) %>%
  mutate(
    redip = recode(redip,
                   AORe = "Early rediploidization",
                   LORe = "Late rediploidization"),
    tissue_combo = factor(tissue_combo, levels = stage_order),
    category = factor(category, levels = c("shared","alignable","exclusive"))
  )

ggplot(plot_data, aes(tissue_combo, n_peaks, fill = category)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = cat_cols) +
  facet_wrap(~ redip, scales = "free_x") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 70, hjust = 1)) +
  labs(x = "Sample type",
       y = "Number of tissue-specific enhancers")

ggsave("enhancer_sample_specific_counts.pdf", width = 25, height = 12)


#########################
## get activity of enhancers across tissues 
############################

dups_act_enh_activity<-active_enh_filter_counts_lore_aore %>%
  filter(n == 2) %>%
  select(c("stage","redip")) %>%
  count(stage,redip) %>%
  mutate(nos = n*2 ) %>% ## account for duplication here
  select(-c(n))

sing_act_enh_activity<- active_enh_filter_counts_sing_lore_aore %>%
  filter(n == 1) %>%
  select(c("stage","redip")) %>%
  count(stage,redip) %>%
  mutate(nos = n)%>%
  select(-c(n))

twoone_act_enh_activity<-active_enh_filter_counts_twoone_lore_aore %>%
  filter(n == 1) %>%
  select(c("stage","redip")) %>%
  count(stage,redip) %>%
  mutate(nos = n)%>%
  select(-c(n))

act_enh_redip_activity<- bind_rows(dups_act_enh_activity,sing_act_enh_activity,twoone_act_enh_activity) %>% group_by(stage,redip) %>%
  summarise(n=sum(nos))


## get proportions for three categories across LORe AORe

dups_conserved_props_active_enh_lore_aore<-active_enh_filter_counts_lore_aore %>%
  filter(n == 2) %>%
  select(c("stage","redip")) %>%
  count(stage,redip) %>%
  mutate(nos = n*2) %>%
  left_join(act_enh_redip_activity, by = c("stage","redip")) %>%
  mutate(perc_conserved = round((nos / n.y) * 100, 2)) %>%
  mutate(Activity="Shared")

sing_conserved_props_active_enh_lore_aore <- active_enh_filter_counts_sing_lore_aore %>%
  filter(n == 1) %>%
  select(c("stage","redip")) %>%
  count(stage,redip) %>%
  mutate(nos = n) %>%
  left_join(act_enh_redip_activity, by = c("stage","redip")) %>%
  mutate(perc_conserved = round((nos / n.y) * 100, 2))%>%
  mutate(Activity="Exclusive")

twoone_conserved_props_active_enh_lore_aore <- active_enh_filter_counts_twoone_lore_aore %>%
  filter(n == 1) %>%
  select(c("stage","redip")) %>%
  count(stage,redip) %>%
  mutate(nos = n) %>%
  left_join(act_enh_redip_activity, by = c("stage","redip")) %>%
  mutate(perc_conserved = round((nos / n.y) * 100, 2)) %>%
  mutate(Activity="Alignable")

#### save these as bed files - shared, exclusive, aligned_active_enhancers_aore_lore_list.bed

list_active_dup_enh_filter_counts_lore_aore <- active_enh_filter_counts_lore_aore %>% filter(n == 2)
write.table(list_active_dup_enh_filter_counts_lore_aore,"shared_active_enhancers_aore_lore_list.txt",
            quote = F, na="NA", row.names = F, col.names = T)


list_active_enh_filter_counts_sing_lore_aore<- active_enh_filter_counts_sing_lore_aore %>% filter(n == 1)
write.table(list_active_enh_filter_counts_sing_lore_aore,"exclusive_active_enhancers_aore_lore_list.txt",
            quote = F, na="NA", row.names = F, col.names = T)

list_active_enh_filter_counts_twoone_lore_aore <- active_enh_filter_counts_twoone_lore_aore %>% filter(n == 1)
write.table(list_active_enh_filter_counts_twoone_lore_aore,"alignable_active_enhancers_aore_lore_list.txt",
            quote = F, na="NA", row.names = F, col.names = T)

### generate  a detailed list

dups_conserved_props_active_enh_lore_aore_list<-active_enh_filter_counts_lore_aore %>%
  filter(n == 2) %>% mutate(chromosome = str_replace(chromosome, "^ssa0*(\\d+)$", "ssa\\1")) %>%
  unite(c(chromosome,start,end),col="Peak_ID", sep="_") %>% select(-c(stage)) %>%
  distinct(Peak_ID,redip) %>%  mutate(Peak_ID = paste0(Peak_ID, ".maf")) %>% left_join(dup_active_enh_detailed,by="Peak_ID")

write.table(dups_conserved_props_active_enh_lore_aore_list,"shared_active_enhancers_aore_lore_detailed.txt",
            quote = F, na="NA", row.names = F, col.names = T)


sing_conserved_props_active_enh_lore_aore_list <- active_enh_filter_counts_sing_lore_aore %>%
  filter(n == 1) %>% mutate(chromosome = str_replace(chromosome, "^ssa0*(\\d+)$", "ssa\\1")) %>%
  unite(c(chromosome,start,end),col="Peak_ID", sep="_") %>% select(-c(stage)) %>%
  distinct(Peak_ID,redip) %>%  mutate(Peak_ID = paste0(Peak_ID, ".maf")) %>% left_join(sing_active_enh_detailed,by="Peak_ID") 

write.table(sing_conserved_props_active_enh_lore_aore_list,"exclusive_active_enhancers_aore_lore_detailed.txt",
            quote = F, na="NA", row.names = F, col.names = T)

twoone_conserved_props_active_enh_lore_aore_list <- active_enh_filter_counts_twoone_lore_aore %>%
  filter(n == 1) %>% mutate(chromosome = str_replace(chromosome, "^ssa0*(\\d+)$", "ssa\\1")) %>%
  unite(c(chromosome,start,end),col="Peak_ID", sep="_") %>% select(-c(stage)) %>%
  distinct(Peak_ID,redip) %>%  mutate(Peak_ID = paste0(Peak_ID, ".maf")) %>% left_join(twoone_active_enh_detailed,by="Peak_ID") 

write.table(twoone_conserved_props_active_enh_lore_aore_list,"aligned_active_enhancers_aore_lore_detailed.txt",
            quote = F, na="NA", row.names = F, col.names = T)

## plot this

summ2_active_enh_lore_aore <- dups_conserved_props_active_enh_lore_aore %>% mutate(cons = "Shared")
summ1_active_enh_lore_aore <- sing_conserved_props_active_enh_lore_aore %>% mutate(cons = "Exclusive")
summ21_active_enh_lore_aore <- twoone_conserved_props_active_enh_lore_aore %>% mutate(cons = "Alignable")
summ12_active_enh_lore_aore <- bind_rows(summ2_active_enh_lore_aore, summ1_active_enh_lore_aore, summ21_active_enh_lore_aore)

summ12_active_enh_lore_aore_actual_counts <-summ12_active_enh_lore_aore

summ12_active_enh_lore_aore_ft<-summ12_active_enh_lore_aore_actual_counts %>% separate(stage,into=c("tissue","maturity","sex"),sep="_",remove=F) %>% replace(is.na(.), "embryo") %>%
  group_by(redip,cons,tissue,maturity) %>% summarise(cons_peaks=round(mean(nos),0), total_peaks = round(mean(n.y),0)) %>% ungroup() %>%
  unite(col="tissue", tissue, maturity, sep ="_")


write.xlsx(summ12_active_enh_lore_aore_ft,"active_enhancer_conservation_actual_numbers_lore_aore_ontogeny_2025.xlsx") #to be used as supplementary data

cons_vec_active_enh_lore_aore<- unique(summ12_active_enh_lore_aore_ft$cons)

result_fishers_active_enh_aore<-lapply(cons_vec_active_enh_lore_aore, function(a1) { ## for AORe regions
  
  t_group<- summ12_active_enh_lore_aore_ft %>% filter(redip=="AORe") %>% filter(cons==a1)
  
  idx = t(combn(seq_along(t_group$tissue),2))
  
  lapply(1:nrow(idx),function(i){
    
    tissue1_idx <- idx[i, 1]
    tissue2_idx <- idx[i, 2]
    
    test = fisher.test(t_group[idx[i,],c("cons_peaks","total_peaks")])
    
    data.frame(
      tissue1 = t_group$tissue[tissue1_idx],
      conserved_peaks_1 = t_group$cons_peaks[tissue1_idx],
      total_number_peaks_1 = t_group$total_peaks[tissue1_idx],
      tissue2 = t_group$tissue[tissue2_idx],
      conserved_peaks_2 = t_group$cons_peaks[tissue2_idx],
      total_number_peaks_2 = t_group$total_peaks[tissue2_idx],
      odds_ratio = as.numeric(test$estimate),
      p = as.numeric(test$p.value)
      
      
    )
  })
  
})

names(result_fishers_active_enh_aore)<-cons_vec_active_enh_lore_aore

##merge all datasets 

# Create an empty list to store merged data frames
merged_dfs_active_enh_AORe <- list()

# Iterate over each element in my_list
for (i in seq_along(result_fishers_active_enh_aore)) {
  # Extract the sublists for each letter
  sublists <- result_fishers_active_enh_aore[[i]]
  
  # Merge the data frames within each sublist
  merged_df <- do.call(bind_rows, sublists)
  
  # Add an additional column with the name of i
  merged_df <- mutate(merged_df, tissue = names(result_fishers_active_enh_aore)[i])
  
  # Store the merged data frame in merged_dfs list
  merged_dfs_active_enh_AORe[[i]] <- merged_df
  
}

# Combine all merged data frames into a single data frame

final_merged_results_fishers_active_enhancers_AORe_conservation <- do.call(bind_rows, merged_dfs_active_enh_AORe)

final_merged_results_fishers_active_enhancers_AORe_conservation <- final_merged_results_fishers_active_enhancers_AORe_conservation %>%
  mutate(
    element = "Active Enhancer",
    redip = "Early rediploidization",
    p_bonferroni = p * n(),  # Bonferroni correction
    p_bonferroni = ifelse(p_bonferroni > 1, 1, p_bonferroni),  # cap at 1
    significant_pre_correction = ifelse(p < 0.05, "yes", "no"),
    significant_post_correction = ifelse(p_bonferroni < 0.05, "yes", "no")
  )

write.xlsx(final_merged_results_fishers_active_enhancers_AORe_conservation,"active_enhancers_AORE_conservation_across_stages_fishers_test_2025.xlsx")

### now repeat it for LORe regions

result_fishers_active_enh_lore<-lapply(cons_vec_active_enh_lore_aore, function(a1) { ## for LORe regions
  
  t_group<- summ12_active_enh_lore_aore_ft %>% filter(redip=="LORe") %>% filter(cons==a1)
  
  idx = t(combn(seq_along(t_group$tissue),2))
  
  lapply(1:nrow(idx),function(i){
    
    tissue1_idx <- idx[i, 1]
    tissue2_idx <- idx[i, 2]
    
    test = fisher.test(t_group[idx[i,],c("cons_peaks","total_peaks")])
    
    data.frame(
      tissue1 = t_group$tissue[tissue1_idx],
      conserved_peaks_1 = t_group$cons_peaks[tissue1_idx],
      total_number_peaks_1 = t_group$total_peaks[tissue1_idx],
      tissue2 = t_group$tissue[tissue2_idx],
      conserved_peaks_2 = t_group$cons_peaks[tissue2_idx],
      total_number_peaks_2 = t_group$total_peaks[tissue2_idx],
      odds_ratio = as.numeric(test$estimate),
      p = as.numeric(test$p.value)
      
      
    )
  })
  
})

names(result_fishers_active_enh_lore)<-cons_vec_active_enh_lore_aore

##merge all datasets 

# Create an empty list to store merged data frames
merged_dfs_active_enh_LORe <- list()

# Iterate over each element in my_list
for (i in seq_along(result_fishers_active_enh_lore)) {
  # Extract the sublists for each letter
  sublists <- result_fishers_active_enh_lore[[i]]
  
  # Merge the data frames within each sublist
  merged_df <- do.call(bind_rows, sublists)
  
  # Add an additional column with the name of i
  merged_df <- mutate(merged_df, tissue = names(result_fishers_active_enh_lore)[i])
  
  # Store the merged data frame in merged_dfs list
  merged_dfs_active_enh_LORe[[i]] <- merged_df
  
}

# Combine all merged data frames into a single data frame

final_merged_results_fishers_active_enhancers_LORe_conservation <- do.call(bind_rows, merged_dfs_active_enh_LORe)

final_merged_results_fishers_active_enhancers_LORe_conservation <- final_merged_results_fishers_active_enhancers_LORe_conservation %>%
  mutate(
    element = "Active Enhancer",
    redip = "Late rediploidization",
    p_bonferroni = p * n(),  # Bonferroni correction
    p_bonferroni = ifelse(p_bonferroni > 1, 1, p_bonferroni),  # cap at 1
    significant_pre_correction = ifelse(p < 0.05, "yes", "no"),
    significant_post_correction = ifelse(p_bonferroni < 0.05, "yes", "no")
  )

write.xlsx(final_merged_results_fishers_active_enhancers_LORe_conservation,"active_enhancers_LORE_conservation_across_stages_fishers_test_2025.xlsx")

############################
## active enhancer conservation plot
############################

library(tidyverse)
library(RColorBrewer)
library(scales)

#------------------------------------------------------------
# 1. Prepare input table
#------------------------------------------------------------

summ12_active_enh_lore_aore_actual_counts <- summ12_active_enh_lore_aore_actual_counts %>%
  separate(stage,
           into = c("tissue", "maturity", "sex"),
           sep = "_",
           remove = FALSE) %>%
  mutate(
    tissue   = replace_na(tissue, "embryo")
  ) %>%
  group_by(cons, tissue, maturity, redip) %>%
  summarise(avg_perc = mean(perc_conserved), .groups = "drop") %>%
  mutate(
    map = if_else(maturity == "embryo", "devmap", "bodymap"),
    tissue_maturity = paste(tissue, maturity, sep = "_")
  )

#------------------------------------------------------------
# 2. Compute proportions within each facet
#------------------------------------------------------------

plot_df <- summ12_active_enh_lore_aore_actual_counts %>%
  group_by(redip, tissue_maturity) %>%
  mutate(prop = avg_perc / sum(avg_perc)) %>%
  ungroup()

#------------------------------------------------------------
# 3. Color palettes
#------------------------------------------------------------

blue_shades <- brewer.pal(3, "Blues")[3:1]
red_shades  <- brewer.pal(3, "Reds")[3:1]

fill_colors <- c(
  "Exclusive" = red_shades[3],
  "Alignable" = red_shades[2],
  "Shared"    = red_shades[1],
  "AORe.Exclusive" = blue_shades[3],
  "AORe.Alignable" = blue_shades[2],
  "AORe.Shared"    = blue_shades[1]
)

facet_labels <- c(
  "AORe" = "Early rediploidization",
  "LORe" = "Late rediploidization"
)

#------------------------------------------------------------
# 4. Build combined color identity
#------------------------------------------------------------

plot_df <- plot_df %>%
  mutate(
    cons_redip = if_else(redip == "AORe",
                         paste0("AORe.", cons),
                         cons)
  )

#------------------------------------------------------------
# 5. Add gap category
#------------------------------------------------------------

stage_order <- c(
  "lateblastulation_embryo","midgastrulation_embryo",
  "earlysomitogenesis_embryo","midsomitogenesis_embryo",
  "latesomitogenesis_embryo","gap",
  "brain_immature","brain_mature",
  "muscle_immature","muscle_mature",
  "liver_immature","liver_mature",
  "ovary_immature","ovary_mature",
  "testis_immature","testis_mature"
)

gap_df <- expand.grid(
  tissue_maturity = "gap",
  redip = unique(plot_df$redip),
  cons_redip = unique(plot_df$cons_redip)
) %>%
  mutate(prop = 0)

plot_df_aug <- bind_rows(plot_df, gap_df) %>%
  mutate(
    tissue_maturity = factor(tissue_maturity, levels = stage_order)
  )

#------------------------------------------------------------
# 6. Final factor ordering (clean + safe)
#------------------------------------------------------------

plot_df_aug <- plot_df_aug %>%
  mutate(
    cons = factor(cons, levels = c("Shared", "Alignable", "Exclusive")),
    cons_redip = factor(cons_redip)
  )

#------------------------------------------------------------
# 7. Plot
#------------------------------------------------------------

p <- ggplot(plot_df_aug,
            aes(x = tissue_maturity, y = prop, fill = cons_redip)) +
  geom_bar(stat = "identity",
           color = "black",
           width = 0.85,
           linewidth = 0.07) +
  scale_fill_manual(values = fill_colors) +
  scale_x_discrete(labels = c(
    "LB", "MG", "ES", "MS", "LS",
    "",
    "Br\n(I)", "Br\n(M)",
    "Mu\n(I)", "Mu\n(M)",
    "Li\n(I)", "Li\n(M)",
    "Ov\n(I)", "Ov\n(M)",
    "Te\n(I)", "Te\n(M)"
  )) +
  facet_wrap(~ redip,
             nrow = 1,
             strip.position = "top",
             labeller = as_labeller(facet_labels)) +
  scale_y_continuous(labels = number_format(accuracy = 0.1)) +
  theme_classic(base_size = 9) +
  theme(
    strip.background = element_rect(fill = "white", color = "black"),
    panel.border = element_rect(color = "black", fill = NA),
    axis.text.x = element_text(size = 6, vjust = 0.5),
    axis.text.y = element_text(size = 6),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 8),
    legend.position = "none"
  ) +
  ylab("Proportion of active enhancers")

#------------------------------------------------------------
# 8. Save output : Figure 4D
#------------------------------------------------------------

ggsave(
  "F4_1and2way_active_enhancer_conservation_devmap_bodymap_redip_2025_V2.tiff",
  plot = p,
  dpi = 600,
  width = 11.9,
  height = 6,
  units = "cm"
)

ggsave(
  "F4_1and2way_active_enhancer_conservation_devmap_bodymap_redip_2025_V2.pdf",
  plot = p,
  device = cairo_pdf,
  width = 11.9,
  height = 6,
  dpi = 600,
  units = "cm"
)

##########################################
#### Manuscript summary statistics
##########################################

#------------------------------------------------------------
# 1. Prepare LORe regions
#------------------------------------------------------------

lore_clean <- salmon_LORe %>%
  mutate(
    chr = as.integer(str_remove(chromosome, "ssa")),
    start_lore = start,
    end_lore = end
  )

#------------------------------------------------------------
# 2. Helper function for AORe / LORe coverage
#------------------------------------------------------------

calculate_lore_overlap <- function(active_df) {
  
  total_len <- sum(active_df$peaksize, na.rm = TRUE)
  
  overlap_df <- active_df %>%
    inner_join(lore_clean, by = "chr") %>%
    filter(
      start_active < end_lore,
      end_active > start_lore
    ) %>%
    mutate(
      overlap_start = pmax(start_active, start_lore),
      overlap_end   = pmin(end_active, end_lore),
      overlap_len   = overlap_end - overlap_start
    )
  
  in_lore_len  <- sum(overlap_df$overlap_len, na.rm = TRUE)
  out_lore_len <- total_len - in_lore_len
  
  tibble(
    total_length = total_len,
    lore_length = in_lore_len,
    aore_length = out_lore_len,
    lore_prop = round(in_lore_len / total_len, 4),
    aore_prop = round(out_lore_len / total_len, 4)
  )
}

#------------------------------------------------------------
# 3. Active promoters
#------------------------------------------------------------

active_prom_clean <- active_prom_2 %>%
  distinct(Origin, chromosome, start, end) %>%
  mutate(
    chr = as.integer(str_remove(chromosome, "Ssal")),
    start_active = start,
    end_active = end,
    peaksize = end - start
  )

promoter_peak_stats <- active_prom_clean %>%
  summarise(
    n_promoters = n(),
    total_bp = sum(peaksize),
    mean_bp = mean(peaksize),
    median_bp = median(peaksize),
    min_bp = min(peaksize),
    max_bp = max(peaksize)
  )

promoter_lore_stats <- calculate_lore_overlap(active_prom_clean)

#------------------------------------------------------------
# 4. Active enhancers
#------------------------------------------------------------

active_enh_clean <- active_enh_2 %>%
  distinct(Origin, chromosome, start, end) %>%
  mutate(
    chr = as.integer(str_remove(chromosome, "Ssal")),
    start_active = start,
    end_active = end,
    peaksize = end - start
  )

enhancer_peak_stats <- active_enh_clean %>%
  summarise(
    n_enhancers = n(),
    total_bp = sum(peaksize),
    mean_bp = mean(peaksize),
    median_bp = median(peaksize),
    min_bp = min(peaksize),
    max_bp = max(peaksize)
  )

enhancer_lore_stats <- calculate_lore_overlap(active_enh_clean)

#------------------------------------------------------------
# 5. Peak size summaries by conservation class
#------------------------------------------------------------

promoter_conservation_sizes <- bind_rows(
  dups_conserved_props_active_pro_lore_aore_list %>%
    mutate(category = "Shared"),
  
  twoone_conserved_props_active_pro_lore_aore_list %>%
    mutate(category = "Alignable"),
  
  sing_conserved_props_active_pro_lore_aore_list %>%
    mutate(category = "Exclusive")
) %>%
  group_by(category, redip) %>%
  summarise(
    n_peaks = n(),
    total_peak_size = sum(peak_size, na.rm = TRUE),
    mean_peak_size = mean(peak_size, na.rm = TRUE),
    median_peak_size = median(peak_size, na.rm = TRUE),
    .groups = "drop"
  )

enhancer_conservation_sizes <- bind_rows(
  dups_conserved_props_active_enh_lore_aore_list %>%
    mutate(category = "Shared"),
  
  twoone_conserved_props_active_enh_lore_aore_list %>%
    mutate(category = "Alignable"),
  
  sing_conserved_props_active_enh_lore_aore_list %>%
    mutate(category = "Exclusive")
) %>%
  group_by(category, redip) %>%
  summarise(
    n_peaks = n(),
    total_peak_size = sum(peak_size, na.rm = TRUE),
    mean_peak_size = mean(peak_size, na.rm = TRUE),
    median_peak_size = median(peak_size, na.rm = TRUE),
    .groups = "drop"
  )

#------------------------------------------------------------
# 6. Quick output for manuscript writing
#------------------------------------------------------------

promoter_peak_stats
promoter_lore_stats

enhancer_peak_stats
enhancer_lore_stats

promoter_conservation_sizes
enhancer_conservation_sizes


###################################################
### JSD vs enhancer conservation corrrelation
##################################################

sample_jsd <- read.xlsx("Diego_jsd.xlsx")


active_enh_vs_jsd_correlation<-summ12_active_enh_lore_aore_actual_counts %>% ungroup() %>% select(tissue_maturity,redip,avg_perc,cons) %>%
  left_join(sample_jsd, by=(c("tissue_maturity","redip"))) %>%
  mutate( cons_redip = paste0(cons, "_", redip))


# Calculate correlation coefficients per 'cons' group
cor_results_jsd <- active_enh_vs_jsd_correlation %>%
  group_by(redip,cons) %>%
  summarise(
    cor_test = list(cor.test(avg_perc, mean_jsd, method = "pearson")),
    .groups = "drop"
  ) %>%
  mutate(
    tidied = map(cor_test, broom::tidy)
  ) %>%
  unnest(tidied) %>%
  select(redip, cons,estimate, p.value, conf.low, conf.high)

cor_results_jsd <- cor_results_jsd %>%
  mutate(p_adj = p.adjust(p.value, method = "BH"))

write.xlsx(cor_results_jsd, "pearson_correlation_active_enh_meanJSD.xlsx")


fill_colors_jsd<- c(
  "Exclusive_LORe" = red_shades[3],
  "Alignable_LORe" = red_shades[2],
  "Shared_LORe" = red_shades[1],
  "Exclusive_AORe" = blue_shades[3],
  "Alignable_AORe" = blue_shades[2],
  "Shared_AORe" = blue_shades[1]
)

# Custom labels for facet strip
facet_labels <- c(
  "AORe" = "Early rediploidization",
  "LORe" = "Late rediploidization"
)

# Plot with ggplot2: facet by 'redip', color by 'cons'
ggplot(active_enh_vs_jsd_correlation, aes(x = avg_perc, y = mean_jsd)) +
  geom_point(
    aes(fill = cons_redip), 
    shape = 21, size = 1.5, stroke = 0.1, color = "grey45", alpha = 0.99
  ) +
  geom_smooth(
    method = "lm", se = FALSE, 
    aes(group = cons, color = cons_redip),
    linetype = "solid", linewidth = 0.5
  ) +
  facet_wrap(~redip, ncol = 1, scales = "free_y", labeller = as_labeller(facet_labels)) +
  scale_fill_manual(values = fill_colors_jsd) +
  scale_color_manual(values = fill_colors_jsd) +
  theme_classic(base_size = 9) +
  theme(
    strip.background = element_rect(fill = "white", color = "black"),
    panel.border = element_rect(color = "black", fill = NA),
    axis.text.x = element_text(vjust = 0.5, size = 7),
    axis.text.y = element_text(size = 7),
    axis.title.x = element_text(size = 8),
    axis.title.y = element_text(size = 8),
    legend.position = "none"
  ) +
  labs(
    x = "Proportion of active enhancers",
    y = "JSD (2:2 ohnologs)",
    color = "Conservation (cons)"
  )

ggsave(paste0("F4_1and2way_active_enhancer_conservation_JSD_correlation_redip_2025_V1.tiff"), dpi = 600, width = 6.2, height = 9.9, units = "cm") 

ggsave(filename = "F4_1and2way_active_enhancer_conservation_JSD_correlation_redip_2025_V1.pdf",  plot = last_plot(), device = cairo_pdf, width = 6.2, height = 9.9, dpi = 600,units = "cm")

###################################################
### JSD vs promoter conservation corrrelation
##################################################

sample_jsd <- read.xlsx("Diego_jsd.xlsx")


active_pro_vs_jsd_correlation<-summ12_p_lore_aore_actual_counts %>% ungroup() %>% select(tissue_maturity,redip,avg_perc,cons) %>%
  left_join(sample_jsd, by=(c("tissue_maturity","redip"))) %>%
  mutate( cons_redip = paste0(cons, "_", redip))


# Calculate correlation coefficients per 'cons' group
cor_results_jsd_promoter <- active_pro_vs_jsd_correlation %>%
  group_by(redip,cons) %>%
  summarise(
    cor_test = list(cor.test(avg_perc, mean_jsd, method = "pearson")),
    .groups = "drop"
  ) %>%
  mutate(
    tidied = map(cor_test, broom::tidy)
  ) %>%
  unnest(tidied) %>%
  select(redip, cons,estimate, p.value, conf.low, conf.high)

cor_results_jsd_promoter <- cor_results_jsd_promoter %>%
  mutate(p_adj = p.adjust(p.value, method = "BH"))

write.xlsx(cor_results_jsd_promoter, "pearson_correlation_active_promoter_meanJSD.xlsx")


fill_colors_jsd<- c(
  "Exclusive_LORe" = red_shades[3],
  "Alignable_LORe" = red_shades[2],
  "Shared_LORe" = red_shades[1],
  "Exclusive_AORe" = blue_shades[3],
  "Alignable_AORe" = blue_shades[2],
  "Shared_AORe" = blue_shades[1]
)

# Custom labels for facet strip
facet_labels <- c(
  "AORe" = "Early rediploidization",
  "LORe" = "Late rediploidization"
)

# Plot with ggplot2: facet by 'redip', color by 'cons'
ggplot(active_pro_vs_jsd_correlation, aes(x = avg_perc, y = mean_jsd)) +
  geom_point(
    aes(fill = cons_redip), 
    shape = 21, size = 1.5, stroke = 0.1, color = "grey45", alpha = 0.99
  ) +
  geom_smooth(
    method = "lm", se = FALSE, 
    aes(group = cons, color = cons_redip),
    linetype = "solid", linewidth = 0.5
  ) +
  facet_wrap(~redip, ncol = 1, scales = "free_y", labeller = as_labeller(facet_labels)) +
  scale_fill_manual(values = fill_colors_jsd) +
  scale_color_manual(values = fill_colors_jsd) +
  theme_classic(base_size = 9) +
  theme(
    strip.background = element_rect(fill = "white", color = "black"),
    panel.border = element_rect(color = "black", fill = NA),
    axis.text.x = element_text(vjust = 0.5, size = 7),
    axis.text.y = element_text(size = 7),
    axis.title.x = element_text(size = 8),
    axis.title.y = element_text(size = 8),
    legend.position = "none"
  ) +
  labs(
    x = "Proportion of active promoters",
    y = "JSD (2:2 ohnologs)",
    color = "Conservation (cons)"
  )

ggsave(paste0("F4_1and2way_active_promoter_conservation_JSD_correlation_redip_2025_V1.tiff"), dpi = 600, width = 6.2, height = 9.9, units = "cm") 

ggsave(filename = "F4_1and2way_active_promoter_conservation_JSD_correlation_redip_2025_V1.pdf",  plot = last_plot(), device = cairo_pdf, width = 6.2, height = 9.9, dpi = 600,units = "cm")



###############################################################################
# CNE ENHANCER OVERLAP ANALYSIS IN LORe AND AORe REGIONS
#
# Purpose:
#   1. Load Pike CNE coordinates and calculate age-specific genome coverage.
#   2. Classify Atlantic salmon orthologous regions as LORe or AORe.
#   3. Retain only high-confidence CNEs from the MAF filtering workflow.
#   4. Quantify active enhancer overlap across:
#        - Shared CNEs (2-way conserved)
#        - Alignable CNEs (2→1 conserved)
#        - Exclusive CNEs (1-way conserved)
#   5. Normalise enhancer coverage by total available CNE sequence length.
###############################################################################

library(tidyverse)
library(fuzzyjoin)
library(openxlsx)
library(tidygenomics)
library(janitor)

###############################################################################
# SECTION 1. LOAD PIKE CNE COORDINATES
###############################################################################

# Pike CNE catalogue generated from previous analyses
pike_CNE <- "./CNEs_pike_2023.03.23_coordinates.txt"

pike_cne_coords <- read.table(
  pike_CNE,
  col.names = c(
    "chromosome",
    "start",
    "end",
    "cne_id",
    "age"
  )
)

# Number of unique CNEs
pike_cne_coords %>%
  distinct(cne_id)

###############################################################################
# SECTION 2. SUMMARISE PIKE CNE COUNTS AND GENOMIC COVERAGE
###############################################################################

# Count CNEs per evolutionary age class
pike_counts <- pike_cne_coords %>%
  count(age)

# Total CNE length per age class
pike_length <- pike_cne_coords %>%
  mutate(length = end - start) %>%
  group_by(age) %>%
  summarise(
    length_sum = sum(length),
    .groups = "drop"
  )

# Overall genomic coverage of all Pike CNEs
pike_overall_length <- pike_cne_coords %>%
  mutate(length = end - start) %>%
  summarise(
    length_sum = sum(length)
  )

# Combined summary table
pike_cne_summary <- pike_counts %>%
  left_join(pike_length, by = "age")

write.xlsx(
  pike_cne_summary,
  "pike_CNE_length_Count_summary.xlsx"
)

###############################################################################
# SECTION 3. LOAD PIKE ENHANCER OVERLAP DATASETS
###############################################################################

# All files generated from enhancer overlap analyses
files <- list.files(
  pattern = "_pike_sequences.xlsx$",
  recursive = TRUE
)

all_files <- lapply(files, read.xlsx)

pike_numbers <- bind_rows(all_files)

###############################################################################
# SECTION 4. MAP ENHANCER OVERLAPS BACK TO PIKE CNE COORDINATES
###############################################################################

AllDat_maf_pike <- bind_rows(all_files) %>%
  genome_left_join(
    pike_cne_coords,
    by = c("chromosome", "start", "end")
  ) %>%
  mutate(
    overlap  = pmin(end.x, end.y) -
      pmax(start.x, start.y),
    cne_size = end.y - start.y
  )

# Create unique genomic region identifier
AllDat_maf_pike <- AllDat_maf_pike %>%
  unite(
    region,
    chromosome.y,
    start.y,
    end.y,
    sep = "_",
    remove = FALSE
  )

rm(all_files)

###############################################################################
# SECTION 5. LOAD ATLANTIC SALMON LORe REGIONS
###############################################################################

# Late rediploidized regions following salmonid WGD
salmon_LORe <- read.table(
  "./Ssal_late_rediploidized_regions.tsv",
  header = TRUE
)

salmon_LORe <- salmon_LORe %>%
  rename(chromosome = chrom)

# Convert chromosome names to salmon genome format
salmon_LORe$chromosome <- sprintf(
  "%02d",
  salmon_LORe$chromosome
)

salmon_LORe$chromosome <- paste0(
  "ssa",
  salmon_LORe$chromosome
)

salmon_LORe$start <- as.numeric(salmon_LORe$start)
salmon_LORe$end   <- as.numeric(salmon_LORe$end)

###############################################################################
# SECTION 6. ASSIGN PIKE CNEs TO LORe OR AORe REGIONS
###############################################################################

# Extract Atlantic salmon coordinates encoded in Origin
pike_lore_aore_list <- AllDat_maf_pike %>%
  rename_with(
    ~ str_remove(., "\\.x$"),
    ends_with(".x")
  ) %>%
  mutate(
    Origin1 = str_remove(
      Origin,
      "\\.maf$"
    )
  ) %>%
  separate(
    Origin1,
    into = c(
      "chromosome_ssa",
      "start_ssa",
      "end_ssa"
    ),
    sep = "_",
    convert = TRUE
  ) %>%
  mutate(
    chromosome_ssa = str_remove(
      chromosome_ssa,
      "^ssa"
    )
  ) %>%
  mutate(
    chromosome_ssa = as.numeric(
      chromosome_ssa
    )
  )

# Standardise chromosome naming
pike_lore_aore_list$chromosome_ssa <- sprintf(
  "%02d",
  pike_lore_aore_list$chromosome_ssa
)

pike_lore_aore_list$chromosome_ssa <- paste0(
  "ssa",
  pike_lore_aore_list$chromosome_ssa
)

###############################################################################
# Fuzzy genomic overlap join against LORe coordinates.
###############################################################################

pike_lore_aore_list <- pike_lore_aore_list %>%
  rename(
    chromsom_el = chromosome.y,
    start_el    = start.y,
    end_el      = end.y
  ) %>%
  fuzzy_left_join(
    salmon_LORe,
    by = c(
      "chromosome_ssa" = "chromosome",
      "start_ssa"      = "start",
      "end_ssa"        = "end"
    ),
    match_fun = list(
      `==`,
      `>=`,
      `<=`
    )
  ) %>%
  mutate(
    redip = if_else(
      !is.na(start.y),
      "LORe",
      "AORe"
    )
  ) %>%
  select(-ends_with(".y")) %>%
  rename_with(
    ~ str_remove(., "\\.x$"),
    ends_with(".x")
  )

# Summary of LORe vs AORe CNEs
pike_lore_aore_list %>%
  distinct(cne_id, age, redip) %>%
  count(age, redip)

###############################################################################
# SECTION 7. LOAD HIGH-CONFIDENCE CNE SETS
#
# These CNEs passed the alignment-quality filtering workflow performed in the
# previous script. Only retained CNEs are used for downstream analyses.
#
# Categories:
#   Shared      = conserved in both salmon duplicates (2:2)
#   Alignable   = conserved in one duplicate pair (2:1)
#   Exclusive   = conserved in a single copy (1:1)
###############################################################################

shared_cnes_filtered <- read.table(
  "./cne_mafs/shared_CNEs_post_filtering.tsv",
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE
) %>%
  select(cne_id) %>%
  distinct(cne_id) %>%
  mutate(retain = "YES")

alignable_cnes_filtered <- read.table(
  "./cne_mafs/alignable_CNEs_post_filtering.tsv",
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE
) %>%
  select(cne_id) %>%
  distinct(cne_id) %>%
  mutate(retain = "YES")

exclusive_cnes_filtered <- read.table(
  "./cne_mafs/exclusive_CNEs_post_filtering.tsv",
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE
) %>%
  select(cne_id) %>%
  distinct(cne_id) %>%
  mutate(retain = "YES")

# Master list of all retained CNEs
valid_cnes <- bind_rows(
  shared_cnes_filtered,
  alignable_cnes_filtered,
  exclusive_cnes_filtered
)

###############################################################################
# SECTION 8. CREATE WEIGHTED CNE REFERENCE SET
#
# Shared CNEs are represented twice in Atlantic salmon due to retention of
# both post-WGD copies.
# 
#
# This weighting preserves duplicated regulatory sequence space when
# calculating counts and genomic coverage.
###############################################################################

pike_lore_aore_list2 <- pike_lore_aore_list %>%
  ungroup() %>%
  left_join(valid_cnes, by = "cne_id") %>%
  filter(retain == "YES") %>%
  mutate(
    weight = ifelse(
      cne_id %in% shared_cnes_filtered$cne_id,
      2,
      1
    )
  )

###############################################################################
# SECTION 9. TOTAL AVAILABLE CNE SEQUENCE SPACE
#
# These values provide denominators for later normalisation of enhancer
# overlap frequencies.
###############################################################################

# Weighted CNE counts
pike_lore_aore_list_counts <- pike_lore_aore_list2 %>%
  distinct(cne_id, age, redip, weight) %>%
  group_by(age, redip) %>%
  summarise(
    count = sum(weight),
    .groups = "drop"
  )

# Weighted genomic coverage
pike_lore_aore_list_length <- pike_lore_aore_list2 %>%
  distinct(
    cne_id,
    age,
    redip,
    cne_size,
    weight
  ) %>%
  group_by(age, redip) %>%
  summarise(
    length_sum = sum(cne_size * weight),
    .groups = "drop"
  )

###############################################################################
# SECTION 10. LOAD ACTIVE ENHANCER OVERLAP DATASETS
#
# Naming convention:
#
#   22 = Shared conservation (2:2)
#   21 = Alignable conservation (2:1)
#   11 = Exclusive conservation (1:1)
###############################################################################

shared_enhancers <- read.table(
  "shared_active_enhancers_aore_lore_list.txt",
  header = TRUE
)

exclusive_enhancers <- read.table(
  "exclusive_active_enhancers_aore_lore_list.txt",
  header = TRUE
)

alignable_enhancers <- read.table(
  "alignable_active_enhancers_aore_lore_list.txt",
  header = TRUE
)

###############################################################
##### overlap with active enhancers
#############################################################

##shared elements
dup_active_enh_cons_redip<-dups_conserved_props_active_enh_lore_aore_list %>% 
  left_join(AllDat_maf_pike,by="Origin")%>% mutate(age=as.character(age)) %>% left_join(valid_CNEs, by="cne_id")

dup_active_enh_combs22_redip <- dup_active_enh_cons_redip %>%
  mutate(age = ifelse(is.na(age), 'NO_ovrlap_CNE', age)) %>%
  filter(!(age == 'NO_ovrlap_CNE')) %>%
  filter(retain == "YES") %>%
  mutate(
    weight = ifelse(cne_id %in% curated_shared_CNEs$cne_id, 2, 1),   # weight duplicates
    overlap_rate = overlap / cne_size
  ) %>%
  filter(overlap_rate == 1) %>%
  group_by(stage, age, redip) %>%
  summarise(length_total = sum(cne_size * weight), .groups = "drop") %>%  # weighted sum
  left_join(pike_lore_aore_list_length, by = c("age", "redip")) %>%
  group_by(age, redip) %>%
  mutate(normalise = (length_total / length_sum) * 100) %>%
  separate(stage, into = c("tissue", "maturity", "sex"), sep = "_", remove = F) %>%
  replace(is.na(.), "embryo") %>%
  group_by(redip, age, tissue, maturity) %>%
  summarise(avg_norm = round(mean(normalise), 2), .groups = "drop")


cne_active_enh_counts_size_22_redip <- dup_active_enh_cons_redip %>%
  filter(!is.na(chromosome.y)) %>%
  mutate(
    overlap_rate = overlap / cne_size,
    weight = ifelse(cne_id %in% curated_shared_CNEs$cne_id, 2, 1)  # weight duplicates
  ) %>%
  filter(overlap_rate == 1, retain == "YES") %>%
  select(redip, age, cne_id, cne_size, weight) %>%
  distinct(cne_id, .keep_all = TRUE) %>%
  group_by(age, redip) %>%
  summarise(
    n    = sum(weight),                    # count with duplication weight
    size = sum(cne_size * weight, na.rm = TRUE),  # length with duplication weight
    .groups = "drop"
  )

sum(cne_active_enh_counts_size_22_redip$size)

###exclusive elements

sing_active_enh_cons_redip<-sing_conserved_props_active_enh_lore_aore_list %>% 
  left_join(AllDat_maf_pike,by="Origin")%>% mutate(age=as.character(age)) %>% left_join(valid_CNEs, by="cne_id")

sing_active_enh_combs11_redip<-sing_active_enh_cons_redip %>% mutate(age = ifelse(is.na(age), 'NO_ovrlap_CNE', age)) %>%  filter(retain=="YES") %>%
  filter(!(age=='NO_ovrlap_CNE')) %>% 
  mutate(overlap_rate=overlap/cne_size) %>% filter(overlap_rate==1) %>%
  group_by(stage,age,redip)%>%
  summarise(length_total=sum(cne_size)) %>% ungroup %>% left_join(pike_lore_aore_list_length, by=c("age","redip")) %>% group_by(age,redip) %>%
  mutate(normalise=(length_total/length_sum)*100) %>% 
  separate(stage,into=c("tissue","maturity","sex"),sep="_",remove=F) %>% replace(is.na(.), "embryo") %>%
  group_by(redip,age,tissue,maturity) %>% summarise(avg_norm=round(mean(normalise),2)) %>% ungroup()

cne_active_enh_counts_11_redip<- sing_active_enh_cons_redip %>% filter(!(is.na(chromosome.y))) %>% mutate(overlap_rate=overlap/cne_size) %>% filter(overlap_rate==1) %>% filter(retain=="YES") %>%
  select(redip,age,cne_id,cne_size) %>% distinct(cne_id, .keep_all = TRUE) %>% 
  group_by(age,redip) %>% 
  summarise(
    n = n(),             # Count of unique CNEs
    size = sum(cne_size, na.rm = TRUE),  # Sum of cne_size
    .groups = "drop"
  )

## alignable elements
twoone_active_enh_cons_redip<-twoone_conserved_props_active_enh_lore_aore_list %>% 
  left_join(AllDat_maf_pike,by="Origin")%>% mutate(age=as.character(age)) %>% left_join(valid_CNEs, by="cne_id")

twoone_active_enh_combs21_redip<-twoone_active_enh_cons_redip %>% mutate(age = ifelse(is.na(age), 'NO_ovrlap_CNE', age)) %>% 
  filter(!(age=='NO_ovrlap_CNE')) %>% filter(retain=="YES") %>%
  mutate(overlap_rate=overlap/cne_size) %>% filter(overlap_rate==1) %>%
  group_by(stage,age,redip)%>%
  summarise(length_total=sum(cne_size)) %>% ungroup %>% left_join(pike_lore_aore_list_length, by=c("age","redip")) %>% group_by(age,redip) %>%
  mutate(normalise=(length_total/length_sum)*100) %>%  
  separate(stage,into=c("tissue","maturity","sex"),sep="_",remove=F) %>% replace(is.na(.), "embryo") %>%
  group_by(redip,age,tissue,maturity) %>% summarise(avg_norm=round(mean(normalise),2)) %>% ungroup()

cne_active_enh_counts_21_redip<- twoone_active_enh_cons_redip %>% filter(!(is.na(chromosome.y))) %>% mutate(overlap_rate=overlap/cne_size) %>% 
  filter(overlap_rate==1) %>% filter(retain=="YES") %>%
  select(redip,age,cne_id,cne_size) %>% distinct(cne_id, .keep_all = TRUE) %>% 
  group_by(age,redip) %>% 
  summarise(
    n = n(),             # Count of unique CNEs
    size = sum(cne_size, na.rm = TRUE),  # Sum of cne_size
    .groups = "drop"
  )

####################################################################################
### merge all the outputs and summarise
####################################################################################

cne_peak_categories_size_length_redip <- cne_active_enh_counts_size_22_redip %>% 
  left_join(cne_active_enh_counts_21_redip, by = c("age","redip"), suffix = c("_shared","_alignable")) %>% 
  left_join(cne_active_enh_counts_11_redip, by = c("age","redip")) %>%
  rename(
    n_exclusive    = n,
    size_exclusive = size
  ) %>%
  mutate(
    total_counts    = n_shared + n_alignable + n_exclusive,       # sum of counts
    length_sum = size_shared + size_alignable + size_exclusive # sum of lengths
  )

write.xlsx(cne_peak_categories_size_length_redip,"cne_peak_Categories_size_length_redip_v2.xlsx")

###############################################################################
# Test AORe vs LORe enrichment of active enhancer-associated CNEs
###############################################################################

# Total genomic space assigned to each rediploidisation category
aore_genome_bp <- 2023914877   # Total AORe sequence (bp)
lore_genome_bp <- 649544050    # Total LORe sequence (bp)

# ---------------------------------------------------------------------------
# Fisher's exact test:
# Is enhancer-associated CNE sequence enriched in AORe or LORe
# relative to available genomic space?
#
# Performed independently for each evolutionary age class.
# ---------------------------------------------------------------------------

results_cne_enrichment <- cne_peak_categories_size_length_redip %>%
  select(age, redip, length_sum) %>%
  pivot_wider(
    names_from  = redip,
    values_from = length_sum,
    values_fill = 0
  ) %>%
  mutate(
    length_AORe = AORe,
    length_LORe = LORe
  ) %>%
  rowwise() %>%
  mutate(
    fisher_test = list(
      tryCatch(
        fisher.test(
          matrix(
            c(
              length_AORe,
              aore_genome_bp - length_AORe,
              length_LORe,
              lore_genome_bp - length_LORe
            ),
            nrow = 2,
            byrow = TRUE
          )
        ),
        error = function(e) NA
      )
    )
  ) %>%
  ungroup() %>%
  mutate(
    p.value = map_dbl(
      fisher_test,
      ~ ifelse(is.list(.x), broom::tidy(.x)$p.value, NA_real_)
    ),
    odds_ratio = map_dbl(
      fisher_test,
      ~ ifelse(is.list(.x), broom::tidy(.x)$estimate, NA_real_)
    )
  ) %>%
  mutate(
    enrichment = case_when(
      p.value < 0.05 & odds_ratio > 1 ~ "AORe enriched",
      p.value < 0.05 & odds_ratio < 1 ~ "LORe enriched",
      TRUE                            ~ "NS"
    ),
    p.adj = p.adjust(p.value, method = "fdr")
  ) %>%
  select(
    age,
    length_AORe,
    length_LORe,
    odds_ratio,
    p.value,
    p.adj,
    enrichment
  )

# Format values for reporting
results_cne_enrichment <- results_cne_enrichment %>%
  mutate(
    odds_ratio = signif(odds_ratio, 4),
    p.value    = signif(p.value, 3),
    p.adj      = signif(p.adj, 3)
  )

write.xlsx(
  results_cne_enrichment,
  "cne_enrichment_genomewide_fishers.xlsx"
)

###############################################################################
# Repeat AORe/LORe enrichment tests separately for
# Shared, Alignable and Exclusive CNE categories
###############################################################################

# Convert category-specific sequence lengths into long format

df_long <- cne_peak_categories_size_length_redip %>%
  select(
    age,
    redip,
    size_shared,
    size_alignable,
    size_exclusive
  ) %>%
  pivot_longer(
    cols      = starts_with("size_"),
    names_to  = "category",
    values_to = "cne_bp"
  ) %>%
  mutate(
    category = gsub("size_", "", category)
  )

# ---------------------------------------------------------------------------
# Fisher's exact test for each:
#
# age × conservation category
# ---------------------------------------------------------------------------

results_cne_age_category <- df_long %>%
  pivot_wider(
    names_from  = redip,
    values_from = cne_bp,
    values_fill = 0
  ) %>%
  rowwise() %>%
  mutate(
    fisher_test = list(
      tryCatch(
        fisher.test(
          matrix(
            c(
              AORe,
              aore_genome_bp - AORe,
              LORe,
              lore_genome_bp - LORe
            ),
            nrow = 2,
            byrow = TRUE
          )
        ),
        error = function(e) NA
      )
    )
  ) %>%
  ungroup() %>%
  mutate(
    p.value = map_dbl(
      fisher_test,
      ~ ifelse(is.list(.x), broom::tidy(.x)$p.value, NA_real_)
    ),
    odds_ratio = map_dbl(
      fisher_test,
      ~ ifelse(is.list(.x), broom::tidy(.x)$estimate, NA_real_)
    )
  ) %>%
  mutate(
    enrichment = case_when(
      p.value < 0.05 & odds_ratio > 1 ~ "AORe enriched",
      p.value < 0.05 & odds_ratio < 1 ~ "LORe enriched",
      TRUE                            ~ "NS"
    ),
    p.adj = p.adjust(p.value, method = "fdr")
  ) %>%
  select(
    age,
    category,
    AORe,
    LORe,
    odds_ratio,
    p.value,
    p.adj,
    enrichment
  ) %>%
  mutate(
    odds_ratio = signif(odds_ratio, 4),
    p.value    = signif(p.value, 3),
    p.adj      = signif(p.adj, 3)
  ) %>%
  arrange(category, age)

###############################################################################
# Build enhancer activity datasets across tissues and ontogeny
###############################################################################

# ---------------------------------------------------------------------------
# Combine activity summaries from:
#
# 2-to-2 conserved enhancers  (Shared)
# 2-to-1 conserved enhancers  (Alignable)
# 1-to-1 conserved enhancers  (Exclusive)
# ---------------------------------------------------------------------------

one_to_one_active_e_redip <- sing_active_enh_combs11_redip %>%
  mutate(type = "1_to_1") %>%
  unite(
    col = "tissue_maturity",
    tissue,
    maturity,
    sep = "_",
    remove = FALSE
  )

two_to_two_active_e_redip <- dup_active_enh_combs22_redip %>%
  mutate(type = "2_to_2") %>%
  unite(
    col = "tissue_maturity",
    tissue,
    maturity,
    sep = "_",
    remove = FALSE
  )

two_to_one_active_e_redip <- twoone_active_enh_combs21_redip %>%
  mutate(type = "2_to_1") %>%
  unite(
    col = "tissue_maturity",
    tissue,
    maturity,
    sep = "_",
    remove = FALSE
  )

# Pool all conservation classes

pooled_combs_active_e_redip <- bind_rows(
  one_to_one_active_e_redip,
  two_to_two_active_e_redip,
  two_to_one_active_e_redip
)

# Standardise ordering for plotting

pooled_combs_active_e_redip$age <- factor(
  pooled_combs_active_e_redip$age,
  levels = c(
    "Euteleostomi",
    "Actinopteri",
    "Neopterygii",
    "Teleostei",
    "Clupeocephala",
    "Euteleosteomorpha",
    "Protacanthopterygii"
  )
)

pooled_combs_active_e_redip$type <- factor(
  pooled_combs_active_e_redip$type,
  levels = c(
    "2_to_2",
    "2_to_1",
    "1_to_1"
  )
)

# Age groups used repeatedly in downstream figures

target <- c(
  "Euteleostomi",
  "Actinopteri",
  "Teleostei",
  "Protacanthopterygii"
)

target2 <- c(
  "Neopterygii",
  "Clupeocephala",
  "Euteleosteomorpha"
)

# Developmental stages

embryos <- c(
  "lateblastulation_embryo",
  "midgastrulation_embryo",
  "earlysomitogenesis_embryo",
  "midsomitogenesis_embryo",
  "latesomitogenesis_embryo"
)

# Adult tissues

tissues <- c(
  "brain_immature",
  "brain_mature",
  "liver_immature",
  "liver_mature",
  "muscle_immature",
  "muscle_mature",
  "ovary_immature",
  "ovary_mature",
  "testis_immature",
  "testis_mature"
)

#################################################################################################
###Fishers test for AORe data across the three conservation categories each at a time
### do fishers test in all combinations for tissues for active enhancers
#################################################################################################

pike_lengths_true<-cne_peak_categories_size_length_redip %>% select(age,redip,length_sum)

######################################
#### prepare shared, exclusive and alignable data first
######################################

dup_active_enh_combs22_redip_fishers <- dup_active_enh_cons_redip %>%
  ungroup() %>%
  mutate(age = ifelse(is.na(age), 'NO_ovrlap_CNE', age)) %>%
  filter(age != 'NO_ovrlap_CNE', retain == "YES") %>%       # keep only retained CNEs
  mutate(
    overlap_rate = overlap / cne_size,
    weight       = ifelse(cne_id %in% curated_shared_CNEs$cne_id, 2, 1)  # 2 for curated shared
  ) %>%
  filter(overlap_rate == 1) %>%
  group_by(stage, age, redip) %>%
  summarise(
    length_total = sum(cne_size * weight),   # weighted sum of CNE lengths
    .groups = "drop"
  ) %>%
  left_join(pike_lengths_true, by = c("age","redip")) %>%
  separate(stage, into = c("tissue","maturity","sex"), sep = "_", remove = FALSE) %>%
  replace(is.na(.), "embryo") %>%
  group_by(age, tissue, maturity, redip) %>%
  summarise(
    avg_CNE  = round(mean(length_total), 0),  # mean of weighted lengths
    avg_age  = round(mean(length_sum), 0),    # total CNE length from reference table
    .groups = "drop"
  ) %>%
  unite(col = "tissue", tissue, maturity, sep = "_") %>%
  mutate(category = "Shared")


sing_active_enh_combs11_redip_fishers<-sing_active_enh_cons_redip %>% ungroup() %>% 
  mutate(age = ifelse(is.na(age), 'NO_ovrlap_CNE', age)) %>% 
  filter(!(age=='NO_ovrlap_CNE')) %>% 
  mutate(overlap_rate=overlap/cne_size) %>% 
  filter(overlap_rate==1) %>%
  group_by(stage,age,redip)%>%
  summarise(length_total=sum(cne_size)) %>% 
  ungroup %>% 
  left_join(pike_lengths_true, by=c("age","redip")) %>% 
  separate(stage,into=c("tissue","maturity","sex"),sep="_",remove=F) %>% 
  replace(is.na(.), "embryo") %>%
  group_by(age,tissue,maturity,redip) %>%
  summarise(avg_CNE=round(mean(length_total),0),avg_age = round(mean(length_sum),0)) %>% 
  ungroup() %>%
  unite(col="tissue", tissue, maturity, sep ="_") %>% 
  mutate(category="Exclusive")

twoone_active_enh_combs21_redip_fishers<-twoone_active_enh_cons_redip %>% 
  ungroup() %>% mutate(age = ifelse(is.na(age), 'NO_ovrlap_CNE', age)) %>% 
  filter(!(age=='NO_ovrlap_CNE')) %>% 
  mutate(overlap_rate=overlap/cne_size) %>% 
  filter(overlap_rate==1) %>%
  group_by(stage,age,redip)%>%
  summarise(length_total=sum(cne_size)) %>% 
  ungroup %>% 
  left_join(pike_lengths_true, by=c("age","redip")) %>% 
  separate(stage,into=c("tissue","maturity","sex"),sep="_",remove=F) %>% 
  replace(is.na(.), "embryo") %>%
  group_by(age,tissue,maturity,redip) %>%
  summarise(avg_CNE=round(mean(length_total),0),avg_age = round(mean(length_sum),0)) %>% 
  ungroup() %>%
  unite(col="tissue", tissue, maturity, sep ="_") %>% 
  mutate(category="Alignable")

CNE_activity_through_ontogeny<-bind_rows(dup_active_enh_combs22_redip_fishers,
                                         twoone_active_enh_combs21_redip_fishers,
                                         sing_active_enh_combs11_redip_fishers)

write.xlsx(CNE_activity_through_ontogeny,"cne_activity_through_ontogeny.xlsx")

##################################################################################
############generate the same activity data but in counts and not length
##################################################################################

pike_counts_true<-cne_peak_categories_size_length_redip %>% select(age,redip,total_counts)

dup_active_enh_combs22_redip_fishers_counts <- dup_active_enh_cons_redip %>%
  ungroup() %>%
  mutate(age = ifelse(is.na(age), 'NO_ovrlap_CNE', age)) %>%
  filter(age != 'NO_ovrlap_CNE', retain == "YES") %>%       # keep only retained CNEs
  mutate(
    overlap_rate = overlap / cne_size,
    weight       = ifelse(cne_id %in% curated_shared_CNEs$cne_id, 2, 1)  # 2 for curated shared
  ) %>%
  filter(overlap_rate == 1) %>%
  group_by(stage, age, redip) %>%
  summarise(
    counts_total = sum(weight),   # weighted counts
    .groups = "drop"
  ) %>%
  left_join(pike_counts_true, by = c("age","redip")) %>%
  separate(stage, into = c("tissue","maturity","sex"), sep = "_", remove = FALSE) %>%
  replace(is.na(.), "embryo") %>%
  group_by(age, tissue, maturity, redip) %>%
  summarise(
    CNE_counts  = round(mean(counts_total), 0),  # mean of weighted counts
    total_CNE_counts  = round(mean(total_counts), 0),    # total CNE counts from reference table
    .groups = "drop"
  ) %>%
  unite(col = "tissue", tissue, maturity, sep = "_") %>%
  mutate(category = "Shared")


sing_active_enh_combs11_redip_fishers_counts<-sing_active_enh_cons_redip %>% 
  ungroup() %>% mutate(age = ifelse(is.na(age), 'NO_ovrlap_CNE', age)) %>% 
  filter(!(age=='NO_ovrlap_CNE')) %>% 
  mutate(overlap_rate=overlap/cne_size) %>% 
  filter(overlap_rate==1) %>%
  group_by(stage,age,redip)%>%
  summarise(cne_counts=n()) %>% ungroup %>% 
  left_join(pike_counts_true, by=c("age","redip")) %>% 
  separate(stage,into=c("tissue","maturity","sex"),sep="_",remove=F) %>% 
  replace(is.na(.), "embryo") %>%
  group_by(age,tissue,maturity,redip) %>%
  summarise(CNE_counts=round(mean(cne_counts),0),total_CNE_counts = round(mean(total_counts),0)) %>% 
  ungroup() %>%
  unite(col="tissue", tissue, maturity, sep ="_") %>% 
  mutate(category="Exclusive")

twoone_active_enh_combs21_redip_fishers_counts<-twoone_active_enh_cons_redip %>% 
  ungroup() %>% mutate(age = ifelse(is.na(age), 'NO_ovrlap_CNE', age)) %>% 
  filter(!(age=='NO_ovrlap_CNE')) %>% 
  mutate(overlap_rate=overlap/cne_size) %>% 
  filter(overlap_rate==1) %>%
  group_by(stage,age,redip)%>%
  summarise(cne_counts=n()) %>% ungroup %>% 
  left_join(pike_counts_true, by=c("age","redip")) %>% 
  separate(stage,into=c("tissue","maturity","sex"),sep="_",remove=F) %>% 
  replace(is.na(.), "embryo") %>%
  group_by(age,tissue,maturity,redip) %>%
  summarise(CNE_counts=round(mean(cne_counts),0),total_CNE_counts = round(mean(total_counts),0)) %>% 
  ungroup() %>%
  unite(col="tissue", tissue, maturity, sep ="_") %>% 
  mutate(category="Alignable")

CNE_counts_activity_through_ontogeny<-bind_rows(dup_active_enh_combs22_redip_fishers_counts,twoone_active_enh_combs21_redip_fishers_counts,sing_active_enh_combs11_redip_fishers_counts)

write.xlsx(CNE_counts_activity_through_ontogeny,"cne_counts_activity_through_ontogeny.xlsx")

overall_cne_activity <- CNE_activity_through_ontogeny %>% left_join(CNE_counts_activity_through_ontogeny, by=(c("age","tissue","redip","category")))

write.xlsx(overall_cne_activity,"cne_counts_lenths_merged_activity_through_ontogeny.xlsx")

################################################################################
#### fishers test for CNEs overlapping duplicated active enhancers
################################################################################

tissue_vec_dup_active_enh<- unique(dup_active_enh_combs22_redip_fishers$tissue)

result_fishers_dup_act_enh_redip<-lapply(tissue_vec_dup_active_enh, function(a1) {
  
  t_group<- dup_active_enh_combs22_redip_fishers %>% filter(tissue==a1) %>% 
    filter(redip=="LORe") ## run in turns for AORe and LORe, make sure you swap them
  
  idx = t(combn(seq_along(t_group$age),2))
  
  lapply(1:nrow(idx),function(i){
    
    age1_idx <- idx[i, 1]
    age2_idx <- idx[i, 2]
    
    test = fisher.test(t_group[idx[i,],c("avg_CNE","avg_age")])
    
    data.frame(
      age1 = t_group$age[age1_idx],
      tissue_specific_CNE_1 = t_group$avg_CNE[age1_idx],
      overall_CNE_length_1 = t_group$avg_age[age1_idx],
      age2 = t_group$age[age2_idx],
      tissue_specific_CNE_2 = t_group$avg_CNE[age2_idx],
      overall_CNE_length_2 = t_group$avg_age[age2_idx],
      odds_ratio = as.numeric(test$estimate),
      p = as.numeric(test$p.value)
      
      
    )
  })
  
})

names(result_fishers_dup_act_enh_redip)<-tissue_vec_dup_active_enh

##merge all datasets 

# Create an empty list to store merged data frames
merged_dfs_dup_act_enh_redip <- list()

# Iterate over each element in my_list
for (i in seq_along(result_fishers_dup_act_enh_redip)) {
  # Extract the sublists for each letter
  sublists <- result_fishers_dup_act_enh_redip[[i]]
  
  # Merge the data frames within each sublist
  merged_df <- do.call(bind_rows, sublists)
  
  # Add an additional column with the name of i
  merged_df <- mutate(merged_df, tissue = names(result_fishers_dup_act_enh_redip)[i])
  
  # Store the merged data frame in merged_dfs list
  merged_dfs_dup_act_enh_redip[[i]] <- merged_df
  
}

# Combine all merged data frames into a single data frame

final_merged_results_fishers_duplicate_all_active_enchancers_AORe <- do.call(bind_rows, merged_dfs_dup_act_enh_redip)

final_merged_results_fishers_duplicate_all_active_enchancers_AORe <- final_merged_results_fishers_duplicate_all_active_enchancers_AORe %>% 
  mutate(
    element = "Active Enhancer",
    redip = "Early rediploidization",
    conservation="Shared",
    p_bonferroni = p * n(),  # Bonferroni correction
    p_bonferroni = ifelse(p_bonferroni > 1, 1, p_bonferroni),  # cap at 1
    significant_pre_correction = ifelse(p < 0.05, "yes", "no"),
    significant_post_correction = ifelse(p_bonferroni < 0.05, "yes", "no")
  )

final_merged_results_fishers_duplicate_all_active_enchancers_LORe <- do.call(bind_rows, merged_dfs_dup_act_enh_redip)

final_merged_results_fishers_duplicate_all_active_enchancers_LORe <- final_merged_results_fishers_duplicate_all_active_enchancers_LORe %>% 
  mutate(
    element = "Active Enhancer",
    redip = "Late rediploidization",
    conservation="Shared",
    p_bonferroni = p * n(),  # Bonferroni correction
    p_bonferroni = ifelse(p_bonferroni > 1, 1, p_bonferroni),  # cap at 1
    significant_pre_correction = ifelse(p < 0.05, "yes", "no"),
    significant_post_correction = ifelse(p_bonferroni < 0.05, "yes", "no")
  )

final_merged_results_fishers_duplicate_all_active_enchancers_redip<-bind_rows(final_merged_results_fishers_duplicate_all_active_enchancers_AORe,
                                                                              final_merged_results_fishers_duplicate_all_active_enchancers_LORe)


write.xlsx(final_merged_results_fishers_duplicate_all_active_enchancers_redip,
           "Shared_active_enhancers_CNE_conservation_across_stages_fishers_test_2025_v2.xlsx")

######################################################################################
### fishers test for CNEs overlapping exclusive active enhancers
##############################################################################################

tissue_vec_sing_active_enh<- sing_active_enh_combs11_redip_fishers %>% filter(redip == "LORe") %>% distinct(tissue) %>% pull(tissue)

result_fishers_sing_act_enh_redip <- lapply(tissue_vec_sing_active_enh, function(a1) {
  
  t_group <- sing_active_enh_combs11_redip_fishers %>% 
    filter(tissue == a1, redip == "LORe") ### make sure to swap LORe and AORe
  
  if (nrow(t_group) < 2) return(NULL)
  
  idx <- t(combn(seq_along(t_group$age), 2))
  
  lapply(1:nrow(idx), function(i) {
    age1_idx <- idx[i, 1]
    age2_idx <- idx[i, 2]
    
    test <- fisher.test(t_group[idx[i, ], c("avg_CNE", "avg_age")])
    
    data.frame(
      age1 = t_group$age[age1_idx],
      tissue_specific_CNE_1 = t_group$avg_CNE[age1_idx],
      overall_CNE_length_1 = t_group$avg_age[age1_idx],
      age2 = t_group$age[age2_idx],
      tissue_specific_CNE_2 = t_group$avg_CNE[age2_idx],
      overall_CNE_length_2 = t_group$avg_age[age2_idx],
      odds_ratio = as.numeric(test$estimate),
      p = as.numeric(test$p.value)
    )
  })
})

names(result_fishers_sing_act_enh_redip)<-tissue_vec_sing_active_enh

##merge all datasets 

# Create an empty list to store merged data frames
merged_dfs_sing_act_enh_redip <- list()

for (i in seq_along(result_fishers_sing_act_enh_redip)) {
  sublists <- result_fishers_sing_act_enh_redip[[i]]
  
  # Only proceed if sublists is not NULL and is a list
  if (!is.null(sublists) && is.list(sublists)) {
    merged_df <- dplyr::bind_rows(sublists)
    merged_df <- dplyr::mutate(merged_df, tissue = names(result_fishers_sing_act_enh_redip)[i])
    merged_dfs_sing_act_enh_redip[[i]] <- merged_df
  }
}

#Combine all merged data frames into a single data frame

final_merged_results_fishers_singleton_all_active_enchancers_AORe <- do.call(bind_rows, merged_dfs_sing_act_enh_redip)

final_merged_results_fishers_singleton_all_active_enchancers_AORe <- final_merged_results_fishers_singleton_all_active_enchancers_AORe %>% 
  mutate(
    element = "Active Enhancer",
    redip = "Early rediploidization",
    conservation="Exclusive",
    p_bonferroni = p * n(),  # Bonferroni correction
    p_bonferroni = ifelse(p_bonferroni > 1, 1, p_bonferroni),  # cap at 1
    significant_pre_correction = ifelse(p < 0.05, "yes", "no"),
    significant_post_correction = ifelse(p_bonferroni < 0.05, "yes", "no")
  )

final_merged_results_fishers_singleton_all_active_enchancers_LORe <- do.call(bind_rows, merged_dfs_sing_act_enh_redip)

final_merged_results_fishers_singleton_all_active_enchancers_LORe <- final_merged_results_fishers_singleton_all_active_enchancers_LORe %>% 
  mutate(
    element = "Active Enhancer",
    redip = "Late rediploidization",
    conservation="Exclusive",
    p_bonferroni = p * n(),  # Bonferroni correction
    p_bonferroni = ifelse(p_bonferroni > 1, 1, p_bonferroni),  # cap at 1
    significant_pre_correction = ifelse(p < 0.05, "yes", "no"),
    significant_post_correction = ifelse(p_bonferroni < 0.05, "yes", "no")
  )

final_merged_results_fishers_singleton_all_active_enchancers_redip<-bind_rows(final_merged_results_fishers_singleton_all_active_enchancers_AORe,
                                                                              final_merged_results_fishers_singleton_all_active_enchancers_LORe)


write.xlsx(final_merged_results_fishers_singleton_all_active_enchancers_redip,
           "Exclusive_active_enhancers_CNE_conservation_across_stages_fishers_test_2025_v2.xlsx")

#####################################################################################################
#### fishers test for CNEs overlapping alignable active enhancers
#####################################################################################################


tissue_vec_twoone_active_enh<- twoone_active_enh_combs21_redip_fishers %>% filter(redip == "LORe") %>% distinct(tissue) %>% pull(tissue)

result_fishers_twoone_act_enh_redip <- lapply(tissue_vec_twoone_active_enh, function(a1) {
  
  t_group <- twoone_active_enh_combs21_redip_fishers %>% 
    filter(tissue == a1, redip == "LORe") ## swap AORe and LORe here
  
  if (nrow(t_group) < 2) return(NULL)
  
  idx <- t(combn(seq_along(t_group$age), 2))
  
  lapply(1:nrow(idx), function(i) {
    age1_idx <- idx[i, 1]
    age2_idx <- idx[i, 2]
    
    test <- fisher.test(t_group[idx[i, ], c("avg_CNE", "avg_age")])
    
    data.frame(
      age1 = t_group$age[age1_idx],
      tissue_specific_CNE_1 = t_group$avg_CNE[age1_idx],
      overall_CNE_length_1 = t_group$avg_age[age1_idx],
      age2 = t_group$age[age2_idx],
      tissue_specific_CNE_2 = t_group$avg_CNE[age2_idx],
      overall_CNE_length_2 = t_group$avg_age[age2_idx],
      odds_ratio = as.numeric(test$estimate),
      p = as.numeric(test$p.value)
    )
  })
})

names(result_fishers_twoone_act_enh_redip)<-tissue_vec_twoone_active_enh

##merge all datasets 

# Create an empty list to store merged data frames
merged_dfs_twoone_act_enh_redip <- list()

for (i in seq_along(result_fishers_twoone_act_enh_redip)) {
  sublists <- result_fishers_twoone_act_enh_redip[[i]]
  
  # Only proceed if sublists is not NULL and is a list
  if (!is.null(sublists) && is.list(sublists)) {
    merged_df <- dplyr::bind_rows(sublists)
    merged_df <- dplyr::mutate(merged_df, tissue = names(result_fishers_twoone_act_enh_redip)[i])
    merged_dfs_twoone_act_enh_redip[[i]] <- merged_df
  }
}

# Combine all merged data frames into a twoone data frame

final_merged_results_fishers_alignable_all_active_enchancers_AORe <- do.call(bind_rows, merged_dfs_twoone_act_enh_redip)
final_merged_results_fishers_alignable_all_active_enchancers_AORe <- final_merged_results_fishers_alignable_all_active_enchancers_AORe %>% 
  mutate(
    element = "Active Enhancer",
    redip = "Early rediploidization",
    conservation="Alignable",
    p_bonferroni = p * n(),  # Bonferroni correction
    p_bonferroni = ifelse(p_bonferroni > 1, 1, p_bonferroni),  # cap at 1
    significant_pre_correction = ifelse(p < 0.05, "yes", "no"),
    significant_post_correction = ifelse(p_bonferroni < 0.05, "yes", "no")
  )

final_merged_results_fishers_alignable_all_active_enchancers_LORe <- do.call(bind_rows, merged_dfs_twoone_act_enh_redip)
final_merged_results_fishers_alignable_all_active_enchancers_LORe <- final_merged_results_fishers_alignable_all_active_enchancers_LORe %>% 
  mutate(
    element = "Active Enhancer",
    redip = "Late rediploidization",
    conservation="Alignable",
    p_bonferroni = p * n(),  # Bonferroni correction
    p_bonferroni = ifelse(p_bonferroni > 1, 1, p_bonferroni),  # cap at 1
    significant_pre_correction = ifelse(p < 0.05, "yes", "no"),
    significant_post_correction = ifelse(p_bonferroni < 0.05, "yes", "no")
  )

final_merged_results_fishers_alignable_all_active_enchancers_redip<-bind_rows(final_merged_results_fishers_alignable_all_active_enchancers_AORe,
                                                                              final_merged_results_fishers_alignable_all_active_enchancers_LORe)

write.xlsx(final_merged_results_fishers_alignable_all_active_enchancers_redip,
           "Alignable_active_enhancers_CNE_conservation_across_stages_fishers_test_2025_v2.xlsx")


########################################################################
###new fishers test combinations
######################################################################

dup_active_enh_combs22_redip_fishers_try2 <- dup_active_enh_cons_redip %>% #shared
  ungroup() %>%
  mutate(age = ifelse(is.na(age), 'NO_ovrlap_CNE', age)) %>%
  filter(age != 'NO_ovrlap_CNE', retain == "YES") %>%       # keep only retained CNEs
  mutate(
    overlap_rate = overlap / cne_size,
    weight       = ifelse(cne_id %in% curated_shared_CNEs$cne_id, 2, 1)  # 2 for curated shared
  ) %>%
  filter(overlap_rate == 1) %>%
  group_by(stage, age, redip) %>%
  summarise(
    length_total = sum(cne_size * weight),   # weighted sum of CNE lengths
    .groups = "drop"
  ) %>%
  separate(stage,into=c("tissue","maturity","sex"),sep="_",remove=F) %>% 
  replace(is.na(.), "embryo") %>%
  group_by(age,tissue,maturity,redip) %>%
  summarise(avg_CNE=round(mean(length_total),0)) %>% ungroup() %>%
  unite(col="tissue", tissue, maturity, sep ="_") %>%
  mutate(category = "shared")


sing_active_enh_combs11_redip_fishers_try2<-sing_active_enh_cons_redip %>%  ## exclusive
  ungroup() %>% mutate(age = ifelse(is.na(age), 'NO_ovrlap_CNE', age)) %>% 
  filter(!(age=='NO_ovrlap_CNE')) %>% 
  mutate(overlap_rate=overlap/cne_size) %>% 
  filter(overlap_rate==1) %>%
  group_by(stage,age,redip)%>%
  summarise(length_total=sum(cne_size)) %>% ungroup %>% 
  separate(stage,into=c("tissue","maturity","sex"),sep="_",remove=F) %>% 
  replace(is.na(.), "embryo") %>%
  group_by(age,tissue,maturity,redip) %>%
  summarise(avg_CNE=round(mean(length_total),0)) %>% 
  ungroup() %>%
  unite(col="tissue", tissue, maturity, sep ="_") %>% 
  mutate(category="exclusive")

twoone_active_enh_combs21_redip_fishers_try2<-twoone_active_enh_cons_redip %>% 
  ungroup() %>% mutate(age = ifelse(is.na(age), 'NO_ovrlap_CNE', age)) %>% 
  filter(!(age=='NO_ovrlap_CNE')) %>% 
  mutate(overlap_rate=overlap/cne_size) %>% 
  filter(overlap_rate==1) %>%
  group_by(stage,age,redip)%>%
  summarise(length_total=sum(cne_size)) %>% 
  ungroup %>% 
  separate(stage,into=c("tissue","maturity","sex"),sep="_",remove=F) %>% 
  replace(is.na(.), "embryo") %>%
  group_by(age,tissue,maturity,redip) %>%
  summarise(avg_CNE=round(mean(length_total),0)) %>% 
  ungroup() %>%
  unite(col="tissue", tissue, maturity, sep ="_") %>% 
  mutate(category="alignable")

##################################################################################################################
### combination 1 - Within a sample type and per age category– test for differences in conservation category 
##(e.g. is there a greater proportion of shared enhancer CNEs vs. alignable enhancer CNEs in LS with actinopteri)
###################################################################################################################

merged_data_fishers_test_try2<-bind_rows(dup_active_enh_combs22_redip_fishers_try2,sing_active_enh_combs11_redip_fishers_try2,
                                         twoone_active_enh_combs21_redip_fishers_try2)

length_group_one<-merged_data_fishers_test_try2 %>% group_by(redip,tissue,category) %>% summarise(length_total=sum(avg_CNE)) 


merged_data_fishers_test_try2<-merged_data_fishers_test_try2 %>% left_join(length_group_one, by=c("redip","tissue","category")) 

# get all unique tissues
tissue_vec <- unique(merged_data_fishers_test_try2$tissue)

# iterate over tissues
result_fishers_cat <- lapply(tissue_vec, function(tiss) {
  
  tiss_df <- merged_data_fishers_test_try2 %>% filter(tissue == tiss,redip == "AORe")
  
  # get unique ages for this tissue
  age_vec <- unique(tiss_df$age)
  
  # iterate over ages
  age_results <- lapply(age_vec, function(a) {
    
    age_df <- tiss_df %>% filter(age == a)
    
    # get all pairwise combinations of categories
    if(nrow(age_df) < 2) return(NULL)
    
    idx <- t(combn(seq_len(nrow(age_df)), 2))
    
    lapply(1:nrow(idx), function(i) {
      cat1_idx <- idx[i, 1]
      cat2_idx <- idx[i, 2]
      
      # build 2x2 contingency table
      mat <- matrix(c(
        age_df$avg_CNE[cat1_idx],
        age_df$length_total[cat1_idx] - age_df$avg_CNE[cat1_idx],
        age_df$avg_CNE[cat2_idx],
        age_df$length_total[cat2_idx] - age_df$avg_CNE[cat2_idx]
      ), nrow = 2, byrow = TRUE)
      
      test <- fisher.test(mat)
      
      data.frame(
        age = a,
        category1 = age_df$category[cat1_idx],
        CNE1 = age_df$avg_CNE[cat1_idx],
        total1 = age_df$length_total[cat1_idx],
        category2 = age_df$category[cat2_idx],
        CNE2 = age_df$avg_CNE[cat2_idx],
        total2 = age_df$length_total[cat2_idx],
        odds_ratio = as.numeric(test$estimate),
        p_value = as.numeric(test$p.value),
        tissue = tiss
      )
    })
  })
  
  # flatten age_results
  dplyr::bind_rows(unlist(age_results, recursive = FALSE))
})


# bind all tissues together for AORe data
fishers_try2_result_AOre <- dplyr::bind_rows(result_fishers_cat)
fishers_try2_result_AOre <- fishers_try2_result_AOre %>% mutate(redip="Early rediploidization")

# bind all tissues together for LORe data
fishers_try2_result_LOre <- dplyr::bind_rows(result_fishers_cat)
fishers_try2_result_LOre <- fishers_try2_result_LOre %>% mutate(redip="Late rediploidization")

fishers_try2_overall_results <- bind_rows(fishers_try2_result_AOre,fishers_try2_result_LOre)

fishers_try2_overall_results<-fishers_try2_overall_results %>% mutate(p_bonferroni = p_value * n(),  # Bonferroni correction
                                                                      p_bonferroni = ifelse(p_bonferroni > 1, 1, p_bonferroni),  # cap at 1
                                                                      significant_pre_correction = ifelse(p_value < 0.05, "yes", "no"),
                                                                      significant_post_correction = ifelse(p_bonferroni < 0.05, "yes", "no"))


write.xlsx(fishers_try2_overall_results,"sampleandage_vs_category_fishers_test_2025_v1.xlsx")


##############################################################################################################
### coombination 2: Within age category within alignment category difference across tissues 
##(e.g. for actinopterygi is the proportion of shared CNE enhancers greater in LS vs. ES)
######################################################################################################################

merged_data_fishers_test_try3<-bind_rows(dup_active_enh_combs22_redip_fishers_try2,sing_active_enh_combs11_redip_fishers_try2,
                                         twoone_active_enh_combs21_redip_fishers_try2)

length_group_two<-merged_data_fishers_test_try3 %>% group_by(redip,tissue,category) %>% summarise(length_total=sum(avg_CNE)) 


merged_data_fishers_test_try3<-merged_data_fishers_test_try3 %>% left_join(length_group_two, by=c("redip","tissue","category"))


# get all unique categories
category_vec <- unique(merged_data_fishers_test_try3$category)

result_fishers_tissue <- lapply(category_vec, function(cat) {
  
  cat_df <- merged_data_fishers_test_try3 %>% filter(category == cat,redip == "LORe")
  
  # unique ages
  age_vec <- unique(cat_df$age)
  
  age_results <- lapply(age_vec, function(a) {
    
    age_df <- cat_df %>% filter(age == a)
    
    # all pairwise combinations of tissues
    if(nrow(age_df) < 2) return(NULL)
    
    idx <- t(combn(seq_len(nrow(age_df)), 2))
    
    lapply(1:nrow(idx), function(i) {
      tissue1_idx <- idx[i, 1]
      tissue2_idx <- idx[i, 2]
      
      # 2x2 contingency table
      mat <- matrix(c(
        age_df$avg_CNE[tissue1_idx],
        age_df$length_total[tissue1_idx] - age_df$avg_CNE[tissue1_idx],
        age_df$avg_CNE[tissue2_idx],
        age_df$length_total[tissue2_idx] - age_df$avg_CNE[tissue2_idx]
      ), nrow = 2, byrow = TRUE)
      
      test <- fisher.test(mat)
      
      data.frame(
        age = a,
        tissue1 = age_df$tissue[tissue1_idx],
        CNE1 = age_df$avg_CNE[tissue1_idx],
        total1 = age_df$length_total[tissue1_idx],
        tissue2 = age_df$tissue[tissue2_idx],
        CNE2 = age_df$avg_CNE[tissue2_idx],
        total2 = age_df$length_total[tissue2_idx],
        odds_ratio = as.numeric(test$estimate),
        p_value = as.numeric(test$p.value),
        category = cat
      )
    })
  })
  
  dplyr::bind_rows(unlist(age_results, recursive = FALSE))
})

# merge all categories
final_result_tissue_AORe <- dplyr::bind_rows(result_fishers_tissue)
final_result_tissue_AORe<-final_result_tissue_AORe %>% mutate(redip="Early rediploidization")

# merge all categories
final_result_tissue_LORe <- dplyr::bind_rows(result_fishers_tissue)
final_result_tissue_LORe<-final_result_tissue_LORe %>% mutate(redip="Late rediploidization")

fishers_results_tissue_overall<-bind_rows(final_result_tissue_AORe,final_result_tissue_LORe) %>% mutate(p_bonferroni = p_value * n(),  # Bonferroni correction
                                                                                                        p_bonferroni = ifelse(p_bonferroni > 1, 1, p_bonferroni),  # cap at 1
                                                                                                        significant_pre_correction = ifelse(p_value < 0.05, "yes", "no"),
                                                                                                        significant_post_correction = ifelse(p_bonferroni < 0.05, "yes", "no"))
write.xlsx(fishers_results_tissue_overall,"categoryandage_vs_sample_fishers_test_2025_v1.xlsx")



#####################################################################################
##### plot data for figure 6 here
#########################################################################################


tissue_order <- c("lateblastulation_embryo", "midgastrulation_embryo", "earlysomitogenesis_embryo",
                  "midsomitogenesis_embryo", "latesomitogenesis_embryo", "brain_immature",
                  "brain_mature", "muscle_immature", "muscle_mature", "liver_immature",
                  "liver_mature", "ovary_immature", "ovary_mature", "testis_immature", "testis_mature")

tissue_labels <- c("LB", "MG", "ES", "MS", "LS", 
                   "Br\n(I)", "Br\n(M)", "Mu\n(I)", "Mu\n(M)",
                   "Li\n(I)", "Li\n(M)", "Ov\n(I)", "Ov\n(M)", "Te\n(I)", "Te\n(M)")
# Prepare the data
df_plot <- pooled_combs_e_active_redip %>%
  mutate(
    type = factor(type, levels = c("2_to_2","2_to_1","1_to_1")),
    age = factor(age, levels = age_levels),
    tissue_maturity = factor(tissue_maturity, levels = tissue_order, labels = tissue_labels)
  )

################################
####### Plot AORe data now
####################################

age_levels <- c("Euteleostomi", "Actinopteri", "Teleostei", "Protacanthopterygii")

fill_colors_aore <- c(
  "Ancient" = "#00441b",  # very dark green
  "Old"     = "#1b7837",  # dark-medium green
  "Mid"     = "#5aae61",  # medium green
  "Young"   = "#a6dba0"   # soft but visible light green
)

names(fill_colors_aore) <- age_levels

###labels
label_df_aore <- df_plot %>%
  filter(redip == "AORe") %>%
  group_by(tissue_maturity, type) %>%
  summarise(label_pos = sum(avg_norm), .groups = "drop") %>%
  mutate(
    label = case_when(
      type == "2_to_2" ~ "S",
      type == "2_to_1" ~ "A",
      type == "1_to_1" ~ "E",
      TRUE ~ NA_character_
    ))

df_plot %>%
  filter(redip == "AORe") %>%
  mutate(type = factor(type, levels = c("1_to_1","2_to_1", "2_to_2"))) %>%
  ggplot(aes(y = type, x = avg_norm, fill = age)) +
  geom_bar(position = "stack", stat = "identity", color = "black", linewidth = 0.07, width = 0.85) +
  facet_wrap(vars(tissue_maturity), ncol = 1, strip.position = "left") +
  theme_classic(base_size = 9) +
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 9),
    axis.ticks.y = element_line(size = 0.2),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    strip.placement = "outside",
    strip.text.y.left = element_text(size = 10, angle = 0),
    strip.background = element_rect(fill = NA, color = "white", size = 0),
    panel.spacing = unit(0.1, "cm"),
    panel.grid = element_blank(),
    plot.margin = margin(t = 2, r = 5, b = 2, l = 5),
    legend.position = "none"
  ) +
  scale_fill_manual(values = fill_colors_aore) +
  scale_x_continuous(
    limits = c(0, 50),
    breaks = seq(0, 50, 10),
    minor_breaks = seq(0, 50, 5),
    expand = expansion(mult = c(0.05, 0.05)))

ggsave("F6_AORe_CNE_stacked_by_type_dodged_by_type_final.tiff", dpi = 600, width = 18, height = 24, units = "cm")

ggsave(filename = "F6_AORe_CNE_stacked_by_type_dodged_by_type_final.pdf",  plot = last_plot(), device = cairo_pdf, width = 18, height = 24, dpi = 600,units = "cm")


###########################################
#plot LORe data - figure 6
#############################################

age_levels <- c("Euteleostomi", "Actinopteri", "Teleostei", "Protacanthopterygii")

fill_colors_lore <- c(
  "Group1" = "#4B0082",  # Indigo (deep purple)
  "Group2" = "#6A0DAD",  # Royal purple
  "Group3" = "#9B30FF",  # Purple (vivid)
  "Group4" = "#CFA8E4"   # Soft violet (still visible)
)

names(fill_colors_lore) <- age_levels

###labels
label_df_lore <- df_plot %>%
  filter(redip == "LORe") %>%
  group_by(tissue_maturity, type) %>%
  summarise(label_pos = sum(avg_norm), .groups = "drop") %>%
  mutate(
    label = case_when(
      type == "2_to_2" ~ "S",
      type == "2_to_1" ~ "A",
      type == "1_to_1" ~ "E",
      TRUE ~ NA_character_
    ))

df_plot %>%
  filter(redip == "LORe") %>%
  mutate(type = factor(type, levels = c("1_to_1","2_to_1","2_to_2"))) %>%
  ggplot(aes(y = type, x = avg_norm, fill = age)) +  # note: y is now type
  geom_bar(position = "stack", stat = "identity", color = "black", linewidth = 0.07, width = 0.85) +
  facet_wrap(vars(tissue_maturity), ncol = 1, strip.position = "left") +
  theme_classic(base_size = 9) +
  theme(
    axis.text.y = element_blank(), 
    axis.text.x = element_text(size = 9),
    axis.ticks.y = element_line(size = 0.2),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    strip.placement = "outside",
    strip.text.y.left = element_text(size = 10, angle = 0),  # << this is the key line
    strip.background = element_rect(fill = NA, color = "white", size = 0),
    panel.spacing = unit(0.1, "cm"),
    panel.grid = element_blank(),
    plot.margin = margin(t = 2, r = 5, b = 2, l = 5),
    legend.position = "none"
  ) +
  scale_x_continuous(limits = c(0, 80),
                     breaks = seq(0, 80, 10),
                     minor_breaks = seq(0, 80, 5),
                     expand = expansion(mult = c(0.05, 0.05))) +
  scale_fill_manual(values = fill_colors_lore)

ggsave("F5_LORe_CNE_stacked_by_type_dodged_by_type_final.tiff", dpi = 600, width = 18, height = 24, units = "cm")

ggsave(filename = "F5_LORe_CNE_stacked_by_type_dodged_by_type_final.pdf",  plot = last_plot(), device = cairo_pdf, width = 18, height = 24, dpi = 600,units = "cm")

###################################################################
#### generate bar plots for numbers now
###################################################################

##AORe data

merged_donut_active_enh_redip <-bind_rows(cne_active_enh_counts_size_22_redip,
                                          cne_active_enh_counts_21_redip,
                                          cne_active_enh_counts_11_redip) %>% 
  group_by(redip,age) %>% 
  summarise(n_sum=sum(as.numeric(n)),n_merged=paste(n, collapse = ':')) 

merged_donut_active_enh_redip<-merged_donut_active_enh_redip %>% dplyr::select(age,n_merged) %>% separate(n_merged, into=c("shared","alignable","exclusive"), sep=":") %>% 
  pivot_longer(!c(age,redip),names_to = "category", values_to = "counts") %>% mutate(counts=as.numeric(counts))

merged_donut_active_enh_redip$age <- factor(merged_donut_active_enh_redip$age,levels=c('Euteleostomi','Actinopteri','Neopterygii','Teleostei',
                                                                                       'Clupeocephala','Euteleosteomorpha','Protacanthopterygii'))

merged_donut_active_enh_redip %>%
  filter(redip == "AORe") %>%
  mutate(
    category = factor(category, levels = c("exclusive", "alignable", "shared")),
    age = factor(age, levels = rev(c(
      'Euteleostomi','Actinopteri','Neopterygii','Teleostei',
      'Clupeocephala','Euteleosteomorpha','Protacanthopterygii')))
  ) %>%
  ggplot(aes(
    fill   = category,
    colour = category,
    x      = counts,
    y      = age
  )) +
  geom_bar(position = "stack", stat = "identity", size = 0.1) +
  scale_y_discrete(
    labels = rev(c('Et','Ac','Ne','Te','Cl','Eu','Pr')),
    position = "left"
  ) +
  scale_x_continuous(
    limits = c(0,6000),
    breaks = seq(0,6000,3000)
  ) +
  # match color order to factor levels
  scale_fill_manual(values = c("exclusive" = "white",
                               "alignable" = "grey50",
                               "shared"    = "black")) +
  scale_colour_manual(values = c("exclusive" = "grey",
                                 "alignable" = "grey50",
                                 "shared"    = "black")) +
  labs(x = "Number of CNEs", y = "CNE age") +
  theme_classic() +
  theme(
    axis.text  = element_text(colour = "black", size = 8),
    axis.ticks = element_line(colour = "black", size = 0.25),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.line    = element_line(size = 0.25),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    legend.position = "none"
  )

ggsave("CNE_props_active_enhancers_AORE_stacked_bar_2025_v2.tiff", dpi = 600, width = 6, height = 5, units = "cm")

ggsave(filename = "CNE_props_active_enhancers_AORE_stacked_bar_2025_v2.pdf",  plot = last_plot(), device = cairo_pdf, width = 6, height = 5, dpi = 600,units = "cm")


# LORe data now

merged_donut_active_enh_redip %>%
  filter(redip == "LORe") %>%
  mutate(
    category = factor(category, levels = c("exclusive", "alignable", "shared")),
    age = factor(age, levels = rev(c(
      'Euteleostomi','Actinopteri','Neopterygii','Teleostei',
      'Clupeocephala','Euteleosteomorpha','Protacanthopterygii')))
  ) %>%
  ggplot(aes(
    fill   = category,
    colour = category,
    x      = counts,
    y      = age
  )) +
  geom_bar(position = "stack", stat = "identity", size = 0.1) +
  scale_y_discrete(
    labels = rev(c('Et','Ac','Ne','Te','Cl','Eu','Pr')),
    position = "left"
  ) +
  scale_x_continuous(
    limits = c(0,3000),
    breaks = seq(0,3000,1500)
  ) +
  # match color order to factor levels
  scale_fill_manual(values = c("exclusive" = "white",
                               "alignable" = "grey50",
                               "shared"    = "black")) +
  scale_colour_manual(values = c("exclusive" = "grey",
                                 "alignable" = "grey50",
                                 "shared"    = "black")) +
  labs(x = "Number of CNEs", y = "CNE age") +
  theme_classic() +
  theme(
    axis.text  = element_text(colour = "black", size = 8),
    axis.ticks = element_line(colour = "black", size = 0.25),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.line    = element_line(size = 0.25),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    legend.position = "none"
  )

ggsave("CNE_props_active_enhancers_LORE_stacked_bar_2025_v2.tiff", dpi = 600, width = 6, height = 5, units = "cm")

ggsave(filename = "CNE_props_active_enhancers_LORE_stacked_bar_2025_v2.pdf",  plot = last_plot(), device = cairo_pdf, width = 6, height = 5, dpi = 600,units = "cm") 






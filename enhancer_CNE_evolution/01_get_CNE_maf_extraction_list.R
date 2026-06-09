############################################################
## CNE–Enhancer list extraction
#######################################################

library(tidyverse)
library(fuzzyjoin)
library(openxlsx)
library(tidygenomics)
library(janitor)
library(stringr)

############################################################
## 1. Load Pike CNE coordinates
############################################################

pike_cne_file <- "./CNEs_pike_2023.03.23_coordinates.txt"

pike_cne_coords <- read.table(
  pike_cne_file,
  header = FALSE,
  col.names = c("chromosome", "start", "end", "cne_id", "age")
)

# Quick sanity checks
pike_cne_coords %>% distinct(cne_id)


############################################################
## 2. Summarise CNE landscape
############################################################

pike_cne_counts <- pike_cne_coords %>%
  count(age, name = "n_cnes")

pike_cne_length_by_age <- pike_cne_coords %>%
  mutate(length = end - start) %>%
  group_by(age) %>%
  summarise(total_length = sum(length), .groups = "drop")

pike_cne_total_length <- pike_cne_coords %>%
  mutate(length = end - start) %>%
  summarise(total_length = sum(length))

pike_cne_summary <- pike_cne_counts %>%
  left_join(pike_cne_length_by_age, by = "age")

###########################################################
## 3. load pike and atac data # outputs from the step 3 of the promoter enhancer conservation step - 03_overlay_maf_with_atac.R
############################################################

maf_files <- list.files(pattern = "_pike_sequences.xlsx$", recursive = TRUE)

maf_list <- lapply(maf_files, read.xlsx)

# Combine all enhancer/MAF datasets
all_maf_pike <- bind_rows(maf_list)

############################################################
## 4. Join CNEs with enhancer/Maf coordinates
############################################################

all_dat_maf_pike <- all_maf_pike %>%
  genome_left_join(
    pike_cne_coords,
    by = c("chromosome", "start", "end")
  ) %>%
  mutate(
    overlap = pmin(end.x, end.y) - pmax(start.x, start.x),
    cne_size = end.y - start.y
  )

# Create a unique region identifier for downstream tracking
all_dat_maf_pike <- all_dat_maf_pike %>%
  unite(
    col = region,
    chromosome.y, start.y, end.y,
    sep = "_",
    remove = FALSE
  )

rm(maf_list)

############################################################
## 5. Load Salmon LORe / AORE regions
############################################################

salmon_lore <- read.table(
  "./Ssal_late_rediploidized_regions.tsv",
  header = TRUE
) %>%
  rename(chromosome = chrom) %>%
  mutate(
    chromosome = sprintf("ssa%02d", as.numeric(chromosome)),
    start = as.numeric(start),
    end = as.numeric(end)
  )

############################################################
## 6. Load enhancer conservation category lists
############################################################

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

############################################################
## 7. Map enhancers to CNE annotations
############################################################

shared_joined <- shared_enhancers %>%
  left_join(all_dat_maf_pike, by = "Origin") %>%
  mutate(age = as.character(age))

exclusive_joined <- exclusive_enhancers %>%
  left_join(all_dat_maf_pike, by = "Origin") %>%
  mutate(age = as.character(age))

alignable_joined <- alignable_enhancers %>%
  left_join(all_dat_maf_pike, by = "Origin") %>%
  mutate(age = as.character(age))

############################################################
## 8. Extract perfectly overlapping CNEs (consistency check)
############################################################

extract_full_overlap_cnes <- function(df, output_file) {
  
  out <- df %>%
    mutate(age = replace_na(age, "NO_CNE_OVERLAP")) %>%
    filter(age != "NO_CNE_OVERLAP") %>%
    mutate(overlap_rate = overlap / cne_size) %>%
    filter(overlap_rate == 1) %>%
    select(Origin, region, chromosome.y, start.y, end.y) %>%
    mutate(
      ssal_chr = str_extract(Origin, "ssa\\d+"),
      ssal_chr = str_replace(
        ssal_chr,
        "ssa(\\d+)",
        function(x) {
          num <- str_match(x, "ssa(\\d+)")[, 2]
          paste0("ssa", str_pad(num, 2, pad = "0"))
        }
      )
    ) %>%
    select(ssal_chr, chromosome.y, start.y, end.y, Origin) %>%
    distinct(chromosome.y, start.y, end.y, .keep_all = TRUE)
  
  write.table(
    out,
    file = output_file,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  
  return(out)
}

############################################################
## 9. Run extraction for each enhancer class
############################################################

shared_cne_extraction_list <- extract_full_overlap_cnes(
  shared_joined,
  "shared_CNE_extraction_list.tsv"
)

exclusive_cne_extraction_list <- extract_full_overlap_cnes(
  exclusive_joined,
  "exclusive_CNE_extraction_list.tsv"
)

alignable_cne_extraction_list <- extract_full_overlap_cnes(
  alignable_joined,
  "alignable_CNE_extraction_list.tsv"
)
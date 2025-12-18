library(tidyverse)

# Make sample design csv files for Atlantic salmon Devmap samples, to run nf-core pipelines

results <- "/path/to/project/Aqua-Faang/nfcore"

# RNA samples
path <- paste(results, "/AtlanticSalmon/DevMap/RNA")
read_csv(paste0(path, "/seq_data/samplesheet/samplesheet.csv"))%>%
  transmute(
    sample = paste0(
      "AtlanticSalmon_RNA_", 
      gsub("[ -]", "", str_to_title(sub(".*, ", "", sample_description))), 
      "_R", 
      str_sub(library_name, 8, 8)
    ),
    fastq_1 = paste0(path, fastq_1), 
    fastq_2 = paste0(path, fastq_2), 
    strandedness = "reverse") %>%
  write_csv(paste0(path, "design.csv"))

# ATAC samples
path <- paste(results, "/AtlanticSalmon/DevMap/ATAC")
read_csv(paste0(path, "/seq_data/samplesheet/samplesheet.csv")) %>%
  transmute(
    group = paste0(
      "AtlanticSalmon_ATAC_", 
      gsub("[ -,]", "", str_to_title(sub(".*, ", "", sample_description)))
    ),
    replicate = str_sub(library_name, 8, 8),
    fastq_1 = paste0(path, fastq_1), 
    fastq_2 = paste0(path, fastq_2)
  ) %>%
  write_csv(paste0(path, "design.csv"))

# ChIP
path <- paste(results, "/AtlanticSalmon/DevMap/ChIP")
read_csv(paste0(path, "/seq_data/samplesheet/samplesheet.csv")) %>%
  transmute(
    group = paste0(
      "AtlanticSalmon_ChIP-", sub(".*_", "", experiment_alias), "_",
      gsub("[ -,]", "", str_to_title(sub(".*, ", "", sample_description)))
    ),
    replicate = sub(".*Rep", "", run_alias),
    fastq_1 = paste0(path, fastq_1), 
    fastq_2 = paste0(path, fastq_2),
    antibody = sub(".*_", "", experiment_alias),
  ) %>%
  mutate(control = sub("ChIP-[^_]*", "ChIP-input", group)) %>%
  mutate(antibody = ifelse(grepl("input", group), "", antibody),
         control = ifelse(grepl("input", group), "", control)) %>%
  write_csv(paste0(path, "design.csv"))

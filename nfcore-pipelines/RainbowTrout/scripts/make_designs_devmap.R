library(tidyverse)

# Make sample design csv files for Atlantic salmon Devmap samples, to run nf-core pipelines

results <- "/mnt/project/Aqua-Faang/nfcore"

# RNA samples
path <- paste(results, "/RainbowTrout/DevMap/RNA")
read_csv(paste0(path, "seq_data/samplesheet/samplesheet.csv")) %>%
  transmute(
    sample = paste0(
      "RainbowTrout_RNA_", 
      gsub("[ -]", "", str_to_title(sub(".*, ", "", sample_description))), 
      "_R", 
      str_sub(library_name, 8, 8)
    ),
    fastq_1 = paste0(path, fastq_1), 
    fastq_2 = paste0(path, fastq_2), 
    strandedness = "reverse") %>%
  write_csv(paste0(path, "design.csv"))

# ATAC samples
path <- paste(results, "/RainbowTrout/DevMap/ATAC")
read_csv(paste0(path, "seq_data/samplesheet/samplesheet.csv")) %>%
  transmute(
    group = paste0(
      "RainbowTrout_ATAC_", 
      gsub("[ -,]", "", str_to_title(sub(".*, ", "", sample_description)))
    ),
    replicate = sub(".*_", "", experiment_alias),
    fastq_1 = paste0(path, fastq_1),
    fastq_2 = paste0(path, fastq_2),
    experiment_alias,
  ) %>%
  arrange(experiment_alias) %>%
  write_csv(paste0(path, "design.csv"))

# ChIP samples
path <- paste(results, "/RainbowTrout/DevMap/ChIP")
read_csv(paste0(path, "seq_data/samplesheet/samplesheet.csv")) %>%
  transmute(
    sample = paste0(
      "RainbowTrout_ChIP-",
      str_to_title(sub(".*_", "", experiment_alias)), "_",
      gsub("[ -]", "", str_to_title(sub(".*, ", "", sample_description)))
    ),
    replicate = sub(".*Rep", "", run_alias),
    fastq_1 = paste0(path, fastq_1), 
    fastq_2 = paste0(path, fastq_2), 
    antibody = ifelse(grepl("Input", sample), "", str_to_title(sub(".*_", "", experiment_alias))),
    control = ifelse(grepl("Input", sample), "", paste0("RainbowTrout_ChIP-Input_", gsub("[ -]", "", str_to_title(sub(".*, ", "", sample_description)))))) %>%
  write_csv(paste0(path, "design.csv"))

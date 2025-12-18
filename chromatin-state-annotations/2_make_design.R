library(tidyverse)

# Results from nfcore pipelines
nfcore <- "/mnt/project/Aqua-Faang/nfcore"

# Atlantic salmon

# BodyMap

# ATAC samples
ATAC_samples <- do.call(
	bind_rows, 
	lapply(list.files(paste0(nfcore, "/AtlanticSalmon/BodyMap/ATAC"), "design.csv", recursive = TRUE, full.names = TRUE), 
				 function(x) mutate(read_csv(x), design = x))
  ) %>%
	mutate(results = sub("design.csv", "results/bwa/mergedLibrary", design)) %>%
	unite("file", results, group, sep = "/", remove = FALSE) %>%
	unite("file", file, replicate, sep = "_R", remove = FALSE) %>%
	transmute(
		cell = sub("AtlanticSalmon_ATAC_", "", group), 
		mark = "ATAC",
		file = paste0(file, ".mLb.clN.sorted.bam"),
	) %>%
  filter(!grepl("Gill", cell)) %>% # Exclude these tissues
	distinct()

# ChIP samples, exluding certain marks and tissue
ChIP_samples <- do.call(
	bind_rows, 
	lapply(list.files(paste0(nfcore, "/BodyMap/ChIP"), "design.csv", recursive = TRUE, full.names = TRUE), 
				 function(x) mutate(read_csv(x), design = x))
	) %>%
	mutate(results = sub("design.csv", "results/bwa/mergedLibrary", design)) %>%
	unite("file", results, group, sep = "/", remove = FALSE) %>%
	unite("file", file, replicate, sep = "_R", remove = FALSE) %>%
	unite("control", results, control, sep = "/", remove = FALSE) %>%
	transmute(
		cell = sub("AtlanticSalmon_ChIP-[^_]*_", "", group), 
		mark = antibody,
		file = paste0(file, ".mLb.clN.sorted.bam"),
		control = paste0(control, "_R1.mLb.clN.sorted.bam")
	) %>%
	filter(!is.na(mark)) %>%
	filter(!grepl("Sperm|Gill", cell)) %>% # Exclude these tissues
	filter(!mark %in% c("CTCF", "H3K36me3", "DMC1", "Input")) %>% # Exclude these marks
	distinct()

# Path to chromHMM results
results = "/mnt/project/Aqua-Faang/chromatin_states/AtlanticSalmon/BodyMap"

# Total design
total_design <- bind_rows(ATAC_samples, ChIP_samples) %>%
  arrange(cell)
total_design %>%
  mutate(file = basename(file), control = basename(control)) %>%
  write_tsv(paste0(results, "/design.tsv"), na = "", col_names = FALSE)

# Create links to bam files
dir.create(paste0(results, "/bams"), recursive = TRUE)
for (x in unique(c(total_design$file, na.omit(total_design$control)))) {
  system(paste("ln -s", x, paste0(path, "/bams/", basename(x))))
}


# DevMap

# ATAC samples
ATAC_samples <- do.call(
  bind_rows, 
  lapply(list.files(paste0(nfcore, "/AtlanticSalmon/DevMap/ATAC"), "design.csv", recursive = TRUE, full.names = TRUE), 
         function(x) mutate(read_csv(x), design = x))
) %>%
  mutate(results = sub("design.csv", "results/bwa/mergedLibrary", design)) %>%
  unite("file", results, group, sep = "/", remove = FALSE) %>%
  unite("file", file, replicate, sep = "_R", remove = FALSE) %>%
  transmute(
    cell = sub("AtlanticSalmon_ATAC_", "", group), 
    mark = "ATAC",
    file = paste0(file, ".mLb.clN.sorted.bam"),
  ) %>%
  filter(!grepl("LateEyed", cell)) %>% # Exclude these samples
  distinct()

# ChIP samples
ChIP_samples <- do.call(
  bind_rows, 
  lapply(list.files(paste0(nfcore, "/AtlanticSalmon/DevMap/ChIP"), "design.csv", recursive = TRUE, full.names = TRUE), 
         function(x) mutate(read_csv(x), design = x))
) %>%
  mutate(results = sub("design.csv", "results/bwa/mergedLibrary", design)) %>%
  unite("file", results, group, sep = "/", remove = FALSE) %>%
  unite("file", file, replicate, sep = "_R", remove = FALSE) %>%
  unite("control", results, control, sep = "/", remove = FALSE) %>%
  transmute(
    cell = sub("AtlanticSalmon_ChIP-[^_]*_", "", group), 
    mark = antibody,
    file = paste0(file, ".mLb.clN.sorted.bam"),
    control = paste0(control, "_R1.mLb.clN.sorted.bam")
  ) %>%
  filter(!is.na(mark)) %>%
  filter(!grepl("LateEyed", cell)) %>% # Exclude these samples
  filter(!mark %in% c("Input")) %>%
  distinct()

# Path to chromHMM results
results = "/mnt/project/Aqua-Faang/chromatin_states/AtlanticSalmon/DevMap"

# Total design
total_design <- bind_rows(ATAC_samples, ChIP_samples) %>%
  arrange(cell)
total_design %>%
  mutate(file = basename(file), control = basename(control)) %>%
  write_tsv(paste0(results, "/design.tsv"), na = "", col_names = FALSE)

# Create links to bam files
dir.create(paste0(results, "/bams"), recursive = TRUE)
for (x in unique(c(total_design$file, na.omit(total_design$control)))) {
	system(paste("ln -s", x, paste0(results, "/bams/", basename(x))))
}

# Rainbow Trout

# Bodymap data

# ATAC samples
ATAC_samples <- do.call(
  bind_rows, 
  lapply(list.files(paste0(nfcore, "/RainbowTrout/BodyMap/ATAC"), "design.csv", recursive = TRUE, full.names = TRUE), 
         function(x) mutate(read_csv(x), design = x))
) %>%
  mutate(results = sub("design.csv", "results/bwa/mergedLibrary", design)) %>%
  unite("file", results, group, sep = "/", remove = FALSE) %>%
  unite("file", file, replicate, sep = "_R", remove = FALSE) %>%
  transmute(
    cell = sub("RainbowTrout_ATAC_", "", group), 
    mark = "ATAC",
    file = paste0(file, ".mLb.clN.sorted.bam"),
  ) %>%
  filter(!grepl("Gill|Gonad", cell)) %>% # Exclude these tissues
  distinct()

# ChIP samples
ChIP_samples <- do.call(
  bind_rows, 
  lapply(list.files(paste0(nfcore, "/RainbowTrout/BodyMap/ChIP"), "design.csv", recursive = TRUE, full.names = TRUE), 
         function(x) mutate(read_csv(x), design = x))
) %>%
  mutate(results = sub("design.csv", "results/bwa/mergedLibrary", design)) %>%
  unite("file", results, group, sep = "/", remove = FALSE) %>%
  unite("file", file, replicate, sep = "_R", remove = FALSE) %>%
  unite("control", results, control, sep = "/", remove = FALSE) %>%
  transmute(
    cell = sub("RainbowTrout_ChIP-[^_]*_", "", group), 
    mark = antibody,
    file = paste0(file, ".mLb.clN.sorted.bam"),
    control = paste0(control, "_R1.mLb.clN.sorted.bam")
  ) %>%
  filter(!is.na(mark)) %>%
  filter(!grepl("Gill|Gonad", cell)) %>% # Exclude these tissues
  filter(!mark %in% c("CTCF", "H3K36me3", "DMC1", "Input")) %>% # Exclude these marks
  distinct()

# Path to chromHMM results
results = "/mnt/project/Aqua-Faang/chromatin_states/RainbowTrout/BodyMap"

# Total design
total_design <- bind_rows(ATAC_samples, ChIP_samples) %>%
  arrange(cell)
total_design %>%
  mutate(file = basename(file), control = basename(control)) %>%
  write_tsv(paste0(results, "/design.tsv"), na = "", col_names = FALSE)

# Create links to bam files
dir.create(paste0(results, "/bams"), recursive = TRUE)
for (x in unique(c(total_design$file, na.omit(total_design$control)))) {
  system(paste("ln -s", x, paste0(path, "/bams/", basename(x))))
}

# DevMap data

# ATAC samples
ATAC_samples <- do.call(
  bind_rows, 
  lapply(list.files(paste0(nfcore, "/RainbowTrout/DevMap/ATAC"), "design.csv", recursive = TRUE, full.names = TRUE), 
         function(x) mutate(read_csv(x), design = x))
) %>%
  mutate(results = sub("design.csv", "results/bwa/mergedLibrary", design)) %>%
  unite("file", results, group, sep = "/", remove = FALSE) %>%
  unite("file", file, replicate, sep = "_R", remove = FALSE) %>%
  transmute(
    cell = sub("RainbowTrout_ATAC_", "", group), 
    mark = "ATAC",
    file = paste0(file, ".mLb.clN.sorted.bam"),
  ) %>%
  filter(!grepl("LateEyed", cell)) %>% # Exclude these samples
  distinct()

# ChIP samples
ChIP_samples <- do.call(
  bind_rows, 
  lapply(list.files(paste0(AF_path, "/RainbowTrout/DevMap/ChIP"), "design.csv", recursive = TRUE, full.names = TRUE), 
         function(x) mutate(read_csv(x), design = x))
) %>%
  mutate(results = sub("design.csv", "results/bwa/mergedLibrary", design)) %>%
  unite("file", results, group, sep = "/", remove = FALSE) %>%
  unite("file", file, replicate, sep = "_R", remove = FALSE) %>%
  unite("control", results, control, sep = "/", remove = FALSE) %>%
  transmute(
    cell = sub("RainbowTrout_ChIP-[^_]*_", "", group), 
    mark = antibody,
    file = paste0(file, ".mLb.clN.sorted.bam"),
    control = paste0(control, "_R1.mLb.clN.sorted.bam")
  ) %>%
  filter(!is.na(mark)) %>%
  filter(!grepl("LateEyed", cell)) %>% # Exclude these samples
  filter(!mark %in% c("Input")) %>%
  distinct()

# Path to chromHMM results
results = "/mnt/project/Aqua-Faang/chromatin_states/RainbowTrout/DevMap"

# Total design
total_design <- bind_rows(ATAC_samples, ChIP_samples) %>%
  arrange(cell)
total_design %>%
  mutate(file = basename(file), control = basename(control)) %>%
  write_tsv(paste0(results, "/design.tsv"), na = "", col_names = FALSE)

# Create links to bam files
dir.create(paste0(results, "/bams"), recursive = TRUE)
for (x in unique(c(total_design$file, na.omit(total_design$control)))) {
  system(paste("ln -s", x, paste0(results, "/bams/", basename(x))))
}

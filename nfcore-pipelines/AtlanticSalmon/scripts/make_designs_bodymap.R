library(tidyverse)

# Make sample design csv files for Atlantic salmon BodyMap samples, to run nf-core pipelines

# Decode function for Aqua-Faang sample id codes
decodeSampleID <- function (id) {
	
	partnerCode <- c(
		"A" = "HCMR",
		"B" = "INRAe",
		"C" = "NMBU",
		"D" = "UEDIN",
		"E" = "UNIABDN",
		"F" = "UniPD",
		"G" = "USC",
		"H" = "WU",
		"I" = "ZIGR",
		"J" = "UoB"
	)
	
	speciesCode <- c(
		"1" = "AtlanticSalmon",
		"2" = "Carp",
		"3" = "RainbowTrout",
		"4" = "SeaBass",
		"5" = "SeaBream",
		"6" = "Turbot"
	)
	
	activityCode <- c(
		"1" = "Immature_Female",
		"2" = "Immature_Male",
		"3" = "Mature_Female",
		"4" = "Mature_Male"
	)
	
	libraryTypeCode <- c(
		"1" = "ATAC",
		"2" = "ChIP",
		"3" = "RNA",
		"4" = "WGS",
		"5" = "sRNA"
	)
	
	antibodyCode <- c(
		"1" = "CTCF",
		"2" = "H3K27ac",
		"3" = "H3K27me3",
		"4" = "H3K4me1",
		"5" = "H3K4me3",
		"6" = "Input",
		"7" = "H3K36me3",
		"8" = "DMC1",
		"0" = NA
	)
	
	tissueCode <- c(
		"1" = "Brain",
		"2" = "Gill",
		"3" = "Gonad",
		"4" = "HeadKidney",
		"5" = "Liver",
		"6" = "Muscle",
		"7" = "SuppDistalIntestine",
		"8" = "Sperm"
	)
	
	partner <- partnerCode[str_sub(id, 1, 1)]
	species <- speciesCode[str_sub(id, 2, 2)]
	library <- libraryTypeCode[str_sub(id, 3, 3)]
	activity <- activityCode[str_sub(id, 4, 4)]
	fishNumber <- str_sub(id, 5, 6)
	antibody <- antibodyCode[str_sub(id, 7, 7)]
	tissue <- tissueCode[str_sub(id, 8, 8)]
	
	fullName <- paste0(
		partner, "_",
		species, "_",
		library, "_",
		activity, "_",
		fishNumber, "_",
		ifelse(!is.na(antibody), paste0(antibody, "_"), ""),
		tissue
	)
	group <- paste0(
		species, "_",
		ifelse(is.na(antibody), paste0(library, "_"), paste0(library, "-", antibody, "_")),
		tissue, "_",
		activity
		
	)
	return (group)
} 

data <- "/mnt/labfiles/data/Aqua-Faang/"
results <- "/mnt/project/Aqua-Faang/nfcore"

# ChIP samples
inPath = paste0(data, "/AtlanticSalmon/ChIP")
outPath = paste0(results, "/AtlanticSalmon/BodyMap/ChIP")
if (!dir.exists(outPath)) {
	dir.create(outPath)
}
design <- tibble(sample = unique(sub("_[12]\\.fq\\.gz$", "", list.files(inPath, ".*fq.gz")))) %>%
	mutate(
		group = decodeSampleID(sample),
		replicate = 1,
		fastq_1 = paste0(inPath, "/", sample, "_1.fq.gz"),
		fastq_2 = paste0(inPath, "/", sample, "_2.fq.gz")
	) %>%
	group_by(group) %>%
	mutate(replicate = 1:n()) %>%
	ungroup() %>%
	mutate(antibody = sapply(group, function(x) sub("_.*", "", sub(".*_ChIP-", "", x)))) %>%
	mutate(control = sapply(group, function(x) sub(sub("_.*", "", sub(".*_ChIP-", "", x)), "Input", x))) %>%
	mutate(antibody = ifelse(grepl("Input", group), "", antibody),
				 control = ifelse(grepl("Input", group), "", control)) %>%
	select(-sample) %>%
  arrange(group)
for (tissue in unique(sapply(design$group, function(x) unlist(str_split(x, "_"))[3]))) { # Split by tissue
	if (!dir.exists(paste0(outPath, "/", tissue))) {
		dir.create(paste0(outPath, "/", tissue))
	}
	design %>%
		filter(grepl(tissue, group)) %>%
		write_csv(paste0(outPath, "/", tissue, "/design.csv"))
}

# ATAC samples
inPath = paste0(data, "/AtlanticSalmon/ATAC")
outPath = paste0(results, "/AtlanticSalmon/BodyMap/ATAC")
if (!dir.exists(outPath)) {
	dir.create(outPath)
}
design <- tibble(sample = unique(sub("_[12]\\.fq\\.gz$", "", list.files(inPath, ".*fq.gz")))) %>%
	mutate(
		group = decodeSampleID(sample),
		replicate = 1,
		fastq_1 = paste0(inPath, "/", sample, "_1.fq.gz"),
		fastq_2 = paste0(inPath, "/", sample, "_2.fq.gz")
	) %>%
	group_by(group) %>%
	mutate(replicate = 1:n()) %>%
	ungroup() %>%
	select(-sample) %>%
  arrange(group)
for (tissue in unique(sapply(design$group, function(x) unlist(str_split(x, "_"))[3]))) { # Split by tissue
	if (!dir.exists(paste0(outPath, "/", tissue))) {
		dir.create(paste0(outPath, "/", tissue))
	}
	design %>%
		filter(grepl(tissue, group)) %>%
		write_csv(paste0(outPath, "/", tissue, "/design.csv"))
}

# RNA samples
inPath = paste0(data, "/AtlanticSalmon/RNA")
outPath = paste0(results, "/AtlanticSalmon/BodyMap/RNA")
if (!dir.exists(outPath)) {
	dir.create(outPath)
}
design <- tibble(sample = unique(sub("_[12]\\.fq\\.gz$", "", list.files(inPath, ".*fq.gz")))) %>%
	mutate(
		group = decodeSampleID(sample),
		replicate = 1,
		fastq_1 = paste0(inPath, "/", sample, "_1.fq.gz"),
		fastq_2 = paste0(inPath, "/", sample, "_2.fq.gz"),
		strandedness = "reverse"
	) %>%
	group_by(group) %>%
	mutate(replicate = 1:n()) %>%
	ungroup() %>%
	transmute(sample = paste0(group, "_R", replicate), fastq_1, fastq_2, strandedness) %>%
  arrange(sample)
for (tissue in unique(sapply(design$sample, function(x) unlist(str_split(x, "_"))[3]))) { # Split by tissue
	if (!dir.exists(paste0(outPath, "/", tissue))) {
		dir.create(paste0(outPath, "/", tissue))
	}
	design %>%
		filter(grepl(tissue, sample)) %>%
		write_csv(paste0(outPath, "/", tissue, "/design.csv"))
}

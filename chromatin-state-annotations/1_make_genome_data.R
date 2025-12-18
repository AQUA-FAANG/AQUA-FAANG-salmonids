library(tidyverse)

# Make Atlantic salmon and rainbowtrout genome size, coords, and anchor files for chromHMM

# Path to chromHMM
chromhmm <- "/mnt/project/Aqua-Faang/chromatin_states/scripts/ChromHMM"

# Atlantic salmon

# Path to genome files
genome <- "/mnt/project/Aqua-Faang/nfcore/AtlanticSalmon/genome"

# Assembly name
assembly <- "Ssal_v3.1"

# Chromosome sizes
chrs <- read_tsv(paste0(genome, "/Salmo_salar-GCA_905237065.2-softmasked.fa.sizes"), col_names = c("seqname", "length")) %>%
  filter(!grepl("CAJNNT", seqname))
chrs %>%
  write_tsv(paste0(chromhmm, "/CHROMSIZES/Ssal_v3.1.txt"), col_names = FALSE)

# Read in gtf
GTF <- read_tsv(paste0(genome, "/Salmo_salar-GCA_905237065.2-softmasked_genes.gtf"),
								col_names = c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"), 
								comment = "#") %>%
	right_join(chrs, by = "seqname")

# Anchor and coord directories 
dir.create(paste0(chromhmm, "/ANCHORFILES/", assembly))
dir.create(paste0(chromhmm, "/COORDS/", assembly))

# TSS
GTF %>%
	filter(feature == "gene") %>%
	transmute(
		seqname, 
		start = ifelse(strand == "+", start, end - 1), 
		end = ifelse(strand == "+", start + 1, end), 
		strand
	) %>%
	write_tsv(paste0(chromhmm, "/COORDS/", assembly, "/TSS.tsv"), col_names = FALSE)
GTF %>%
	filter(feature == "gene") %>%
	transmute(
		seqname, 
		start = ifelse(strand == "+", start, end), 
		strand
	) %>%
	write_tsv(paste0(chromhmm, "/ANCHORFILES/", assembly, "/TSS.tsv"), col_names = FALSE)

# TES
GTF %>%
	filter(feature == "gene") %>%
	transmute(
		seqname, 
		start = ifelse(strand == "+", end, start - 1), 
		end = ifelse(strand == "+", end + 1, start + 1), 
		strand	
	) %>%
	write_tsv(paste0(chromhmm, "/COORDS/", assembly, "/TES.tsv"), col_names = FALSE)
GTF %>%
	filter(feature == "gene") %>%
	transmute(
		seqname, 
		start = ifelse(strand == "+", end, start), 
		strand	
	) %>%
	write_tsv(paste0(chromhmm, "/ANCHORFILES/", assembly, "/TES.tsv"), col_names = FALSE)

# Gene
GTF %>%
	filter(feature == "gene") %>%
	transmute(seqname, start, end, strand) %>%
	write_tsv(paste0(chromhmm, "/COORDS/", assembly, "/gene.tsv"), col_names = FALSE)

# Exon
GTF %>%
	filter(feature == "exon") %>%
	transmute(seqname, start, end, strand) %>%
	write_tsv(paste0(chromhmm, "/COORDS/", assembly, "/exon.tsv"), col_names = FALSE)

# 2kb upstream of TSS
GTF %>%
	filter(feature == "gene") %>%
	transmute(
		seqname, 
		start = ifelse(strand == "+", start - 2000, end), 
		end = ifelse(strand == "+", start + 2000, end + 2000),
		strand,
		length) %>%
	mutate(start = ifelse(start < 0, 0, start)) %>% 
	mutate(end = ifelse(end > length, length, end)) %>%
	filter((end - start) > 0) %>%
	select(-length) %>%
	write_tsv(paste0(chromhmm, "/COORDS/", assembly, "/TSS2kb.tsv"), col_names = FALSE)

# 2kb downstream of TES
GTF %>%
	filter(feature == "gene") %>%
	transmute(
		seqname, 
		start = ifelse(strand == "+", end, start - 2000),
		end = ifelse(strand == "+", end + 2000, start + 2000),
		strand,
		length) %>%
	mutate(start = ifelse(start < 0, 0, start)) %>% 
	mutate(end = ifelse(end > length, length, end)) %>%
	filter((end - start) > 0) %>%
	select(-length) %>%
	write_tsv(paste0(chromhmm, "/COORDS/", assembly, "/TES2kb.tsv"), col_names = FALSE)

# Rainbow trout

# Path to genome files
genome <- "/mnt/project/Aqua-Faang/seq_results/RainbowTrout/genome/"

# Assembly name
assembly <- "Omyk"

# Chromosome sizes
chrs <- read_tsv(paste0(genome, "/Oncorhynchus_mykiss-GCA_013265735.3-softmasked.fa.sizes"), col_names = c("seqname", "length")) %>%
  filter(!grepl("JAAXML", seqname))
chrs %>%
  write_tsv(paste0(chromhmm, "/CHROMSIZES/Omyk.txt"), col_names = FALSE)

# Read in gtf
GTF <- read_tsv(paste0(genome, "/Oncorhynchus_mykiss-GCA_013265735.3-2020_12-genes.gtf"),
                col_names = c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"), 
                comment = "#") %>%
  right_join(chrs, by = "seqname")

# Anchor and coord directories 
dir.create(paste0(chromhmm, "/ANCHORFILES/", assembly))
dir.create(paste0(chromhmm, "/COORDS/", assembly))

# TSS
GTF %>%
  filter(feature == "gene") %>%
  transmute(
    seqname, 
    start = ifelse(strand == "+", start, end - 1), 
    end = ifelse(strand == "+", start + 1, end), 
    strand
  ) %>%
  write_tsv(paste0(chromhmm, "/COORDS/", assembly, "/TSS.tsv"), col_names = FALSE)
GTF %>%
  filter(feature == "gene") %>%
  transmute(
    seqname, 
    start = ifelse(strand == "+", start, end), 
    strand
  ) %>%
  write_tsv(paste0(chromhmm, "/ANCHORFILES/", assembly, "/TSS.tsv"), col_names = FALSE)

# TES
GTF %>%
  filter(feature == "gene") %>%
  transmute(
    seqname, 
    start = ifelse(strand == "+", end, start - 1), 
    end = ifelse(strand == "+", end + 1, start + 1), 
    strand	
  ) %>%
  write_tsv(paste0(chromhmm, "/COORDS/", assembly, "/TES.tsv"), col_names = FALSE)
GTF %>%
  filter(feature == "gene") %>%
  transmute(
    seqname, 
    start = ifelse(strand == "+", end, start), 
    strand	
  ) %>%
  write_tsv(paste0(chromhmm, "/ANCHORFILES/", assembly, "/TES.tsv"), col_names = FALSE)

# Gene
GTF %>%
  filter(feature == "gene") %>%
  transmute(seqname, start, end, strand) %>%
  write_tsv(paste0(chromhmm, "/COORDS/", assembly, "/gene.tsv"), col_names = FALSE)

# Exon
GTF %>%
  filter(feature == "exon") %>%
  transmute(seqname, start, end, strand) %>%
  write_tsv(paste0(chromhmm, "/COORDS/", assembly, "/exon.tsv"), col_names = FALSE)

# 2kb upstream of TSS
GTF %>%
  filter(feature == "gene") %>%
  transmute(
    seqname, 
    start = ifelse(strand == "+", start - 2000, end), 
    end = ifelse(strand == "+", start + 2000, end + 2000),
    strand,
    length) %>%
  mutate(start = ifelse(start < 0, 0, start)) %>% 
  mutate(end = ifelse(end > length, length, end)) %>%
  filter((end - start) > 0) %>%
  select(-length) %>%
  write_tsv(paste0(chromhmm, "/COORDS/", assembly, "/TSS2kb.tsv"), col_names = FALSE)

# 2kb downstream of TES
GTF %>%
  filter(feature == "gene") %>%
  transmute(
    seqname, 
    start = ifelse(strand == "+", end, start - 2000),
    end = ifelse(strand == "+", end + 2000, start + 2000),
    strand,
    length) %>%
  mutate(start = ifelse(start < 0, 0, start)) %>% 
  mutate(end = ifelse(end > length, length, end)) %>%
  filter((end - start) > 0) %>%
  select(-length) %>%
  write_tsv(paste0(chromhmm, "/COORDS/", assembly, "/TES2kb.tsv"), col_names = FALSE)

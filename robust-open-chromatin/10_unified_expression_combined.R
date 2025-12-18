library(tidyverse)
library(edgeR)

# For each species
for (species in c("AtlanticSalmon", "RainbowTrout")) {
  
  # Read in robust ATAC peaks
  peaks <- read_tsv(paste0("unified_annotated_peaks/", species, "_unified_peaks.bed"), 
                    col_names = c("chr", "start", "end", "name", "score", "strand", "count"), 
                    show_col_types = FALSE)
  
  # Read in assay counts against robust ATAC peaks
  counts <- list.files("unified_annotated_peaks/assay_expression/counts", paste0(species, ".*counts.txt$"), recursive = TRUE, full.names = TRUE) %>%
    lapply(., function (file) {
      read_tsv(file, show_col_types = FALSE, comment = "#") %>%
        select(7:ncol(.)) %>%
        setNames(., sub(".*bams/", "", sub(".mLb.clN.sorted.bam", "", colnames(.))))
    }) %>%
    bind_cols()
  
  # Read in assay read library sizes
  library_sizes <- list.files("unified_annotated_peaks/assay_expression/counts", paste0(species, ".*counts.txt.summary$"), recursive = TRUE, full.names = TRUE) %>%
    lapply(., function (file) {
      read_tsv(file, show_col_types = FALSE, comment = "#") %>%
        select(2:ncol(.)) %>%
        setNames(., sub(".*bams/", "", sub(".mLb.clN.sorted.bam", "", colnames(.))))
    }) %>%
    bind_cols() %>%
    colSums()
  
  # Normalise by map
  for (map in c("DevMap", "BodyMap")) {
    # Normalise by assay type
    for (assay in c("ATAC", "ChIP-H3K4me1", "ChIP-H3K4me3", "ChIP-H3K27ac", "ChIP-H3K27me3")) {
      
      # Select samples
      samples <- colnames(counts)[grepl(assay, colnames(counts))]
      if (map == "DevMap") {
        samples <- samples[!grepl("(Female)|(Male)", samples)]  
      } else {
        samples <- samples[grepl("(Female)|(Male)", samples)]
      }
      
      # Calculate normalised CPM values
      dge <- DGEList(counts = counts[, samples], lib.size = library_sizes[samples])
      dge <- calcNormFactors(dge)
      cpms <- cpm(dge, normalized.lib.sizes = TRUE)
      cpm_table <- peaks %>%
        bind_cols(cpms)
      
      # Make ordered table
      if (map == "DevMap") {
        cpm_table %>%
          select(chr, start, end, name, score, strand,
                 contains("LateBlastulation"),
                 contains("MidGastrulation"),
                 contains("EarlySomitogenesis"),
                 contains("MidSomitogenesis"),
                 contains("LateSomitogenesis"))
      } else {
        cpm_table %>%
          select(chr, start, end, name, score, strand,
                 contains("Brain_Immature_Female"), contains("Brain_Immature_Male"), contains("Brain_Mature_Female"), contains("Brain_Mature_Male"),
                 contains("Gonad_Immature_Female"), contains("Gonad_Immature_Male"), contains("Gonad_Mature_Female"), contains("Gonad_Mature_Male"),
                 contains("Liver_Immature_Female"), contains("Liver_Immature_Male"), contains("Liver_Mature_Female"), contains("Liver_Mature_Male"),
                 contains("Muscle_Immature_Female"), contains("Muscle_Immature_Male"), contains("Muscle_Mature_Female"), contains("Muscle_Mature_Male"))
      }
      
      # Write to TSV
      cpm_table %>%
        write_tsv(paste0("unified_annotated_peaks/assay_expression/", species, "_unified_peaks_", map, "_", assay, "_CPM.tsv"))
    }
  }
}

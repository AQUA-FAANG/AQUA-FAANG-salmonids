library(tidyverse)

# Create output directories
setwd("/mnt/project/Aqua-Faang/TFBS")

# Categories
promoters <- c("Active_TSS")
enhancers <- c("Active_enhancer")
DevMap <- c("LateBlastulation", "MidGastrulation", "EarlySomitogenesis", "MidSomitogenesis", "LateSomitogenesis")
BodyMap <- c("Brain", "Gonad", "Liver", "Muscle")

for (species in c("AtlanticSalmon", "RainbowTrout")) {
  for (map in c("DevMap", "BodyMap")) {
    # Robust ATAC peak assay CPM
    counts <- read_tsv(paste0("/mnt/project/Aqua-Faang/robust_ATAC_peaks/unified_annotated_peaks/assay_expression/", 
                               species, "_unified_peaks_", map, "_ATAC_CPM.tsv"), 
                       show_col_types = FALSE)
    
    # Element counts
    for (element in c("promoters", "enhancers")) {
      if (element == "promoters") {
        elements <- promoters
      } else {
        elements <- enhancers
      }
      
      if (map == "DevMap") {
        samples <- DevMap
      } else {
        samples <- BodyMap
        if (species == "RainbowTrout") {
          samples <- BodyMap[-2]
        }
      }
      
      select_rows <- sapply(samples, function(x) paste0(".*", x, ".*:", elements, ".*")) %>%
        sapply(., function(y) grepl(y, counts$name, perl = TRUE)) %>%
        rowSums() > 0
      
      # Count table
      count_table <- counts %>%
        filter(select_rows) %>%
        mutate(loc = paste0(chr, ":", start, "-", end)) %>%
        select(loc, contains("ATAC")) %>%
        gather(sample, cpm, 2:ncol(.)) %>%
        mutate(sample = sub(".*ATAC_", "", sub("_R.", "", sample))) %>%
        group_by(loc, sample) %>%
        mutate(cpm = mean(cpm)) %>%
        distinct() %>%
        spread(sample, cpm) %>%
        data.frame(., row.names = 1)
      if (map == "DevMap") {
        count_table <- count_table[, samples]
      }
      
      # Log2 transform and mean-center scale CPMs
      count_table <- log2((count_table + 1)) %>% 
        t() %>% scale(center = TRUE) %>% t() %>%
        as_tibble(rownames = "loc")
    
      # Write TSV
      count_table %>%
        write_tsv(paste0(paste(species, map, element, sep = "_"), ".tsv"))
    }
  }
}

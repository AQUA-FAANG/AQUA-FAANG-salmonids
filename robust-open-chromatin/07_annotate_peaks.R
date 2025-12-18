library(tidyverse)

state_annotations <- read_tsv("../chromatin_states/annotations/state_annotations.tsv") %>%
  mutate(rgb = sapply(colour, function (x) paste(col2rgb(x), collapse = ","))) %>%
  select(annotation, rgb)

list.files("add_chromatin_states", "bed", recursive = TRUE, full.names = TRUE) %>%
  lapply(., function (file) {
    
    print(file)
    
    out <- sub("add_chromatin_states", "annotated_peaks", sub(".bed", "", file))
    
    if(!dir.exists(dirname(out))) {
      dir.create(dirname(out), recursive = TRUE)
    }
    
    data <- file %>% 
      read_tsv(col_select = c(1:6,11,8,15), col_names = FALSE, show_col_types = FALSE) %>%
      setNames(., c("chrom", "chromStart", "chromEnd", "name", "score", "strand", "signalValue", "peak", "annotation"))
    
    data$reps <- data$name %>% 
      sapply(., function(x) {
        paste("Replicates ", paste(sort(unique(unlist(str_split(x, "[_,]")))), collapse = ","))
      })
    
    data <- data %>%
      left_join(
        state_annotations %>%
          rowid_to_column(var = "rank"),
        by = "annotation"
      ) %>%
      distinct() %>%
      arrange(rank) %>%
      select(-rank, -rgb)
    
    data <- data %>%
      group_by(across(c(-annotation))) %>%
      mutate(annotation = paste(annotation, collapse = ", ")) %>%
      ungroup() %>%
      distinct() %>%
      mutate(name = sub(",.*", "", annotation)) %>%
      left_join(state_annotations, by = c("name" = "annotation"))
    
    data %>%
      arrange(chrom, chromStart) %>%
      transmute(
        chrom, 
        chromStart, 
        chromEnd,
        name,
        score,
        strand = ".",
        signalValue,
        pValue = -1,
        qValue = -1,
        peak = peak - chromStart,
        reps,
        annotation,
        rgb
      ) %>%
      write_tsv(paste0(out, ".narrowPeak"), col_names = FALSE)
    
    data %>%
      arrange(chrom, chromStart) %>%
      transmute(
        chrom, 
        chromStart, 
        chromEnd,
        name,
        score,
        strand = ".",
        thickStart = peak,
        thickEnd = peak + 1,
        rgb,
        blockCount = ".",
        blockSizes = ".",
        blockStarts = ".",
        reps,
        annotation
      ) %>%
      write_tsv(paste0(out, ".bed"), col_names = FALSE)
  })








       
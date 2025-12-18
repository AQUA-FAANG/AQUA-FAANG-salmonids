
library(tidyverse)

# Read in curated state annotations for selected models of species / maps
state_annotations <- read_tsv("annotation_reports/state_annotations.tsv") %>%
  mutate(rgb = sapply(colour, function (x) paste(col2rgb(x), collapse = ",")))

# Atlantic salmon - DevMap
data <- state_annotations %>%
  select(state = AS_DM_m12, annotation, rgb) %>%
  filter(!is.na(state)) %>%
  separate_rows(state, sep = "&") %>%
  mutate(Estate = paste0("E", state))
data %>%
  select(Estate, annotation) %>%
  write_tsv("AtlanticSalmon/DevMap/state_annotation.tsv", col_names = FALSE)
data %>%
  select(state, rgb) %>%
  write_tsv("AtlanticSalmon/DevMap/annotation_colours.tsv", col_names = FALSE)

# Atlantic salmon - BodyMap
data <- state_annotations %>%
  select(state = AS_BM_m12, annotation, rgb) %>%
  filter(!is.na(state)) %>%
  separate_rows(state, sep = "&") %>%
  mutate(Estate = paste0("E", state))
data %>%
  select(Estate, annotation) %>%
  write_tsv("AtlanticSalmon/BodyMap/state_annotation.tsv", col_names = FALSE)
data %>%
  select(state, rgb) %>%
  write_tsv("AtlanticSalmon/BodyMap/annotation_colours.tsv", col_names = FALSE)  

# Rainbow trout - DevMap
data <- state_annotations %>%
  select(state = RT_DM_m12, annotation, rgb) %>%
  filter(!is.na(state)) %>%
  separate_rows(state, sep = "&") %>%
  mutate(Estate = paste0("E", state))
data %>%
  select(Estate, annotation) %>%
  write_tsv("RainbowTrout/DevMap/state_annotation.tsv", col_names = FALSE)
data %>%
  select(state, rgb) %>%
  write_tsv("RainbowTrout/DevMap/annotation_colours.tsv", col_names = FALSE)

# Rainbow trout - BodyMap - Brain
data <- state_annotations %>%
  select(state = RT_BM_m10_Brain, annotation, rgb) %>%
  filter(!is.na(state)) %>%
  separate_rows(state, sep = "&") %>%
  mutate(Estate = paste0("E", state))
data %>%
  select(Estate, annotation) %>%
  write_tsv("RainbowTrout/BodyMap/Brain/state_annotation.tsv", col_names = FALSE)
data %>%
  select(state, rgb) %>%
  write_tsv("RainbowTrout/BodyMap/Brain/annotation_colours.tsv", col_names = FALSE)

# Rainbow trout - BodyMap - Liver
data <- state_annotations %>%
  select(state = RT_BM_m10_Liver, annotation, rgb) %>%
  filter(!is.na(state)) %>%
  separate_rows(state, sep = "&") %>%
  mutate(Estate = paste0("E", state))
data %>%
  select(Estate, annotation) %>%
  write_tsv("RainbowTrout/BodyMap/Liver/state_annotation.tsv", col_names = FALSE)
data %>%
  select(state, rgb) %>%
  write_tsv("RainbowTrout/BodyMap/Liver/annotation_colours.tsv", col_names = FALSE)

# Rainbow trout - BodyMap - Muscle
data <- state_annotations %>%
  select(state = RT_BM_m10_Muscle, annotation, rgb) %>%
  filter(!is.na(state)) %>%
  separate_rows(state, sep = "&") %>%
  mutate(Estate = paste0("E", state))
data %>%
  select(Estate, annotation) %>%
  write_tsv("RainbowTrout/BodyMap/Muscle/state_annotation.tsv", col_names = FALSE)
data %>%
  select(state, rgb) %>%
  write_tsv("RainbowTrout/BodyMap/Muscle/annotation_colours.tsv", col_names = FALSE)


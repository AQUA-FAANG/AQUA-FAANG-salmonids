library(fuzzyjoin)
library(dplyr)
library(readr)

OmykUniPeak <- 
  read_tsv("data/unified_peaks_liftover/OmykUniPeak.bed",
           col_names = c("chr","start","end","peakID_Omyk"), col_types = "ciic")
SsalUniPeak_liftover_to_Omyk_A <- 
  read_tsv("data/unified_peaks_liftover/SsalUniPeak_liftover_to_Omyk_A.bed",
           col_names = c("chr","start","end","peakID_Ssal"), col_types = "ciic")

UniPeak_Ssal_overlapping_Omyk <-
  genome_inner_join( OmykUniPeak, SsalUniPeak_liftover_to_Omyk_A ) %>% 
  distinct(peakID_Ssal,peakID_Omyk)

write_tsv(UniPeak_Ssal_overlapping_Omyk,"data/unified_peaks_liftover/UniPeak_Ssal_overlapping_Omyk.tsv")

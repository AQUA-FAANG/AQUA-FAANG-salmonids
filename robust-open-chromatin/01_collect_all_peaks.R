library(tidyverse)

atac_peaks <- list.files(
  paste0("../../seq_results/", 
         c("AtlanticSalmon/BodyMap/ATAC/Brain",
           "AtlanticSalmon/BodyMap/ATAC/Gonad",
           "AtlanticSalmon/BodyMap/ATAC/Liver",
           "AtlanticSalmon/BodyMap/ATAC/Muscle",
           "AtlanticSalmon/DevMap/ATAC",
           "RainbowTrout/BodyMap/ATAC/Brain",
           "RainbowTrout/BodyMap/ATAC/Gonad",
           "RainbowTrout/BodyMap/ATAC/Liver",
           "RainbowTrout/BodyMap/ATAC/Muscle",
           "RainbowTrout/DevMap/ATAC"),
         "/results/bwa/mergedLibrary/macs/narrowPeak"),
  pattern = ".narrowPeak$", full.names = TRUE)

system("rm -fr ATAC_narrowPeak")

atac_peaks %>%
  sapply(., function (file) {
    out <- sub("../../seq_results/", "ATAC_narrowPeak/", 
               sub("ATAC.*narrowPeak/", "",
                   sub(".mLb.clN_peaks", "", file)))
    if (!dir.exists(dirname(out))) {
      dir.create(dirname(out), recursive = TRUE)
    }
    file.link(file, out)
  })

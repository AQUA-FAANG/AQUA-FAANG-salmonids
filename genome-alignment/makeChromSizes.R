
# fasta files for each genome
spc2fasta=list(
  Ssal="/mnt/project/ELIXIR/salmobase/datafiles/genomes/AtlanticSalmon/Ssal_v3.1/sequence_Ensembl/Salmo_salar.Ssal_v3.1.dna_sm.toplevel.fa.gz",
  Omyk="/mnt/project/ELIXIR/salmobase/datafiles/genomes/RainbowTrout/USDA_OmykA_1.1/sequence_Ensembl/Oncorhynchus_mykiss.USDA_OmykA_1.1.dna_sm.toplevel.fa.gz",
  Eluc="data/Eluc_genome/Esox_lucius.Eluc_v4.dna_sm.toplevel.fa"
)

library(tidyverse)

lapply(spc2fasta, function(fasta){
  fai_file <- paste0(fasta,".fai")
  # read first two columns of fai file
  read_tsv(fai_file, col_names = c("chr","size"), col_types = "cc---")
}) %>%
  bind_rows(.id="genome") %>%
  write_tsv("chromSizes.tsv", col_names=F)



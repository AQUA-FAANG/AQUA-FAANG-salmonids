library(tidyverse)

out_dir = "data/blocks"
syntenyRegionXL = "Salmonid_Synteny_for_alignments_2023.04.12.xlsx"

spc2fasta=c(
  Ssal="/mnt/project/ELIXIR/salmobase/datafiles/genomes/AtlanticSalmon/Ssal_v3.1/sequence_Ensembl/Salmo_salar.Ssal_v3.1.dna_sm.toplevel.fa.gz",
  Omyk="/mnt/project/ELIXIR/salmobase/datafiles/genomes/RainbowTrout/USDA_OmykA_1.1/sequence_Ensembl/Oncorhynchus_mykiss.USDA_OmykA_1.1.dna_sm.toplevel.fa.gz",
  Eluc="data/Eluc_genome/Esox_lucius.Eluc_v4.dna_sm.toplevel.fa"
)


seqFile <- 
  c("(((Ssal_A,Omyk_A),(Ssal_B,Omyk_B)),Eluc);",
    "",
    "Ssal_A Ssal_A.fa",
    "Ssal_B Ssal_B.fa",
    "Omyk_A Omyk_A.fa",
    "Omyk_B Omyk_B.fa",
    "Eluc Eluc.fa"
  )


extractBlockSequenes <- function(tbl){
  blockID = tbl$blockID[1]
  block_dir = file.path(out_dir,blockID)
  dir.create(block_dir,recursive = T,showWarnings = F)

  tbl %>% 
    group_by(genome) %>% 
    mutate(region=str_glue("{chr}:{start}-{end}")) %>% 
    summarise( cmd = str_glue("samtools faidx {spc2fasta[spc[1]]} {paste(region,collapse=' ')}") ) %>% 
    mutate(cmd = paste(cmd, str_glue("> {block_dir}/{genome}.fa"))) %>% 
    with(cmd) %>% 
    write_lines(file.path(block_dir,"extractSequences.sh"))
  
  write_lines(seqFile,file.path(block_dir,"seqFile.txt"))
}

readxl::read_xlsx(syntenyRegionXL) %>%
  mutate(blockID=`Ssal_A-chr`) %>% 
  pivot_longer(cols = -blockID,
               names_to = c("genome",".value"),
               names_pattern = "(.*)-(.*)") %>% 
  mutate(spc=sub("_[AB]$","",genome)) %>%
  mutate(chr=sub("(omy|ssa)0?","",chr)) %>%
  na.omit() %>% 
  split(., f = .$blockID) %>% 
  sapply( extractBlockSequenes )
  
  
  
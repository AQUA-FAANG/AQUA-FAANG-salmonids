# split and convert bed to block coordinates

# input: 
# - bed file with Omyk (or Eluc for that matter) coordinates
# - genome (e.g. Omyk_A)
# - synteny block table
# - output prefix
# output:
# - One bed file per salmon block with coordinates relative to the blocks


args <- commandArgs(trailingOnly = TRUE)
if(length(args) != 4){
  stop("Usage: Rscript split_bed_to_blocks.R <bed_file> <syntenyRegion.xlsx> <input_genome> <output_prefix>")
}

bed_file <- args[1]
syntenyRegionXL <- args[2]
input_genome <- args[3]
output_prefix <- args[4]

# bed_file = "~/Downloads/RainbowTrout_consensus_peaks.bed"
# syntenyRegionXL = "~/Downloads/Salmonid_Synteny_for_alignments_2023.04.12.xlsx"
# input_genome = "Omyk_A"
# output_prefix = "~/Downloads/RainbowTrout_consensus_peaks/"


library(tidyverse)
library(fuzzyjoin)



# read the synteny block table
blockTbl <- 
  readxl::read_xlsx(syntenyRegionXL) %>%
  mutate(blockID=`Ssal_A-chr`) %>% 
  pivot_longer(cols = -blockID,
               names_to = c("genome",".value"),
               names_pattern = "(.*)-(.*)") %>% 
  mutate(spc=sub("_[AB]$","",genome)) %>%
  mutate(chr=sub("(omy|ssa)0?","",chr)) %>%
  na.omit()


bed_tbl <- read_tsv(bed_file,col_names = F,
                    col_types = cols(X2=col_integer(), X3=col_integer(), .default = col_character()))

#
# 0-based:   0 1 2 3 4
#            v v v v v
# sequence:   A C G T
#             ^ ^ ^ ^
# 1-based:    1 2 3 4 
#
# E.g. The subsequence "CG" has 1-based: start=2, end=3 and 0-based: start=1, end=3 
# genome_join uses 1-based, blocks are extracted using samtools which uses 1-based
# bed file is 0-based

# Note regarding the 0-based/1-based coordinates of blocks in the synteny table:
# The coordinates are used to extract the subsequence using "samtools faidx" which
# expects 1-based coordinates but tolerates being given 0.
# This means that when converting coordinates, the block coordinates should be
# treated as 1-based and any values of 0 needs to corrected to 1.
# However, since the block coordinates is also used in the chromosome name they 
# have to remain the same in that context, i.e. don't convert 0's to 1's.
# (All this could have been avoided if the )

bed_tbl %>% 
  mutate(X2=X2+1) %>% # make start coordinates 1-based for the genome join
  # join the blocks that match coordinates of given genome
  genome_inner_join(filter(blockTbl,genome==input_genome) %>% select(-genome,-spc), 
                    by=c("X1"="chr","X2"="start","X3"="end")) %>% 
  mutate(X2=X2-1) %>% # return start coordinates to 0-based 
  # set chr = block_chr:start-end
  mutate(X1=paste0(chr,":",start,"-",end)) %>% 
  # replace block start 0's to 1's (see above note)
  mutate(start = pmax(1,start)) %>% 
# convert start and end by subtracting block start. Also keep the coordinates within the block
  mutate(X2=pmax(0,X2-(start-1))) %>% 
  mutate(X3=pmin(end-start+1,X3-(start-1))) %>% 
  select(-chr,-start,-end) %>% 
  # split by blockID
  split(., f = .$blockID) %>% 
  lapply(function(bed_block){
    # save each block to file
    fileName = paste0(output_prefix,bed_block$blockID[1],".bed")
    write_tsv(bed_block[,-which(colnames(bed_block)=="blockID")], file=fileName,col_names=F)
  }) %>% 
  invisible()
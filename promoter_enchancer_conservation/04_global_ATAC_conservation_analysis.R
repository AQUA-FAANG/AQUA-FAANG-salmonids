#### load libraries
library(tidyverse)
library(openxlsx)
library(fuzzyjoin)
library(tidygenomics)

########################################################
# 1. Load chromosome-level outputs
# 2. Merge all chromosomes
# 3. Remove problematic alignments
# 4. Generate final curated dataset
# 5. Generate conservation summaries
# 6. Export manuscript tables
# 7. Generate salmon duplicate peak dataset
# 8. Generate salmon singleton peak dataset
# 9. Generate alignable peak dataset
###############################################################

###############################################################################
# Load chromosome-specific conservation results
###############################################################################

files<-list.files(pattern="_salmon_trout_atac_maf_final.txt$",recursive = TRUE)

All_files <- lapply(files,function(i){
  read.table(i,header = T)
})

all_atac_data<- bind_rows(All_files)

all_atac_data<-all_atac_data %>% unite(col = "atac_loc", chromosome.y,start.y,end.y, sep = "_",remove = FALSE,na.rm = TRUE)

###############################################################################
# Remove problematic alignments
#
# 1. Duplicate alignments
# 2. Nested peak mappings
# 3. Same peak mapping to duplicated regions
# 4. Multiple peak assignments
###############################################################################


duplicate_alignments<-all_atac_data %>% group_by(Origin,species) %>% mutate(dupe = n()>1) %>% ungroup %>% filter(dupe==TRUE) %>% distinct(Origin)

all_atac_data<-all_atac_data[!(all_atac_data$Origin %in% duplicate_alignments$Origin),]

nested_peaks<-all_atac_data %>%
  filter(!is.na(chromosome.y)) %>% 
  filter(grepl("Ssal*", species)) %>%
  group_by(atac_loc) %>%
  filter(n()>2) 

nested_to_remove<-nested_peaks %>% distinct(Origin) ##9250 peaks nested

curated_no_nested<-all_atac_data[!(all_atac_data$Origin %in% nested_to_remove$Origin),]
curated_no_nested<-curated_no_nested %>% mutate_all(na_if,"")


duplicated_peak_mappings<-curated_no_nested %>% filter(!is.na(chromosome.y)) %>% 
  filter(grepl("Ssal*", species)) %>%
  select(atac_loc,Origin) %>% 
  group_by(atac_loc) %>%
  filter(n()==2) %>%
  summarise(y=paste(Origin,collapse = ",")) %>%
  separate(col=y, into=c("o1","o2"),sep = ",") %>%
  separate(col=o1, into=c("z1","z2"),sep = "_",remove = FALSE) %>%
  separate(col=o2, into=c("z3","z4"),sep = "_",remove = FALSE) %>%
  select(-c("z2","z4")) %>%
  mutate(dups=ifelse(z1==z3,"dups","nodups")) %>%
  filter(dups=="dups")

remove_odds<- rbind(duplicated_peak_mappings %>% select(o1),
                    duplicated_peak_mappings %>% select(o2) %>% rename(o1 = o2))

curated_no_duplicates<-curated_no_nested[!(curated_no_nested$Origin %in% remove_odds$o1),]
curated_no_duplicates<-curated_no_duplicates %>% mutate_all(na_if,"")

curated_no_duplicates %>% distinct(Origin)

overlapping_peak_regions<-curated_no_duplicates %>%
  filter(!is.na(chromosome.y)) %>% 
  filter(grepl("Ssal*", species)) %>%
  group_by(atac_loc) %>%
  filter(n()>1) %>%
  ungroup()%>%
  arrange(atac_loc,Origin) %>%
  distinct(atac_loc,.keep_all = TRUE) %>%
  distinct(Origin)

curated_final<-curated_no_duplicates[!(curated_no_duplicates$Origin %in% overlapping_peak_regions$Origin),]
curated_final<-curated_final %>% mutate_all(na_if,"")

curated_final %>% distinct(Origin)## to get the counts of disctinct alignment regions

###############################################################################
# Export final curated conservation dataset
###############################################################################
write.table(curated_final,"Final_merged_ATAC_peaks_curated_2024.txt",quote = F,na = "NA",col.names = T,row.names = F)

###############################################################################
# Global conservation patterns
#
# Summarize:
#   - Peak conservation
#   - Sequence conservation
#   - Species combinations
###############################################################################

peak_chr_combs<-curated_final %>% distinct(Origin,.keep_all=T) %>% count(chr2sp_simple,peak2sp_simple)

write.xlsx(peak_chr_combs,"global_salmon_peaks_chromosomes_4way_combs_2024.xlsx")

peaks_combination<-curated_final %>% select(Origin,peak2sp) %>% 
  distinct(Origin,.keep_all = TRUE) %>% count(peak2sp,name="counts") %>%
  group_by(test = map_chr(strsplit(as.character(peak2sp),","), ~paste(sort(as.character(.x)),collapse = ","))) %>% 
  summarise(peak2sp=peak2sp, count=sum(counts),combs=map_chr(strsplit(as.character(peak2sp),","), ~length(as.character(.x)))) %>% 
  arrange(test) %>% group_by(test) %>% slice(1)%>% ungroup %>% select(-c(test))

write.xlsx(peaks_combination,"global_peaks_comb_2024.xlsx")

seq_combination<-curated_final %>% select(Origin,chr2sp) %>% 
  distinct(Origin,.keep_all = TRUE) %>% count(chr2sp,name="counts") %>%
  group_by(test = map_chr(strsplit(as.character(chr2sp),","), ~paste(sort(as.character(.x)),collapse = ","))) %>% 
  summarise(chr2sp=chr2sp, count=sum(counts),combs=map_chr(strsplit(as.character(chr2sp),","), ~length(as.character(.x)))) %>% 
  arrange(test) %>% group_by(test) %>% slice(1)%>% ungroup %>% select(-c(test))

write.xlsx(seq_combination,"global_seq_comb_2024.xlsx")

simple_peak_comb <- curated_final %>% select(Origin,peak2sp_simple) %>% 
  distinct(Origin,.keep_all = TRUE) %>% count(peak2sp_simple,name="counts") 

write.xlsx(simple_peak_comb,"global_simple_peak_comb_2024.xlsx")

simple_seq_comb <- curated_final %>% select(Origin,chr2sp_simple) %>% 
  distinct(Origin,.keep_all = TRUE) %>% count(chr2sp_simple,name="counts") 

write.xlsx(simple_seq_comb,"global_simple_seq_comb_2024.xlsx")


###############################################################################
# Four-way conserved regulatory elements
#
# Conserved in:
#   Ssal_A
#   Ssal_B
#   Omyk_A
#   Omyk_B
###############################################################################

## check for names here

conserved_4way<- curated_final %>% filter(peak_count==4) %>% select(1,c(18:38))

summary_4vs4<-conserved_4way %>% group_by(Origin) %>% summarise(across(.fns = sum)) %>%
  pivot_longer(cols=c(2:22),
               names_to='stage',
               values_to='peaks_open') %>%
  filter(peaks_open==4)


summ<-summary_4vs4 %>% count(stage) %>% left_join(atac_numbers, by="stage") %>% mutate(perc_conserved=round((n/total)*100,2))

write.xlsx(summ,"global_atac_conservation_4way_2024.xlsx")

###############################################################################
# Salmon duplicated regulatory elements
#
# Peak present on both salmon homeologous chromosomes - Shared peaks
###############################################################################

# use peak summit data and filter for peak overlap

salmon_summit_files<-list.files("../../atac_peak_summits/salmon/",pattern='summits',recursive = TRUE,full.names=TRUE)

All_summits <- lapply(salmon_summit_files,function(i){
  read.table(i,header = F,col.names = c("chromosome.y","overlap_start","overlap_end"))
})

salmon_summits<- bind_rows(All_summits)
salmon_summits$chromosome.y <- sprintf("%02d",salmon_summits$chromosome.y)
salmon_summits$chromosome.y = paste0("Ssal",salmon_summits$chromosome.y)

salmon_summits<-salmon_summits %>% mutate(chromosome=paste("Ssal",chromosome,sep = ""))

####### filter for salmon

salmon_atac<-curated_final %>% filter(grepl("Ssal*", species)) %>% 
  select(-c(sequence_count,peak_count,chr2sp,peak2sp,chr2sp_simple,peak2sp_simple,Species_order)) %>%
  group_by(Origin) %>%
  mutate(sequence_count=length(Origin)) %>%
  mutate(peak_count=sum(!is.na(chromosome.y))) %>% dplyr::rename(testis_immature_male="gonad_immature_male",
                                                                 ovary_immature_female="gonad_immature_female",
                                                                 ovary_mature_female="gonad_mature_female",
                                                                 testis_mature_male="gonad_mature_male") %>% 
  ungroup %>% mutate(species=as.character(species))


duplicated_peak_candidates<- salmon_atac %>% filter(peak_count==2) %>% filter(sequence_count==2) %>% 
  mutate(overlap_start=ifelse(species=="Ssal_B",pmax(start.x,start.y),NA),overlap_end=ifelse(species=="Ssal_B",pmin(end.x,end.y),NA)) %>% 
  mutate(chromosome.y=as.character(chromosome.y)) 


####get overlaps start and end for ssalA

files<-list.files(pattern="salmon_trout_atac_maf_raw.xlsx$",recursive = TRUE)

All_raw_files <- lapply(files,function(i){
  read.xlsx(i,colNames = T)
})

#####bind rows
All_raw_Data<- bind_rows(All_raw_files)

All_raw_Data<-All_raw_Data %>% unite(col = "atac_loc", chromosome.y,start.y,end.y, sep = "_",remove = FALSE,na.rm = TRUE)

salmon_overlap_atac<-All_raw_Data %>% filter(grepl("Ssal*", species)) %>% 
  group_by(Origin) %>%
  mutate(sequence_count=length(Origin)) %>%
  mutate(peak_count=sum(!is.na(chromosome.y))) %>% dplyr::rename(testis_immature_male="gonad_immature_male",
                                                                 ovary_immature_female="gonad_immature_female",
                                                                 ovary_mature_female="gonad_mature_female",
                                                                 testis_mature_male="gonad_mature_male") %>% ungroup

two_peaks_overlap<- salmon_overlap_atac %>% filter(peak_count==2) %>% filter(sequence_count==2) %>%
  mutate(overlap_s=pmax(start.x,start.y), overlap_e=pmin(end.x,end.y)) %>% mutate(chromosome.y=as.character(chromosome.y)) %>%
  select(c("atac_loc","overlap_s","overlap_e"))

resecued_to_add_later<-duplicated_peak_candidates %>% left_join(two_peaks_overlap,by="atac_loc") %>%
  unite(col = "a", overlap_start,overlap_end, sep = "_",remove = FALSE) %>%
  unite(col = "b", overlap_s,overlap_e, sep = "_",remove = FALSE) %>%
  group_by(Origin,species) %>%
  mutate(o_size=overlap_e - overlap_s) %>% 
  ungroup() %>%
  mutate(overlap_start=ifelse(is.na(overlap_start),overlap_s,overlap_start)) %>%
  mutate(overlap_end=ifelse(is.na(overlap_end),overlap_e,overlap_end)) %>%
  select(-c("a","b","overlap_s","overlap_e","sequence_count","peak_count","o_size")) %>%
  filter(overlap_start==start.y & overlap_end==end.y) %>%
  group_by(Origin,species) %>%
  distinct(Origin,.keep_all = TRUE) %>% ungroup %>%
  group_by(Origin) %>%
  mutate(sequence_count=length(Origin)) %>%
  mutate(peak_count=sum(!is.na(chromosome.y))) %>% ungroup %>%
  unite(col = "addit", Origin, species, sep = "_",remove = FALSE)


duplicated_peak_filtered<-duplicated_peak_candidates %>% left_join(two_peaks_overlap,by="atac_loc") %>%
  unite(col = "a", overlap_start,overlap_end, sep = "_",remove = FALSE) %>%
  unite(col = "b", overlap_s,overlap_e, sep = "_",remove = FALSE) %>%
  group_by(Origin,species) %>%
  mutate(o_size=overlap_e - overlap_s) %>% 
  ungroup() %>%
  mutate(overlap_start=ifelse(is.na(overlap_start),overlap_s,overlap_start)) %>%
  mutate(overlap_end=ifelse(is.na(overlap_end),overlap_e,overlap_end)) %>%
  select(-c("a","b","overlap_s","overlap_e","sequence_count","peak_count","o_size")) %>%
  filter(!(overlap_start==start.y & overlap_end==end.y)) %>%
  group_by(Origin,species) %>%
  distinct(Origin,.keep_all = TRUE) %>% ungroup %>%
  group_by(Origin) %>%
  mutate(sequence_count=length(Origin)) %>%
  mutate(peak_count=sum(!is.na(chromosome.y))) %>% ungroup


#pullout 1 sequence and 1 peak combos to add back to the table above

to_add<-duplicated_peak_filtered %>% filter(sequence_count==1 & peak_count==1) %>% select(Origin,species) %>% mutate(spp=ifelse(species=="Ssal_B","Ssal_A","Ssal_B")) %>%
  select(-c(species)) %>%  unite(col = "addition", Origin, spp, sep = "_")

rescue_origins<-to_add$addition

final_rescued<-resecued_to_add_later %>% filter(addit %in% rescue_origins) %>% select(-c(addit))

duplicated_peak_final<-bind_rows(duplicated_peak_filtered,final_rescued) %>% select(-c("sequence_count","peak_count")) %>% group_by(Origin) %>%
  mutate(sequence_count=length(Origin)) %>%
  mutate(peak_count=sum(!is.na(chromosome.y))) %>% ungroup

###overlap summit data on the peak data

two_two_summits<-duplicated_peak_final %>% fuzzyjoin::genome_join(salmon_summits, by = (c("chromosome.y","overlap_start","overlap_end")),mode="left") %>%
  group_by(Origin,species) %>%
  distinct(Origin,species,.keep_all = TRUE) %>% filter(!(is.na(chromosome.y.y))) %>% group_by(Origin) %>%
  mutate(count=length(Origin)) %>% ungroup %>% filter(count==2) %>% mutate(Origin=as.character(Origin)) %>% select(-c("sequence_count","peak_count"))

summary_2peaks_summit<-two_two_summits %>% select(1,c(15:35)) %>% group_by(Origin) %>% summarise(across(.fns = sum)) %>%
  pivot_longer(cols=c(2:22),
               names_to='stage',
               values_to='peaks_open') %>%
  filter(peaks_open==2) 

summary_2peaks_summit %>% distinct(Origin)

summ2<-summary_2peaks_summit %>% count(stage) %>% left_join(atac_numbers, by="stage") %>% mutate(perc_conserved=round((n/total)*100,2))

write.xlsx(summ2,"salmon_summit_atac_conservation_2way_summary_2024.xlsx")

final_two_peaks_summit<-summary_2peaks_summit %>% distinct(Origin) %>% left_join(duplicated_peak_candidates, by="Origin")

write.xlsx(final_two_peaks_summit,"salmon_atac_conservation_2to2_summit_cutoff_2024.xlsx")
write.table(final_two_peaks_summit,"salmon_atac_conservation_2to2_summit_cutoff_2024.txt",quote = F,na = "NA",col.names = T,row.names = F)

###to share for publication
final_two_peaks_summit_2025_share<-final_two_peaks_summit %>% rename("Peak_ID" = "Origin", "Alignment_chromosome" = "chromosome.x",  "Alignment_start" = "start.x",  "Alignment_end" = "end.x",
                                                                     "ATAC_ID" = "atac_loc", "ATAC_chromosome" = "chromosome.y",  "ATAC_start" = "start.y",  "ATAC_end" = "end.y") %>%
  select(-c(species_b))


write.table(final_two_peaks_summit_2025_share,"Final_Shared_ATAC_peaks_curated_2025_manuscript.txt", quote = F,na = "NA",col.names = T,row.names = F)

################ get simplified representations in bed format

duplicated_peaks_bed<-final_two_peaks_summit %>% select(Origin,chromosome.y,start.y,end.y)

duplicated_peaks_bed$ID <- cumsum(!duplicated(duplicated_peaks_bed$Origin))

duplicated_peaks_bed <- duplicated_peaks_bed %>% select(-c(Origin))

write.table(duplicated_peaks_bed,"salmon_duplicated_peaks_2024.bed",quote = F,na = "NA",col.names = T,row.names = F)

###############################################################################
# Duplicate-specific regulatory elements - Exclusive peaks
#
# One sequence
# One peak
###############################################################################

exclusive_peaks<- salmon_atac %>% filter(sequence_count==1) %>% filter(peak_count==1) %>% select(1,c(15:35))

one_peak_1<- salmon_atac %>% filter(sequence_count==1) %>% filter(peak_count==1)

write.xlsx(one_peak_1,"salmon_atac_conservation_1to1_summit_cutoff_2024.xlsx")

write.table(one_peak_1,"salmon_atac_conservation_1to1_summit_cutoff_2024.txt",quote = F,na = "NA",col.names = T,row.names = F)

### save a version for manuscript

one_peak_2025_share<-one_peak_1 %>% rename("Peak_ID" = "Origin", "Alignment_chromosome" = "chromosome.x",  "Alignment_start" = "start.x",  "Alignment_end" = "end.x",
                                           "ATAC_ID" = "atac_loc", "ATAC_chromosome" = "chromosome.y",  "ATAC_start" = "start.y",  "ATAC_end" = "end.y") %>%
  select(-c(species_b))


write.table(one_peak_2025_share,"Final_exclusive_ATAC_peaks_curated_2025_manuscript.txt", quote = F,na = "NA",col.names = T,row.names = F)


### get a bed file

singleton_peaks.bed<-one_peak_1 %>% select(chromosome.y,start.y,end.y)

write.table(singleton_peaks.bed,"salmon_singleton_peaks_2024.bed",quote = F,na = "NA",col.names = T,row.names = F)

summary_1peak<-exclusive_peaks %>% group_by(Origin) %>% summarise(across(.fns = sum)) %>%
  pivot_longer(cols=c(2:22),
               names_to='stage',
               values_to='peaks_open') %>%
  filter(peaks_open==1)

write.table(summary_1peak,"salmon_one_way_conservation_data_2024.txt",quote = F,na = "NA",col.names = T,row.names = F)

summ1<-summary_1peak %>% count(stage) %>% left_join(atac_numbers, by="stage") %>% mutate(perc_conserved=round((n/total)*100,2))

write.xlsx(summ1,"salmon_atac_conservation_1seq_1peak_2024.xlsx") 


###############################################################################
# Alignable but non-conserved elements
#
# Two aligned sequences
# One ATAC peak
###############################################################################

alignable_single_peak<- salmon_atac %>% filter(sequence_count==2) %>% filter(peak_count==1) %>% select(1,c(15:35))

twoone_peak_1<- salmon_atac %>% filter(sequence_count==2) %>% filter(peak_count==1)

write.xlsx(twoone_peak_1,"salmon_atac_conservation_2to1_2024.xlsx")

write.table(twoone_peak_1,"salmon_atac_conservation_2to1_2024.txt",quote = F,na = "NA",col.names = T,row.names = F)

twoone_peak_2025_share<-twoone_peak_1 %>% rename("Peak_ID" = "Origin", "Alignment_chromosome" = "chromosome.x",  "Alignment_start" = "start.x",  "Alignment_end" = "end.x",
                                                 "ATAC_ID" = "atac_loc", "ATAC_chromosome" = "chromosome.y",  "ATAC_start" = "start.y",  "ATAC_end" = "end.y") %>%
  select(-c(species_b))


write.table(twoone_peak_2025_share,"Final_aligned_ATAC_peaks_curated_2025_manuscript.txt", quote = F,na = "NA",col.names = T,row.names = F)

summary_21peak<-alignable_single_peak %>% replace(is.na(.), 0) %>% group_by(Origin) %>% summarise(across(.fns = sum)) %>%
  pivot_longer(cols=c(2:22),
               names_to='stage',
               values_to='peaks_open') %>%
  filter(peaks_open==1)

write.table(summary_21peak,"salmon_twoseq_onepeak_conservation_data_2024.txt",quote = F,na = "NA",col.names = T,row.names = F)

summ21<-summary_21peak %>% count(stage) %>% left_join(atac_numbers, by="stage") %>% mutate(perc_conserved=round((n/total)*100,2))

write.xlsx(summ21,"salmon_atac_conservation_2seq_1peak_2024.xlsx")





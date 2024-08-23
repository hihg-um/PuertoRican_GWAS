# SPDX-License-Identifier: GPL-2.0

library(data.table)
library(tidyverse)
library(gtools)

args <- commandArgs(trailingOnly = TRUE)
dir <- args[1]
file1 <- args[2]
file2 <- args[3]


### If using one data file, use this

df <- fread(file1)
ref <- fread(file2)

### Merge study and HGDP reference plink files

ref$V2 <- paste0("chr",ref$V1,":",ref$V4)
df$V2 <- paste0("chr",df$V1,":",df$V4)

stot.rs <- inner_join(df,ref, by='V2',suffix = c('_1','_2'))
dup <- stot.rs %>% filter(duplicated(V2) | duplicated(V2,fromLast = T))
stot.rs <- stot.rs %>% filter(!V2 %in% dup$V2)

### Remove all rows with missing allele status, D allele, and I allele

stot1 <- stot.rs %>% mutate(V5_1 = replace(V5_1, V5_1 == '.', NA))
stot1 <- stot1 %>% mutate(V6_1 = replace(V6_1, V6_1 == '.', NA))
stot1 <- stot1 %>% mutate(V5_2 = replace(V5_2, V5_2 == '.', NA))
stot1 <- stot1 %>% mutate(V6_2 = replace(V6_2, V6_2 == '.', NA))

stot1 <- stot1 %>% mutate(V5_1 = replace(V5_1, V5_1 == 'I', NA))
stot1 <- stot1 %>% mutate(V6_1 = replace(V6_1, V6_1 == 'I', NA))
stot1 <- stot1 %>% mutate(V5_2 = replace(V5_2, V5_2 == 'I', NA))
stot1 <- stot1 %>% mutate(V6_2 = replace(V6_2, V6_2 == 'I', NA))

stot1 <- stot1 %>% mutate(V5_1 = replace(V5_1, V5_1 == 'D', NA))
stot1 <- stot1 %>% mutate(V6_1 = replace(V6_1, V6_1 == 'D', NA))
stot1 <- stot1 %>% mutate(V5_2 = replace(V5_2, V5_2 == 'D', NA))
stot1 <- stot1 %>% mutate(V6_2 = replace(V6_2, V6_2 == 'D', NA))

stot1 <- stot1 %>% mutate(V5_1 = replace(V5_1, V5_1 == '0', NA))
stot1 <- stot1 %>% mutate(V6_1 = replace(V6_1, V6_1 == '0', NA))
stot1 <- stot1 %>% mutate(V5_2 = replace(V5_2, V5_2 == '0', NA))
stot1 <- stot1 %>% mutate(V6_2 = replace(V6_2, V6_2 == '0', NA))

stot2 <- na.omit(stot1)
#table(stot.rs$V5_1)
#table(stot2$V5_1)


### Remove duplicate snpids

dup <- stot2 %>% filter(duplicated(V2) | duplicated(V2,fromLast = T))
stot2 <- stot2 %>% filter(!V2 %in% dup$V2)

stot2$gt_1 <- paste0(stot2$V5_1,"_",stot2$V6_1)
stot2$gt_2 <- paste0(stot2$V5_2,"_",stot2$V6_2)



### Filter for matching or flipped genotypes

gt_ac_tg <- filter(stot2, (gt_1 == 'A_C' | gt_1 == 'C_A' | gt_1 == 'T_G' | gt_1 == 'G_T') 
                        & (gt_2 == "A_C" | gt_2 == "C_A" | gt_2 == 'T_G' | gt_2 == 'G_T'))
             
gt_ag_ct <- filter(stot2, (gt_1 == 'A_G' | gt_1 == 'G_A' | gt_1 == 'T_C' | gt_1 == 'C_T') 
                        & (gt_2 == "A_G" | gt_2 == "G_A" | gt_2 == 'T_C' | gt_2 == 'C_T'))
             

#gt_ac <- filter(stot2, (gt_1 == 'A_C' | gt_1 == 'C_A') & (gt_2 == "A_C" | gt_2 == "C_A")) 
#gt_ag <- filter(stot2, (gt_1 == 'A_G' | gt_1 == 'G_A') & (gt_2 == "A_G" | gt_2 == "G_A")) 
#gt_at <- filter(stot2, (gt_1 == 'A_T' | gt_1 == 'T_A') & (gt_2 == "A_T" | gt_2 == "T_A")) 
#gt_cg <- filter(stot2, (gt_1 == 'C_G' | gt_1 == 'G_C') & (gt_2 == "C_G" | gt_2 == "G_C")) 
#gt_ct <- filter(stot2, (gt_1 == 'C_T' | gt_1 == 'T_C') & (gt_2 == "C_T" | gt_2 == "T_C")) 
#gt_tg <- filter(stot2, (gt_1 == 'T_G' | gt_1 == 'G_T') & (gt_2 == "T_G" | gt_2 == "G_T")) 
#gt_merge <- rbind(gt_ac,gt_ag,gt_ct,gt_tg)

gt_merge_wflip <- rbind(gt_ac_tg,gt_ag_ct) %>% arrange(V1_1,V4_1)
gt_merge_wflip$flip <- gt_merge_wflip$gt_1==gt_merge_wflip$gt_2

knitr::kable(table(gt_merge_wflip$gt_1,gt_merge_wflip$gt_2))

gt_merge_wflip$id1 <- paste0(gt_merge_wflip$V2,":",gt_merge_wflip$V5_1,":",gt_merge_wflip$V6_1)
gt_merge_wflip$id2 <- paste0(gt_merge_wflip$V2,":",gt_merge_wflip$V5_2,":",gt_merge_wflip$V6_2)

keep_rsid <- gt_merge_wflip[,c('V2')]
keep_rsid1 <- gt_merge_wflip[,c('id1')]
keep_rsid2 <- gt_merge_wflip[,c('id2')]

keep_rsid_vcf <- gt_merge_wflip$V2

fwrite(keep_rsid,paste0(dir,"keep_rsid.txt"),col.names = FALSE, row.names = FALSE, quote = FALSE, sep = " ")

#fwrite(keep_rsid1,paste0(dir,"keep_rsid_df.txt"),col.names = FALSE, row.names = FALSE, quote = FALSE, sep = " ")
#fwrite(keep_rsid2,paste0(dir,"keep_rsid_ref.txt"),col.names = FALSE, row.names = FALSE, quote = FALSE, sep = " ")


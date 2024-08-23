# SPDX-License-Identifier: GPL-2.0

##Required libraries
library(GWASTools)
library(GENESIS)
library(GWASdata)
library(gdsfmt)
library(data.table)
library(tidyr)
library(tidyverse)
library(dplyr)


#---------------------------------
rm(list = ls())
showfile.gds(closeall=TRUE)

#Set working directory
dir <- "" ### directory
grm_path <- ""  ### grm matrix path

setwd(dir)

#---------------------------------
### Import covariate file that also contain phenotypes and edit the format
covar <- fread(paste0(dir,"saige_pheno_file.tsv"))

covar$sex[covar$sex == '0'] <- 'M'
covar$sex[covar$sex == '1'] <- 'F'

sample.id_covar = covar$id

pheno = covar$dx

covar_sub = data.frame(age = covar$age,
                       sex = covar$sex,
                       apoe = covar$apoe)


grm_matrix = readRDS(grm_path)

#-------------------------------------
study = "MAIN_pc1_4"

scanAnnot <- ScanAnnotationDataFrame(data.frame(scanID = sample.id_covar, 
                                                covar_sub, pheno, stringsAsFactors=FALSE))

# fit the null mixed model
nullmod <- fitNullModel(scanAnnot, outcome = "pheno", 
                        covars = c("sex","age","pc1","pc2","pc3","pc4"),
                        cov.mat = grm_matrix)


anc_pop = c("afr","amr","eur")   ## Ancestries used in the Phasing step

chr_list = paste0("chr",c(1:22))
files_list = list.files()   ##RFmix2 outputs (msp files) should be inside of the working directory

for (nn in chr_list) {
  print(nn)
  
  rfmix = fread(files_list[grep(paste0("AFR.AMR.EUR.",nn,".LA.msp"),files_list)])  ## RFmix2 outputs (msp files) for each chromosome
  colnames(rfmix)[1] = "chr"

  #afr = 1; amr = 2; eur = 3 (Make sure about the order)
  
  anc_pop_val = c(0,1,2)+1
  
  #snp.chromosome = parse_number(rfmix$chr)
  snp.chromosome =as.numeric(rfmix$chr)
  
  snp.id = as.character(paste0(rfmix$spos,"_",rfmix$epos))
  snp.position = as.numeric(rfmix$spos)
  
  for (i in anc_pop){
    
    colIDs = colnames(rfmix)[7:ncol(rfmix)]
    substr_colIDs <- sapply(strsplit(colIDs,"[.]"),function(x) x[[1]])
    matched.id <- sapply(substr_colIDs, function(x) any(sapply(sample.id_covar, grepl, x = x))) %>% which() %>% unname() + 6 
    temp = rfmix %>% select(all_of(matched.id)) %>% as.matrix() + 1
    
    
    criteria = which(temp %in% anc_pop_val[-which(anc_pop == i)])
    
    temp[criteria] = 0
    temp[which(temp != 0)] = 1
    
    col.ind = c(1:ncol(temp)) %% 2
    temp_tot = temp[,col.ind == 0] + temp[,col.ind == 1]
    
    colnames(temp_tot) = sapply(strsplit(colnames(temp_tot),"[.]"), function(x) x[[1]])
    temp_tot_t = t(temp_tot)
    assign(paste0("dosage_",i),temp_tot)
  }
  
  #----------------------------------
  # Make gds file
  
  gfile = createfn.gds("admix_map.gds")
  add.gdsn(gfile, "snp.id",snp.id )
  add.gdsn(gfile, "sample.id", sample.id_covar)
  add.gdsn(gfile, "snp.chromosome", snp.chromosome)
  add.gdsn(gfile, "snp.position", snp.position)
  
  
  for (i in anc_pop) {
    add.gdsn(gfile, paste0("dosage_",i), get(paste0("dosage_",i)))  
  }
  
  
  # read in GDS data
  genoDataList <- list()
  for (anc in anc_pop){
    gdsr <- GdsGenotypeReader(gfile, genotypeVar=paste0("dosage_", anc))
    genoDataList[[anc]] <- GenotypeData(gdsr, scanAnnot=scanAnnot)
  }
  
  # run the association test
  print(paste0(nn,"--association test"))
  anc_pop_run = c(anc_pop,"afr_amr","afr_eur","amr_eur")
  run_comb = list(1,2,3,c(1,2),c(1,3),c(2,3))
  # iterators
  for (i in 1:4) {
    genoIterators <- lapply(genoDataList[run_comb[[i]]], GenotypeBlockIterator)
    myassoc <- admixMap(genoIterators, nullmod)
    
    if( nn == "chr1"){
      assign(paste0("myassoc_", anc_pop_run[i]), myassoc)    
    }else{
      temp = rbind(get(paste0("myassoc_",anc_pop_run[i])), myassoc)
      assign(paste0("myassoc_", anc_pop_run[i]), temp)
    }
    
  }
  
  closefn.gds(gfile)
  showfile.gds(closeall=TRUE)
  
}

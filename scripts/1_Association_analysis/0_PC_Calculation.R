##Required libraries
library(data.table)
library(tidyverse)
library(gdsfmt)
library(SNPRelate)
library(GWASTools)
library(GENESIS)

#---------------------------------
rm(list = ls())
showfile.gds(closeall=TRUE)

#-------------------------------------------------------------------------------
# STEP 0---PREP
#-------------------------------------------------------------------------------
#Set working directory
setwd("")

dir <- "" ### directory
plink <- ""  ### plink file prefix


#Import covariate file that  contains phenotypes and edit the format

covar = fread("Covariates_Pheno.txt")
colnames(covar)[1] = "ID"

#Convert plink files to gds format

snpgdsBED2GDS(bed.fn = paste0(dir,plink,".bed"),
              bim.fn = paste0(dir,plink,".bim"),
              fam.fn = paste0(dir,plink,".fam"),
              out.gdsfn = paste0(dir,plink,".gds"),
              family=FALSE,snpfirstdim=NA, compress.annotation="ZIP_RA.max", compress.geno="",
              option=NULL, cvt.chr=c("int", "char"), cvt.snpid=c("auto", "int"),
              verbose=TRUE)

showfile.gds(closeall=TRUE)


#-------------------------------------------------------------------------------
# STEP 1---calc kinship matrix
#-------------------------------------------------------------------------------

showfile.gds(closeall=TRUE)

gds = snpgdsOpen(paste0(dir,plink,".gds"))

#pruned gds file
snpset <- snpgdsLDpruning(gds, method="corr", slide.max.bp=10e6,maf=0.05, 
                          ld.threshold=sqrt(0.2), verbose=FALSE)

pruned <- unlist(snpset, use.names=FALSE)
length(pruned)

snpgdsClose(gds)
showfile.gds(closeall=TRUE)

gds = snpgdsOpen(paste0(dir,plink,".gds"))

ibd_KING = snpgdsIBDKING(gds, sample.id=NULL, snp.id=pruned,autosome.only = TRUE, maf=NaN, type="KING-robust", family.id=NULL, missing.rate = NaN)
data_ibdKING <- snpgdsIBDSelection(ibd_KING)

#Extract matrix of kinship coefficients
ibd_KING_matrix = ibd_KING$kinship
dimnames(ibd_KING_matrix) = list(ibd_KING$sample.id, ibd_KING$sample.id)

snpgdsClose(gds)
showfile.gds(closeall=TRUE)

#-------------------------------------------------------------------------------
# STEP 2---Principal Components Analysis
#-------------------------------------------------------------------------------

# PCAir calculation and GRM calculation
showfile.gds(closeall=TRUE)

mygeno <- GdsGenotypeReader(paste0(dir,plink,".gds"))
mygenoData <- GenotypeData(mygeno)

mypcair <- pcair(mygenoData,
                 kinobj = ibd_KING_matrix, kin.thresh=2^(-11/2),
                 divobj = ibd_KING_matrix, div.thresh=-2^(-11/2),
                 snp.include = pruned,
                 sample.include = NULL,
                 unrel.set = NULL)
pcair_sample_sets <- pcairPartition(kinobj = ibd_KING_matrix, divobj = ibd_KING_matrix)
mygenoData <- GenotypeBlockIterator(mygenoData, snpInclude = pruned)

### format file
  pcs = mypcair$vectors
  pcs = data.frame(pcs)
  colnames(pcs) = paste0("pc", c(1:ncol(pcs)))
  pcs$ID = row.names(pcs)

  covar_pcs <- inner_join(covar, pcs, by=c("ID"))
  covar_pcs$Pop <- "STUDY"
  pcs3 <- covar_pcs %>% select(pc1, pc2, pc3, pc4, pc5, pc6, pc7, pc8, pc9, pc10, Pop)
  

### Logistic Regression on PCA results

cov.glm <- glm(dx_strict ~ age + sex + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10,family = 'binomial', data = covar_pcs) ; summary(cov.glm)

a <- summary(cov.glm)
b <- as.data.frame(a$coefficients)
fwrite(b,paste0(dir,"pc_log_regression.txt"),col.names = T,row.names = T,quote = F,sep = "\t")

a <- summary(cov.glm.apoe)
b <- as.data.frame(a$coefficients)
fwrite(b,paste0(dir,"pc_log_regression_apoe.txt"),col.names = T,row.names = T,quote = F,sep = "\t")

saige_pheno_file <- covar_pcs %>% select(id,sex,age,apoe4,dx,pc1,pc2,pc3,pc4,pc5,pc6,pc7,pc8,pc9,pc10)
fwrite(saige_pheno_file,paste0(dir,"saige_pheno_file.tsv"),col.names = T,row.names = F,quote = F,sep = "\t")


### Log regression done to determine PCs to use for GRM creation

#Step 2a---Related and unrelated sets will be defined in pcair step above, output objects from pcair are input for pcrelate, GRM is also created at pcrelate step.
#PCRealts calculation => PCs and GRM
pcrelate_pcairKING <- pcrelate(mygenoData, pcs = mypcair$vectors[,1:3],
                               scale = "overall", training.set = mypcair$unrels,
                               sample.include = NULL)

grm_matrix <- pcrelateToMatrix(pcrelate_pcairKING)


#-------------------------------------------------------------------------------
# STEP 3---GRM  and relatedness adjusted PC are saved to output RDS files to be used in null model in another *R.  Above steps do not need to be run over and over when you have multiple models, etc.
#-------------------------------------------------------------------------------

saveRDS(grm_matrix, file = paste0(dir,plink,"grm",".rds"))
saveRDS(pcs, file = paste0(dir,plink,"pcs",".rds"))
saveRDS(mypcair, file  = paste0(dir,plink,"pcair",".rds"))

pruned <- as.data.frame(pruned)
fwrite(pruned,paste0(dir,"pruned_variant_list.txt"))


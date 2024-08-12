##Required libraries
library(CARMA)
library(data.table)
library(tidyverse)
library(dplyr)
library(R.utils)

#---------------------------------
rm(list = ls())

#-------------------------------------------------------------------------------
# STEP 0---PREP
#-------------------------------------------------------------------------------

#Set working directory
setwd("")

#---------------------------------
### Import the required input files and edit the format
Z <- fread("Region1.Zscores.txt")   ### Z score for each variant in the locus
annot <- read.table("Region1.AnnotationMatrix.txt", header=T)   ### CADD scores for each variant in the locus
ld = read.csv("Region1.ld",sep="\t",header=F)   ###Ld matrix for the locus

ld.mat = matrix(as.vector(data.matrix(ld)), nrow=nrow(Z), ncol=nrow(Z))
 
z.list<-list()
ld.list<-list()
lambda.list<-list()
annot.list<-list()

z.list[[1]]<-Z$Z
lambda.list[[1]]<-1
ld.list[[1]] <- ld.mat
annot.list[[1]]<-as.matrix(cbind(1, annot))

#-------------------------------------------------------------------------------
# Run CARMA for each suggestive significant region
#-------------------------------------------------------------------------------
Region1.CARMA.results<-CARMA(z.list,ld.list,w.list=annot.list,lambda.list=lambda.list,outlier.switch=F,effect.size.prior="Cauchy")


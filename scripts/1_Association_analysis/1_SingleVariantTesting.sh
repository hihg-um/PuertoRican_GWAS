#!/bin/bash
# SPDX-License-Identifier: GPL-2.0

if [ "$#" -ne 6 ]
then
      echo "usage: $0 <Plink> <VCF> <VCF_index> <Covar_File> <OUT_Prefix> <OUT_DIR>"
      exit 1
fi

PLINK=$1
VCF=$2
INDEX=$3
COVAR=$4
PREFIX=$5
OUT_DIR=$6

THREADS=16

mkdir -p $OUT_DIR
cd $OUT_DIR

#------------------------------------------------------------------------------------------------
# Step 0: Pruning
#------------------------------------------------------------------------------------------------

plink --bfile $PLINK --maf 0.05 --indep-pairwise 50 5 0.2 
plink --bfile $PLINK --extract plink.prune.in --make-bed --out $PLINK.pruned 

#------------------------------------------------------------------------------------------------
# Step 1: fitting the null logistic mixed model - MODEL 1 (adjusted for AGE, SEX, and PCs)
#------------------------------------------------------------------------------------------------

saige step1_fitNULLGLMM.R     \
    --plinkFile=$PLINK.pruned  \
    --phenoFile=$COVAR \
    --phenoCol=dx \
    --covarColList=age,sex,pc1,pc2,pc3,pc4 \
	--qCovarColList=sex \
    --sampleIDColinphenoFile=id \
    --traitType=binary        \
    --outputPrefix=$PREFIX.MODEL1_NULLGLMMR_LOCO \
    --nThreads=$THREADS \
	--LOCO=TRUE \
    --IsOverwriteVarianceRatioFile TRUE


#------------------------------------------------------------------------------------------------
#Step 2: performing single-variant association tests - MODEL 1
#------------------------------------------------------------------------------------------------

for CHR in {1..22}; do \
saige step2_SPAtests.R \
	--vcfFile=$VCF \
	--vcfFileIndex=$INDEX \
	--vcfField=GT \
	--chrom=chr$CHR \
	--minMAC=0.5 \
	--GMMATmodelFile=$PREFIX.MODEL1_NULLGLMMR_LOCO.rda \
	--varianceRatioFile=$PREFIX.MODEL1_NULLGLMMR_LOCO.varianceRatio.txt \
	--SAIGEOutputFile=${PREFIX}_chr$CHR.SAIGE.MODEL1.results.txt \
	--LOCO=TRUE \
	--is_output_moreDetails=TRUE \
	--is_Firth_beta=TRUE \
	--pCutoffforFirth=0.05 \
; done



#--------------------------------------------------------------------------------------------------------
# Step 1: fitting the null logistic mixed model - MODEL 2 (adjusted for APOE4 dosage, AGE, SEX, and PCs)
#--------------------------------------------------------------------------------------------------------

saige step1_fitNULLGLMM.R     \
    --plinkFile=$PLINK.pruned  \
    --phenoFile=$COVAR \
    --phenoCol=dx \
    --covarColList=age,sex,pc1,pc2,pc3,pc4,apoe4 \
	--qCovarColList=sex \
    --sampleIDColinphenoFile=id \
    --traitType=binary        \
    --outputPrefix=$PREFIX.MODEL2_NULLGLMMR_LOCO \
    --nThreads=$THREADS \
	--LOCO=TRUE \
    --IsOverwriteVarianceRatioFile TRUE


#--------------------------------------------------------------------------------------------------------
#Step 2: performing single-variant association tests - MODEL 2
#--------------------------------------------------------------------------------------------------------


for CHR in {1..22}; do \
saige step2_SPAtests.R \
	--vcfFile=$VCF \
	--vcfFileIndex=$INDEX \
	--vcfField=GT \
	--chrom=chr$CHR \
	--minMAC=0.5 \
	--GMMATmodelFile=$PREFIX.MODEL2_NULLGLMMR_LOCO.rda \
	--varianceRatioFile=$PREFIX.MODEL2_NULLGLMMR_LOCO.varianceRatio.txt \
	--SAIGEOutputFile=${PREFIX}_chr$CHR.SAIGE.MODEL2.results.txt \
	--LOCO=TRUE \
	--is_output_moreDetails=TRUE \
	--is_Firth_beta=TRUE \
	--pCutoffforFirth=0.05 \
; done

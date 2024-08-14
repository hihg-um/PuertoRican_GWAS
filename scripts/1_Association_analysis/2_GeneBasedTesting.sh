#!/bin/bash

if [ "$#" -ne 7 ]
then
      echo "usage: $0 <Plink_Pruned> <VCF> <VCF_index> <AnnotationFile_for_groupTest> <Covar_File> <OUT_Prefix> <OUT_DIR>"
      exit 1
fi

PLINK=$1
VCF=$2
INDEX=$3
ANNO=$4
COVAR=$5
PREFIX=$6
OUT_DIR=$7

THREADS=16

mkdir -p $OUT_DIR
cd $OUT_DIR

#------------------------------------------------------------------------------------------------
# Step 0: create a sparse GRM
#------------------------------------------------------------------------------------------------

saige createSparseGRM.R       \
	--plinkFile=$PLINK \
	--nThreads=$THREADS  \
	--outputPrefix=sparseGRM_extract    \
	--numRandomMarkerforSparseKin=2000      \
	--relatednessCutoff=0.125


#------------------------------------------------------------------------------------------------
# Step 1: fitting the null logistic mixed model - MODEL 1 (adjusted for AGE, SEX, and PCs)
#------------------------------------------------------------------------------------------------

saige step1_fitNULLGLMM.R    \
    --plinkFile=$PLINK \
	--sparseGRMFile=sparseGRM_extract_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx  \
	--sparseGRMSampleIDFile=sparseGRM_extract_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt \
	--phenoFile=$COVAR \
	--phenoCol=dx \
	--covarColList=age,sex,pc1,pc2,pc3,pc4 \
	--qCovarColList=sex  \
	--sampleIDColinphenoFile=id \
	--traitType=binary \
	--outputPrefix=$PREFIX.MODEL1_SAIGEGENE \
	--outputPrefix_varRatio=$PREFIX.MODEL1_SAIGEGENE \
	--useSparseGRMtoFitNULL=FALSE \
	--useSparseGRMforVarRatio=TRUE \
	--isCateVarianceRatio=TRUE \
	--IsOverwriteVarianceRatioFile=TRUE \
	--nThreads=$THREADS

#------------------------------------------------------------------------------------------------
# Step 2: performing the gene-based association tests - MODEL 1
#------------------------------------------------------------------------------------------------

for CHR in {1..22}; do \
saige step2_SPAtests.R \
	--vcfFile=$VCF \
	--vcfFileIndex=$INDEX \
	--vcfField=GT \
	--SAIGEOutputFile=${PREFIX}_chr$CHR.SAIGEGENE.MODEL1.results.txt \
	--chrom=chr$CHR \
	--LOCO=TRUE \
	--minMAF=0 \
	--minMAC=0.5 \
	--maxMissing=0.80 \
	--GMMATmodelFile=$PREFIX.MODEL1_SAIGEGENE.rda \
	--varianceRatioFile=$PREFIX.MODEL1_SAIGEGENE.varianceRatio.txt \
	--sparseGRMFile=sparseGRM_extract_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx \
	--sparseGRMSampleIDFile=sparseGRM_extract_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt \
	--groupFile=$ANNO \
	--annotation_in_groupTest=cadd_0-10:cadd_10-20:cadd_20+,cadd_10-20:cadd_20+,cadd_20+ \
	--maxMAF_in_groupTest=0.01 \
	--is_output_markerList_in_groupTest=TRUE \
	--is_output_moreDetails=TRUE 	\
	--is_Firth_beta=TRUE \
	--pCutoffforFirth=0.05 \
; done


#--------------------------------------------------------------------------------------------------------
# Step 1: fitting the null logistic mixed model - MODEL 2 (adjusted for APOE4 dosage, AGE, SEX, and PCs)
#--------------------------------------------------------------------------------------------------------

saige step1_fitNULLGLMM.R    \
    --plinkFile=$PLINK \
	--sparseGRMFile=sparseGRM_extract_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx  \
	--sparseGRMSampleIDFile=sparseGRM_extract_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt \
	--phenoFile=$COVAR \
	--phenoCol=dx \
	--covarColList=age,sex,pc1,pc2,pc3,pc4,apoe4 \
	--qCovarColList=sex  \
	--sampleIDColinphenoFile=id \
	--traitType=binary \
	--outputPrefix=$PREFIX.MODEL2_SAIGEGENE \
	--outputPrefix_varRatio=$PREFIX.MODEL2_SAIGEGENE \
	--useSparseGRMtoFitNULL=FALSE \
	--useSparseGRMforVarRatio=TRUE \
	--isCateVarianceRatio=TRUE \
	--IsOverwriteVarianceRatioFile=TRUE \
	--nThreads=$THREADS


#------------------------------------------------------------------------------------------------
# Step 2: performing the gene-based association tests - MODEL 2
#------------------------------------------------------------------------------------------------

for CHR in {1..22}; do \
saige step2_SPAtests.R \
	--vcfFile=$VCF \
	--vcfFileIndex=$INDEX \
	--vcfField=GT \
	--SAIGEOutputFile=${PREFIX}_chr$CHR.SAIGEGENE.MODEL2.results.txt \
	--chrom=chr$CHR \
	--LOCO=TRUE \
	--minMAF=0 \
	--minMAC=0.5 \
	--maxMissing=0.80 \
	--GMMATmodelFile=$PREFIX.MODEL2_SAIGEGENE.rda \
	--varianceRatioFile=$PREFIX.MODEL2_SAIGEGENE.varianceRatio.txt \
	--sparseGRMFile=sparseGRM_extract_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx \
	--sparseGRMSampleIDFile=sparseGRM_extract_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt \
	--groupFile=$ANNO \
	--annotation_in_groupTest=cadd_0-10:cadd_10-20:cadd_20+,cadd_10-20:cadd_20+,cadd_20+ \
	--maxMAF_in_groupTest=0.01 \
	--is_output_markerList_in_groupTest=TRUE \
	--is_output_moreDetails=TRUE 	\
	--is_Firth_beta=TRUE \
	--pCutoffforFirth=0.05 \
; done

#!/bin/bash

if [ "$#" -ne 3 ]
then
      echo "usage: $0 <Plink> <Reference_Plink>  <OUT_Prefix> <OUT_DIR> <variant_overlaps_R>"
      exit 1
fi

PLINK=$1
REF_PLINK=$2
PREFIX=$3
OUT_DIR=$4
variant_overlaps_R=$5

THREADS=16

mkdir -p $OUT_DIR
cd $OUT_DIR


#------------------------------------------------------------------------------------------------
# Define overlappint variants between target and reference datasets
#------------------------------------------------------------------------------------------------

Rscript $variant_overlaps_R $OUT_DIR ${PLINK}.bim ${REF_PLINK}.bim

#-----------------------------------------------------------------------------------------------
# Keep overlapping variants and merge datasets
#-----------------------------------------------------------------------------------------------

plink --bfile $PLINK --extract keep_rsid.txt --make-bed --out file_ov
plink --bfile $REF_PLINK --extract keep_rsid.txt --make-bed --out ref_ov

plink --bfile file_ov --bmerge ref_ov --make-bed --out $PREFIX
plink --bfile ref_ov --flip ${PREFIX}-merge.missnp --make-bed --out flip
plink --bfile file_ov --bmerge flip --make-bed --out $PREFIX

#-----------------------------------------------------------------------------------------------
# Keep overlapping variants and merge datasets
#-----------------------------------------------------------------------------------------------

plink2 --bfile ${PREFIX} \
	--set-all-var-ids "chr"@:# \
	--rm-dup exclude-all \
	--snps-only just-acgt \
	--geno 0.03 \
	--hwe 1e-6 \
	--maf 0.05 \
	--make-bed --out ${PREFIX}_rmdup_bi_snp_hwe_maf_filters

#-----------------------------------------------------------------------------------------------
# Calculate ROHs
#-----------------------------------------------------------------------------------------------

plink --bfile ${PREFIX}_rmdup_bi_snp_hwe_maf_filters \
	--homozyg-snp 50 \
	--homozyg-kb 300 \
	--homozyg-density 50 \
	--homozyg-gap 1000 \
	--homozyg-window-snp 50 \
	--homozyg-window-het 1 \
	--homozyg-window-missing 1 \
	--homozyg-window-threshold 0.05 \
	--out ${PREFIX}_rmdup_bi_snp_hwe_maf_filters_5hetperwin \

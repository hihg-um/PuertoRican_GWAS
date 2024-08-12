#!/bin/bash

if [ "$#" -ne 3 ]
then
      echo "usage: $0 <Target_VCF> <Base_Reference> <Pheno_for_Plink> <Covar_File> <APOE_range_file> <OUT_Prefix> <OUT_DIR>"
      exit 1
fi

VCF=$1
BASE=$2
PHENO=$3
COVAR=$4
APOEREG=$5
PREFIX=$6
OUT_DIR=$7

THREADS=16

mkdir -p $OUT_DIR
cd $OUT_DIR

#-------------------------------------------------------------------------------
# convert VCF to PLINK
#-------------------------------------------------------------------------------
plink2 \
	--vcf ${VCF} \
	--set-all-var-ids 'chr'@:# \
	--pheno ${PHENO} --mpheno 2 \
	--maf 0.01 \
	--geno 0.01 \
	--hwe 1e-6 \
	--make-bed \
	--out ${PREFIX}


#-------------------------------------------------------------------------------
# run PRSice
#-------------------------------------------------------------------------------

PRSice \
	--a1 A1 \
	--a2 A2 \
	--base ${BASE} \
	--beta \
	--binary-target T \
	--bp BP \
	--chr CHR \
	--clump-kb 250kb \
	--clump-p 1.000000 \
	--clump-r2 0.100000 \
	--cov ${COVAR} \
	--cov-col age,sex,pc1,pc2,pc3,pc4 \
	--cov-factor sex \
	--ignore-fid \
	--interval 5e-05 \
	--lower 5e-08 \
	--maf 0.01 \
	--num-auto 22 \
	--out ${PREFIX} \
	--print-snp \
	--proxy 0.8 \
	--pvalue P \
	--seed 809 \
	--snp SNP \
	--stat BETA \
	--target ${PREFIX} \
	--thread ${THREADS} \
	--upper 0.5 \
	--x-range ${APOEREG}

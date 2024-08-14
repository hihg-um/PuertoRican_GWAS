#!/bin/bash

if [ "$#" -ne 4 ]
then
      echo "usage: $0 <VCF> <REGION> <OUT_Prefix> <OUT_DIR>"
      exit 1
fi

VCF=$1
REGION=$2
PREFIX=$3
OUT_DIR=$4

mkdir -p $OUT_DIR
cd $OUT_DIR


#------------------------------------------------------------------------------------------------
# VCF to LD matrix for each suggestive significant region
#------------------------------------------------------------------------------------------------

plink2 --vcf  $VCF --make-bed --keep-allele-order --set-all-var-ids chr@:# --extract $REGION --out $PREFIX

plink --bfile $REGION --a1-allele $REGION 2 1 '#' --r square --out $PREFIX

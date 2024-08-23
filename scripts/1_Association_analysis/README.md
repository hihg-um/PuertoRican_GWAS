---

## Description
### Input files
- WGS VCF and index files
- Plink files for the same WGS dataset (bed, bim fam)
- Covariates (this file should include diagnosis, age, sex and apoe4 dosages)
- Annotation file for groupTest (for gene-based testing)

---

### Execution

1. Principal Components Analysis (PCA)  [scripts/1_Association_analysis/0_PC_Calculation.R](0_PC_Calculation.R)
2. Single variant testing [scripts/1_Association_analysis/1_SingleVariantTesting.sh](1_SingleVariantTesting.sh)
	- The "saige_pheno_file.tsv" file generated in the PCA step can be used as a covar file.
3. Gene-based testing [scripts/1_Association_analysis/2_GeneBasedTesting.sh](2_GeneBasedTesting.sh)
	- Pruned plink files generated in the single variant testing step can be used as an input plink file.
	- The "saige_pheno_file.tsv" file generated in the PCA step can be used as a covar file.
	- An annotation file containing a list of variants at different CADD thresholds (CADD0-10, CADD10-20, CADD_20+) ​​is required for each set to be tested.

---

## Dependencies
### Tools
- PLINK 1.9 and PLINK 2.0 (both versions are required)
- SAIGE
### R packages
- GENESIS
- GWASTools
- SNPRelate
- gdsfmt
- tidyverse
- table

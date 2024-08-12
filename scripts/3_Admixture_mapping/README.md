---

## Description
### Input files
- Genetic Relationship Matrix (GRM)
- Covariates (this file should include diagnosis, age, apoe4 dosages, and principal components)
- MSP files (RFmix2 output)

---

### Execution

1. Admixture mapping  [scripts/3_Admixture_mapping/Admixture_Mapping.sh](Admixture_Mapping.R)
	- GRM RDS files generated in the PCA step can be used as GRM input files.
	- The "saige_pheno_file.tsv" file generated in the PCA step can be used as a covar file.
	- The script should be revised according to the ancestry order in the MSP files.

---

## Dependencies
### R packages
- GENESIS
- GWASTools
- GWASdata
- gdsfmt
- dplyr
- data.table
- tidyr
- tidyverse

---

## Description
### Input files
- WGS VCF files (target dataset)
- Effect sizes from summary statistics from the reference study (base dataset)
- Phenotypes (this file should include diagnosis)
- Covariates (this file should include age, apoe4 dosages, and principal components)
- APOE region file (Chromosomal location required to exclude APOE region from PRS)

---

### Execution

1. Polygenic Risk Score calculation  [scripts/5_PRS/PRS.sh](PRS.sh)
	- APOE region file format: CHR START END (eg.19 43907789 45908821)
	- The "saige_pheno_file.tsv" file generated in the PCA step can be used as a covar file.

---

## Dependencies
### Tools
- PLINK 2.0
- PRSice
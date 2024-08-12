---

## Description
### Input files
- WGS VCF files
- Region file of suggested suggestive significant loci
- Z-scores file for suggested suggestive significant loci
- Annotation matrix

---

### Execution

1. LD matrix generation for suggestive significant loci  [scripts/2_Fine-mapping/1_VCFtoLDmatrix.sh](1_VCFtoLDmatrix.sh)
	- Region file format: 1st column = variant id & 2nd column = reference allele
2. Fine mapping [scripts/2_Fine-mapping/2_CARMA.R](2_CARMA.R)
	- LD matrix files generated in the single variant testing step can be used as an input ld files.
	- Z-scores file format: 1st column = variant id & 2nd column = Z-values (The Z score can be calculated by dividing the BETA by the StandartError values ​​in the summary statistics file obtained as a result of the single variant testing.)
	- Annotation matrix file format: 1st column = 0 & 2nd column = CADD PHRED values
	- Since each line in the input files will represent a variant, their order must be the same.

---

## Dependencies
### Tools
- PLINK 1.9
- PLINK 2.0
### R packages
- CARMA
- data.table
- R.utils
- tidyverse
- dplyr
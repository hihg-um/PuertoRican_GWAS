---

## Description
### Input files
- Plink files for the WGS dataset (bed, bim fam)
- Plink files for the reference population dataset (combined reference panels of the Human Genome Diversity Project (HGDP) and 1000 Genomes Phase 3)
- variant_overlaps.R (R script to find overlapping variants between two datasets)

---

### Execution

1. Runs of homozygosity (ROH) analysis  [scripts/4_ROH/ROH.sh](ROH.sh)
	- The variant_overlaps.R script should be given as input. [scripts/4_ROH/variant_overlaps.R](variant_overlaps.R)

---

## Dependencies
### Tools
- PLINK 1.9 and PLINK 2.0 (both versions are required)
### R packages
- data.table
- tidyverse
- gtools

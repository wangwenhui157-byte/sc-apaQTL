# sc-apaQTL

## Description
Alternative polyadenylation (APA) is a pervasive post-transcriptional mechanism regulating gene expression and disease susceptibility, yet its genetic architecture across cell types and developmental states remains incompletely defined. To address this, we provides a high-resolution, cell-type-resolved atlas of APA genetic regulation in human immune cells and establishes APA as a critical molecular layer linking genetic variation to immune function and autoimmune disease susceptibility.

## 01.genotype.sh
The shell script `01.genotype.sh` was used to perform quality control and preprocessing of genotype data.

## 02.PASTA.sh
The shell script `02.PASTA.sh` was used to quantify alternative polyadenylation (APA) events and generate polyA site usage matrices.

## 03.apaQTL_calculate.R
The R script `03.apaQTL_calculate.R` was used to identify apaQTLs for each cell type. The R package Matrix eQTL was utlilized to test for association between genotype and ployA residual.

## 04.apaQTL_MASH.R
The R script `04.apaQTL_MASH.R` utilized the mash method to elucidate the heterogeneity of apaQTL effect sizes across different cell types.

## 05.Dynamic_PAS_identify.R
The R script `05.Dynamic_PAS_identify.R` was used to identify polyA sites with dynamic APA usage along the B cell maturation pseudotime.

## 06.Dynamic_apaQTL_identify.R
The R script `06.Dynamic_apaQTL_identify.R` was used to identify apaQTLs with pseudotime-dependent effects, revealing dynamic genetic regulation of APA during B cell maturation.

## 07.LDSC.sh
The shell script `07.LDSC.sh` was used to compute partitioned heritability.

## 08.COLOC.R
The R script `08.COLOC.R` was used for colocalization analysis of apaQTLs with GWAS associated autoimmune diseases.

## 09.SMR.sh
The shell script `09.SMR.sh` was used to perform summary-data-based Mendelian randomization (SMR) analysis.

## 10.MOLOC.R
The R script `10.MOLOC.R` was used to perform multi-trait colocalization analysis integrating apaQTL, pQTL, and GWAS data.

## 11.Mediation.R
The R script `11.Mediation.R` was used to perform mediation analysis integrating apaQTL, pQTL, and GWAS data to evaluate potential causal pathways linking genetic variants, APA regulation, protein abundance, and disease risk.

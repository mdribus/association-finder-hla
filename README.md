# association-finder-hla-pretend
This repository contains scripts to perform HLA disease association studies using a disease association pipeline developed in R called AssociationFinderHLA. AssociationFinderHLA was designed to run on the command line. It can identify allele, haplotype, and genotype associations. 
## Scripts in AssociationFinderHLA:
### 1_find_associations_final.R:
This script identifies the associations and calculates summary statistics for each analysis category. It produces output in the /Multi and /Summary subdirectories.
### 2_fdr_correction_final.R: 
Since multiple comparisons are being made in each analysis category, correction for false discovery rate is required. This script produces output in the /FDR and /Summary subdirectories.
### 3_factor_analysis_final.R: 
Factor analysis is performed to organize associated HLA variants into groups based on the underlying structure in the data. The output goes to the /Multi subdirectory.
### 4_master_tables_final.R: 
This script organizes the associations and factor analysis group assignments into tables for each population. It produces output in the /Masters subdirectory.
### 5_final_tables_final.R: 
This script writes an Excel file with the association results for each population sorted by p-value and and factor analysis group assignment. It produces output in the /Masters subdirectory. 
## Other required files:
### HLA typing files:
These files contain the IDs, disease indicators, and HLA typing for the cases and controls. They should be named with the prefix "HLA_input_" followed by the disease, population, imputation replicate number, and the file extension ".txt". This is how the a file should be named if it contains HLA typing for cases with severe aplastic anemia (SAA) in the Asian & Pacific Islander population and is the first imputation replicate: 
**HLA_input_SAA_API_1.txt**

In the HLA typing files, the IDs occupy the first column, **ID**, and the disease indicators occupy the second column, **Disease_Ind**. A value of "1" in the Disease_Ind column identifies a case, while a value of "0" identifies a control. The other columns in the file contain HLA typing for individual alleles, haplotypes, and genotypes. This is how the input data should look for a case with the HLA typing *A\*11:01\~C\*12:03\~B\*38:01\+A\*24:02\~C\*08:02\~B\*14:02* and a control with the HLA typing *A\*01:01\~C\*07:01\~B\*08:01\+A\*01:01\~C\*07:01\~B\*08:01*:
|    ID    | Disease_Ind |  A1  |  C1  |  B1  |  A2  |  C2  |  B2  | Hap_A1_C1 | Hap_A2_C2 |    Geno_A_C_B    |
| -------- | ----------- | ---- | ---- | ---- | ---- | ---- | ---- | --------- | --------- | ---------------- |
|24537890|1|A\*11:01|C\*12:03|B\*38:01|A\*24:02|C\*08:02|B\*14:02|A\*11:01\~C\*12:03|A\*24:02\~C\*08:02|A\*11:01\~C\*12:03\~B\*38:01\+A\*24:02\~C\*08:02\~B\*14:02|
|43654356|0|A\*01:01|C\*07:01|B\*08:01|A\*01:01|C\*07:01|B\*08:01|A\*01:01~C\*07:01|A\*01:01~C\*07:01|A\*01:01\~C\*07:01\~B\*08:01\+A\*01:01\~C\*07:01\~B\*08:01|

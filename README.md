# rg-matix
Calculate pairwise genetic correlations between a set of GWAS summary statistics  
Version 2.0.0  

Input: "munged" GWAS summary statistics, see https://github.com/bulik/ldsc  
NOTE: script further requires Python version >=2.7 & <3; preferably 2.7.9  

Output:  
3x tab-delmited text files; N rows by (N+1) columns:  
  - rg.tab:	symmetric matrix of genetic correlations  
  - se.tab:	symmetric matrix of SE's for the genetic correlations  
  - cti.tab:	symmetric matrix of cross-trait-intercepts  

1x tab-delimited text file containing all estimated values reported by LDSC; (N+1) rows (including header) by 12 columns  

How to use:  
 1) copy this script and all sumstats to the working directory  
 2) run the script using the following command:  
 ```
 bash script.make_rg_mat.sh
 ```  

NOTE: the script assumes no files ending on ".sumstats.gz" are present except for the GWASs for which the genetic correlations are to be calculated

Additional flags:  
Restrict the diagonal of both rg.tab and cti.tab to all 1's, and restrict the diagonal of se.tab to 0's
```
-c
```
Save LDSC log-files
```
-s
```
Add prefix to output files
```
-n [value]
```
Change the number of parallel runs of LDSC to start; default 5
```
-p [value]
```
Supply the location of LDSC.py (see https://github.com/bulik/ldsc)
```
-l [value]
```
Indicate which files to use as reference LD scores (see https://github.com/bulik/ldsc)
```
-r [value]
```
Indicate which files to use as regression weights (see https://github.com/bulik/ldsc)
```
-w [value]
```

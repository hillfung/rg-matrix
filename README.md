# rg-matrix
Create a matrix of genetic correlations between a set of GWAS summary statistics  
Current version: 2.2.0  

## Input: ##
  - "munged" GWAS summary statistics, see https://github.com/bulik/ldsc  
  - Python version >=2.7 & <3; preferably 2.7.9  

## Output: ##
- 4x tab-delmited text files; N rows by (N+1) columns:  
  - rg.tab:  symmetric matrix of genetic correlations  
  - rg_se.tab:  symmetric matrix of SE's for the genetic correlations  
  - rg_p.tab: symmetric matrix of P-values for the genetic correlations
  - cti.tab: symmetric matrix of cross-trait-intercepts  
  - cti_se.tab: symmetric matrix of SE's for the cross-trait-intercepts
- 1x tab-delimited text file containing all estimated values reported by LDSC; (N+1) rows (including header) by 12 columns  

## How to use: ##
 Run the script using the following command: ```bash make_rg_mat.sh```

The script requires at least the three following flags:

Supply the location of LDSC.py (see https://github.com/bulik/ldsc)  
`-l [value]`

Indicate which files to use as reference LD scores (see https://github.com/bulik/ldsc)  
`-r [value]`

Indicate which files to use as regression weights (see https://github.com/bulik/ldsc)  
`-w [value]`

## Additional flags: ##
Restrict the diagonal of both rg.tab and cti.tab to all 1's, and restrict the diagonal of rg_se.tab and cti_se.tab to 0's  
`-c`

Save LDSC log-files  
`-s`

Add prefix to output files  
`-n [value]`

Change the number of parallel runs of LDSC to start; default 5  
`-p [value]`

Input directory; defaults to working directory  
`-i [value]`


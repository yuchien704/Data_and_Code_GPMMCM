# GPMMCM

## Supplementary Material for "Flexible clustering via Gaussian parsimonious mixture models with censored and missing values" by Wan-Lun Wang, Victor H. Lachos, Yu-Chien Chen, and Tsung-I Lin

### Author responsible for the code
For questions, comments or remarks about the code please contact responsible author, Tsung-I Lin (tilin@nchu.edu.tw).

### Configurations
The code was written/evaluated in R with the following software versions:
R version 4.2.1 (2022-06-23 ucrt)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 22621)

Matrix products: default

locale:
[1] LC_COLLATE=Chinese (Traditional)_Taiwan.utf8  LC_CTYPE=Chinese (Traditional)_Taiwan.utf8   
[3] LC_MONETARY=Chinese (Traditional)_Taiwan.utf8 LC_NUMERIC=C                                 
[5] LC_TIME=Chinese (Traditional)_Taiwan.utf8   

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods   base  

other attached packages:
 [1] RColorBrewer_1.1-3 mixture_2.0.5      lattice_0.20-45    mixtools_2.0.0     MomTrunc_6.0       mclust_6.0.0      
 [7] tmvtnorm_1.5       gmm_1.7            sandwich_3.0-2     Matrix_1.5-3       mvtnorm_1.1-3     

loaded via a namespace (and not attached):
 [1] tlrmvnmvt_1.1.2   deSolve_1.34      zoo_1.8-11        tidyselect_1.2.0  kernlab_0.9-31    purrr_0.3.5      
 [7] splines_4.2.1     colorspace_2.0-3  vctrs_0.5.0       generics_0.1.3    htmltools_0.5.3   viridisLite_0.4.1
[13] survival_3.3-1    utf8_1.2.2        plotly_4.10.1     rlang_1.0.6       pillar_1.8.1      glue_1.6.2       
[19] DBI_1.1.3         segmented_1.6-2   lifecycle_1.0.3   munsell_0.5.0     gtable_0.3.1      elliptic_1.4-0   
[25] contfrac_1.1-12   htmlwidgets_1.5.4 fastmap_1.1.0     fansi_1.0.3       Rcpp_1.0.10       scales_1.2.1     
[31] hypergeo_1.2-13   jsonlite_1.8.3    ggplot2_3.4.0     digest_0.6.30     dplyr_1.0.10      grid_4.2.1       
[37] cli_3.4.0         tools_4.2.1       magrittr_2.0.3    lazyeval_0.2.2    tibble_3.1.8      tidyr_1.2.1      
[43] pkgconfig_2.0.3   MASS_7.3-57       data.table_1.14.4 assertthat_0.2.1  httr_1.4.4        rstudioapi_0.14  
[49] R6_2.5.1          nlme_3.1-157      compiler_4.2.1   


### Descriptions of the codes 
Please copy the files to the "current working directory" of the R package.
The 'getwd()' function shall determine an absolute pathname of the "current working directory".

Before running all of the codes, one needs to install the following R packages:

    install.packages("mvtnorm")  Version: 1.1-3
    install.packages("tmvtnorm") Version: 1.5
    install.packages("mclust") Version： 6.0.0
    install.packages("MomTrunc")   Version： 6.0
    install.packages("mixtools") Version： 2.0.0
    install.packages("RColorBrewer")   Version： 1.1-3

R codes for the implementation of our methodology are provided.

#### Subfolder: ./Function ####
'./Function/'
       contains 
       
       (1) the subfolder 'Both_gmixcm', in this folder contains '.R' file code for 14 different structures that can deal with the data including both interval-censored values and missing values;
       (2) 'gen_cen_na_new.R' main script for generating the data with both interval-censored values and missing values, and the matrix "cen" indicates the corresponding censored and missing positions; 
       (3) 'F-G.r' main script for solving closed-form solutions of matrix D in the "EVE" and "VVE" structure by the FG algorithm; and
       (4) 'MM1_alg.R' main script for using MM algorithm to solve closed-form solutions of matrix D in the "EVE" and "VVE" structure.

'./Function/Both_gmixcm'
	subfolder collects functions for maximum likelihood (ML) estimation for 14 Gaussian parsimonious mixture models with censored and missing data (GPMM-CM), including
	
    	(1) 'EEE_GMIXCMB.R' that can perform the ECM algorithm for fitting the GPMM-CM with EEE structure; 
    	(2) 'EEV_GMIXCMB.R' that can perform the ECM algorithm for fitting the GPMM-CM with EEV structure; 
    	(3) 'EVE_GMIXCMB.R' that can perform the ECM algorithm for fitting the GPMM-CM with EVE structure; 
    	(4) 'VEE_GMIXCMB.R' that can perform the ECM algorithm for fitting the GPMM-CM with VEE structure; 
    	(5) 'VEV_GMIXCMB.R' that can perform the ECM algorithm for fitting the GPMM-CM with VEV structure; 
    	(6) 'VVE_GMIXCMB.R' that can perform the ECM algorithm for fitting the GPMM-CM with VVE structure; 
    	(7) 'EVV_GMIXCMB.R' that can perform the ECM algorithm for fitting the GPMM-CM with EVV structure; 
    	(8) 'VVV_GMIXCMB.R' that can perform the ECM algorithm for fitting the GPMM-CM with VVV structure; 
    	(9) 'VVI_GMIXCMB.R' that can perform the ECM algorithm for fitting the GPMM-CM with VVI structure; 
    	(10) 'EVI_GMIXCMB.R' that can perform the ECM algorithm for fitting the GPMM-CM with EVI structure; 
    	(11) 'VEI_GMIXCMB.R' that can perform the ECM algorithm for fitting the GPMM-CM with VEI structure; 
    	(12) 'EEI_GMIXCMB.R' that can perform the ECM algorithm for fitting the GPMM-CM with EEI structure; 
    	(13) 'VII_GMIXCMB.R' that can perform the ECM algorithm for fitting the GPMM-CM with VII structure; 
    	(14) 'EII_GMIXCMB.R' that can perform the ECM algorithm for fitting the GPMM-CM with EII structure; 
 
   	(15) 'EVE_GMIXCMB_MM.R' that can perform the MM algorithm for fitting the GPMM-CM with EVE structure; 
    	(16) 'VVE_GMIXCMB_MM.R' that can perform the MM algorithm for fitting the GPMM-CM with EVE structure;
    	
###### Note for 14 structures function ######
When fitting the data, sometimes the function cannot work because one of the subgroups has too many missing values. 
We can utilize the 'try()' in R to overcome the situation. 
Therefore, we need other functions to implement the procedure. 

    	(17) 'VVV_GMIXCMB_1.R' that carries out ECM-based ML estimation coupled with 'try()' for re-generating initial clustering to fit the GPMM-CM with VVV structure;
    	(18) 'VVE_GMIXCMB_1.R' that carries out ECM-based ML estimation coupled with 'try()' for re-generating initial clustering to fit the GPMM-CM with VVE structure;
    	(19) 'VEV_GMIXCMB_1.R' that carries out ECM-based ML estimation coupled with 'try()' for re-generating initial clustering to fit the GPMM-CM with VEV structure;
    	(20) 'EVV_GMIXCMB_1.R' that carries out ECM-based ML estimation coupled with 'try()' for re-generating initial clustering to fit the GPMM-CM with EVV structure;
    	(21) 'EVE_GMIXCMB_1.R' that carries out ECM-based ML estimation coupled with 'try()' for re-generating initial clustering to fit the GPMM-CM with EVE structure;
    	(22) 'EEV_GMIXCMB_1.R' that carries out ECM-based ML estimation coupled with 'try()' for re-generating initial clustering to fit the GPMM-CM with EEV structure;
    	(23) 'EEE_GMIXCMB_1.R' that carries out ECM-based ML estimation coupled with 'try()' for re-generating initial clustering to fit the GPMM-CM with EEE structure; and
    	(24) 'EEI_GMIXCMB_1.R' that carries out ECM-based ML estimation coupled with 'try()' for re-generating initial clustering to fit the GPMM-CM with EEI structure.


#### Subfolder: ./Code ####
'./Code/'
       contains 
       
       (1) 'fig1.R' main script for reproducing drawing Figure 1 that shows 14 GPMMs with three groups in two dimensions;
       (2) 'fig2.R' main script for reproducing Figure 2; 
       (3) 'fig3.4.R' main script for reproducing Figure 3 and Figure 4;
       (4) 'fig5.R' main script for reproducing Figure 5;
       (5) 'fig6.R' main script for reproducing Figure 6; 
       (6) 'table2.3.R' main script for reproducing Tables 2 and 3;
       (7) 'table4.R' main script for Tables 4; 

###### Note for Section 5 - Illustrative examples - faithful data:
Because the 'faithful_code.R' and 'faithful_code_01.R' code takes a huge amount of time to run the ECM procedure for fitting the GPMM-CM model, we record these intermediate results in 'faithful.RData' and 'faithful_01.RData' so that one can use the R codes 'fig2.R' and 'fig3.4.R' to obtain the final results immediately.
To reproduce the results presented in Figure 2, just load 'faithful_01.RData' file in the './Data/' and then run the script 'fig2.R' in the subfolder './Code/';  
To reproduce the results presented in Figures 3 and 4, just load 'faithful.RData' file in the './Data/' and then run the script 'fig3.4.R' in the subfolder './Code/'

       (8) 'faithful_code.R' main script for fitting the GPMM-CM models with 14 parsimonious structures and g = 1-5 to the faithful data including both interval-censored values and missing values;
       (9) 'faithful_code_01.R' main script for fitting the GPMM-CM models with 14 parsimonious structures and g = 1-5 to the faithful data;

###### Note for Section 5 - Illustrative examples - hawks data:
Because the 'hawks_code.R' code takes a huge amount of time to run the ECM procedure for fitting the GPMM-CM model, we record these intermediate results in 'hawks.RData' so that one can use the R codes 'fig5.R', and 'table2.3.R' to obtain the final results immediately.
To reproduce the results presented in Figures 5 and Tables 2 and 3, just load 'hawks.RData' file in the './Data/' and then run the script 'fig5.R', and 'table2.3.R' in the subfolder './Code/';  

       (10) 'hawks_code.R' main script for fitting the GPMM-CM models with 14 parsimonious structures and g = 3 to the standardized hawks data including both interval-censored values and missing values;

###### Note for Section 5 - Illustrative examples - simulated data:
Because the 'sim_VVE_code_code.R' code takes a huge amount of time to run the ECM procedure for fitting the GPMM-CM model, we record every 20 replications result in 'CM_1.RData', 'CM_2.RData', 'CM_3.RData', 'CM_4.RData', and 'CM_5.RData' respectively so that one can use the R codes 'table7.R' to obtain the final results immediately.
To reproduce the results presented in Table 4, just load 'CM_1.RData', 'CM_2.RData', 'CM_3.RData', 'CM_4.RData', and 'CM_5.RData' files in the './Data/' and then run the script 'table4.R'in the subfolder './Code/';  

       (11) 'sim_VVE_code.R' main script for fitting the GPMM-CM models with 14 parsimonious structures and g = 1 to 5 to simulated datasets including both interval-censored values and missing values.

#### Subfolder: ./Data ####
'./Data/'
       contains 
       
       (1) 'CM_1.RData', 'CM_2.RData', 'CM_3.RData', 'CM_4.RData', and 'CM_5.RData' collecting the frequencies and average BIC values obtained by fitting GPMM-CM models with 14 parsimonious structures to the 20 simulated datasets;
       (2) 'faithful.RData' collecting the fitting results of the 14 GPMM-CM models to the old faithful data with interval-censored and missing values, where g ranging from 1 to 5;
       (3) 'faithful_01.RData' collecting the fitting results of the 14 GPMM-CM models to the old faithful data, where g ranging from 1 to 5;
       (4) 'hawks.RData' collecting the fitting and classification results of the 3-component GPMM-CM with 14 parsimonious structures to the hawks data with missing values and synthetic censored values; 
       (5) The subfolder 'source', which contains the hawks dataset used in Section 5.

'./Data/source'
	subfolder contains
	
    	(1) 'Hawks.csv' for the hawks dataset usued in Section 5.3; 
    
#### Subfolder: ./Results ####
'./Results/'
       contains 
       
       (1) 'fig1.eps': an illustration of elliptical contours and scatter plots for 14 GPMMs with three groups in two dimensions;
       (2) 'fig2.eps'; perspective plot and bivariate density contours obtained from the best-fitted (EEE, 3) model; 
       (3) 'fig3.eps': comparison of BIC values for 14 GPMM-CMs with g = 1 − 5 fitted to the old faithful data involving censored and missing values;
       (4) 'fig4.eps': the histograms and scatter plot overlaid with contour curves obtained from the best fitted (EVI, 3) model for the old faithful data involving censored and missing values;
       (5) 'fig5.eps': the BIC curve for 14 GPMM-CMs with g = 3 fitted to the hawks data; 
       (6) 'fig6.eps': the scatter-histogram plot of one simulated case with 450 random samples containing missing and censored values whose true values of parameters are specified in Section 5.3;

       (7) 'Table2.csv': an overview of 5 variables of the hawks data with three species;
       (8) 'Table3.csv': the results of the ARI and CCR of the 14 clustering models with fitting the hawks data into three clusters; and
       (9) 'Table4.csv': performance comparisons obtained from fitting GPMM-CM models with 14 parsimonious structures to the simulated dataset over 100 replications.

# Additional Remark 
One can directly run each "source(.)" described in 'master.r' file in the seperate R session to obtain the results.

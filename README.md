# GPMMCM

## Supplementary Material for "Flexible clustering via Gaussian parsimonious mixture models with censored and missing values" by Wan-Lun Wang, Victor H. Lachos, Yu-Chien Chen, and Tsung-I Lin

### Author responsible for the code
For questions, comments or remarks about the code please contact responsible author, Tsung-I Lin (tilin@nchu.edu.tw).

### Configurations
The code was written/evaluated in R with the following software versions:
R version 4.3.2 (2023-10-31 ucrt)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 11 x64 (build 22631)

Matrix products: default

locale:
[1] LC_COLLATE=Chinese (Traditional)_Taiwan.utf8  LC_CTYPE=Chinese (Traditional)_Taiwan.utf8   
[3] LC_MONETARY=Chinese (Traditional)_Taiwan.utf8 LC_NUMERIC=C                                 
[5] LC_TIME=Chinese (Traditional)_Taiwan.utf8   

attached base packages:
[1] stats4    grid      stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] mi_1.1             mice_3.16.0        Amelia_1.8.1       Rcpp_1.0.11        kpodclustr_1.1     Stat2Data_2.0.0    CensMFM_3.1       
 [8] RColorBrewer_1.1-3 mixtools_2.0.0     MomTrunc_6.0       mclust_6.0.1       tmvtnorm_1.5       gmm_1.8            sandwich_3.0-2    
[15] Matrix_1.6-1.1     mvtnorm_1.2-4      GGally_2.2.1       ggplot2_3.4.4      VIM_6.2.2          colorspace_2.1-0  

loaded via a namespace (and not attached):
 [1] gridExtra_2.3     rlang_1.1.2       magrittr_2.0.3    e1071_1.7-14      compiler_4.3.2    vctrs_0.6.4       pkgconfig_2.0.3  
 [8] shape_1.4.6       crayon_1.5.2      fastmap_1.1.1     arm_1.13-1        backports_1.4.1   labeling_0.4.3    utf8_1.2.4       
[15] deSolve_1.40      nloptr_2.0.3      purrr_1.0.2       glmnet_4.1-8      jomo_2.7-6        jsonlite_1.8.7    progress_1.2.3   
[22] pan_1.9           parallel_4.3.2    broom_1.0.5       prettyunits_1.2.0 R6_2.5.1          vcd_1.4-12        ranger_0.16.0    
[29] rpart_4.1.21      car_3.1-2         boot_1.3-28.1     lmtest_0.9-40     iterators_1.0.14  zoo_1.8-12        splines_4.3.2    
[36] nnet_7.3-19       tidyselect_1.2.0  abind_1.4-5       codetools_0.2-19  lattice_0.21-9    tibble_3.2.1      plyr_1.8.9       
[43] withr_2.5.2       coda_0.19-4       foreign_0.8-85    survival_3.5-7    ggstats_0.7.0     proxy_0.4-27      kernlab_0.9-32   
[50] pillar_1.9.0      carData_3.0-5     foreach_1.5.2     plotly_4.10.3     generics_0.1.3    sp_2.1-3          hms_1.1.3        
[57] elliptic_1.4-0    munsell_0.5.0     scales_1.2.1      minqa_1.2.6       laeken_0.5.3      class_7.3-22      glue_1.6.2       
[64] lazyeval_0.2.2    tools_4.3.2       robustbase_0.99-2 data.table_1.14.8 lme4_1.1-35.1     contfrac_1.1-12   hypergeo_1.2-13  
[71] tidyr_1.3.0       nlme_3.1-163      cli_3.6.1         fansi_1.0.5       segmented_1.6-4   viridisLite_0.4.2 dplyr_1.1.4      
[78] gtable_0.3.4      DEoptimR_1.1-3    digest_0.6.33     htmlwidgets_1.6.3 farver_2.1.1      htmltools_0.5.7   lifecycle_1.0.4  
[85] tlrmvnmvt_1.1.2   httr_1.4.7        mitml_0.4-5       MASS_7.3-60      

### Descriptions of the codes 
Please copy the files to the "current working directory" of the R package.
The 'getwd()' function shall determine an absolute pathname of the "current working directory".

Before running all of the codes, one needs to install the following R packages:

    install.packages("mvtnorm")  Version: 1.2-4
    install.packages("tmvtnorm") Version: 1.5
    install.packages("mclust") Version: 6.0.1
    install.packages("MomTrunc") Version: 6.0
    install.packages("mixtools") Version: 2.0.0
    install.packages("RColorBrewer") Version: 1.1-3
    install.packages("VIM") Version: 6.2.2
    install.packages("GGally") Version: 2.2.1
    install.packages("ggplot2") Version: 3.4.4
    install.packages("CensMFM") Version: 3.1
    install.packages("Stat2Data"): 2.0.0
    install.packages("kpodclustr"): 1.1
    install.packages("Amelia"): 1.8.1
    install.packages("mice"): 3.16.0
    install.packages("mi"): 1.1
    
R codes for the implementation of our methodology are provided.

#### Subfolder: ./Function ####
'./Function/'
       contains 
       
       (1) the subfolder 'fn', in this folder contains '.R' file code for 14 different structures that can deal with the data including both interval-censored values and missing values;
       (2) 'gen_cen_na_new.R' main script for generating the data with both interval-censored values and missing values, and the matrix "cen" indicates the corresponding censored and missing positions; and
       (3) 'F-G.r' main script for solving closed-form solutions of matrix D in the "EVE" and "VVE" structure by the FG algorithm; 

'./Function/fn'
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
    	
#### Subfolder: ./Code ####
'./Code/'
       contains 
       
       (1) 'fig1.R' main script for reproducing drawing Figure 1 that shows 14 GPMMs with three groups in two dimensions;
       (2) 'fig2.R' main script for reproducing Figure 2; 
       (3) 'fig3.4.R' main script for reproducing Figure 3 and Figure 4;
       (4) 'fig5.R' main script for reproducing Figure 5;
       (5) 'Fig.S1.R' main script for reproducing Figure S1;
       (6) 'Fig.S2(a).R' main script for reproducing Figure S2(a); 
       (7) 'Fig.S2(b).R' main script for reproducing Figure S2(b);
       (8) 'table1.2.R' main script for reproducing Tables 1 and 2;
       (9) 'table3.R' main script for Tables 3; 

###### Note for Section 5 - Illustrative examples - faithful data:
Because the 'faithful_code.R' and 'faithful_code_01.R' code takes a huge amount of time to run the ECM procedure for fitting the GPMM-CM model, we record these intermediate results in 'faithful.RData' and 'faithful_01.RData' so that one can use the R codes 'fig2.R' and 'fig3.4.R' to obtain the final results immediately.
To reproduce the results presented in Figure 2, just load 'faithful_01.RData' file in the './Data/' and then run the script 'fig2.R' in the subfolder './Code/';  
To reproduce the results presented in Figures 3 and 4, just load 'faithful.RData' file in the './Data/' and then run the script 'fig3.4.R' in the subfolder './Code/'

       (10) 'faithful_code.R' main script for fitting the GPMM-CM models with 14 parsimonious structures and g = 1-5 to the faithful data including both interval-censored values and missing values;
       (11) 'faithful_code_01.R' main script for fitting the GPMM-CM models with 14 parsimonious structures and g = 1-5 to the faithful data;

###### Note for Section 5 - Illustrative examples - hawks data:
Because the 'fit_hawkdata.R' code takes a huge amount of time to run the ECM procedure for fitting the GPMM-CM model, we record these intermediate results in 'hawks_new.RData' so that one can use the R codes 'fig5.R', and 'table1.2.R' to obtain the final results immediately.
To reproduce the results presented in Figures 5 and Tables 1 and 2, just load 'hawks_new.RData' file in the './Data/' and then run the script 'fig5.R', and 'table1.2.R' in the subfolder './Code/';  

       (12) 'fit_hawkdata.R' main script for fitting the GPMM-CM models with 14 parsimonious structures and g = 3 to the hawks data including both interval-censored values and missing values;
       

#### Subfolder: ./Data ####
'./Data/'
       contains 
       
       (1) 'faithful.RData' collecting the fitting results of the 14 GPMM-CM models to the old faithful data with interval-censored and missing values, where g ranging from 1 to 5;
       (2) 'faithful_01.RData' collecting the fitting results of the 14 GPMM-CM models to the old faithful data, where g ranging from 1 to 5;
       (3) 'hawks_new.RData' collecting the fitting and classification results of the 3-component GPMM-CM with 14 parsimonious structures to the hawks data with missing values and synthetic censored values; 
       (4) The subfolder 'source', which contains the DMS dataset used in Section 5.

'./Data/source'
	subfolder contains
	
    	(1) 'DMSdata.csv' for the DMS dataset used in Section 5.3.

#### Subfolder: ./Results ####
'./Results/'
       contains 
       
       (1) 'fig1.eps': an illustration of elliptical contours and scatter plots for 14 GPMMs with three groups in two dimensions;
       (2) 'fig2.eps'; perspective plot and bivariate density contours obtained from the best-fitted (EEE, 3) model; 
       (3) 'fig3.eps': comparison of BIC values for 14 GPMM-CMs with g = 1 − 5 fitted to the old faithful data involving censored and missing values;
       (4) 'fig4.eps': the histograms and scatter plot overlaid with contour curves obtained from the best fitted (EVI, 3) model for the old faithful data involving censored and missing values;
       (5) 'fig5.eps': the BIC curve for 14 GPMM-CMs with g = 3 fitted to the hawks data; 
       (6) 'Fig.S1.eps': Pairwise scatter and correlation plots for the 5 attributes of the DMS data;
       (7) 'Fig.S2(a).eps': Scatter-histogram plots of one simulation case with poorly separated mixture samples containing censored and missing values;
       (8) 'Fig.S2(b).eps': Scatter-histogram plots of one simulation case with well separated mixture samples containing censored and missing values;
       (9) 'Table1.csv': an overview of 5 variables of the hawks data with three species;
       (10) 'Table2.csv': the results of the ARI and CCR of the 14 clustering models with fitting the hawks data into three clusters; and 
       (11) Table3.csv': an overview of 5 variables of the DMS data with two California Current System.
       
# Additional Remark 
One can directly run each "source(.)" described in 'master.r' file in the seperate R session to obtain the results.

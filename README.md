# nc_covid19_bavaria

## Nowcast
This is the public code repository for the manuscript:

"Felix Günther, Andreas Bender, Katharina Katz, Helmut Küchenhoff, Michael Höhle: Nowcasting the COVID-19 pandemic in Bavaria"

To be published in Biometrical Journal.
A preprint of the manuscript is available here:
https://www.medrxiv.org/content/10.1101/2020.06.26.20140210v2

Code was mainly written by Günther, Bender, Höhle. We thank Titus Laska for initial code on the R(t) estimation.

### Content

#### Folder data_public:
Contains data to reproduce the reported analyses. The following files are included:
  - /main_analysis/synth_dat.csv: A synthetic dataset mimicking closely the original person-specific data on SARS-CoV-2 cases in Bavaria. It contains the columns 'rep_date_reg' (Reporting date at regional health authority, LGL), 'disease start' (symptom onset date if available, NA otherwise), 'age' (age of the person), 'sympt_expl' (dummy on reported COVID-19 related symptoms), and 'no_sympt_expl' (dummy indicating that an individual does explicitly has no symptoms at time-point of reporting). There are several individuals with no information on symptoms, i.e., neither known symptoms, nor explicitly symptom-free. The artificial data was created the following way:
We utilize for each of the 29246 cases the official reporting date at LGL and sampled an artifical age by adding/substracting a random number to the reported age. Based on an persons age, the reporting week and weekday, we sampled the expected reporting delay from the original Weibull GAMLSS imputation model and derived the artificial disease onset date for each case. We removed the sampled onset date for all persons without onset in the original data.

  - /main_analysis/onsets_true_retrospec_07-31.csv: contains the retrospective 'true' number of individuals with disease onset on a given day ('ntInf_true') from February, 24 until April, 06; based on reported and imputed disease onsets from person-specific data available at July, 31.
  
  - /evaluation/dat_mod_synthetic.RData: RData object containing synthetic data based for the retrospective evaluation of different nowcast models on synthetic data. The code to create the dataset is given in the code file /code_public/2_evaluation/est_synthetic_data/1_make_data.R. 15051 'person-specific' observations with reporting date (column 'rep_date'), weekday of reporting ('rep_date_weekday') and disease onset date ('disease_start')
  
  - /evaluation/epidemic_curve_synthetic_data.csv: Contains the epidemic curve (number of individuals with disease onset on a given day) used to create the synthetic data set dat_mod_synthetic.RData in the R-script /code_public/2_evaluation/est_synthetic_data/1_make_data.R.
  
  

In case of questions, feel free to contact felix.guenther(at)stat.uni-muenchen.de.

### Package and software versions

The computations were performed on a device with the following software versions installed and loaded (R-call to sessionInfo):
R version 3.6.3 (2020-02-29)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 18.04.5 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.7.1
LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.7.1

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8   
 [6] LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] parallel  splines   stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] future.apply_1.4.0       scoringRules_1.0.0       Rcpp_1.0.4.6             scales_1.1.0             furrr_0.1.0             
 [6] ggpubr_0.2.5             magrittr_1.5             surveillance_1.18.0.9000 xtable_1.8-4             sp_1.4-1                
[11] rstan_2.19.3             StanHeaders_2.19.2       gamlss_5.1-6             nlme_3.1-149             gamlss.dist_5.1-6       
[16] gamlss.data_5.1-4        future_1.16.0            broom_0.5.5              lubridate_1.7.4          nleqslv_3.3.2           
[21] R0_1.2-6                 MASS_7.3-52              data.table_1.12.8        forcats_0.5.0            stringr_1.4.0           
[26] dplyr_0.8.5              purrr_0.3.4              readr_1.3.1              tidyr_1.0.2              tibble_3.0.1            
[31] ggplot2_3.3.0            tidyverse_1.3.0         

loaded via a namespace (and not attached):
 [1] matrixStats_0.55.0    fs_1.3.2              httr_1.4.1            tools_3.6.3           backports_1.1.8       R6_2.4.1             
 [7] rpart_4.1-15          DBI_1.1.0             mgcv_1.8-33           colorspace_1.4-1      withr_2.2.0           tidyselect_1.0.0     
[13] gridExtra_2.3         prettyunits_1.1.1     processx_3.4.3        compiler_3.6.3        cli_2.0.2             rvest_0.3.5          
[19] xml2_1.3.1            spatstat.data_1.4-3   callr_3.4.3           spatstat_1.63-3       goftest_1.2-2         digest_0.6.25        
[25] spatstat.utils_1.17-0 pkgconfig_2.0.3       dbplyr_1.4.2          rlang_0.4.7           readxl_1.3.1          rstudioapi_0.11      
[31] generics_0.0.2        jsonlite_1.6.1        inline_0.3.15         loo_2.2.0             Matrix_1.2-18         munsell_0.5.0        
[37] fansi_0.4.1           abind_1.4-5           polyCub_0.7.1         lifecycle_0.2.0       stringi_1.4.6         pkgbuild_1.1.0       
[43] grid_3.6.3            listenv_0.8.0         crayon_1.3.4          deldir_0.1-25         lattice_0.20-41       haven_2.2.0          
[49] tensor_1.5            hms_0.5.3             knitr_1.28            ps_1.3.3              pillar_1.4.4          ggsignif_0.6.0       
[55] codetools_0.2-16      stats4_3.6.3          reprex_0.3.0          glue_1.4.1            modelr_0.1.6          vctrs_0.2.4          
[61] cellranger_1.1.0      gtable_0.3.0          polyclip_1.10-0       assertthat_0.2.1      xfun_0.12             survival_3.1-12      
[67] globals_0.12.5        ellipsis_0.3.1    
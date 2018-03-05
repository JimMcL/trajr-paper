# Trajr-paper - figures, analysis and examples for the `trajr` paper

This project contains data and R code used to generate figures and perform calculations used in the paper "Trajr: an R package for characterisation of animal trajectories" which describes the [`trajr`](https://cran.r-project.org/package=trajr) R package for analysing animal trajectories. As such, it provides small, self-contained examples of use of the package. To regenerate all figures and reports, run `source("R/generate-figures-reports.R")`. Figures are generated in the `figures` directory, and text reports in the `reports` directory. 

To reproduce the clearwing moth example and generate figure 4, it is first necessary to download the moth, bee and wasp trajectory files into the directory `data/clearwing-moths`. The files are accessable at (https://doi.org/10.5061/dryad.682dc).


## R files

* `R/generate-figures-reports.R` generates all figures and text files. It calls functions in `R/whale-analysis.R` and `R/clearwing-analysis.R` as required.

* `R/whale-analysis.R` defines functions which implement the sample analysis for whale trajectories.

* `R/clearwing-analysis.R` defines functions which implement the sample analysis for clearwing moth locomotor mimicry.


## Environment

```
Session info --------------------------------------------------------------
 setting  value                       
 version  R version 3.4.2 (2017-09-28)
 system   x86_64, mingw32             
 ui       RStudio (1.1.383)           
 language (EN)                        
 collate  English_United States.1252  
 tz       Australia/Sydney            
 date     2018-01-13                  

Packages ------------------------------------------------------------------
 package     * version    date       source                             
 base        * 3.4.2      2017-09-28 local                              
 codetools     0.2-15     2016-10-05 CRAN (R 3.4.2)                     
 colorspace    1.3-2      2016-12-14 CRAN (R 3.4.2)                     
 compiler      3.4.2      2017-09-28 local                              
 datasets    * 3.4.2      2017-09-28 local                              
 devtools      1.13.4     2017-11-09 CRAN (R 3.4.3)                     
 digest        0.6.12     2017-01-27 CRAN (R 3.4.2)                     
 geosphere     1.5-7      2017-11-05 CRAN (R 3.4.2)                     
 ggmap       * 2.6.1      2016-01-23 CRAN (R 3.4.3)                     
 ggplot2     * 2.2.1      2016-12-30 CRAN (R 3.4.2)                     
 graphics    * 3.4.2      2017-09-28 local                              
 grDevices   * 3.4.2      2017-09-28 local                              
 grid          3.4.2      2017-09-28 local                              
 gtable        0.2.0      2016-02-26 CRAN (R 3.4.2)                     
 jpeg          0.1-8      2014-01-23 CRAN (R 3.4.1)                     
 labeling      0.3        2014-08-23 CRAN (R 3.4.1)                     
 lattice       0.20-35    2017-03-25 CRAN (R 3.4.2)                     
 lazyeval      0.2.1      2017-10-29 CRAN (R 3.4.2)                     
 magrittr      1.5        2014-11-22 CRAN (R 3.4.2)                     
 mapproj       1.2-5      2017-06-08 CRAN (R 3.4.3)                     
 maps          3.2.0      2017-06-08 CRAN (R 3.4.3)                     
 MASS          7.3-47     2017-02-26 CRAN (R 3.4.2)                     
 memoise       1.1.0      2017-04-21 CRAN (R 3.4.2)                     
 methods     * 3.4.2      2017-09-28 local                              
 munsell       0.4.3      2016-02-13 CRAN (R 3.4.2)                     
 plyr          1.8.4      2016-06-08 CRAN (R 3.4.2)                     
 png           0.1-7      2013-12-03 CRAN (R 3.4.1)                     
 proto         1.0.0      2016-10-29 CRAN (R 3.4.2)                     
 Rcpp          0.12.13    2017-09-28 CRAN (R 3.4.2)                     
 reshape2      1.4.2      2016-10-22 CRAN (R 3.4.2)                     
 rgdal       * 1.2-15     2017-10-30 CRAN (R 3.4.2)                     
 RgoogleMaps   1.4.1      2016-09-18 CRAN (R 3.4.2)                     
 rjson         0.2.15     2014-11-03 CRAN (R 3.4.1)                     
 rlang         0.1.4.9000 2017-11-18 Github (tidyverse/rlang@1e54041)   
 rstudioapi    0.7.0-9000 2017-11-18 Github (rstudio/rstudioapi@335f257)
 scales        0.5.0      2017-08-24 CRAN (R 3.4.2)                     
 signal        0.7-6      2015-07-30 CRAN (R 3.4.2)                     
 sp          * 1.2-5      2017-06-29 CRAN (R 3.4.2)                     
 stats       * 3.4.2      2017-09-28 local                              
 stringi       1.1.6      2017-11-17 CRAN (R 3.4.2)                     
 stringr       1.2.0      2017-02-18 CRAN (R 3.4.2)                     
 tibble        1.3.4      2017-08-22 CRAN (R 3.4.2)                     
 tools         3.4.2      2017-09-28 local                              
 trajr       * 1.0.0      2018-01-12 local                              
 utils       * 3.4.2      2017-09-28 local                              
 withr         2.1.0      2017-11-01 CRAN (R 3.4.2)                     
 yaml          2.1.14     2016-11-12 CRAN (R 3.4.2)               
 ```

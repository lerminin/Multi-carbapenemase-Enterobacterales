
This repository accompanies the manuscript: Global genomic epidemiology of multi-carbapenemase Enterobacterales (available as a preprint on bioRxiv: [https://doi.org/10.64898/2026.10.07.757280](https://doi.org/10.64898/2026.10.07.757280))

Scripts for plots and statistical analysis in R are provided as a Quarto file `code.qmd` with the corresponding rendered HTML file `code.html`. Data files required to run the code are included in the `data/` folder. The HTML file can be viewed in GitHub by clicking the `code.html.md` file. 

Note that several figures included in the manuscript (Figure 2 legend, Figure 3) were made or modified with external programs (Inkscape, Cytoscape) and are not included here.

This code was run with the following R and package versions (output of `sessionInfo()`) and takes less than 15 mins to run:

```R
R version 4.6.0 (2026-04-24 ucrt)
Platform: x86_64-w64-mingw32/x64
Running under: Windows 11 x64 (build 26200)

Matrix products: default
  LAPACK version 3.12.1

locale:
[1] LC_COLLATE=English_United Kingdom.utf8  LC_CTYPE=English_United Kingdom.utf8    LC_MONETARY=English_United Kingdom.utf8
[4] LC_NUMERIC=C                            LC_TIME=English_United Kingdom.utf8    

time zone: Europe/London
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] goftest_1.2-3      rstatix_1.1.0      Matrix_1.7-5       logistf_1.26.1     treemapify_2.6.0   networkD3_0.4.1    gridpattern_1.4.2 
 [8] ggpattern_1.3.1    ComplexUpset_1.3.3 patchwork_1.3.2    tibble_3.3.1       ggplot2_4.0.3      readr_2.2.0        here_1.0.2        
[15] stringr_1.6.0      tidyr_1.3.2        dplyr_1.2.1       

loaded via a namespace (and not attached):
 [1] tidyselect_1.2.1       farver_2.1.2           S7_0.2.2               fastmap_1.2.0          digest_0.6.39          rpart_4.1.27          
 [7] lifecycle_1.0.5        survival_3.8-6         magrittr_2.0.5         compiler_4.6.0         rlang_1.2.0            tools_4.6.0           
[13] igraph_2.3.1           htmlwidgets_1.6.4      RColorBrewer_1.1-3     abind_1.4-8            withr_3.0.2            purrr_1.2.2           
[19] nnet_7.3-21            grid_4.6.0             jomo_2.7-6             colorspace_2.1-2       mice_3.19.0            data.tree_1.2.0       
[25] scales_1.4.0           iterators_1.0.14       MASS_7.3-65            cli_3.6.6              reformulas_0.4.4       generics_0.1.4        
[31] tzdb_0.5.0             minqa_1.2.8            cachem_1.1.0           operator.tools_1.6.3.1 splines_4.6.0          vctrs_0.7.3           
[37] boot_1.3-32            glmnet_5.0             carData_3.0-6          jsonlite_2.0.0         car_3.1-5              hms_1.1.4             
[43] mitml_0.4-5            Formula_1.2-5          systemfonts_1.3.2      foreach_1.5.2          glue_1.8.1             nloptr_2.2.1          
[49] pan_2.0                codetools_0.2-20       stringi_1.8.7          gtable_0.3.6           shape_1.4.6.1          lme4_2.0-6            
[55] pillar_1.11.1          htmltools_0.5.9        ggfittext_0.10.3       R6_2.6.1               textshaping_1.0.5      Rdpack_2.6.6          
[61] formula.tools_1.7.1    rprojroot_2.1.1        lattice_0.22-9         rbibutils_2.4.1        backports_1.5.1        memoise_2.0.1         
[67] broom_1.0.13           Rcpp_1.1.1-1.1         svglite_2.2.2          nlme_3.1-169           mgcv_1.9-4             pkgconfig_2.0.3  
```

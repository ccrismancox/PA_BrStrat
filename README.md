# Replication Archive for "Detecting and Correcting for Separation in Strategic Choice Models"

### Casey Crisman-Cox, Olga Gasparyan, and Curtis S. Signorino

### August 15, 2022

This archive contains the replication code and data for "Detecting and Correcting for Separation in Strategic Choice Models." All computations were performed on a Ubuntu 20.04.4 computer in R 4.1.3. The original computer has 64 GB memory and dual Intel Xeon Silver processors (40 cores total). On that machine the complete archive is runs in about 45 minutes.

## R packages and session info

The session information is printed below

    R version 4.1.3 (2022-03-10)
    Platform: x86_64-pc-linux-gnu (64-bit)
    Running under: Ubuntu 20.04.4 LTS

    Matrix products: default
    BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.9.0
    LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.9.0

    locale:
     [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
     [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
     [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
     [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
     [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

    attached base packages:
    [1] stats     graphics  grDevices utils     datasets  methods   base     

    other attached packages:
     [1] mc2d_0.1-21          mvtnorm_1.1-3        stringr_1.4.0       
     [4] gridExtra_2.3        ggplot2_3.3.6        readstata13_0.10.0  
     [7] numDeriv_2016.8-1.1  knitr_1.33           matrixStats_0.60.1  
    [10] detectseparation_0.2 brglm_0.7.2          profileModel_0.6.1  
    [13] doRNG_1.8.2          rngtools_1.5.2       foreach_1.5.2       
    [16] games2_0.1.1         MASS_7.3-55          Formula_1.2-4       
    [19] maxLik_1.5-2         miscTools_0.6-26     devtools_2.4.4      
    [22] usethis_2.1.6       

    loaded via a namespace (and not attached):
     [1] pkgload_1.3.0            shiny_1.6.0              assertthat_0.2.1        
     [4] ROI.plugin.lpsolve_1.0-1 remotes_2.4.2            slam_0.1-50             
     [7] sessioninfo_1.2.2        pillar_1.8.0             lattice_0.20-45         
    [10] glue_1.6.2               digest_0.6.29            promises_1.2.0.1        
    [13] colorspace_2.0-3         sandwich_3.0-2           htmltools_0.5.2         
    [16] httpuv_1.6.2             pkgconfig_2.0.3          purrr_0.3.4             
    [19] xtable_1.8-4             scales_1.2.0             processx_3.5.2          
    [22] later_1.3.0              tibble_3.1.8             generics_0.1.3          
    [25] ellipsis_0.3.2           withr_2.5.0              cachem_1.0.6            
    [28] cli_3.3.0                magrittr_2.0.3           crayon_1.5.1            
    [31] mime_0.12                memoise_2.0.1            ps_1.6.0                
    [34] fs_1.5.2                 fansi_1.0.3              pkgbuild_1.3.1          
    [37] profvis_0.3.7            tools_4.1.3              registry_0.5-1          
    [40] prettyunits_1.1.1        lifecycle_1.0.1          ROI_1.0-0               
    [43] munsell_0.5.0            callr_3.7.0              compiler_4.1.3          
    [46] rlang_1.0.4              grid_4.1.3               iterators_1.0.14        
    [49] htmlwidgets_1.5.3        miniUI_0.1.1.1           gtable_0.3.0            
    [52] codetools_0.2-18         DBI_1.1.1                curl_4.3.2              
    [55] R6_2.5.1                 lpSolveAPI_5.5.2.0-17.8  zoo_1.8-10              
    [58] dplyr_1.0.7              fastmap_1.1.0            utf8_1.2.2              
    [61] stringi_1.7.8            parallel_4.1.3           Rcpp_1.0.9              
    [64] vctrs_0.4.1              tidyselect_1.1.1         xfun_0.32               
    [67] urlchecker_1.0.1

## Archive contents

In this section, we describe the files contained in this replication archive and what they contain.

### Main level

The main level contains three files and three folders. The files are:

-   `README.md` This file the in text format
-   `README.pdf` This file in pdf format
-   `master.r` An R code file that replicates the paper. This file changes the working directory to the `code` folder and runs all of the R scripts.

The three folders are described below

### Data

The folder `data` contains one file.

-   `huth.dta` The replication data from Signorino and Tarar (2006) and originally from Huth (1998). The variable descriptions and full citations can be found in the Online Appendix.

### Code

The folder `code` contains 25 files: 13 R scripts and 12 rdata output files.

-   `packages.r` This file installs the package versions used in the analysis.
-   `extraFunctions.r` This file contains helper functions for the simulations and replication
-   `mainSimulation.R` This files runs the main simulation reported in the text. It reproduces Tables 1, 2, and B.1 in files `tables_and_figures/Table1.md`, `tables_and_figures/Table2.md`, and `tables_and_figures/TableB1.md`, respectively. Its raw output is saved in the file `table1.rdata`.
-   `SignorinoAndTararReplication.R` This file replicates the Signorino and Tarar (2006) study. This file reproduces Tables 3 and 4 and Figure 3 in files `tables_and_figures/Table3.md`, `tables_and_figures/Table4.md`, and `tables_and_figures/figure3.pdf`, respectively. Its main output is saved in the files `sigTarar_output.rdata` and `signorinorarar_brfit.rdata`.
-   `TableB2.R` This file runs the simulation in Appendix B.2. It reproduces Table B.2 in the file `tables_and_figures/TableB2.md`. Its raw output is saved in the file `tableB2.rdata`
-   `TableB3.R` This file runs the first simulation in Appendix B.3. It reproduces Table B.3 in the file `tables_and_figures/TableB3.md`. Its raw output is saved in the file `tableB3.rdata`
-   `TableB4.R` This file runs the second simulation in Appendix B.3. It reproduces Table B.4 in the file `tables_and_figures/TableB4.md`. Its raw output is saved in the file `tableB4.rdata`
-   `TableB5.R` This file runs the simulation in Appendix B.4. It reproduces Table B.5 in the file `tables_and_figures/TableB5.md`. Its raw output is saved in the file `tableB5.rdata`
-   `TableB6.R` This file runs first the simulation in Appendix B.4.1. It reproduces Table B.6 in the file `tables_and_figures/TableB6.md`. Its raw output is saved in the file `tableB6.rdata`
-   `TableB7.R` This file runs second the simulation in Appendix B.4.1. It reproduces Table B.7 in the file `tables_and_figures/TableB7.md`. Its raw output is saved in the file `tableB7.rdata`
-   `TableB8.R` This file runs first the simulation in Appendix B.5. It reproduces Table B.8 in the file `tables_and_figures/TableB8.md`. Its raw output is saved in the file `tableB8.rdata`
-   `TableB9.R` This file runs second the simulation in Appendix B.5. It reproduces Table B.9 in the file `tables_and_figures/TableB9.md`. Its raw output is saved in the file `tableB9.rdata`
-   `Leblang.R` This file replicates results from Leblang (2003) using the data packaged within the `games2` package. It reproduces Tables D.1-2 in files `tables_and_figures/TableD1.md` and `tables_and_figures/TableD2.md`, respectively. The raw output is saved in the file `Leblang_output.rdata`

### Table and figures

The folder `tables_and_figures` contains 16 files. The generation and contents of these files are described in the Code section above

## Running the code

Files may be run individually from the `code` folder or the file `master.r` can be run from the main folder.

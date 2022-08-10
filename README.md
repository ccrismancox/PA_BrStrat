---
output:
  html_document: default
  pdf_document: default
---
# Replication Arhive for "Detecting and Correcting for Separation in Strategic Choice Models" 

## Note to the PA replication team
Dear replication team, 

Thank you for taking the time to reproduce the results in our paper. We have one note about differences we found in preparing the replication archive.  Specifically, some of the values in Table D.2 (Online Appendix) have changed. We attribute this to slight adjustments in coding and software routines. The changes are minor and lead to no changes in the text or any conclusions. The values in `tableD2.md` are correct and are what we ask be published. Thank you.

## R packages and session info

The session information is printed below 

```
R version 4.1.3 (2022-03-10)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 20.04.4 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.9.0
LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.9.0

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8       
 [4] LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
[10] LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] parallel  stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] xtable_1.8-4         knitr_1.33           matrixStats_0.60.1   detectseparation_0.2
 [5] brglm_0.7.2          profileModel_0.6.1   doRNG_1.8.2          rngtools_1.5.2      
 [9] doParallel_1.0.16    iterators_1.0.13     foreach_1.5.1        games2_0.1.1        
[13] MASS_7.3-55          Formula_1.2-4        maxLik_1.5-2         miscTools_0.6-26    
[17] devtools_2.4.2       usethis_2.0.1       

loaded via a namespace (and not attached):
 [1] pkgload_1.2.1            carData_3.0-4            ROI.plugin.lpsolve_1.0-1 cellranger_1.1.0        
 [5] remotes_2.4.0            slam_0.1-48              sessioninfo_1.1.1        numDeriv_2016.8-1.1     
 [9] pillar_1.6.2             lattice_0.20-45          glue_1.6.2               digest_0.6.29           
[13] sandwich_3.0-2           pkgconfig_2.0.3          haven_2.4.3              purrr_0.3.4             
[17] processx_3.5.2           openxlsx_4.2.4           rio_0.5.27               tibble_3.1.4            
[21] generics_0.1.3           car_3.0-11               ellipsis_0.3.2           cachem_1.0.6            
[25] withr_2.4.2              cli_3.0.1                magrittr_2.0.3           crayon_1.4.1            
[29] readxl_1.3.1             memoise_2.0.0            ps_1.6.0                 fs_1.5.0                
[33] fansi_0.5.0              forcats_0.5.1            foreign_0.8-82           pkgbuild_1.2.0          
[37] tools_4.1.3              registry_0.5-1           data.table_1.14.0        prettyunits_1.1.1       
[41] hms_1.1.0                lifecycle_1.0.0          ROI_1.0-0                stringr_1.4.0           
[45] zip_2.2.0                callr_3.7.0              compiler_4.1.3           rlang_0.4.11            
[49] grid_4.1.3               rstudioapi_0.13          testthat_3.0.4           codetools_0.2-18        
[53] abind_1.4-5              curl_4.3.2               R6_2.5.1                 lpSolveAPI_5.5.2.0-17.7 
[57] zoo_1.8-10               fastmap_1.1.0            utf8_1.2.2               rprojroot_2.0.2         
[61] desc_1.3.0               stringi_1.7.8            Rcpp_1.0.7               vctrs_0.3.8             
[65] xfun_0.25   
```

## Archive contents
In this section, we describe the files contained in this replication archive and what they contain.

### Main level
The main level contains three files and three folders. The files are:

- `README.md` This file the in text format
- `README.pdf` This file in pdf format
- `master.r` An R code file that replicates the paper. This file changes the working directory to the `code` folder and runs all of the R scripts.

The three folders are described below

### Data
The folder `data` contains one file.

- `huth.dta` The replication data from Huth (1998). These data were graciously provided by Curt Signorino who used them in Signorino and Tarar (2006). The variable descriptions and full citations can be found in the Online Appendix.

### Code

The folder `code` contains 25 files.

- `packages.r` This file installs the package versions used in the analysis.
- `extraFunctions.r` This file contains helper functions for the simulations and replication
- `mainSimulation.R` This files runs the main simulation reported in the text. It reproduces Tables 1, 2, and B.1 in files `tables_and_figures/Table1.md`, `tables_and_figures/Table2.md`, and `tables_and_figures/TableB1.md`, respectively. Its raw output is saved in the file `table1.rdata`.
- `SignorinoAndTararReplication.R` This file replicates the Signorino and Tarar (2006) study. This file reproduces Tables 3 and 4 and Figure 3  in files `tables_and_figures/Table3.md`, `tables_and_figures/Table4.md`, and `tables_and_figures/figure3.pdf`, respectively. Its main output is saved in the files `sigTarar_output.rdata` and `signorinorarar_brfit.rdata`.
- `TableB2.R` This file runs the simulation in Appendix B.2.  It reproduces Table B.2 in the file `tables_and_figures/TableB2.md. Its raw output is saved in the file `tableB2.rdata`
- `TableB3.R` This file runs the first simulation in Appendix B.3.  It reproduces Table B.3 in the file `tables_and_figures/TableB3.md. Its raw output is saved in the file `tableB3.rdata`
- `TableB4.R` This file runs the second simulation in Appendix B.3.  It reproduces Table B.4 in the file `tables_and_figures/TableB4.md. Its raw output is saved in the file `tableB4.rdata`
- `TableB5.R` This file runs the simulation in Appendix B.4.  It reproduces Table B.5 in the file `tables_and_figures/TableB5.md. Its raw output is saved in the file `tableB5.rdata`
- `TableB6.R` This file runs first the simulation in Appendix B.4.1.  It reproduces Table B.6 in the file `tables_and_figures/TableB6.md. Its raw output is saved in the file `tableB6.rdata`
- `TableB7.R` This file runs second the simulation in Appendix B.4.1.  It reproduces Table B.7 in the file `tables_and_figures/TableB7.md. Its raw output is saved in the file `tableB7.rdata`
- `TableB8.R` This file runs first the simulation in Appendix B.5.  It reproduces Table B.8 in the file `tables_and_figures/TableB8.md. Its raw output is saved in the file `tableB8.rdata`
- `TableB9.R` This file runs second the simulation in Appendix B.5.  It reproduces Table B.9 in the file `tables_and_figures/TableB9.md. Its raw output is saved in the file `tableB9.rdata`
- `Leblang.R` This file replicates results from Leblang (2003) using the data packaged within the `games2` package. It reproduces Tables D.1-2 in files `tables_and_figures/TableD1.md` and `tables_and_figures/TableD2.md`, respectively. The raw output is saved in the file `Leblang_output.rdata`


### Table and figures
The folder `tables_and_figures` contains 16 files. The generation and contents of these files are described in the Code section above

## Running the code
Files may be run individually from the `code` folder or the file `master.r` can be run from the main folder.

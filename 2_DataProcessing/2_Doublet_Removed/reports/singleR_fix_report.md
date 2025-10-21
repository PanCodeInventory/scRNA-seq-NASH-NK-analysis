# SingleR 空 data 层修复报告

- 输入对象：/home/harry/NASH/scRNA-seq/Files/Doublet_Removed/RDS/nk.integrated.doublet_scored.rds
- 选用 assay：RNA
- 处理策略：从 counts 构建 SCE，使用 scater::logNormCounts 生成 logcounts，然后以 logcounts 作为 SingleR 的输入。
- 注释列：singleR.label（必有），singleR.pruned（若可用）
- 输出对象：/home/harry/NASH/scRNA-seq/Files/Doublet_Removed/RDS/nk.integrated.singleR_annotated.rds

## 标签统计（Top10）
```
      label count
3       ILC  9854
5  NK cells  9153
1   B cells    67
6       NKT    33
2        DC     5
4 Monocytes     4
```

## R 会话信息
```
R version 4.4.1 (2024-06-14)
Platform: x86_64-pc-linux-gnu
Running under: Ubuntu 20.04.6 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/liblapack.so.3;  LAPACK version 3.9.0

locale:
 [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
 [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
 [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
[10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   

time zone: Asia/Shanghai
tzcode source: system (glibc)

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods  
[8] base     

other attached packages:
 [1] scater_1.32.1               ggplot2_3.5.1              
 [3] scuttle_1.14.0              celldex_1.14.0             
 [5] SingleR_2.6.0               SingleCellExperiment_1.26.0
 [7] SummarizedExperiment_1.34.0 Biobase_2.64.0             
 [9] GenomicRanges_1.56.1        GenomeInfoDb_1.40.1        
[11] IRanges_2.38.1              S4Vectors_0.42.1           
[13] BiocGenerics_0.50.0         MatrixGenerics_1.16.0      
[15] matrixStats_1.3.0           Seurat_5.1.0               
[17] SeuratObject_5.0.2          sp_2.1-4                   

loaded via a namespace (and not attached):
  [1] RcppAnnoy_0.0.22          splines_4.4.1            
  [3] later_1.3.2               filelock_1.0.3           
  [5] tibble_3.2.1              polyclip_1.10-7          
  [7] fastDummies_1.7.4         httr2_1.0.3              
  [9] lifecycle_1.0.4           globals_0.16.3           
 [11] lattice_0.22-6            MASS_7.3-65              
 [13] alabaster.base_1.4.2      magrittr_2.0.3           
 [15] plotly_4.10.4             yaml_2.3.10              
 [17] httpuv_1.6.15             sctransform_0.4.1        
 [19] spam_2.10-0               spatstat.sparse_3.1-0    
 [21] reticulate_1.38.0         cowplot_1.1.3            
 [23] pbapply_1.7-2             DBI_1.2.3                
 [25] RColorBrewer_1.1-3        abind_1.4-5              
 [27] zlibbioc_1.50.0           Rtsne_0.17               
 [29] purrr_1.0.2               rappdirs_0.3.3           
 [31] GenomeInfoDbData_1.2.12   ggrepel_0.9.5            
 [33] irlba_2.3.5.1             listenv_0.9.1            
 [35] spatstat.utils_3.1-0      goftest_1.2-3            
 [37] RSpectra_0.16-2           spatstat.random_3.3-1    
 [39] fitdistrplus_1.2-1        parallelly_1.38.0        
 [41] DelayedMatrixStats_1.26.0 leiden_0.4.3.1           
 [43] codetools_0.2-20          DelayedArray_0.30.1      
 [45] tidyselect_1.2.1          farver_2.1.2             
 [47] UCSC.utils_1.0.0          viridis_0.6.5            
 [49] ScaledMatrix_1.12.0       BiocFileCache_2.12.0     
 [51] spatstat.explore_3.3-2    jsonlite_1.8.8           
 [53] BiocNeighbors_1.22.0      progressr_0.14.0         
 [55] ggridges_0.5.6            survival_3.8-3           
 [57] systemfonts_1.1.0         tools_4.4.1              
 [59] ragg_1.3.2                ica_1.0-3                
 [61] Rcpp_1.0.13               glue_1.7.0               
 [63] gridExtra_2.3             SparseArray_1.4.8        
 [65] HDF5Array_1.32.1          gypsum_1.0.1             
 [67] dplyr_1.1.4               withr_3.0.1              
 [69] BiocManager_1.30.25       fastmap_1.2.0            
 [71] rhdf5filters_1.16.0       fansi_1.0.6              
 [73] digest_0.6.37             rsvd_1.0.5               
 [75] R6_2.5.1                  mime_0.12                
 [77] textshaping_0.4.0         colorspace_2.1-1         
 [79] scattermore_1.2           tensor_1.5               
 [81] spatstat.data_3.1-2       RSQLite_2.3.7            
 [83] utf8_1.2.4                tidyr_1.3.1              
 [85] generics_0.1.3            data.table_1.16.0        
 [87] httr_1.4.7                htmlwidgets_1.6.4        
 [89] S4Arrays_1.4.1            uwot_0.2.2               
 [91] pkgconfig_2.0.3           gtable_0.3.5             
 [93] blob_1.2.4                lmtest_0.9-40            
 [95] XVector_0.44.0            htmltools_0.5.8.1        
 [97] dotCall64_1.1-1           alabaster.matrix_1.4.2   
 [99] scales_1.3.0              png_0.1-8                
[101] spatstat.univar_3.0-0     reshape2_1.4.4           
[103] nlme_3.1-166              curl_5.2.2               
[105] rhdf5_2.48.0              zoo_1.8-12               
[107] cachem_1.1.0              stringr_1.5.1            
[109] BiocVersion_3.19.1        KernSmooth_2.23-26       
[111] vipor_0.4.7               parallel_4.4.1           
[113] miniUI_0.1.1.1            AnnotationDbi_1.66.0     
[115] alabaster.schemas_1.4.0   pillar_1.9.0             
[117] grid_4.4.1                vctrs_0.6.5              
[119] RANN_2.6.2                promises_1.3.0           
[121] BiocSingular_1.20.0       dbplyr_2.5.0             
[123] beachmat_2.20.0           xtable_1.8-4             
[125] cluster_2.1.8.1           beeswarm_0.4.0           
[127] cli_3.6.3                 compiler_4.4.1           
[129] rlang_1.1.4               crayon_1.5.3             
[131] future.apply_1.11.2       labeling_0.4.3           
[133] ggbeeswarm_0.7.2          plyr_1.8.9               
[135] stringi_1.8.4             alabaster.se_1.4.1       
[137] viridisLite_0.4.2         deldir_2.0-4             
[139] BiocParallel_1.38.0       munsell_0.5.1            
[141] Biostrings_2.72.1         lazyeval_0.2.2           
[143] spatstat.geom_3.3-2       Matrix_1.7-0             
[145] ExperimentHub_2.12.0      RcppHNSW_0.6.0           
[147] patchwork_1.2.0           sparseMatrixStats_1.16.0 
[149] bit64_4.0.5               future_1.34.0            
[151] Rhdf5lib_1.26.0           KEGGREST_1.44.1          
[153] shiny_1.9.1               alabaster.ranges_1.4.2   
[155] AnnotationHub_3.12.0      ROCR_1.0-11              
[157] igraph_2.0.3              memoise_2.0.1            
[159] bit_4.0.5                
```

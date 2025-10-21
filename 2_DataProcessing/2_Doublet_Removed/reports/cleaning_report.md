# NK细胞去双胞与去污染报告

- 输入对象：/home/harry/NASH/scRNA-seq/Files/UMAP/Results/rds/nk.integrated.rds
- 备份对象：/home/harry/NASH/scRNA-seq/Files/Doublet_Removed/RDS/nk.integrated.backup.rds
- 双胞评分对象：/home/harry/NASH/scRNA-seq/Files/Doublet_Removed/RDS/nk.integrated.doublet_scored.rds
- 清理后对象：/home/harry/NASH/scRNA-seq/Files/Doublet_Removed/RDS/nk.integrated.filtered.rds

## 参数与阈值
- UCell NK分位阈值：P60 (0.6107)
- UCell Δ(NK - max(nonNK)) 阈值：0.050
- 簇级非NK占比剔除阈值：0.70
- 聚类分辨率：0.5；PCA/UMAP维度：1:30

## 统计
- 初始细胞数：19438
- 去双胞后细胞数：19126（移除 312，1.61%）
- 去污染后细胞数：19040（再移除 86，0.45%）

## 簇级非NK占比
```
# A tibble: 9 × 2
  seurat_clusters frac_nonNK
  <fct>                <dbl>
1 0                  0.00240
2 1                  0      
3 2                  0      
4 3                  0      
5 4                  0      
6 5                  0.00141
7 6                  0      
8 7                  0.411  
9 8                  0      
```

## 已保存图件
- DoubletScore_Density.png
- DoubletScore_Density_bySample.png（若按样本分组）
- Cluster_nonNK_fraction.png
- UMAP_filtered_by_SingleR.png
- UMAP_filtered_clusters_by_Timepoint.png（若含timepoint）
- NK_UCell_by_SingleR_label.png
- Cluster_UCell_signature_means.png

## R会话信息
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
 [1] tidyselect_1.2.1            tidyr_1.3.1                
 [3] readr_2.1.5                 UCell_2.8.0                
 [5] scater_1.32.1               scuttle_1.14.0             
 [7] celldex_1.14.0              SingleR_2.6.0              
 [9] scDblFinder_1.18.0          SingleCellExperiment_1.26.0
[11] SummarizedExperiment_1.34.0 Biobase_2.64.0             
[13] GenomicRanges_1.56.1        GenomeInfoDb_1.40.1        
[15] IRanges_2.38.1              S4Vectors_0.42.1           
[17] BiocGenerics_0.50.0         MatrixGenerics_1.16.0      
[19] matrixStats_1.3.0           patchwork_1.2.0            
[21] ggplot2_3.5.1               dplyr_1.1.4                
[23] Seurat_5.1.0                SeuratObject_5.0.2         
[25] sp_2.1-4                   

loaded via a namespace (and not attached):
  [1] spatstat.sparse_3.1-0     bitops_1.0-8             
  [3] httr_1.4.7                RColorBrewer_1.1-3       
  [5] tools_4.4.1               sctransform_0.4.1        
  [7] alabaster.base_1.4.2      utf8_1.2.4               
  [9] R6_2.5.1                  HDF5Array_1.32.1         
 [11] lazyeval_0.2.2            uwot_0.2.2               
 [13] rhdf5filters_1.16.0       withr_3.0.1              
 [15] gridExtra_2.3             progressr_0.14.0         
 [17] textshaping_0.4.0         cli_3.6.3                
 [19] spatstat.explore_3.3-2    fastDummies_1.7.4        
 [21] labeling_0.4.3            alabaster.se_1.4.1       
 [23] spatstat.data_3.1-2       ggridges_0.5.6           
 [25] pbapply_1.7-2             systemfonts_1.1.0        
 [27] Rsamtools_2.20.0          parallelly_1.38.0        
 [29] limma_3.60.4              RSQLite_2.3.7            
 [31] generics_0.1.3            BiocIO_1.14.0            
 [33] vroom_1.6.5               ica_1.0-3                
 [35] spatstat.random_3.3-1     Matrix_1.7-0             
 [37] ggbeeswarm_0.7.2          fansi_1.0.6              
 [39] abind_1.4-5               lifecycle_1.0.4          
 [41] yaml_2.3.10               edgeR_4.2.1              
 [43] rhdf5_2.48.0              SparseArray_1.4.8        
 [45] BiocFileCache_2.12.0      Rtsne_0.17               
 [47] grid_4.4.1                blob_1.2.4               
 [49] promises_1.3.0            dqrng_0.4.1              
 [51] ExperimentHub_2.12.0      crayon_1.5.3             
 [53] miniUI_0.1.1.1            lattice_0.22-6           
 [55] beachmat_2.20.0           cowplot_1.1.3            
 [57] KEGGREST_1.44.1           pillar_1.9.0             
 [59] metapod_1.12.0            rjson_0.2.22             
 [61] xgboost_1.7.8.1           future.apply_1.11.2      
 [63] codetools_0.2-20          leiden_0.4.3.1           
 [65] glue_1.7.0                spatstat.univar_3.0-0    
 [67] data.table_1.16.0         vctrs_0.6.5              
 [69] png_0.1-8                 gypsum_1.0.1             
 [71] spam_2.10-0               gtable_0.3.5             
 [73] cachem_1.1.0              S4Arrays_1.4.1           
 [75] mime_0.12                 survival_3.8-3           
 [77] statmod_1.5.0             bluster_1.14.0           
 [79] fitdistrplus_1.2-1        ROCR_1.0-11              
 [81] nlme_3.1-166              bit64_4.0.5              
 [83] alabaster.ranges_1.4.2    filelock_1.0.3           
 [85] RcppAnnoy_0.0.22          irlba_2.3.5.1            
 [87] vipor_0.4.7               KernSmooth_2.23-26       
 [89] colorspace_2.1-1          DBI_1.2.3                
 [91] bit_4.0.5                 compiler_4.4.1           
 [93] curl_5.2.2                httr2_1.0.3              
 [95] BiocNeighbors_1.22.0      DelayedArray_0.30.1      
 [97] plotly_4.10.4             rtracklayer_1.64.0       
 [99] scales_1.3.0              lmtest_0.9-40            
[101] rappdirs_0.3.3            stringr_1.5.1            
[103] digest_0.6.37             goftest_1.2-3            
[105] spatstat.utils_3.1-0      alabaster.matrix_1.4.2   
[107] XVector_0.44.0            htmltools_0.5.8.1        
[109] pkgconfig_2.0.3           sparseMatrixStats_1.16.0 
[111] dbplyr_2.5.0              fastmap_1.2.0            
[113] rlang_1.1.4               htmlwidgets_1.6.4        
[115] UCSC.utils_1.0.0          shiny_1.9.1              
[117] DelayedMatrixStats_1.26.0 farver_2.1.2             
[119] zoo_1.8-12                jsonlite_1.8.8           
[121] BiocParallel_1.38.0       BiocSingular_1.20.0      
[123] RCurl_1.98-1.16           magrittr_2.0.3           
[125] GenomeInfoDbData_1.2.12   dotCall64_1.1-1          
[127] Rhdf5lib_1.26.0           munsell_0.5.1            
[129] Rcpp_1.0.13               viridis_0.6.5            
[131] reticulate_1.38.0         stringi_1.8.4            
[133] alabaster.schemas_1.4.0   zlibbioc_1.50.0          
[135] MASS_7.3-65               AnnotationHub_3.12.0     
[137] plyr_1.8.9                parallel_4.4.1           
[139] listenv_0.9.1             ggrepel_0.9.5            
[141] deldir_2.0-4              Biostrings_2.72.1        
[143] splines_4.4.1             tensor_1.5               
[145] hms_1.1.3                 locfit_1.5-9.10          
[147] igraph_2.0.3              spatstat.geom_3.3-2      
[149] RcppHNSW_0.6.0            reshape2_1.4.4           
[151] ScaledMatrix_1.12.0       BiocVersion_3.19.1       
[153] XML_3.99-0.17             scran_1.32.0             
[155] BiocManager_1.30.25       tzdb_0.4.0               
[157] httpuv_1.6.15             RANN_2.6.2               
[159] purrr_1.0.2               polyclip_1.10-7          
[161] future_1.34.0             scattermore_1.2          
[163] rsvd_1.0.5                xtable_1.8-4             
[165] restfulr_0.0.15           RSpectra_0.16-2          
[167] later_1.3.2               ragg_1.3.2               
[169] viridisLite_0.4.2         tibble_3.2.1             
[171] memoise_2.0.1             beeswarm_0.4.0           
[173] AnnotationDbi_1.66.0      GenomicAlignments_1.40.0 
[175] cluster_2.1.8.1           globals_0.16.3           
```

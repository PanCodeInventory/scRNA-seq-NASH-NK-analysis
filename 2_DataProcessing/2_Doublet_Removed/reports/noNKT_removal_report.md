# NKT 剔除报告

- 输入对象：/home/harry/NASH/scRNA-seq/2_DataProcessing/2_Doublet_Removed/RDS/nk.integrated.singleR_annotated.rds
- 标签列：singleR.pruned

## 规则
- 保留：NK/ILC（匹配 nk / nk cell / natural killer / ilc / innate lymphoid）
- 排除：NKT（匹配 nkt / nk t / natural killer t）

## 统计
- 剔除前：19126 细胞
- 剔除后：19007 细胞（移除 119，0.62%）

## 标签计数（Top10 剔除前/后）
```
剔除前（Top10）：

      label count_before
3       ILC         9854
5  NK cells         9153
1   B cells           67
6       NKT           33
2        DC            5
4 Monocytes            4

剔除后（Top10）：

     label count_after
1      ILC        9854
2 NK cells        9153

```

## 产物
- 新 RDS：/home/harry/NASH/scRNA-seq/2_DataProcessing/2_Doublet_Removed/RDS/nk.integrated.noNKT.rds
- 移除细胞 ID：removed_NKT_cell_ids.csv
- 标签计数对比图：NKT_removal_label_counts_before_after.png

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
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] readr_2.1.5        patchwork_1.2.0    ggplot2_3.5.1      dplyr_1.1.4       
[5] Seurat_5.1.0       SeuratObject_5.0.2 sp_2.1-4          

loaded via a namespace (and not attached):
  [1] deldir_2.0-4           pbapply_1.7-2          gridExtra_2.3         
  [4] rlang_1.1.4            magrittr_2.0.3         RcppAnnoy_0.0.22      
  [7] matrixStats_1.3.0      ggridges_0.5.6         compiler_4.4.1        
 [10] spatstat.geom_3.3-2    systemfonts_1.1.0      png_0.1-8             
 [13] vctrs_0.6.5            reshape2_1.4.4         stringr_1.5.1         
 [16] crayon_1.5.3           pkgconfig_2.0.3        fastmap_1.2.0         
 [19] labeling_0.4.3         utf8_1.2.4             promises_1.3.0        
 [22] tzdb_0.4.0             ragg_1.3.2             bit_4.0.5             
 [25] purrr_1.0.2            jsonlite_1.8.8         goftest_1.2-3         
 [28] later_1.3.2            spatstat.utils_3.1-0   irlba_2.3.5.1         
 [31] parallel_4.4.1         cluster_2.1.8.1        R6_2.5.1              
 [34] ica_1.0-3              stringi_1.8.4          RColorBrewer_1.1-3    
 [37] spatstat.data_3.1-2    reticulate_1.38.0      parallelly_1.38.0     
 [40] spatstat.univar_3.0-0  lmtest_0.9-40          scattermore_1.2       
 [43] Rcpp_1.0.13            tensor_1.5             future.apply_1.11.2   
 [46] zoo_1.8-12             sctransform_0.4.1      httpuv_1.6.15         
 [49] Matrix_1.7-0           splines_4.4.1          igraph_2.0.3          
 [52] tidyselect_1.2.1       abind_1.4-5            spatstat.random_3.3-1 
 [55] codetools_0.2-20       miniUI_0.1.1.1         spatstat.explore_3.3-2
 [58] listenv_0.9.1          lattice_0.22-6         tibble_3.2.1          
 [61] plyr_1.8.9             withr_3.0.1            shiny_1.9.1           
 [64] ROCR_1.0-11            Rtsne_0.17             future_1.34.0         
 [67] fastDummies_1.7.4      survival_3.8-3         polyclip_1.10-7       
 [70] fitdistrplus_1.2-1     pillar_1.9.0           KernSmooth_2.23-26    
 [73] plotly_4.10.4          generics_0.1.3         vroom_1.6.5           
 [76] RcppHNSW_0.6.0         hms_1.1.3              munsell_0.5.1         
 [79] scales_1.3.0           globals_0.16.3         xtable_1.8-4          
 [82] glue_1.7.0             lazyeval_0.2.2         tools_4.4.1           
 [85] data.table_1.16.0      RSpectra_0.16-2        RANN_2.6.2            
 [88] leiden_0.4.3.1         dotCall64_1.1-1        cowplot_1.1.3         
 [91] grid_4.4.1             tidyr_1.3.1            colorspace_2.1-1      
 [94] nlme_3.1-166           cli_3.6.3              spatstat.sparse_3.1-0 
 [97] textshaping_0.4.0      spam_2.10-0            fansi_1.0.6           
[100] viridisLite_0.4.2      uwot_0.2.2             gtable_0.3.5          
[103] digest_0.6.37          progressr_0.14.0       ggrepel_0.9.5         
[106] farver_2.1.2           htmlwidgets_1.6.4      htmltools_0.5.8.1     
[109] lifecycle_1.0.4        httr_1.4.7             mime_0.12             
[112] bit64_4.0.5            MASS_7.3-65           
```

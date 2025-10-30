#!/usr/bin/env Rscript

#' 简化UMAP可视化：SingleR细胞类型 vs Cluster分簇对比
#' 
#' 功能：生成左侧SingleR标注、右侧Cluster分簇的UMAP组合图
#' 输入：nk.integrated.v1.rds
#' 输出：组合UMAP图 (PNG/PDF)

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
})

# 设置路径
rds_path <- "2_DataProcessing/RDS/nk.integrated.v1.rds"
output_dir <- "3_Analysis/1_ClusterAnalysis/plots"

# 创建输出目录
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# 主函数
main <- function() {
  message("=== 开始UMAP可视化 ===")
  
  # 加载数据
  message("加载RDS文件...")
  if (!file.exists(rds_path)) {
    stop(sprintf("RDS文件不存在: %s", rds_path))
  }
  
  obj <- readRDS(rds_path)
  message(sprintf("成功加载: %d个细胞, %d个基因", ncol(obj), nrow(obj)))
  
  # 检查SingleR注释列
  singleR_col <- NULL
  if ("singleR.pruned" %in% colnames(obj@meta.data)) {
    singleR_col <- "singleR.pruned"
    message("使用singleR.pruned列作为SingleR注释")
  } else if ("singleR.label" %in% colnames(obj@meta.data)) {
    singleR_col <- "singleR.label" 
    message("使用singleR.label列作为SingleR注释")
  } else {
    message("警告：未找到SingleR注释列，将使用默认细胞类型")
    # 创建默认细胞类型（基于簇）
    obj$cell_type <- paste0("Cluster_", obj$seurat_clusters)
    singleR_col <- "cell_type"
  }
  
  # 确保UMAP存在
  if (!"umap" %in% names(obj@reductions)) {
    stop("错误：RDS文件中未找到UMAP降维结果")
  }
  
  message("生成UMAP图...")
  
  # 筛选细胞：只保留NK和ILC细胞类型（剔除NKT）
  message("筛选NK和ILC细胞（剔除NKT）...")
  nk_ilc_cells <- obj@meta.data %>% 
    filter(grepl("NK|ILC", .data[[singleR_col]], ignore.case = TRUE) & 
           !grepl("NKT", .data[[singleR_col]], ignore.case = TRUE)) %>%
    rownames()
  
  if (length(nk_ilc_cells) > 0) {
    obj_filtered <- subset(obj, cells = nk_ilc_cells)
    message(sprintf("筛选后细胞数: %d (%.1f%%)", 
                   length(nk_ilc_cells), 
                   100 * length(nk_ilc_cells) / ncol(obj)))
  } else {
    message("警告：未找到NK/ILC细胞，使用原始数据")
    obj_filtered <- obj
  }
  
  # 筛选簇：移除簇6
  message("移除簇6...")
  if ("6" %in% unique(obj_filtered$seurat_clusters)) {
    obj_filtered <- subset(obj_filtered, subset = seurat_clusters != "6")
    message("已移除簇6的细胞")
  }
  
  # SingleR细胞类型UMAP（左侧）
  p_singler <- DimPlot(
    obj_filtered,
    reduction = "umap",
    group.by = singleR_col,
    label = TRUE,
    label.size = 4,
    repel = TRUE
  ) + 
    ggtitle("SingleR Cell Type Annotation (NK/ILC only)") +
    theme(
      plot.title = element_text(size = 12, hjust = 0.5),
      legend.position = "right"
    )
  
  # Cluster分簇UMAP（右侧）
  p_clusters <- DimPlot(
    obj_filtered,
    reduction = "umap", 
    group.by = "seurat_clusters",
    label = TRUE,
    label.size = 4,
    repel = TRUE
  ) + 
    ggtitle("Seurat Cluster Assignment (Cluster 6 excluded)") +
    theme(
      plot.title = element_text(size = 12, hjust = 0.5),
      legend.position = "right"
    )
  
  # 组合图
  message("创建组合图...")
  combined_plot <- p_singler | p_clusters
  
  # 保存结果
  message("保存图件...")
  
  # PNG格式
  png_path <- file.path(output_dir, "UMAP_SingleR_vs_Clusters.png")
  ggsave(png_path, plot = combined_plot, width = 16, height = 8, dpi = 300, bg = "white")
  message(sprintf("PNG图件保存: %s", png_path))
  
  # PDF格式
  pdf_path <- file.path(output_dir, "UMAP_SingleR_vs_Clusters.pdf")
  ggsave(pdf_path, plot = combined_plot, width = 16, height = 8)
  message(sprintf("PDF图件保存: %s", pdf_path))
  
  # 保存数据摘要
  summary_path <- file.path("3_Analysis/1_ClusterAnalysis/data", "umap_visualization_summary.txt")
  dir.create(dirname(summary_path), showWarnings = FALSE, recursive = TRUE)
  
  summary_text <- c(
    sprintf("UMAP可视化摘要 - %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    sprintf("输入文件: %s", basename(rds_path)),
    sprintf("细胞总数: %d", ncol(obj)),
    sprintf("基因总数: %d", nrow(obj)),
    sprintf("SingleR注释列: %s", singleR_col),
    sprintf("簇数: %d", length(unique(obj$seurat_clusters))),
    sprintf("细胞类型数: %d", length(unique(obj@meta.data[[singleR_col]]))),
    sprintf("输出图件: UMAP_SingleR_vs_Clusters.png/pdf")
  )
  
  writeLines(summary_text, summary_path)
  message(sprintf("摘要文件保存: %s", summary_path))
  
  message("=== UMAP可视化完成 ===")
  message(sprintf("细胞类型数: %d", length(unique(obj@meta.data[[singleR_col]]))))
  message(sprintf("簇数: %d", length(unique(obj$seurat_clusters))))
}

# 运行主函数
tryCatch({
  main()
}, error = function(e) {
  message(sprintf("错误: %s", e$message))
  quit(status = 1)
})

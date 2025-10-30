#!/usr/bin/env Rscript

#' 根据SingleR标注筛选RDS文件
#' 
#' 功能：只保留"singleR.label"列中标记为"NK cells"和"ILC"的细胞
#' 输入：原始RDS文件
#' 输出：筛选后的v2版本RDS文件

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

# 设置路径
input_rds_path <- "2_DataProcessing/RDS/nk.integrated.v1.rds"
output_dir <- "2_DataProcessing/RDS"

# 主函数
main <- function() {
  message("=== 开始RDS文件筛选 ===")
  
  # 加载数据
  message("加载RDS文件...")
  if (!file.exists(input_rds_path)) {
    stop(sprintf("RDS文件不存在: %s", input_rds_path))
  }
  
  obj <- readRDS(input_rds_path)
  message(sprintf("原始数据: %d个细胞, %d个基因", ncol(obj), nrow(obj)))
  
  # 检查SingleR注释列
  singleR_col <- NULL
  if ("singleR.label" %in% colnames(obj@meta.data)) {
    singleR_col <- "singleR.label"
    message("使用singleR.label列作为SingleR注释")
  } else if ("singleR.pruned" %in% colnames(obj@meta.data)) {
    singleR_col <- "singleR.pruned"
    message("使用singleR.pruned列作为SingleR注释")
  } else {
    stop("错误：未找到SingleR注释列（singleR.label 或 singleR.pruned）")
  }
  
  # 查看原始细胞类型分布
  message("原始细胞类型分布:")
  cell_type_summary <- table(obj@meta.data[[singleR_col]])
  print(cell_type_summary)
  
  # 筛选细胞：只保留"NK cells"和"ILC"细胞
  message("筛选NK cells和ILC细胞...")
  target_cells <- obj@meta.data %>% 
    filter(.data[[singleR_col]] %in% c("NK cells", "ILC")) %>%
    rownames()
  
  if (length(target_cells) == 0) {
    stop("错误：未找到符合条件的细胞（NK cells 或 ILC）")
  }
  
  # 创建筛选后的对象
  obj_filtered <- subset(obj, cells = target_cells)
  
  message(sprintf("筛选后: %d个细胞 (%.1f%%保留)", 
                  length(target_cells), 
                  100 * length(target_cells) / ncol(obj)))
  
  # 查看筛选后的细胞类型分布
  message("筛选后细胞类型分布:")
  filtered_summary <- table(obj_filtered@meta.data[[singleR_col]])
  print(filtered_summary)
  
  # 生成输出文件名
  base_name <- tools::file_path_sans_ext(basename(input_rds_path))
  output_filename <- paste0(base_name, ".v2.rds")
  output_path <- file.path(output_dir, output_filename)
  
  # 保存筛选后的RDS文件
  message(sprintf("保存筛选后的RDS文件: %s", output_path))
  saveRDS(obj_filtered, file = output_path)
  
  # 保存筛选摘要
  summary_path <- file.path(output_dir, paste0(base_name, "_filtering_summary.txt"))
  summary_text <- c(
    sprintf("RDS文件筛选摘要 - %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    sprintf("输入文件: %s", basename(input_rds_path)),
    sprintf("输出文件: %s", output_filename),
    sprintf("原始细胞数: %d", ncol(obj)),
    sprintf("筛选后细胞数: %d", length(target_cells)),
    sprintf("保留比例: %.1f%%", 100 * length(target_cells) / ncol(obj)),
    sprintf("筛选条件: %s列中值为'NK cells'或'ILC'", singleR_col),
    "",
    "原始细胞类型分布:",
    paste(names(cell_type_summary), cell_type_summary, sep = ": "),
    "",
    "筛选后细胞类型分布:",
    paste(names(filtered_summary), filtered_summary, sep = ": ")
  )
  
  writeLines(summary_text, summary_path)
  message(sprintf("筛选摘要保存: %s", summary_path))
  
  message("=== RDS文件筛选完成 ===")
  message(sprintf("新文件路径: %s", output_path))
}

# 运行主函数
tryCatch({
  main()
}, error = function(e) {
  message(sprintf("错误: %s", e$message))
  quit(status = 1)
})

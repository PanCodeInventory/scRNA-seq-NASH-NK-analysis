#!/usr/bin/env Rscript

#' 清理RDS文件metadata中的冗余列
#' 
#' 功能：删除指定的metadata列，包括orig.ident和integrated_snn_res系列列
#' 输入：原始RDS文件
#' 输出：清理后的v3版本RDS文件 (nk.integrated.v3.rds)

suppressPackageStartupMessages({
  library(Seurat)
})

# 设置路径
input_rds_path <- "2_DataProcessing/RDS/nk.integrated.v1.rds"
output_dir <- "2_DataProcessing/RDS"

# 要删除的列
columns_to_remove <- c(
  "orig.ident",
  "integrated_snn_res.0.6", 
  "integrated_snn_res.0.5",
  "integrated_snn_res.0.4", 
  "integrated_snn_res.0.7"
)

# 主函数
main <- function() {
  message("=== 开始清理RDS metadata ===")
  
  # 加载数据
  message("加载RDS文件...")
  if (!file.exists(input_rds_path)) {
    stop(sprintf("RDS文件不存在: %s", input_rds_path))
  }
  
  obj <- readRDS(input_rds_path)
  message(sprintf("原始数据: %d个细胞, %d个基因", ncol(obj), nrow(obj)))
  
  # 检查当前metadata列
  current_meta_cols <- colnames(obj@meta.data)
  message(sprintf("当前metadata列数: %d", length(current_meta_cols)))
  message("当前metadata列:")
  print(current_meta_cols)
  
  # 找出需要删除的列（存在的列）
  cols_to_delete <- intersect(columns_to_remove, current_meta_cols)
  
  if (length(cols_to_delete) > 0) {
    message(sprintf("找到%d个需要删除的列:", length(cols_to_delete)))
    print(cols_to_delete)
    
    # 删除指定的列
    obj@meta.data <- obj@meta.data[, !(colnames(obj@meta.data) %in% cols_to_delete)]
    
    message("已删除指定的metadata列")
  } else {
    message("警告：未找到需要删除的列")
  }
  
  # 检查删除后的metadata列
  new_meta_cols <- colnames(obj@meta.data)
  message(sprintf("清理后metadata列数: %d", length(new_meta_cols)))
  message("清理后metadata列:")
  print(new_meta_cols)
  
  # 保存为v3版本
  output_filename <- "nk.integrated.v3.rds"
  output_path <- file.path(output_dir, output_filename)
  
  # 保存清理后的RDS文件
  message(sprintf("保存清理后的RDS文件: %s", output_path))
  saveRDS(obj, file = output_path)
  
  # 保存清理摘要
  summary_path <- file.path(output_dir, "nk.integrated_metadata_cleaning_summary.txt")
  summary_text <- c(
    sprintf("RDS metadata清理摘要 - %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    sprintf("输入文件: %s", basename(input_rds_path)),
    sprintf("输出文件: %s", output_filename),
    sprintf("原始metadata列数: %d", length(current_meta_cols)),
    sprintf("清理后metadata列数: %d", length(new_meta_cols)),
    sprintf("删除的列数: %d", length(cols_to_delete)),
    "",
    "删除的列:",
    paste("-", cols_to_delete),
    "",
    "原始metadata列:",
    paste("-", current_meta_cols),
    "",
    "清理后metadata列:",
    paste("-", new_meta_cols)
  )
  
  writeLines(summary_text, summary_path)
  message(sprintf("清理摘要保存: %s", summary_path))
  
  message("=== RDS metadata清理完成 ===")
  message(sprintf("新文件路径: %s", output_path))
}

# 运行主函数
tryCatch({
  main()
}, error = function(e) {
  message(sprintf("错误: %s", e$message))
  quit(status = 1)
})

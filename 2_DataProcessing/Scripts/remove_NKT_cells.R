# 基于 SingleR 注释从对象中剔除 NKT 细胞，并保存新 RDS + 对比图 + 报告
# 输入优先：2_Filter/2_Doublet_Removed/RDS/nk.integrated.singleR_annotated.rds
# 回退：Files/Doublet_Removed/RDS/nk.integrated.singleR_annotated.rds
# 输出（统一到 2_Filter/2_Doublet_Removed/* 路径下）：
#   - RDS：nk.integrated.noNKT.rds
#   - CSV：removed_NKT_cell_ids.csv
#   - 图件：NKT_removal_label_counts_before_after.png
#   - 报告：noNKT_removal_report.md
#
# 用法：
#   Rscript 2_DataProcessing/Scripts/remove_NKT_cells.R

suppressPackageStartupMessages({
  library(utils)
})

options(stringsAsFactors = FALSE)
set.seed(1234)
message(sprintf("remove_NKT_cells START: %s", as.character(Sys.time())))
cat("remove_NKT_cells START\n"); flush.console()

# 路径设置
root_dir <- "/home/harry/NASH/scRNA-seq"
in_rds_candidates <- c(
  file.path(root_dir, "2_DataProcessing", "2_Doublet_Removed", "RDS", "nk.integrated.singleR_annotated.rds"),
  file.path(root_dir, "2_Filter", "2_Doublet_Removed", "RDS", "nk.integrated.singleR_annotated.rds"),
  file.path(root_dir, "Files", "Doublet_Removed", "RDS", "nk.integrated.singleR_annotated.rds")
)
out_rds_dir <- file.path(root_dir, "2_DataProcessing", "2_Doublet_Removed", "RDS")
out_plots_dir <- file.path(root_dir, "2_DataProcessing", "2_Doublet_Removed", "plots")
out_reports_dir <- file.path(root_dir, "2_DataProcessing", "2_Doublet_Removed", "reports")

dirs <- c(out_rds_dir, out_plots_dir, out_reports_dir)
for (d in dirs) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# 依赖安装/加载工具
ensure_pkg <- function(pkg, bioc = FALSE) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(sprintf("安装缺失依赖：%s (%s)", pkg, if (bioc) "Bioconductor" else "CRAN"))
    if (bioc) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager", repos = "https://cloud.r-project.org")
      }
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    } else {
      install.packages(pkg, repos = "https://cloud.r-project.org")
    }
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

# 加载依赖
ensure_pkg("Seurat", bioc = FALSE)
ensure_pkg("SeuratObject", bioc = FALSE)
ensure_pkg("dplyr", bioc = FALSE)
ensure_pkg("ggplot2", bioc = FALSE)
ensure_pkg("patchwork", bioc = FALSE)
ensure_pkg("readr", bioc = FALSE)

# 读取对象（优先 2_DataProcessing/2_Doublet_Removed，其次 2_Filter 与 Files）
in_rds <- NULL
for (p in in_rds_candidates) {
  if (file.exists(p)) { in_rds <- p; break }
}
if (is.null(in_rds)) {
  stop(sprintf("未找到输入对象：%s", paste(in_rds_candidates, collapse = " | ")))
}
seu <- readRDS(in_rds)
message(sprintf("读取对象成功：%s", in_rds))
message(sprintf("Assays：%s", paste(Seurat::Assays(seu), collapse = ", ")))
message(sprintf("细胞数：%d", ncol(seu)))

# 选择标签列：优先 singleR.pruned，其次 singleR.label
lab_col <- NULL
if ("singleR.pruned" %in% colnames(seu@meta.data)) {
  lab_col <- "singleR.pruned"
} else if ("singleR.label" %in% colnames(seu@meta.data)) {
  lab_col <- "singleR.label"
}
if (is.null(lab_col)) stop("对象中未找到 SingleR 注释列（singleR.pruned/singleR.label）")

labels_raw <- seu@meta.data[[lab_col]]
if (is.null(labels_raw)) stop(sprintf("元数据列 '%s' 为空", lab_col))
label_lower <- tolower(as.character(labels_raw))

# 定义 NKT 排除与 NK/ILC 保留的规则
# - NKT 排除：匹配 'nkt'、'nk t'、'natural killer t'
is_nkt <- grepl("\\bnkt\\b|\\bnk\\s*t\\b|natural\\s*killer\\s*t", label_lower)
# - NK/ILC 保留：匹配 'nk'（或 'nk cell'、'natural killer (cell)'）、'ilc'、'innate lymphoid'
is_nk_ilc <- grepl("\\bnk\\b|\\bnk\\s*cell\\b|natural\\s*killer(\\s*cell)?\\b|\\bilc\\b|innate\\s*lymphoid", label_lower)

# 保留集合：NK/ILC 且 非 NKT
keep_cells <- is_nk_ilc & (!is_nkt)

n_before <- ncol(seu)
n_keep <- sum(keep_cells, na.rm = TRUE)
n_remove <- n_before - n_keep
pct_remove <- if (n_before > 0) 100 * n_remove / n_before else NA_real_

message(sprintf("NKT 剔除：%d -> %d 细胞（移除 %d，%.2f%%）", n_before, n_keep, n_remove, pct_remove))

# 生成标签计数对比（剔除前/后）
df_before <- as.data.frame(table(labels_raw))
colnames(df_before) <- c("label", "count_before")

labels_after <- labels_raw[keep_cells]
df_after <- as.data.frame(table(labels_after))
colnames(df_after) <- c("label", "count_after")

df_comp <- dplyr::full_join(df_before, df_after, by = "label") %>%
  dplyr::mutate(
    count_before = ifelse(is.na(count_before), 0L, count_before),
    count_after  = ifelse(is.na(count_after), 0L, count_after),
    removed = count_before - count_after
  ) %>%
  dplyr::arrange(dplyr::desc(count_before))

# 保存移除的细胞 ID 列表
removed_ids <- colnames(seu)[!keep_cells]
readr::write_csv(data.frame(cell_id = removed_ids), file.path(out_reports_dir, "removed_NKT_cell_ids.csv"))

# 子集保留并保存新对象
seu_noNKT <- subset(seu, cells = colnames(seu)[keep_cells])
out_rds <- file.path(out_rds_dir, "nk.integrated.noNKT.rds")
saveRDS(seu_noNKT, out_rds)
message(sprintf("已保存无 NKT 对象：%s", out_rds))

# 生成对比图：剔除前/后标签柱状图
p_before <- ggplot2::ggplot(df_before, ggplot2::aes(x = reorder(label, -count_before), y = count_before)) +
  ggplot2::geom_col(fill = "#6baed6") +
  ggplot2::coord_flip() +
  ggplot2::theme_bw() +
  ggplot2::labs(title = sprintf("SingleR labels before (n=%d)", n_before), x = "label", y = "count")

p_after <- ggplot2::ggplot(df_after, ggplot2::aes(x = reorder(label, -count_after), y = count_after)) +
  ggplot2::geom_col(fill = "#3182bd") +
  ggplot2::coord_flip() +
  ggplot2::theme_bw() +
  ggplot2::labs(title = sprintf("SingleR labels after NKT removal (n=%d)", n_keep), x = "label", y = "count")

p_comp <- p_before + p_after + patchwork::plot_annotation(title = "NKT removal: label counts before vs after")
ggplot2::ggsave(filename = file.path(out_plots_dir, "NKT_removal_label_counts_before_after.png"), plot = p_comp, width = 12, height = 10, dpi = 300)

# 报告写入
report_path <- file.path(out_reports_dir, "noNKT_removal_report.md")
report_lines <- c(
  "# NKT 剔除报告",
  "",
  sprintf("- 输入对象：%s", in_rds),
  sprintf("- 标签列：%s", lab_col),
  "",
  "## 规则",
  "- 保留：NK/ILC（匹配 nk / nk cell / natural killer / ilc / innate lymphoid）",
  "- 排除：NKT（匹配 nkt / nk t / natural killer t）",
  "",
  "## 统计",
  sprintf("- 剔除前：%d 细胞", n_before),
  sprintf("- 剔除后：%d 细胞（移除 %d，%.2f%%）", n_keep, n_remove, pct_remove),
  "",
  "## 标签计数（Top10 剔除前/后）",
  "```\n剔除前（Top10）：\n",
  paste(utils::capture.output(print(utils::head(df_before[order(-df_before$count_before), ], 10))), collapse = "\n"),
  "\n剔除后（Top10）：\n",
  paste(utils::capture.output(print(utils::head(df_after[order(-df_after$count_after), ], 10))), collapse = "\n"),
  "\n```",
  "",
  "## 产物",
  sprintf("- 新 RDS：%s", out_rds),
  "- 移除细胞 ID：removed_NKT_cell_ids.csv",
  "- 标签计数对比图：NKT_removal_label_counts_before_after.png",
  "",
  "## R 会话信息",
  paste0("```\n", paste(utils::capture.output(utils::sessionInfo()), collapse = "\n"), "\n```")
)
writeLines(report_lines, con = report_path)
message(sprintf("报告已生成：%s", report_path))

message("remove_NKT_cells 全流程完成。")

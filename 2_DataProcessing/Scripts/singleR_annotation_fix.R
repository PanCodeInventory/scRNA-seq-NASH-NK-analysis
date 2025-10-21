# 针对 SingleR 报警“Layer 'data' is empty”的特异性修复脚本
# 目标：在 data 层为空时，为 SingleR 构造含 logcounts 的 SCE 输入，避免报警并完成自动注释
# 输入优先：Files/Doublet_Removed/RDS/nk.integrated.doublet_scored.rds（若不存在则回退到 Files/UMAP/Results/rds/nk.integrated.rds）
# 输出：
#   - Files/Doublet_Removed/RDS/nk.integrated.singleR_annotated.rds（写入 singleR.label / singleR.pruned）
#   - Files/Doublet_Removed/reports/singleR_fix_report.md（记录处理过程与会话信息）
#   - Files/Doublet_Removed/plots/SingleR_label_barplot.png（注释标签柱状图）
#
# 用法：
#   Rscript Files/UMAP/scripts/singleR_annotation_fix.R
#
# 关键策略：
#   - 从 Seurat 对象中稳健提取 counts（RNA/SCT/integrated），避免依赖空的 data 层
#   - 使用 scater::logNormCounts 计算 logcounts，保证 SingleR 输入满足要求
#   - 若存在 scDblFinder.class，则先保留 singlet 再注释
#   - 输出标签统计与报告

suppressPackageStartupMessages({
  library(utils)
})

options(stringsAsFactors = FALSE)
set.seed(1234)
message(sprintf("SingleR FIX SCRIPT START: %s", as.character(Sys.time())))
cat("SingleR FIX SCRIPT START\n"); flush.console()

# 路径设置
root_dir <- "/home/harry/NASH/scRNA-seq"
in_rds_primary <- file.path(root_dir, "Files/Doublet_Removed", "RDS", "nk.integrated.doublet_scored.rds")
in_rds_fallback <- file.path(root_dir, "Files/UMAP/Results/rds", "nk.integrated.rds")
out_rds_dir <- file.path(root_dir, "Files/Doublet_Removed", "RDS")
out_plots_dir <- file.path(root_dir, "Files/Doublet_Removed", "plots")
out_reports_dir <- file.path(root_dir, "Files/Doublet_Removed", "reports")

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
ensure_pkg("SingleCellExperiment", bioc = TRUE)
ensure_pkg("SingleR", bioc = TRUE)
ensure_pkg("celldex", bioc = TRUE)
ensure_pkg("scater", bioc = TRUE)
ensure_pkg("ggplot2", bioc = FALSE)

# 读取对象（优先 doublet_scored）
in_rds <- if (file.exists(in_rds_primary)) in_rds_primary else in_rds_fallback
if (!file.exists(in_rds)) stop(sprintf("未找到输入对象：%s", in_rds))
seu <- readRDS(in_rds)
message(sprintf("读取对象成功：%s", in_rds))
message(sprintf("Assays 可用：%s", paste(Seurat::Assays(seu), collapse = ", ")))

# 若存在双胞标签，先保留 singlet
if ("scDblFinder.class" %in% colnames(seu@meta.data)) {
  n_before <- ncol(seu)
  seu <- subset(seu, subset = scDblFinder.class == "singlet")
  message(sprintf("保留 singlet：%d -> %d 细胞（移除 %d，%.2f%%）",
                  n_before, ncol(seu), n_before - ncol(seu),
                  100 * (n_before - ncol(seu)) / n_before))
}

# 选择 assay（优先 RNA，其次 SCT，再次 integrated）
assays_avail <- tryCatch(Seurat::Assays(seu), error = function(e) character(0))
if (length(assays_avail) == 0) stop("对象不包含任何 assay。")
chosen_assay <- if ("RNA" %in% assays_avail) "RNA" else if ("SCT" %in% assays_avail) "SCT" else assays_avail[[1]]
Seurat::DefaultAssay(seu) <- chosen_assay
message(sprintf("选用 assay：%s", chosen_assay))

# 稳健提取 counts（优先 layer=counts，若无则 slot=counts；必要时尝试其他 assay）
get_counts <- function(obj, assay_name) {
  # v5 推荐 layer 接口
  m <- tryCatch(SeuratObject::GetAssayData(object = obj[[assay_name]], layer = "counts"), error = function(e) NULL)
  if (!is.null(m) && !is.null(dim(m)) && all(dim(m) > 0)) return(m)
  # 回退 slot 接口（兼容 v4/v5）
  m <- tryCatch(Seurat::GetAssayData(obj[[assay_name]], slot = "counts"), error = function(e) NULL)
  if (!is.null(m) && !is.null(dim(m)) && all(dim(m) > 0)) return(m)
  return(NULL)
}

counts <- get_counts(seu, chosen_assay)
if (is.null(counts) && chosen_assay != "RNA") counts <- get_counts(seu, "RNA")
if (is.null(counts) && chosen_assay != "SCT") counts <- get_counts(seu, "SCT")
if (is.null(counts) && chosen_assay != "integrated") counts <- get_counts(seu, "integrated")

if (is.null(counts)) {
  stop("无法从 Seurat 对象提取非空的 counts（已尝试 RNA/SCT/integrated 的 layer/slot）。")
}

# 与 Seurat 对象细胞对齐
cells <- colnames(seu)
common <- intersect(cells, colnames(counts))
if (length(common) == 0) stop("counts 与 Seurat 对象细胞名不重叠。")
counts <- counts[, common, drop = FALSE]
cd <- seu@meta.data[common, , drop = FALSE]

# 构建 SCE 并计算 logcounts（避免 data 层为空问题）
sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts), colData = cd)
# 使用 scater 进行库大小标准化并生成 logcounts
sce <- scater::logNormCounts(sce)  # 在 assays 中创建 'logcounts'
message("已使用 scater::logNormCounts 生成 logcounts。")

# 参考集
ref <- celldex::ImmGenData()
# SingleR 自动注释（显式指定使用 logcounts）
pred <- SingleR::SingleR(
  test = sce,
  ref = ref,
  labels = ref$label.main,
  assay.type.test = "logcounts",
  assay.type.ref = "logcounts"
)

# 写回 Seurat meta.data
seu <- subset(seu, cells = common)  # 对齐 singleR 预测的细胞集
seu$singleR.label <- pred$labels
if (!is.null(pred$pruned.labels)) {
  seu$singleR.pruned <- pred$pruned.labels
}
message("SingleR 注释已写入 Seurat 对象。")

# 保存对象
out_rds <- file.path(out_rds_dir, "nk.integrated.singleR_annotated.rds")
saveRDS(seu, out_rds)
message(sprintf("已保存注释对象：%s", out_rds))

# 可视化：标签柱状图
lab_col <- if ("singleR.pruned" %in% colnames(seu@meta.data)) "singleR.pruned" else "singleR.label"
df_lab <- as.data.frame(table(seu@meta.data[[lab_col]]))
colnames(df_lab) <- c("label", "count")
p_bar <- ggplot2::ggplot(df_lab, aes(x = reorder(label, -count), y = count)) +
  ggplot2::geom_col(fill = "#3182bd") +
  ggplot2::coord_flip() +
  ggplot2::theme_bw() +
  ggplot2::labs(title = "SingleR annotations (counts)", x = "label", y = "count")
ggplot2::ggsave(filename = file.path(out_plots_dir, "SingleR_label_barplot.png"), plot = p_bar, width = 10, height = 8, dpi = 300)

# 报告写入
report_path <- file.path(out_reports_dir, "singleR_fix_report.md")
report_lines <- c(
  "# SingleR 空 data 层修复报告",
  "",
  sprintf("- 输入对象：%s", in_rds),
  sprintf("- 选用 assay：%s", chosen_assay),
  "- 处理策略：从 counts 构建 SCE，使用 scater::logNormCounts 生成 logcounts，然后以 logcounts 作为 SingleR 的输入。",
  "- 注释列：singleR.label（必有），singleR.pruned（若可用）",
  sprintf("- 输出对象：%s", out_rds),
  "",
  "## 标签统计（Top10）",
  paste0("```\n", paste(utils::capture.output(print(utils::head(df_lab[order(-df_lab$count), ], 10))), collapse = "\n"), "\n```"),
  "",
  "## R 会话信息",
  paste0("```\n", paste(utils::capture.output(utils::sessionInfo()), collapse = "\n"), "\n```")
)
writeLines(report_lines, con = report_path)
message(sprintf("报告已生成：%s", report_path))

message("SingleR FIX 全流程完成。")

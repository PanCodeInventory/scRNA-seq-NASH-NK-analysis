suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

options(stringsAsFactors = FALSE)
set.seed(1234)

# 输入/输出路径（独立于原始结果）
rds_in  <- "Files/UMAP/Results/rds/nk.integrated.rds"
out_dir <- "Files/UMAP/Results/plots_equalized"

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

# 读取既有整合对象（不重跑任何降维/聚类）
message("Loading integrated object from: ", rds_in)
nk.integrated <- readRDS(rds_in)

# 分组列：timepoint（已确认）
group_col <- "timepoint"

# 固定显示顺序（与现有脚本保持一致）
# 注意：Seurat 对 [[ ]] 的方法可能返回 data.frame，这里直接操作 meta.data 列以确保是向量
nk.integrated@meta.data[[group_col]] <- factor(
  nk.integrated@meta.data[[group_col]],
  levels = c("0W_NCD", "1W_MCD", "2W_MCD", "6W_MCD")
)

# 统计每组细胞数并取最小值（等量降采样目标数）
cells <- colnames(nk.integrated)
meta  <- nk.integrated@meta.data
meta$.cell <- cells

count_df <- meta %>% count(.data[[group_col]], name = "n")
min_n <- min(count_df$n)
message(sprintf("Equalizing by '%s': sampling n=%d per group", group_col, min_n))

# 分组等量抽样，得到细胞ID（可重复性由 set.seed 保证）
sampled_cells <- meta %>%
  group_by(.data[[group_col]]) %>%
  slice_sample(n = min_n, replace = FALSE) %>%
  ungroup() %>%
  pull(.cell)

# 创建仅用于可视化的子集对象（不改变原对象与分析结果）
nk.vis <- subset(nk.integrated, cells = sampled_cells)

# 图1：按 timepoint 着色（核心交付，视觉均一化）
p_time <- DimPlot(
  nk.vis,
  reduction = "umap",
  group.by = group_col,
  pt.size = 0.6
) + theme_bw() + theme(legend.position = "right")

ggsave(
  filename = file.path(out_dir, "UMAP_equalized_by_timepoint.png"),
  plot = p_time,
  width = 10, height = 8, dpi = 500
)
ggsave(
  filename = file.path(out_dir, "UMAP_equalized_by_timepoint.pdf"),
  plot = p_time,
  width = 10, height = 8, device = "pdf"
)

# 图2：按 seurat_clusters 着色（可选对照）
p_clu <- DimPlot(
  nk.vis,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.6
) + theme_bw()

ggsave(
  filename = file.path(out_dir, "UMAP_equalized_by_cluster.png"),
  plot = p_clu,
  width = 10, height = 8, dpi = 500
)
ggsave(
  filename = file.path(out_dir, "UMAP_equalized_by_cluster.pdf"),
  plot = p_clu,
  width = 10, height = 8, device = "pdf"
)

# 分面图：与原图一致的布局（按 cluster 着色，按 timepoint 分面），但使用等量降采样后的 nk.vis
p_split_eq <- DimPlot(
  nk.vis,
  reduction = "umap",
  group.by = "seurat_clusters",
  split.by = "timepoint",
  label = TRUE,
  repel = TRUE,
  ncol = 2
) + theme(strip.text.x = element_text(size = 12))

ggsave(
  filename = file.path(out_dir, "UMAP_equalized_Clusters_by_Timepoint.png"),
  plot = p_split_eq,
  width = 12, height = 10, dpi = 500
)
ggsave(
  filename = file.path(out_dir, "UMAP_equalized_Clusters_by_Timepoint.pdf"),
  plot = p_split_eq,
  width = 12, height = 10, device = "pdf"
)

message("Equalized UMAP plots saved to: ", out_dir)

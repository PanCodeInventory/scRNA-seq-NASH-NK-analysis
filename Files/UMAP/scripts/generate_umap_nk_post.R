# 仅运行新增部分：从 RDS 加载整合对象，导出按时间点的独立 UMAP 图与簇组成趋势图
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

options(stringsAsFactors = FALSE)

# 路径
plots_dir <- "/home/harry/NASH/scRNA-seq/Files/UMAP/Results/plots"
rds_dir   <- "/home/harry/NASH/scRNA-seq/Files/UMAP/Results/rds"
if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

# 加载已保存的整合对象
rds_path <- file.path(rds_dir, "nk.integrated.rds")
message("[post] Loading integrated object from: ", rds_path)
if (!file.exists(rds_path)) {
  stop("[post] RDS 未找到：", rds_path, "\n请先运行 generate_umap_nk.R 以生成 nk.integrated.rds")
}
nk.integrated <- readRDS(rds_path)

# 规范 timepoint 因子顺序
desired_levels <- c("0W_NCD","1W_MCD","2W_MCD","6W_MCD")
if (!all(desired_levels %in% unique(nk.integrated$timepoint))) {
  warning("[post] RDS 中的 timepoint 与预期不一致，实际包含：",
          paste(sort(unique(as.character(nk.integrated$timepoint))), collapse = ", "))
}
nk.integrated$timepoint <- factor(nk.integrated$timepoint, levels = desired_levels)

# 1) 按时间点分别导出独立 UMAP 图（PNG/PDF/SVG）
message("[post] Exporting per-timepoint UMAPs...")
tps <- levels(nk.integrated$timepoint)
for (tp in tps) {
  if (is.na(tp)) next
  obj_tp <- subset(nk.integrated, subset = timepoint == tp)
  p_tp <- DimPlot(
    obj_tp,
    reduction = "umap",
    group.by = "seurat_clusters",
    label = TRUE,
    repel = TRUE
  ) + ggtitle(paste("NK Cells UMAP -", tp))
  tp_safe <- gsub("[^A-Za-z0-9_-]", "_", tp)
  out_tp <- file.path(plots_dir, paste0("UMAP_NK_Clusters_", tp_safe))
  ggsave(paste0(out_tp, ".png"), p_tp, width = 10, height = 8, dpi = 500)
  ggsave(paste0(out_tp, ".pdf"), p_tp, width = 10, height = 8, device = "pdf")
  ggsave(paste0(out_tp, ".svg"), p_tp, width = 10, height = 8, device = "svg")
}
message("[post] Per-timepoint UMAPs saved to: ", plots_dir)

# 2) 不同簇在不同样品（时间点）之间的变化趋势图（折线）
message("[post] Computing cluster composition trends...")
cluster_df <- nk.integrated@meta.data %>%
  dplyr::select(seurat_clusters, timepoint) %>%
  dplyr::mutate(
    seurat_clusters = as.factor(seurat_clusters),
    timepoint = factor(timepoint, levels = desired_levels)
  ) %>%
  dplyr::count(timepoint, seurat_clusters, name = "n") %>%
  dplyr::group_by(timepoint) %>%
  dplyr::mutate(percent = 100 * n / sum(n)) %>%
  dplyr::ungroup()

# 导出簇组成比例趋势（CSV：长表与宽表）
data_dir <- file.path(dirname(plots_dir), "data")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

# 长表：timepoint × seurat_clusters × n × percent
cluster_trends_long <- cluster_df %>%
  dplyr::arrange(timepoint, seurat_clusters)
utils::write.csv(cluster_trends_long,
                 file = file.path(data_dir, "NK_Cluster_Composition_Trends_long.csv"),
                 row.names = FALSE)

# 宽表：每行一个簇，每列一个时间点，值为百分比
cluster_trends_wide <- cluster_df %>%
  dplyr::select(seurat_clusters, timepoint, percent) %>%
  tidyr::pivot_wider(names_from = timepoint, values_from = percent) %>%
  dplyr::arrange(as.numeric(as.character(seurat_clusters)))
utils::write.csv(cluster_trends_wide,
                 file = file.path(data_dir, "NK_Cluster_Composition_Trends_wide.csv"),
                 row.names = FALSE)

message("[post] Cluster composition trend CSVs saved to: ", data_dir)

p_trend <- ggplot(cluster_df, aes(x = timepoint, y = percent, group = seurat_clusters, color = seurat_clusters)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  labs(title = "NK Cluster Composition Trends Across Timepoints",
       x = "Timepoint", y = "Percentage of cells") +
  theme_bw() +
  theme(legend.position = "right")

trend_out <- file.path(plots_dir, "UMAP_NK_Cluster_Trends")
ggsave(paste0(trend_out, ".png"), p_trend, width = 12, height = 8, dpi = 500)
ggsave(paste0(trend_out, ".pdf"), p_trend, width = 12, height = 8, device = "pdf")
ggsave(paste0(trend_out, ".svg"), p_trend, width = 12, height = 8, device = "svg")

message("[post] Trend plots saved to: ", plots_dir)

# 新增：堆叠柱状图（按 timepoint 汇总各簇百分比）
p_bar_stacked <- ggplot(cluster_df, aes(x = timepoint, y = percent, fill = seurat_clusters)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "NK Cluster Composition (Stacked Bar)",
       x = "Timepoint", y = "Percentage of cells") +
  theme_bw() +
  theme(legend.position = "right")

bar_out <- file.path(plots_dir, "UMAP_NK_Cluster_Composition_StackedBar")
ggsave(paste0(bar_out, ".png"), p_bar_stacked, width = 12, height = 8, dpi = 500)
ggsave(paste0(bar_out, ".pdf"), p_bar_stacked, width = 12, height = 8, device = "pdf")
ggsave(paste0(bar_out, ".svg"), p_bar_stacked, width = 12, height = 8, device = "svg")

# 新增：饼图（每个 timepoint 一个），在切片中标注比例，并聚合为一张分面图
for (tp in desired_levels) {
  df_tp <- dplyr::filter(cluster_df, timepoint == tp)
  p_pie <- ggplot(df_tp, aes(x = "", y = percent, fill = seurat_clusters)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar(theta = "y") +
    theme_void() +
    labs(title = paste("Cluster Composition -", tp)) +
    theme(legend.position = "right") +
    geom_label(
              data = dplyr::filter(df_tp, as.character(seurat_clusters) %in% c("0","2")),
              aes(label = paste0(round(percent, 1), "%")),
              position = position_stack(vjust = 0.5),
              size = 3,
              label.size = 0.25,
              color = "black",
              fill = "white"
            )
  tp_safe <- gsub("[^A-Za-z0-9_-]", "_", tp)
  pie_out <- file.path(plots_dir, paste0("UMAP_NK_Cluster_Composition_Pie_", tp_safe))
  ggsave(paste0(pie_out, ".png"), p_pie, width = 8, height = 8, dpi = 500)
  ggsave(paste0(pie_out, ".pdf"), p_pie, width = 8, height = 8, device = "pdf")
  ggsave(paste0(pie_out, ".svg"), p_pie, width = 8, height = 8, device = "svg")
}

# 聚合饼图（Demo 风格：高亮 0/2 黑色描边，其余 50% 透明；仅高亮内置比例标签；横向一行）
cluster_df2 <- cluster_df %>%
  dplyr::mutate(is_highlight = as.character(seurat_clusters) %in% c("0","2"))

p_pie_combined_demo <- ggplot(cluster_df2, aes(x = "", y = percent, fill = seurat_clusters)) +
  # 底层：按高亮控制透明度
  geom_bar(stat = "identity", width = 1, aes(alpha = is_highlight), color = NA) +
  scale_alpha_manual(values = c("FALSE" = 0.5, "TRUE" = 1)) +
  guides(alpha = "none") +
  # 高亮描边：仅对 0/2 加黑色边框
  geom_bar(stat = "identity", width = 1, fill = NA, color = "black", linewidth = 0.6,
           data = dplyr::filter(cluster_df2, is_highlight)) +
  coord_polar(theta = "y") +
  facet_wrap(~ timepoint, ncol = 4) +
  theme_void() +
  labs(title = "Cluster Composition by Timepoint (Pie, Demo Style)") +
  theme(legend.position = "right") +
  # 高亮内置比例标签
  geom_label(
    data = dplyr::filter(cluster_df2, is_highlight),
    aes(label = paste0(round(percent, 1), "%")),
    position = position_stack(vjust = 0.5),
    size = 3,
    label.size = 0.25,
    color = "black",
    fill = "white"
  )

pie_demo_out <- file.path(plots_dir, "UMAP_NK_Cluster_Composition_Pie_Combined_DemoStyle")
ggsave(paste0(pie_demo_out, ".png"), p_pie_combined_demo, width = 12, height = 10, dpi = 500)
ggsave(paste0(pie_demo_out, ".pdf"), p_pie_combined_demo, width = 12, height = 10, device = "pdf")
ggsave(paste0(pie_demo_out, ".svg"), p_pie_combined_demo, width = 12, height = 10, device = "svg")

message("[post] Composition stacked bar and pie plots saved to: ", plots_dir)
message("[post] Done.")

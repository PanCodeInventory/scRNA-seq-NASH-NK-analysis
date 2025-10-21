# 为NASH各时间点NK细胞生成UMAP分割图（按cluster着色，按timepoint分面）
# 环境：R + Seurat + dplyr + ggplot2 + patchwork
# 输入：四个已质控的10x目录（包含 barcodes.tsv.gz、features.tsv.gz、matrix.mtx.gz）
# 输出：关键图 p_split 的 PNG/SVG 两种格式到 Results/plots；本脚本保存于 Results/scripts

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(future)
})

options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 8 * 1024^3) # 提高future导出对象大小上限到8GB
future::plan("sequential")                   # 顺序执行，避免全局导出体积过大
set.seed(1234)

# 路径设置
base_dir <- "/home/harry/NASH/scRNA-seq/Files/Filter Files"
plots_dir <- "/home/harry/NASH/scRNA-seq/Files/UMAP/Results/plots"
scripts_dir <- "/home/harry/NASH/scRNA-seq/Files/UMAP/Results/scripts"
rds_dir <- "/home/harry/NASH/scRNA-seq/Files/UMAP/Results/rds"

# 确保输出目录存在
if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(scripts_dir)) dir.create(scripts_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(rds_dir)) dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)

# 样本与时间点映射（按文档要求）
sample_info <- data.frame(
  dir = file.path(base_dir, c("NCD_NK1.1", "MCD-1W_NK1.1", "MCD-2W_NK1.1", "MCD-6W_NK1.1")),
  timepoint = c("0W_NCD", "1W_MCD", "2W_MCD", "6W_MCD"),
  stringsAsFactors = FALSE
)

# 安全读取10x数据：自动判断features/genes文件的列数，选择合适的gene.column（避免“undefined columns selected”）
read10x_safe <- function(data_dir) {
  features_path <- file.path(data_dir, "features.tsv.gz")
  genes_path <- file.path(data_dir, "genes.tsv.gz")
  target_path <- if (file.exists(features_path)) features_path else if (file.exists(genes_path)) genes_path else NA
  if (is.na(target_path)) {
    stop(sprintf("目录缺少features.tsv.gz/genes.tsv.gz：%s", data_dir))
  }
  # 读取首行判断列数
  first_line <- tryCatch({
    con <- gzfile(target_path, "rt")
    on.exit(close(con), add = TRUE)
    readLines(con, n = 1)
  }, error = function(e) "")
  ncols <- if (length(first_line) == 1 && nzchar(first_line)) length(strsplit(first_line, "\t", fixed = TRUE)[[1]]) else 1
  gene_col <- if (ncols >= 2) 2 else 1
  message(sprintf("Read10X: %s using gene.column=%d (detected %d columns)", basename(target_path), gene_col, ncols))
  Read10X(data.dir = data_dir, gene.column = gene_col)
}

# 读取10x数据并创建Seurat对象（直接添加timepoint元数据；按要求跳过QC）
message("[1/5] 读取样本并创建Seurat对象...")
obj_list <- lapply(seq_len(nrow(sample_info)), function(i) {
  data_dir <- sample_info$dir[i]
  if (!dir.exists(data_dir)) {
    stop(sprintf("找不到数据目录：%s", data_dir))
  }
  counts <- read10x_safe(data_dir)
  obj <- CreateSeuratObject(counts = counts, project = "NK", min.cells = 0, min.features = 0)
  obj$timepoint <- sample_info$timepoint[i]
  obj
})
names(obj_list) <- sample_info$timepoint

# 使用SCT工作流进行整合（锚点法）
message("[2/5] SCTransform 归一化...")
obj_list <- lapply(obj_list, function(x) SCTransform(x, vst.flavor = "v2", verbose = FALSE))

message("[3/5] 特征选择与整合锚点计算...")
features <- SelectIntegrationFeatures(object.list = obj_list, nfeatures = 3000)
obj_list <- PrepSCTIntegration(object.list = obj_list, anchor.features = features, verbose = FALSE)
anchors <- FindIntegrationAnchors(object.list = obj_list, normalization.method = "SCT",
                                  anchor.features = features, verbose = FALSE)

message("[4/5] 数据整合...")
nk.integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT", verbose = FALSE)
DefaultAssay(nk.integrated) <- "integrated"

# PCA、邻居图、聚类、UMAP
message("[5/5] PCA / 邻居 / 聚类 / UMAP...")
nk.integrated <- RunPCA(nk.integrated, verbose = FALSE)
# 可选：ElbowPlot 用于选择PC，若需要保存可取消注释下两行
# ggsave(file.path(plots_dir, "ElbowPlot.pdf"), ElbowPlot(nk.integrated), width = 6, height = 4, device = "pdf")
ggsave(file.path(plots_dir, "ElbowPlot.png"), ElbowPlot(nk.integrated), width = 6, height = 4, dpi = 300)

nk.integrated <- FindNeighbors(nk.integrated, dims = 1:10)
nk.integrated <- FindClusters(nk.integrated, resolution = 0.3)
nk.integrated <- RunUMAP(nk.integrated, dims = 1:10)
nk.integrated$timepoint <- factor(nk.integrated$timepoint, levels = c("0W_NCD","1W_MCD","2W_MCD","6W_MCD"))

# 保存整合对象为 RDS，便于后续仅运行新增代码
saveRDS(nk.integrated, file = file.path(rds_dir, "nk.integrated.rds"))

# 总览图（打印显示，按需要可保存）
p_overview <- DimPlot(
  nk.integrated,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  label.size = 5,
  repel = TRUE
) + ggtitle("All NK Cells Clustered")
# print(p_overview) # 交互式环境中可取消注释以查看

# 关键分割图：按时间点分割（split.by = timepoint），细胞按cluster着色
p_split <- DimPlot(
  nk.integrated,
  reduction = "umap",
  group.by = "seurat_clusters",
  split.by = "timepoint",
  label = TRUE,
  repel = TRUE,
  ncol = 2
) + theme(strip.text.x = element_text(size = 12))

# 保存出版级别图件（PNG/PDF/SVG）
output_filename <- file.path(plots_dir, "UMAP_NK_Clusters_by_Timepoint")
plot_width <- 12
plot_height <- 10

ggsave(
  filename = paste0(output_filename, ".png"),
  plot = p_split,
  width = plot_width,
  height = plot_height,
  dpi = 500
)


ggsave(
  filename = paste0(output_filename, ".svg"),
  plot = p_split,
  width = plot_width,
  height = plot_height,
  device = "svg"
)

message("All plots have been saved in PNG and SVG formats at: ", plots_dir)

# 新增：按时间点分别导出独立 UMAP 图（PNG/PDF/SVG）
for (tp in c("0W_NCD","1W_MCD","2W_MCD","6W_MCD")) {
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
    ggsave(paste0(out_tp, ".svg"), p_tp, width = 10, height = 8, device = "svg")
}

# 新增：不同簇在不同样品（时间点）之间的变化趋势图（折线）
cluster_df <- nk.integrated@meta.data %>%
  dplyr::select(seurat_clusters, timepoint) %>%
  dplyr::mutate(
    seurat_clusters = as.factor(seurat_clusters),
    timepoint = factor(timepoint, levels = c("0W_NCD","1W_MCD","2W_MCD","6W_MCD"))
  ) %>%
  dplyr::count(timepoint, seurat_clusters, name = "n") %>%
  dplyr::group_by(timepoint) %>%
  dplyr::mutate(percent = 100 * n / sum(n)) %>%
  dplyr::ungroup()

p_trend <- ggplot(cluster_df, aes(x = timepoint, y = percent, group = seurat_clusters, color = seurat_clusters)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(title = "NK Cluster Composition Trends Across Timepoints",
       x = "Timepoint", y = "Percentage of cells") +
  theme_bw() +
  theme(legend.position = "right")

trend_out <- file.path(plots_dir, "UMAP_NK_Cluster_Trends")
ggsave(paste0(trend_out, ".png"), p_trend, width = 12, height = 8, dpi = 500)
ggsave(paste0(trend_out, ".svg"), p_trend, width = 12, height = 8, device = "svg")

# 新增：堆叠柱状图（按 timepoint 汇总各簇百分比）
p_bar_stacked <- ggplot(cluster_df, aes(x = timepoint, y = percent, fill = seurat_clusters)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "NK Cluster Composition (Stacked Bar)",
       x = "Timepoint", y = "Percentage of cells") +
  theme_bw() +
  theme(legend.position = "right")

bar_out <- file.path(plots_dir, "UMAP_NK_Cluster_Composition_StackedBar")
ggsave(paste0(bar_out, ".png"), p_bar_stacked, width = 12, height = 8, dpi = 500)
ggsave(paste0(bar_out, ".svg"), p_bar_stacked, width = 12, height = 8, device = "svg")

message("Additional per-timepoint UMAP, trend, and stacked bar plots have been saved at: ", plots_dir)

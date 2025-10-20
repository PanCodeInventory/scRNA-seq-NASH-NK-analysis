# 专用于探索最佳 dims 与 resolution 的调参脚本（与正式产出隔离）
# 环境：R + Seurat + dplyr + ggplot2 (+ 可选 cluster)
# 输入：已整合对象 RDS（默认：Files/Filter Files/RDS/nk.integrated.rds）
# 输出：全部写入 Files/UMAP/Results/tuning 下的 plots/、data/、logs/
# 用法：直接运行（根据配置修改网格大小以控制时间/内存）

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  # 可选：用于 silhouette 指标；未安装将自动跳过该指标
  # install.packages("cluster") 如需启用
})

options(stringsAsFactors = FALSE)
set.seed(1234)

# ---------------------------
# 1) 配置
# ---------------------------
# 对象路径（默认使用现有整合对象，避免重复整合）
rds_path <- "Files/Filter Files/RDS/nk.integrated.rds"

# 搜索网格（建议先小网格验证，再逐步加密）
dims_grid <- c(10, 15, 20, 25, 30)          # 主成分个数候选
res_grid  <- seq(0.2, 1.2, by = 0.1)        # 分辨率候选

# 输出目录（完全隔离）
tune_base  <- "Files/UMAP/Results/tuning"
tune_plots <- file.path(tune_base, "plots")
tune_data  <- file.path(tune_base, "data")
tune_logs  <- file.path(tune_base, "logs")
dir.create(tune_plots, showWarnings = FALSE, recursive = TRUE)
dir.create(tune_data,  showWarnings = FALSE, recursive = TRUE)
dir.create(tune_logs,  showWarnings = FALSE, recursive = TRUE)

# 其它参数
timepoint_levels <- c("0W_NCD","1W_MCD","2W_MCD","6W_MCD")
label_size <- 4
dpi_high   <- 500

# 工具函数：分辨率字符串用于文件名（0.6 -> 0p6）
safe_res <- function(x) gsub("\\.", "p", as.character(x))

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

# ---------------------------
# 2) 读取对象与基本准备
# ---------------------------
if (!file.exists(rds_path)) {
  stop(sprintf("找不到整合对象 RDS：%s\n请先在正式脚本中生成 nk.integrated.rds", rds_path))
}
nk <- readRDS(rds_path)

# 确保 timepoint 因子顺序
if (!"timepoint" %in% colnames(nk@meta.data)) {
  stop("对象缺少元数据列 'timepoint'")
}
nk$timepoint <- factor(nk$timepoint, levels = timepoint_levels)

# ---------------------------
# 3) 调参主循环：对每个 dims 只运行一次 Neighbors/UMAP
#    然后对每个 res 运行 FindClusters 并记录指标
# ---------------------------
metrics_list <- list()
has_cluster_pkg <- requireNamespace("cluster", quietly = TRUE)

message(sprintf("开始调参：dims_grid=%s; res_grid=[%.1f..%.1f] by 0.1",
                paste(dims_grid, collapse=","), min(res_grid), max(res_grid)))
message(sprintf("silhouette 指标：%s", if (has_cluster_pkg) "启用" else "未安装 'cluster' 包，跳过"))

for (dims in dims_grid) {
  message(sprintf("[dims=%d] 运行 FindNeighbors / RunUMAP", dims))
  nk <- FindNeighbors(nk, dims = 1:dims, verbose = FALSE)
  nk <- RunUMAP(nk, dims = 1:dims, verbose = FALSE)

  for (res in res_grid) {
    message(sprintf("  └─ FindClusters(res=%.1f)", res))
    nk <- FindClusters(nk, resolution = res, verbose = FALSE)

    # Seurat 会在 meta.data 生成 integrated_snn_res.{res}
    col_res <- paste0("integrated_snn_res.", as.character(res))
    if (!col_res %in% colnames(nk@meta.data)) {
      warning(sprintf("未找到簇列：%s（可能是版本/命名差异）", col_res))
      next
    }

    # 3.1 基础统计（全局构成）
    comp <- nk@meta.data %>%
      dplyr::count(!!rlang::sym(col_res), name = "n") %>%
      dplyr::mutate(frac = n / sum(n))
    n_clusters <- nrow(comp)
    min_frac   <- if (n_clusters > 0) min(comp$frac) else NA_real_
    q95_frac   <- if (n_clusters > 0) as.numeric(quantile(comp$frac, 0.95)) else NA_real_

    # 3.2 轮廓系数（在 PCA 空间；需要 cluster 包；至少两类）
    median_sil <- NA_real_
    pct_sil_neg <- NA_real_
    if (has_cluster_pkg && n_clusters >= 2) {
      emb <- tryCatch({
        Embeddings(nk, "pca")[, 1:dims, drop = FALSE]
      }, error = function(e) NULL)
      cls <- as.integer(factor(nk@meta.data[[col_res]]))
      if (!is.null(emb) && ncol(emb) == dims) {
        # 距离矩阵可用欧氏距离
        dist_mat <- tryCatch(stats::dist(emb), error = function(e) NULL)
        if (!is.null(dist_mat)) {
          sil <- tryCatch(cluster::silhouette(cls, dist_mat), error = function(e) NULL)
          if (!is.null(sil)) {
            sw <- sil[, "sil_width"]
            median_sil  <- stats::median(sw, na.rm = TRUE)
            pct_sil_neg <- mean(sw < 0, na.rm = TRUE) * 100
          }
        }
      }
    }

    # 3.3 记录
    metrics_list[[length(metrics_list) + 1]] <- data.frame(
      dims = dims,
      res  = res,
      n_clusters = n_clusters,
      min_cluster_frac = min_frac,
      q95_cluster_frac = q95_frac,
      median_silhouette = median_sil,
      pct_silhouette_negative = pct_sil_neg,
      stringsAsFactors = FALSE
    )
  }
}

df_metrics <- dplyr::bind_rows(metrics_list)

# 写出指标表
metrics_csv <- file.path(tune_data, "nk_tuning_metrics.csv")
utils::write.csv(df_metrics, metrics_csv, row.names = FALSE)
message("已保存指标表：", metrics_csv)

# ---------------------------
# 4) 可视化：热力图（median_silhouette / n_clusters）
# ---------------------------
if (nrow(df_metrics) > 0) {
  df_plot <- df_metrics %>%
    dplyr::mutate(res_lab = res, dims = as.factor(dims))

  # 4.1 median silhouette 热力图
  p_hm_sil <- ggplot(df_plot, aes(x = res_lab, y = dims, fill = median_silhouette)) +
    geom_tile(color = "white") +
    scale_fill_viridis_c(option = "C", na.value = "#f0f0f0") +
    labs(title = "Median Silhouette across dims x res",
         x = "resolution", y = "dims") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(tune_plots, "heatmap_median_silhouette_dims_vs_res.png"), p_hm_sil, width = 10, height = 6, dpi = dpi_high)
  ggsave(file.path(tune_plots, "heatmap_median_silhouette_dims_vs_res.pdf"), p_hm_sil, width = 10, height = 6, device = "pdf")

  # 4.2 n_clusters 热力图
  p_hm_nc <- ggplot(df_plot, aes(x = res_lab, y = dims, fill = n_clusters)) +
    geom_tile(color = "white") +
    scale_fill_viridis_c(option = "B", na.value = "#f0f0f0") +
    labs(title = "Number of Clusters across dims x res",
         x = "resolution", y = "dims") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(tune_plots, "heatmap_n_clusters_dims_vs_res.png"), p_hm_nc, width = 10, height = 6, dpi = dpi_high)
  ggsave(file.path(tune_plots, "heatmap_n_clusters_dims_vs_res.pdf"), p_hm_nc, width = 10, height = 6, device = "pdf")
}

# ---------------------------
# 5) 选择候选组合并出少量 UMAP 分面图
#    规则示例：簇数介于 6-20、最小簇 >= 0.5%、按 silhouette 降序与负轮廓占比升序
# ---------------------------
select_candidates <- function(df, per_dims_top = 1, global_top = 6) {
  df2 <- df %>%
    dplyr::filter(!is.na(median_silhouette)) %>%
    dplyr::filter(n_clusters >= 6, n_clusters <= 20) %>%
    dplyr::filter(is.na(min_cluster_frac) | min_cluster_frac >= 0.005) %>%
    dplyr::arrange(dplyr::desc(median_silhouette),
                   dplyr::across(pct_silhouette_negative, ~ ifelse(is.na(.x), Inf, .x)),
                   dims)

  best_per_dims <- df2 %>%
    dplyr::group_by(dims) %>%
    dplyr::slice_head(n = per_dims_top) %>%
    dplyr::ungroup()

  best_global <- df2 %>%
    dplyr::slice_head(n = global_top)

  list(best_per_dims = best_per_dims, best_global = best_global)
}

cands <- select_candidates(df_metrics, per_dims_top = 1, global_top = 6)

best_per_dims_csv <- file.path(tune_data, "nk_tuning_best_per_dims.csv")
utils::write.csv(cands$best_per_dims, best_per_dims_csv, row.names = FALSE)
best_global_csv <- file.path(tune_data, "nk_tuning_top_candidates.csv")
utils::write.csv(cands$best_global, best_global_csv, row.names = FALSE)
message("已保存候选组合：", best_per_dims_csv, " ; ", best_global_csv)

# 准备需要出图的 (dims,res) 组合（去重合并）
pairs_to_plot <- unique(rbind(
  cands$best_per_dims %>% dplyr::select(dims, res),
  cands$best_global %>% dplyr::select(dims, res)
))
if (nrow(pairs_to_plot) > 0) {
  message(sprintf("计划出图 %d 个候选组合的 UMAP 分面图", nrow(pairs_to_plot)))
}

# ---------------------------
# 6) 为候选组合出 UMAP 分面图（按 timepoint 分面）
#    每个 dims 切换时重跑 Neighbors/UMAP；每个 res 重跑 FindClusters
# ---------------------------
last_dims <- NA_integer_
for (i in seq_len(nrow(pairs_to_plot))) {
  dims_i <- pairs_to_plot$dims[i]
  res_i  <- pairs_to_plot$res[i]

  if (is.na(last_dims) || dims_i != last_dims) {
    message(sprintf("[绘图] 重新计算 UMAP：dims=%d", dims_i))
    nk <- FindNeighbors(nk, dims = 1:dims_i, verbose = FALSE)
    nk <- RunUMAP(nk, dims = 1:dims_i, verbose = FALSE)
    last_dims <- dims_i
  }

  message(sprintf("[绘图] FindClusters(res=%.1f) 并生成分面图", res_i))
  nk <- FindClusters(nk, resolution = res_i, verbose = FALSE)
  col_res <- paste0("integrated_snn_res.", as.character(res_i))

  # 分面 UMAP（按 timepoint）
  p_split <- DimPlot(
    nk,
    reduction = "umap",
    group.by = col_res,
    split.by = "timepoint",
    label = TRUE,
    repel = TRUE,
    label.size = label_size,
    ncol = 2
  ) + theme(strip.text.x = element_text(size = 12)) +
      ggtitle(sprintf("UMAP tuning | dims=%d, res=%.1f", dims_i, res_i))

  fname_core <- sprintf("UMAP_tuning_dims%d_res%s_byTimepoint", dims_i, safe_res(res_i))
  ggsave(file.path(tune_plots, paste0(fname_core, ".png")), p_split, width = 12, height = 10, dpi = dpi_high)
  ggsave(file.path(tune_plots, paste0(fname_core, ".pdf")), p_split, width = 12, height = 10, device = "pdf")
}

# ---------------------------
# 7) 日志与运行配置快照
# ---------------------------
# 运行配置（简单文本以避免额外依赖）
cfg_file <- file.path(tune_logs, paste0("run_config_", timestamp, ".txt"))
cfg_lines <- c(
  sprintf("rds_path: %s", rds_path),
  sprintf("dims_grid: %s", paste(dims_grid, collapse=",")),
  sprintf("res_grid: from %.1f to %.1f by 0.1", min(res_grid), max(res_grid)),
  sprintf("timepoint_levels: %s", paste(timepoint_levels, collapse=",")),
  sprintf("seed: %d", 1234),
  sprintf("metrics_csv: %s", metrics_csv),
  sprintf("best_per_dims_csv: %s", best_per_dims_csv),
  sprintf("best_global_csv: %s", best_global_csv)
)
writeLines(cfg_lines, cfg_file)
message("已保存运行配置：", cfg_file)

# sessionInfo
sess_file <- file.path(tune_logs, paste0("sessionInfo_", timestamp, ".txt"))
utils::capture.output(utils::sessionInfo(), file = sess_file)
message("已保存 sessionInfo：", sess_file)

message("调参流程完成。所有产物已写入：", tune_base)

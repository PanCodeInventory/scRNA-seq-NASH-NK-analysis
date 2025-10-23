#!/usr/bin/env Rscript

# 基于"无Cluster6"对象的 dims × resolution 调参与最终降维/聚类/UMAP
# 输入：2_DataProcessing/RDS/nk.integrated.singleR_annotated.noCluster6.rds
# 产出：
#   - 指标表/候选组合：2_DataProcessing/3_UMAP-Tuning/data/*
#   - 热力图/候选 UMAP：2_DataProcessing/3_UMAP-Tuning/plots/*
#   - 最终选择参数快照：2_DataProcessing/3_UMAP-Tuning/logs/selected_params.txt
#   - 最终对象与图：2_DataProcessing/RDS/nk.integrated.singleR_annotated.noCluster6.tuned.rds
#                   2_DataProcessing/3_UMAP-Tuning/plots/UMAP_noCluster6_final_byTimepoint.png
#
# 用法：
#   Rscript 2_DataProcessing/Scripts/tune_noCluster6_dims_resolution.R

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(utils)
})

options(stringsAsFactors = FALSE)
set.seed(1234)

# ---------------------------
# 1) 配置
# ---------------------------
rds_path <- "2_DataProcessing/RDS/nk.integrated.singleR_annotated.noCluster6.rds"

dims_grid <- c(10, 15, 20, 25, 30)          # 主成分个数候选
res_grid  <- seq(0.2, 1.2, by = 0.1)        # 分辨率候选

tune_base  <- "2_DataProcessing/3_UMAP-Tuning"
tune_plots <- file.path(tune_base, "plots")
tune_data  <- file.path(tune_base, "data")
tune_logs  <- file.path(tune_base, "logs")
dir.create(tune_plots, showWarnings = FALSE, recursive = TRUE)
dir.create(tune_data,  showWarnings = FALSE, recursive = TRUE)
dir.create(tune_logs,  showWarnings = FALSE, recursive = TRUE)

timepoint_levels <- c("0W_NCD","1W_MCD","2W_MCD","6W_MCD")
label_size <- 4
dpi_high   <- 500
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

safe_res <- function(x) gsub("\\.", "p", as.character(x))

# ---------------------------
# 2) 读取对象与确保 PCA
# ---------------------------
if (!file.exists(rds_path)) {
  stop(sprintf("找不到输入 RDS：%s", rds_path))
}
nk <- readRDS(rds_path)

# 确保 timepoint 因子顺序（如存在）
if ("timepoint" %in% colnames(nk@meta.data)) {
  nk$timepoint <- factor(nk$timepoint, levels = timepoint_levels)
}

# 选择标签列
if ("SingleR.pruned.labels" %in% colnames(nk@meta.data)) {
  label_col <- "SingleR.pruned.labels"
} else if ("SingleR.labels" %in% colnames(nk@meta.data)) {
  label_col <- "SingleR.labels"
} else if ("singleR.pruned.labels" %in% colnames(nk@meta.data)) {
  label_col <- "singleR.pruned.labels"
} else if ("singleR.labels" %in% colnames(nk@meta.data)) {
  label_col <- "singleR.labels"
} else {
  label_col <- NULL
}

# 确保 PCA 存在且至少覆盖 max(dims_grid)
ensure_pca <- function(obj, target_dims = max(dims_grid)) {
  has_pca <- "pca" %in% names(obj@reductions)
  if (has_pca) {
    emb <- tryCatch(Embeddings(obj, "pca"), error = function(e) NULL)
    if (!is.null(emb) && ncol(emb) >= target_dims) {
      message(sprintf("检测到现有 PCA（%d PCs），跳过重算。", ncol(emb)))
      return(obj)
    }
  }
  # 选择用于降维的 assay：优先 integrated 有可变基因，然后 SCT，最后 RNA
  assays_avail <- tryCatch(Assays(obj), error = function(e) character(0))
  dr_assay <- if ("integrated" %in% assays_avail) "integrated" 
              else if ("SCT" %in% assays_avail) "SCT" 
              else "RNA"
  
  if (dr_assay == "integrated") {
    vf <- tryCatch(VariableFeatures(obj), error = function(e) character(0))
    if (length(vf) == 0) dr_assay <- "SCT"
  }
  if (dr_assay == "SCT") {
    vf <- tryCatch(VariableFeatures(obj), error = function(e) character(0))
    if (length(vf) == 0) dr_assay <- "RNA"
  }
  
  DefaultAssay(obj) <- dr_assay
  if (dr_assay == "RNA") {
    message("PCA 保障：使用 RNA 路线（Normalize/FindVariableFeatures/Scale/RunPCA）")
    obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
    if (length(VariableFeatures(obj)) == 0) {
      obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
    }
    feats <- VariableFeatures(obj)
  } else {
    message(sprintf("PCA 保障：使用 %s 路线（VariableFeatures/Scale/RunPCA）", dr_assay))
    feats <- VariableFeatures(obj)
  }
  obj <- ScaleData(obj, features = feats, verbose = FALSE)
  obj <- RunPCA(obj, features = feats, npcs = max(50, target_dims), verbose = FALSE)
  return(obj)
}

nk <- ensure_pca(nk, target_dims = max(dims_grid))

# ---------------------------
# 3) 调参主循环：每个 dims 运行一次 Neighbors/UMAP，然后对每个 res 运行 FindClusters
# ---------------------------
metrics_list <- list()
has_cluster_pkg <- requireNamespace("cluster", quietly = TRUE)

message(sprintf("开始调参：dims_grid=%s; res_grid=[%.1f..%.1f] by 0.1",
                paste(dims_grid, collapse=","), min(res_grid), max(res_grid)))
message(sprintf("silhouette 指标：%s", if (has_cluster_pkg) "启用" else "未安装 'cluster' 包，跳过"))

# 选择可用的图名称（优先使用 assay 相关的 SNN 图；若无则取第一个）
pick_graph_name <- function(obj) {
  g <- names(obj@graphs)
  if (length(g) == 0) stop("未找到邻居图（graphs 为空），请检查 FindNeighbors 是否成功。")
  if ("integrated_snn" %in% g) return("integrated_snn")
  if ("SCT_snn" %in% g) return("SCT_snn")
  if ("RNA_snn" %in% g) return("RNA_snn")
  return(g[[1]])
}
nk_dims_cached <- nk
last_dims <- NA_integer_

for (dims in dims_grid) {
  message(sprintf("[dims=%d] 运行 FindNeighbors / RunUMAP", dims))
  if (is.na(last_dims) || dims != last_dims) {
    nk_dims_cached <- FindNeighbors(nk_dims_cached, dims = 1:dims, verbose = FALSE)
    nk_dims_cached <- RunUMAP(nk_dims_cached, dims = 1:dims, verbose = FALSE)
    last_dims <- dims
  }

  for (res in res_grid) {
    message(sprintf("  └─ FindClusters(res=%.1f)", res))
    graph_nm <- pick_graph_name(nk_dims_cached)
    nk_tmp <- FindClusters(nk_dims_cached, graph.name = graph_nm, resolution = res, verbose = FALSE)

    # 使用通用 seurat_clusters 列
    if (!"seurat_clusters" %in% colnames(nk_tmp@meta.data)) {
      warning("未找到 'seurat_clusters' 列，跳过该组合")
      next
    }

    comp <- nk_tmp@meta.data %>%
      dplyr::count(seurat_clusters, name = "n") %>%
      dplyr::mutate(frac = n / sum(n))
    n_clusters <- nrow(comp)
    min_frac   <- if (n_clusters > 0) min(comp$frac) else NA_real_
    q95_frac   <- if (n_clusters > 0) as.numeric(quantile(comp$frac, 0.95)) else NA_real_

    # 轮廓系数（在 PCA 空间）
    median_sil <- NA_real_
    pct_sil_neg <- NA_real_
    if (has_cluster_pkg && n_clusters >= 2) {
      emb <- tryCatch(Embeddings(nk_tmp, "pca")[, 1:dims, drop = FALSE], error = function(e) NULL)
      cls <- as.integer(factor(nk_tmp$seurat_clusters))
      if (!is.null(emb) && ncol(emb) == dims) {
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
metrics_csv <- file.path(tune_data, "nk_noCluster6_tuning_metrics.csv")
utils::write.csv(df_metrics, metrics_csv, row.names = FALSE)
message("已保存指标表：", metrics_csv)

# ---------------------------
# 4) 可视化：热力图（median_silhouette / n_clusters）
# ---------------------------
if (nrow(df_metrics) > 0) {
  df_plot <- df_metrics %>%
    dplyr::mutate(res_lab = res, dims = as.factor(dims))

  p_hm_sil <- ggplot(df_plot, aes(x = res_lab, y = dims, fill = median_silhouette)) +
    geom_tile(color = "white") +
    scale_fill_viridis_c(option = "C", na.value = "#f0f0f0") +
    labs(title = "Median Silhouette across dims x res (noCluster6)",
         x = "resolution", y = "dims") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(tune_plots, "heatmap_median_silhouette_dims_vs_res_noCluster6.png"), p_hm_sil, width = 10, height = 6, dpi = dpi_high)
  ggsave(file.path(tune_plots, "heatmap_median_silhouette_dims_vs_res_noCluster6.pdf"), p_hm_sil, width = 10, height = 6, device = "pdf")

  p_hm_nc <- ggplot(df_plot, aes(x = res_lab, y = dims, fill = n_clusters)) +
    geom_tile(color = "white") +
    scale_fill_viridis_c(option = "B", na.value = "#f0f0f0") +
    labs(title = "Number of Clusters across dims x res (noCluster6)",
         x = "resolution", y = "dims") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(tune_plots, "heatmap_n_clusters_dims_vs_res_noCluster6.png"), p_hm_nc, width = 10, height = 6, dpi = dpi_high)
  ggsave(file.path(tune_plots, "heatmap_n_clusters_dims_vs_res_noCluster6.pdf"), p_hm_nc, width = 10, height = 6, device = "pdf")
}

# ---------------------------
# 5) 选择候选组合与出少量 UMAP 分面
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

best_per_dims_csv <- file.path(tune_data, "nk_noCluster6_tuning_best_per_dims.csv")
utils::write.csv(cands$best_per_dims, best_per_dims_csv, row.names = FALSE)
best_global_csv <- file.path(tune_data, "nk_noCluster6_tuning_top_candidates.csv")
utils::write.csv(cands$best_global, best_global_csv, row.names = FALSE)
message("已保存候选组合：", best_per_dims_csv, " ; ", best_global_csv)

pairs_to_plot <- unique(rbind(
  if (nrow(cands$best_per_dims) > 0) cands$best_per_dims %>% dplyr::select(dims, res) else data.frame(dims=integer(0), res=numeric(0)),
  if (nrow(cands$best_global) > 0)  cands$best_global %>% dplyr::select(dims, res) else data.frame(dims=integer(0), res=numeric(0))
))
if (nrow(pairs_to_plot) > 0) {
  message(sprintf("计划出图 %d 个候选组合的 UMAP 分面图", nrow(pairs_to_plot)))
}

last_dims_umap <- NA_integer_
for (i in seq_len(nrow(pairs_to_plot))) {
  dims_i <- pairs_to_plot$dims[i]
  res_i  <- pairs_to_plot$res[i]

  if (is.na(last_dims_umap) || dims_i != last_dims_umap) {
    message(sprintf("[绘图] 重新计算 UMAP：dims=%d", dims_i))
    nk <- FindNeighbors(nk, dims = 1:dims_i, verbose = FALSE)
    nk <- RunUMAP(nk, dims = 1:dims_i, verbose = FALSE)
    last_dims_umap <- dims_i
  }

  message(sprintf("[绘图] FindClusters(res=%.1f) 并生成分面图", res_i))
  graph_nm <- pick_graph_name(nk)
  nk <- FindClusters(nk, graph.name = graph_nm, resolution = res_i, verbose = FALSE)

  # 分面 UMAP（按 timepoint）
  if ("timepoint" %in% colnames(nk@meta.data)) {
    p_split <- DimPlot(
      nk,
      reduction = "umap",
      group.by = "seurat_clusters",
      split.by = "timepoint",
      label = TRUE,
      repel = TRUE,
      label.size = label_size,
      ncol = 2
    ) + theme(strip.text.x = element_text(size = 12)) +
        ggtitle(sprintf("UMAP tuning (noCluster6) | dims=%d, res=%.1f", dims_i, res_i))

    fname_core <- sprintf("UMAP_noCluster6_tuning_dims%d_res%s_byTimepoint", dims_i, safe_res(res_i))
    ggsave(file.path(tune_plots, paste0(fname_core, ".png")), p_split, width = 12, height = 10, dpi = dpi_high)
    ggsave(file.path(tune_plots, paste0(fname_core, ".pdf")), p_split, width = 12, height = 10, device = "pdf")
  }

  # 标签 UMAP（SingleR）
  if (!is.null(label_col) && label_col %in% colnames(nk@meta.data)) {
    p_lab <- DimPlot(
      nk,
      reduction = "umap",
      group.by = label_col,
      label = FALSE
    ) + ggtitle(sprintf("UMAP tuning (noCluster6) | dims=%d, res=%.1f | colored by %s", dims_i, res_i, label_col))
    fname_lab <- sprintf("UMAP_noCluster6_tuning_dims%d_res%s_bySingleR", dims_i, safe_res(res_i))
    ggsave(file.path(tune_plots, paste0(fname_lab, ".png")), p_lab, width = 10, height = 8, dpi = dpi_high)
  }
}

# ---------------------------
# 6) 选择最终参数并生成最终对象/图
# ---------------------------
dims_best <- NA_integer_
res_best  <- NA_real_

if (nrow(cands$best_global) > 0) {
  dims_best <- cands$best_global$dims[1]
  res_best  <- cands$best_global$res[1]
  message(sprintf("最终选择（来自 best_global Top1）：dims=%d, res=%.1f", dims_best, res_best))
} else if (nrow(cands$best_per_dims) > 0) {
  dims_best <- cands$best_per_dims$dims[1]
  res_best  <- cands$best_per_dims$res[1]
  message(sprintf("最终选择（来自 best_per_dims Top1）：dims=%d, res=%.1f", dims_best, res_best))
} else {
  # 回退策略
  dims_best <- 20
  res_best  <- 0.4
  warning(sprintf("未找到有效候选组合，使用回退参数：dims=%d, res=%.1f", dims_best, res_best))
}

# 保存最终参数
sel_file <- file.path(tune_logs, "selected_params_noCluster6.txt")
writeLines(c(
  sprintf("dims_best: %d", dims_best),
  sprintf("res_best: %.1f", res_best),
  sprintf("timestamp: %s", timestamp)
), sel_file)
message("已保存选择参数：", sel_file)

# 用最终参数重跑并保存对象与图
nk <- FindNeighbors(nk, dims = 1:dims_best, verbose = FALSE)
nk <- RunUMAP(nk, dims = 1:dims_best, verbose = FALSE)
graph_nm <- pick_graph_name(nk)
nk <- FindClusters(nk, graph.name = graph_nm, resolution = res_best, verbose = FALSE)

# 最终对象保存
out_rds_final <- "2_DataProcessing/RDS/nk.integrated.singleR_annotated.noCluster6.tuned.rds"
saveRDS(nk, out_rds_final)
message("已保存最终对象：", out_rds_final)

# 最终图：按 timepoint 分面
if ("timepoint" %in% colnames(nk@meta.data)) {
  p_final_tp <- DimPlot(
    nk,
    reduction = "umap",
    group.by = "seurat_clusters",
    split.by = "timepoint",
    label = TRUE,
    repel = TRUE,
    label.size = label_size,
    ncol = 2
  ) + theme(strip.text.x = element_text(size = 12)) +
      ggtitle(sprintf("UMAP (noCluster6 final) | dims=%d, res=%.1f | by timepoint", dims_best, res_best))
  ggsave(file.path(tune_plots, "UMAP_noCluster6_final_byTimepoint.png"), p_final_tp, width = 12, height = 10, dpi = dpi_high)
}

# 最终图：按 SingleR 注释
if (!is.null(label_col) && label_col %in% colnames(nk@meta.data)) {
  p_final_lab <- DimPlot(
    nk,
    reduction = "umap",
    group.by = label_col,
    label = FALSE
  ) + ggtitle(sprintf("UMAP (noCluster6 final) | dims=%d, res=%.1f | colored by %s", dims_best, res_best, label_col))
  ggsave(file.path(tune_plots, "UMAP_noCluster6_final_bySingleR.png"), p_final_lab, width = 10, height = 8, dpi = dpi_high)
}

# ---------------------------
# 7) 运行配置快照
# ---------------------------
cfg_file <- file.path(tune_logs, paste0("run_config_noCluster6_", timestamp, ".txt"))
cfg_lines <- c(
  sprintf("rds_path: %s", rds_path),
  sprintf("dims_grid: %s", paste(dims_grid, collapse=",")),
  sprintf("res_grid: from %.1f to %.1f by 0.1", min(res_grid), max(res_grid)),
  sprintf("timepoint_levels: %s", paste(timepoint_levels, collapse=",")),
  sprintf("seed: %d", 1234),
  sprintf("metrics_csv: %s", metrics_csv),
  sprintf("best_per_dims_csv: %s", best_per_dims_csv),
  sprintf("best_global_csv: %s", best_global_csv),
  sprintf("selected_params_file: %s", sel_file),
  sprintf("out_rds_final: %s", out_rds_final)
)
writeLines(cfg_lines, cfg_file)
message("已保存运行配置：", cfg_file)

# sessionInfo
sess_file <- file.path(tune_logs, paste0("sessionInfo_noCluster6_", timestamp, ".txt"))
utils::capture.output(utils::sessionInfo(), file = sess_file)
message("已保存 sessionInfo：", sess_file)

message("noCluster6 调参与最终降维流程完成。所有产物已写入：", tune_base)
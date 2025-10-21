# 去双胞与去除非NK污染的整合处理脚本
# 输入：Files/UMAP/Results/rds/nk.integrated.rds（整合后的Seurat对象）
# 输出：
#   - 清理后对象：Files/Doublet_Removed/RDS/nk.integrated.filtered.rds
#   - 双胞评分中间对象：Files/Doublet_Removed/RDS/nk.integrated.doublet_scored.rds
#   - 可视化图：Files/Doublet_Removed/plots/ 下若干PNG/SVG
#   - 报告：Files/Doublet_Removed/reports/cleaning_report.md
#
# 步骤：
#   1) 备份与依赖检查
#   2) scDblFinder 去双胞（按样本/时间点分组）
#   3) SingleR 自动注释 + UCell 签名评分
#   4) 细胞级与簇级污染剔除
#   5) 重跑降维/聚类/UMAP
#   6) 输出图件与报告
#
# 备注：脚本内含自动安装缺失依赖（CRAN/BiocManager），如需离线环境请事先安装相关包。

suppressPackageStartupMessages({
  # 基本包
  library(utils)
})

options(stringsAsFactors = FALSE)
set.seed(1234)
# 立即输出启动日志，便于诊断“日志空白”的问题
message(sprintf("SCRIPT START: %s", as.character(Sys.time())))
cat("SCRIPT START\n"); flush.console()

# 路径
root_dir <- "/home/harry/NASH/scRNA-seq"
in_rds <- file.path(root_dir, "Files/UMAP/Results/rds", "nk.integrated.rds")
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
ensure_pkg("dplyr", bioc = FALSE)
ensure_pkg("ggplot2", bioc = FALSE)
ensure_pkg("patchwork", bioc = FALSE)
ensure_pkg("SingleCellExperiment", bioc = TRUE)
ensure_pkg("scDblFinder", bioc = TRUE)
ensure_pkg("SingleR", bioc = TRUE)
ensure_pkg("celldex", bioc = TRUE)
ensure_pkg("scater", bioc = TRUE)
ensure_pkg("UCell", bioc = FALSE)
# 额外依赖（用于写CSV与整形数据）
ensure_pkg("readr", bioc = FALSE)
ensure_pkg("tidyr", bioc = FALSE)
ensure_pkg("tidyselect", bioc = FALSE)

# 读取对象与基本检查
if (!file.exists(in_rds)) {
  stop(sprintf("未找到整合对象：%s", in_rds))
}
seu <- readRDS(in_rds)
message("读取整合对象成功")
# 基本对象快照
message(sprintf("细胞数: %d, 基因数(当前默认assay): %d", ncol(seu), nrow(seu)))
message(sprintf("Assays: %s", paste(Seurat::Assays(seu), collapse = ", ")))
if (!("timepoint" %in% colnames(seu@meta.data))) {
  message("警告：meta.data 中缺少 'timepoint' 列")
}
flush.console()

# 备份
backup_path <- file.path(out_rds_dir, "nk.integrated.backup.rds")
saveRDS(seu, backup_path)
message("备份已保存：", backup_path)

# 确保有 RNA assay 与 counts（健壮性增强，避免 %in% 参数为 NULL）
assays_avail <- tryCatch(Seurat::Assays(seu), error = function(e) character(0))
message("可用 assays：", paste(assays_avail, collapse = ", "))
if (length(assays_avail) == 0) {
  stop("对象不包含任何 assay，无法继续。")
}
if ("RNA" %in% assays_avail) {
  Seurat::DefaultAssay(seu) <- "RNA"
} else {
  # 回退到首个可用 assay
  Seurat::DefaultAssay(seu) <- assays_avail[[1]]
  warning("对象不包含 RNA assay，改用可用 assay: ", Seurat::DefaultAssay(seu), "；scDblFinder 将基于该表达层运行。")
}
# 跳过 counts 检查，直接交由 scDblFinder 基于 SCE 对象处理（避免日志解析异常与 Seurat v5 接口差异）
message(sprintf("Default assay: %s", Seurat::DefaultAssay(seu)))

# 确定样本分组列（优先 orig.ident，其次 timepoint）
sample_col <- NULL
if ("orig.ident" %in% colnames(seu@meta.data) && length(unique(seu$orig.ident)) > 1) {
  sample_col <- "orig.ident"
} else if ("timepoint" %in% colnames(seu@meta.data) && length(unique(seu$timepoint)) > 1) {
  sample_col <- "timepoint"
} else {
  message("未找到可用的样本分组列（orig.ident/timepoint），scDblFinder 将不分组运行。")
}

# 工具函数：稳健构建 SingleCellExperiment（兼容 Seurat v4/v5，多层 assay）
to_sce <- function(seu, assay = Seurat::DefaultAssay(seu)) {
  assays_avail <- tryCatch(Seurat::Assays(seu), error = function(e) character(0))
  # layer 获取（v5 推荐）
  get_layer <- function(assay_name, layer_name) {
    if (!(assay_name %in% assays_avail)) return(NULL)
    m <- tryCatch(SeuratObject::GetAssayData(object = seu[[assay_name]], layer = layer_name), error = function(e) NULL)
    if (is.null(m)) return(NULL)
    if (is.null(dim(m)) || any(dim(m) == 0)) return(NULL)
    return(m)
  }
  # slot 获取（v4/v5 兼容，尽管已弃用但更健壮）
  get_slot <- function(assay_name, slot_name) {
    if (!(assay_name %in% assays_avail)) return(NULL)
    m <- tryCatch(Seurat::GetAssayData(seu[[assay_name]], slot = slot_name), error = function(e) NULL)
    if (is.null(m)) return(NULL)
    if (is.null(dim(m)) || any(dim(m) == 0)) return(NULL)
    return(m)
  }
  counts <- NULL
  logcounts <- NULL
  chosen_assay <- assay
  # 第一轮：优先 layer
  layer_candidates <- list(
    list(a = "RNA", l = "counts", t = "counts"),
    list(a = "RNA", l = "data",   t = "logcounts"),
    list(a = "SCT", l = "counts", t = "counts"),
    list(a = "SCT", l = "data",   t = "logcounts"),
    list(a = "integrated", l = "data", t = "logcounts")
  )
  for (cand in layer_candidates) {
    m <- get_layer(cand$a, cand$l)
    if (!is.null(m)) {
      if (cand$t == "counts" && is.null(counts)) {
        counts <- m; chosen_assay <- cand$a
      }
      if (cand$t == "logcounts" && is.null(logcounts)) {
        logcounts <- m; chosen_assay <- cand$a
      }
    }
    if (!is.null(counts) && !is.null(logcounts)) break
  }
  # 第二轮：layer 失败时回退到 slot
  if (is.null(counts)) {
    for (cand in list(
      list(a = "RNA", s = "counts"),
      list(a = "SCT", s = "counts"),
      list(a = "integrated", s = "counts")
    )) {
      m <- get_slot(cand$a, cand$s)
      if (!is.null(m)) { counts <- m; chosen_assay <- cand$a; break }
    }
  }
  if (is.null(logcounts)) {
    for (cand in list(
      list(a = "RNA", s = "data"),
      list(a = "SCT", s = "data"),
      list(a = "integrated", s = "data")
    )) {
      m <- get_slot(cand$a, cand$s)
      if (!is.null(m)) { logcounts <- m; chosen_assay <- cand$a; break }
    }
  }
  if (is.null(counts) && is.null(logcounts)) {
    stop("无法找到非空的 counts 或 data（已尝试 layer 与 slot 的 RNA/SCT/integrated）。")
  }
  # 与 Seurat 对象细胞对齐并按相同顺序
  cells <- colnames(seu)
  src_cols <- if (!is.null(counts)) colnames(counts) else colnames(logcounts)
  common <- intersect(cells, src_cols)
  if (length(common) == 0) stop("表达矩阵与 Seurat 对象细胞名不重叠。")
  if (!is.null(counts)) counts <- counts[, common, drop = FALSE]
  if (!is.null(logcounts)) logcounts <- logcounts[, common, drop = FALSE]
  cd <- seu@meta.data[common, , drop = FALSE]
  assays_list <- list()
  if (!is.null(counts)) assays_list$counts <- counts
  if (!is.null(logcounts)) assays_list$logcounts <- logcounts
  # 额外健壮性：确保 assay 与 colData 维度一致
  if (length(assays_list) == 0 || (nrow(cd) != if (!is.null(counts)) ncol(counts) else ncol(logcounts))) {
    stop("构建 SCE 时维度不一致：colData 行数与 assay 列数不匹配。")
  }
  sce <- SingleCellExperiment::SingleCellExperiment(assays = assays_list, colData = cd)
  rn <- if (!is.null(rownames(seu[[chosen_assay]]))) rownames(seu[[chosen_assay]]) else rownames(assays_list[[1]])
  rownames(sce) <- rn
  return(sce)
}

# 去双胞：scDblFinder
message("[1/5] 去双胞（scDblFinder）...")
sce <- to_sce(seu, assay = Seurat::DefaultAssay(seu))
if (!is.null(sample_col)) {
  message("scDblFinder 按分组列运行：", sample_col)
  dbl <- scDblFinder::scDblFinder(sce, samples = SingleCellExperiment::colData(sce)[[sample_col]])
} else {
  dbl <- scDblFinder::scDblFinder(sce)
}
# 写回元数据（健壮性增强：确保行名为细胞ID，避免 as.data.frame 丢失行名导致 match 报错）
cdf <- SingleCellExperiment::colData(dbl)
meta <- as.data.frame(cdf)
# 强制设置行名为 colData 的行名（细胞ID）
rownames(meta) <- rownames(cdf)

# 按整合对象的细胞ID对齐
idx <- match(colnames(seu), rownames(meta))
if (any(is.na(idx))) {
  # 回退：仅保留共有的细胞并对齐
  common <- intersect(colnames(seu), rownames(meta))
  if (length(common) == 0) {
    stop("scDblFinder 输出与 Seurat 对象的细胞名不重叠，无法写回元数据。")
  }
  meta <- meta[common, , drop = FALSE]
  seu <- subset(seu, cells = common)
  idx <- match(colnames(seu), rownames(meta))
}
meta <- meta[idx, , drop = FALSE]

seu$scDblFinder.class <- meta$scDblFinder.class
seu$scDblFinder.score <- meta$scDblFinder.score
if ("scDblFinder.highconf" %in% colnames(meta)) {
  seu$scDblFinder.highconf <- meta$scDblFinder.highconf
}

# 双胞统计与密度图
doublet_tab <- table(seu$scDblFinder.class)
doublet_df <- as.data.frame(doublet_tab)
readr::write_csv(doublet_df, file.path(out_reports_dir, "doublet_counts.csv"))

score_df <- data.frame(
  score = seu$scDblFinder.score,
  sample = if (!is.null(sample_col)) as.factor(seu@meta.data[[sample_col]]) else factor("All")
)
p_score <- ggplot2::ggplot(score_df, aes(x = score)) +
  ggplot2::geom_density(fill = "#2c7fb8", alpha = 0.4) +
  ggplot2::theme_bw() +
  ggplot2::labs(title = "scDblFinder Score Density", x = "score", y = "density")
ggplot2::ggsave(file.path(out_plots_dir, "DoubletScore_Density.png"), p_score, width = 8, height = 5, dpi = 300)
if (!is.null(sample_col)) {
  p_score_facet <- ggplot2::ggplot(score_df, aes(x = score)) +
    ggplot2::geom_density(fill = "#7fcdbb", alpha = 0.4) +
    ggplot2::facet_wrap(~ sample, scales = "free_y") +
    ggplot2::theme_bw() +
    ggplot2::labs(title = "scDblFinder Score Density by Sample", x = "score", y = "density")
  ggplot2::ggsave(file.path(out_plots_dir, "DoubletScore_Density_bySample.png"), p_score_facet, width = 10, height = 6, dpi = 300)
}

# 保存“仅去双胞（打标签但未剔除）”对象
saveRDS(seu, file.path(out_rds_dir, "nk.integrated.doublet_scored.rds"))

# 剔除标记为 doublet 的细胞
n_before <- ncol(seu)
seu <- subset(seu, subset = scDblFinder.class == "singlet")
n_after_doublet <- ncol(seu)
message(sprintf("去双胞：%d -> %d 细胞（移除 %d，%.2f%%）",
                n_before, n_after_doublet, n_before - n_after_doublet,
                100 * (n_before - n_after_doublet) / n_before))

# 自动注释（SingleR，使用 logcounts 修复 data 层空问题）
message("[2/5] 自动注释（SingleR，logcounts 修复）...")
# 稳健提取 counts（优先 layer=counts，若无则 slot=counts）
get_counts <- function(obj, assay_name) {
  m <- tryCatch(SeuratObject::GetAssayData(object = obj[[assay_name]], layer = "counts"), error = function(e) NULL)
  if (!is.null(m) && !is.null(dim(m)) && all(dim(m) > 0)) return(m)
  m <- tryCatch(Seurat::GetAssayData(obj[[assay_name]], slot = "counts"), error = function(e) NULL)
  if (!is.null(m) && !is.null(dim(m)) && all(dim(m) > 0)) return(m)
  return(NULL)
}
assays_avail <- tryCatch(Seurat::Assays(seu), error = function(e) character(0))
if (length(assays_avail) == 0) stop("SingleR：对象不包含任何 assay。")
chosen_assay <- if ("RNA" %in% assays_avail) "RNA" else if ("SCT" %in% assays_avail) "SCT" else assays_avail[[1]]
Seurat::DefaultAssay(seu) <- chosen_assay
counts <- get_counts(seu, chosen_assay)
if (is.null(counts) && chosen_assay != "RNA") counts <- get_counts(seu, "RNA")
if (is.null(counts) && chosen_assay != "SCT") counts <- get_counts(seu, "SCT")
if (is.null(counts) && chosen_assay != "integrated") counts <- get_counts(seu, "integrated")
if (is.null(counts)) stop("SingleR：无法提取非空 counts。")

# 与 Seurat 对象细胞对齐
cells <- colnames(seu)
common <- intersect(cells, colnames(counts))
if (length(common) == 0) stop("SingleR：counts 与 Seurat 对象细胞名不重叠。")
counts <- counts[, common, drop = FALSE]
cd <- seu@meta.data[common, , drop = FALSE]

# 构建 SCE 并计算 logcounts（避免 data 层为空问题）
sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts), colData = cd)
sce <- scater::logNormCounts(sce)  # 在 assays 中创建 'logcounts'
ref <- celldex::ImmGenData()
pred <- SingleR::SingleR(
  test = sce,
  ref = ref,
  labels = ref$label.main,
  assay.type.test = "logcounts",
  assay.type.ref = "logcounts"
)

# 写回注释（对齐到 common 细胞集）
seu <- subset(seu, cells = common)
seu$singleR.label <- pred$labels
if (!is.null(pred$pruned.labels)) {
  seu$singleR.pruned <- pred$pruned.labels
}

# UCell 签名评分
message("[3/5] 签名评分（UCell）...")
nk_genes <- c("Ncr1","Klrb1c","Nkg7","Prf1","Gzmb","Klrk1","Eomes","Tbx21")
t_genes  <- c("Cd3d","Cd3e","Cd3g","Trac","Cd4","Cd8a","Cd8b1","Lck")
b_genes  <- c("Ms4a1","Cd79a","Cd79b","Pax5","Cd74","Igkc")
myeloid  <- c("Lyz2","Lst1","Itgam","Itgax","Adgre1","Mrc1","Cx3cr1","S100a8","S100a9")
dc_genes <- c("Flt3","Xcr1","Siglech","Itgax")
plasma   <- c("Jchain","Xbp1","Mzb1","Igj")
endo     <- c("Pecam1","Kdr","Klf2")
fibro    <- c("Col1a1","Dcn")
hep      <- c("Alb","Ttr")

sign_list <- list(
  NK = nk_genes,
  T = t_genes,
  B = b_genes,
  Myeloid = myeloid,
  DC = dc_genes,
  Plasma = plasma,
  Endothelium = endo,
  Fibroblast = fibro,
  Hepatocyte = hep
)
# 注意：UCell 将创建 <name>_UCell 列
seu <- UCell::AddModuleScore_UCell(seu, features = sign_list)

# 细胞级污染判定
message("[4/5] 去污染（细胞级 + 簇级规则）...")
lab_col <- if ("singleR.pruned" %in% colnames(seu@meta.data)) {
  "singleR.pruned"
} else if ("singleR.label" %in% colnames(seu@meta.data)) {
  "singleR.label"
} else {
  NULL
}
label_lower <- if (!is.null(lab_col)) tolower(seu@meta.data[[lab_col]]) else rep("", ncol(seu))
is_nk_like <- grepl("nk|natural\\s*killer|ilc", label_lower)

# 计算非NK最高分
other_scores <- c("T_UCell","B_UCell","Myeloid_UCell","DC_UCell","Plasma_UCell","Endothelium_UCell","Fibroblast_UCell","Hepatocyte_UCell")
missing_others <- setdiff(other_scores, colnames(seu@meta.data))
if (length(missing_others) > 0) {
  warning("缺失以下签名列：", paste(missing_others, collapse = ", "))
}
safe_in <- function(x, table) {
  if (is.null(x) || is.null(table)) return(rep(FALSE, length(x)))
  x %in% table
}
get_col <- function(x) {
  v <- tryCatch(seu@meta.data[[x]], error = function(e) NULL)
  if (is.null(v)) {
    rep(NA_real_, ncol(seu))
  } else {
    as.numeric(v)
  }
}
max_nonNK <- do.call(pmax, c(lapply(other_scores, get_col), list(na.rm = TRUE)))
nk_score <- if ("NK_UCell" %in% colnames(seu@meta.data)) seu$NK_UCell else rep(NA_real_, ncol(seu))

qnk <- stats::quantile(nk_score, probs = 0.6, na.rm = TRUE)
delta_thresh <- 0.05
keep_cell <- (nk_score >= qnk & (nk_score - max_nonNK >= delta_thresh)) | is_nk_like
seu$contaminant_flag <- ifelse(keep_cell, "keep", "contaminant_cell_level")

# 如无既有聚类，先生成临时聚类用于簇级判定
if (!("seurat_clusters" %in% colnames(seu@meta.data))) {
  # 选择用于降维的 assay：优先 SCT，但若为多模型且缺少可变基因则回退 RNA
  dr_assay <- if ("SCT" %in% Seurat::Assays(seu)) "SCT" else "RNA"
  if (dr_assay == "SCT") {
    vf_curr <- tryCatch(Seurat::VariableFeatures(seu), error = function(e) character(0))
    if (length(vf_curr) == 0) {
      dr_assay <- "RNA"
    }
  }
  Seurat::DefaultAssay(seu) <- dr_assay
  if (dr_assay == "RNA") {
    seu <- Seurat::NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
    if (length(Seurat::VariableFeatures(seu)) == 0) {
      seu <- Seurat::FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
    }
    feats <- Seurat::VariableFeatures(seu)
  } else {
    # SCT：直接使用现有的可变基因，避免修改多模型的 VariableFeatures 设置
    feats <- Seurat::VariableFeatures(seu)
  }
  seu <- Seurat::ScaleData(seu, features = feats, verbose = FALSE)
  seu <- Seurat::RunPCA(seu, features = feats, verbose = FALSE)
  seu <- Seurat::FindNeighbors(seu, dims = 1:30)
  seu <- Seurat::FindClusters(seu, resolution = 0.5)
}

# 簇级污染判定：若某簇非NK注释占比 >= 阈值，则整体剔除
cluster_thresh <- 0.7
df_cluster <- seu@meta.data %>%
  dplyr::select(seurat_clusters) %>%
  dplyr::mutate(is_nk_like = is_nk_like) %>%
  dplyr::group_by(seurat_clusters) %>%
  dplyr::summarise(frac_nonNK = mean(!is_nk_like)) %>%
  dplyr::ungroup()

clusters_to_drop <- df_cluster$seurat_clusters[df_cluster$frac_nonNK >= cluster_thresh]
clu_char <- as.character(seu$seurat_clusters)
drop_char <- if (length(clusters_to_drop) > 0) as.character(clusters_to_drop) else character(0)
seu$cluster_drop_flag <- ifelse(safe_in(clu_char, drop_char), "drop_cluster", "keep_cluster")

# 综合剔除（细胞级或簇级标记为污染）
final_keep <- (seu$contaminant_flag == "keep") & (seu$cluster_drop_flag == "keep_cluster")
n_after_contam <- sum(final_keep, na.rm = TRUE)
removed_contam <- ncol(seu) - n_after_contam

# 统计汇总
removal_stats <- data.frame(
  step = c("initial", "after_doublet", "after_contaminant"),
  n_cells = c(n_before, n_after_doublet, n_after_contam)
)
readr::write_csv(removal_stats, file.path(out_reports_dir, "removal_stats.csv"))

# 可视化：簇级非NK占比条形图
p_cluster_nonNK <- ggplot2::ggplot(df_cluster, aes(x = seurat_clusters, y = frac_nonNK, fill = frac_nonNK >= cluster_thresh)) +
  ggplot2::geom_col() +
  ggplot2::scale_fill_manual(values = c("#2c7fb8", "#cb181d")) +
  ggplot2::theme_bw() +
  ggplot2::labs(title = "Fraction of non-NK annotation per cluster", x = "cluster", y = "fraction (non-NK)")
ggplot2::ggsave(file.path(out_plots_dir, "Cluster_nonNK_fraction.png"), p_cluster_nonNK, width = 10, height = 5, dpi = 300)

# 生成最终保留对象
seu_filtered <- subset(seu, cells = colnames(seu)[final_keep])

# 重跑降维/聚类/UMAP（最终）
# 最终降维/聚类：优先使用 SCT（若具备可变基因），否则使用 RNA 并计算 HVG
dr_assay2 <- if ("SCT" %in% Seurat::Assays(seu_filtered)) "SCT" else "RNA"
if (dr_assay2 == "SCT") {
  vf2 <- tryCatch(Seurat::VariableFeatures(seu_filtered), error = function(e) character(0))
  if (length(vf2) == 0) {
    dr_assay2 <- "RNA"
  }
}
Seurat::DefaultAssay(seu_filtered) <- dr_assay2
if (dr_assay2 == "RNA") {
  seu_filtered <- Seurat::NormalizeData(seu_filtered, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
  if (length(Seurat::VariableFeatures(seu_filtered)) == 0) {
    seu_filtered <- Seurat::FindVariableFeatures(seu_filtered, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  }
  feats2 <- Seurat::VariableFeatures(seu_filtered)
} else {
  feats2 <- Seurat::VariableFeatures(seu_filtered)
}
seu_filtered <- Seurat::ScaleData(seu_filtered, features = feats2, verbose = FALSE)
seu_filtered <- Seurat::RunPCA(seu_filtered, features = feats2, verbose = FALSE)
seu_filtered <- Seurat::FindNeighbors(seu_filtered, dims = 1:30)
seu_filtered <- Seurat::FindClusters(seu_filtered, resolution = 0.5)
seu_filtered <- Seurat::RunUMAP(seu_filtered, dims = 1:30)
if ("timepoint" %in% colnames(seu_filtered@meta.data)) {
  seu_filtered$timepoint <- factor(seu_filtered$timepoint, levels = c("0W_NCD","1W_MCD","2W_MCD","6W_MCD"))
}

# UMAP 可视化：按 SingleR 注释
p_umap_singleR <- Seurat::DimPlot(
  seu_filtered, reduction = "umap", group.by = if ("singleR.pruned" %in% colnames(seu_filtered@meta.data)) "singleR.pruned" else "singleR.label",
  label = FALSE
) + ggplot2::ggtitle("UMAP colored by SingleR annotation (filtered)")
ggplot2::ggsave(file.path(out_plots_dir, "UMAP_filtered_by_SingleR.png"), p_umap_singleR, width = 10, height = 8, dpi = 500)

# UMAP 分面：按 timepoint
if ("timepoint" %in% colnames(seu_filtered@meta.data)) {
  p_umap_tp <- Seurat::DimPlot(
    seu_filtered, reduction = "umap", group.by = "seurat_clusters", split.by = "timepoint", label = TRUE, repel = TRUE, ncol = 2
  ) + ggplot2::theme(strip.text.x = ggplot2::element_text(size = 12))
  ggplot2::ggsave(file.path(out_plots_dir, "UMAP_filtered_clusters_by_Timepoint.png"), p_umap_tp, width = 12, height = 10, dpi = 500)
}

# NK_UCell 小提琴图（按 SingleR 注释）
nk_ucell <- if ("NK_UCell" %in% colnames(seu_filtered@meta.data)) seu_filtered$NK_UCell else NULL
if (!is.null(nk_ucell)) {
  lab_col <- if ("singleR.pruned" %in% colnames(seu_filtered@meta.data)) "singleR.pruned" else "singleR.label"
  df_violin <- data.frame(label = seu_filtered@meta.data[[lab_col]], NK_UCell = nk_ucell)
  p_violin <- ggplot2::ggplot(df_violin, aes(x = label, y = NK_UCell)) +
    ggplot2::geom_violin(fill = "#9ecae1") +
    ggplot2::geom_boxplot(width = 0.1, outlier.size = 0.5) +
    ggplot2::theme_bw() +
    ggplot2::labs(title = "NK_UCell score across SingleR labels", x = "label", y = "NK_UCell")
  ggplot2::ggsave(file.path(out_plots_dir, "NK_UCell_by_SingleR_label.png"), p_violin, width = 12, height = 6, dpi = 300)
}

# 簇级签名热图（UCell 平均分）
cluster_scores <- seu_filtered@meta.data %>%
  dplyr::select(seurat_clusters, tidyselect::any_of(c("NK_UCell", other_scores))) %>%
  dplyr::group_by(seurat_clusters) %>%
  dplyr::summarise(dplyr::across(everything(), ~ mean(.x, na.rm = TRUE)))
cluster_scores_long <- tidyr::pivot_longer(cluster_scores, cols = -seurat_clusters, names_to = "signature", values_to = "mean_score")
p_heat <- ggplot2::ggplot(cluster_scores_long, aes(x = signature, y = seurat_clusters, fill = mean_score)) +
  ggplot2::geom_tile() +
  ggplot2::scale_fill_viridis_c() +
  ggplot2::theme_bw() +
  ggplot2::labs(title = "Cluster-level UCell signature means", x = "signature", y = "cluster")
ggplot2::ggsave(file.path(out_plots_dir, "Cluster_UCell_signature_means.png"), p_heat, width = 10, height = 8, dpi = 300)

# 保存对象
out_filtered_path <- file.path(out_rds_dir, "nk.integrated.filtered.rds")
saveRDS(seu_filtered, out_filtered_path)
message("清理后对象已保存：", out_filtered_path)

# 报告写入
report_path <- file.path(out_reports_dir, "cleaning_report.md")
report_lines <- c(
  "# NK细胞去双胞与去污染报告",
  "",
  sprintf("- 输入对象：%s", in_rds),
  sprintf("- 备份对象：%s", backup_path),
  sprintf("- 双胞评分对象：%s", file.path(out_rds_dir, "nk.integrated.doublet_scored.rds")),
  sprintf("- 清理后对象：%s", out_filtered_path),
  "",
  "## 参数与阈值",
  sprintf("- UCell NK分位阈值：P60 (%.4f)", as.numeric(qnk)),
  sprintf("- UCell Δ(NK - max(nonNK)) 阈值：%.3f", delta_thresh),
  sprintf("- 簇级非NK占比剔除阈值：%.2f", cluster_thresh),
  sprintf("- 聚类分辨率：0.5；PCA/UMAP维度：1:30"),
  "",
  "## 统计",
  sprintf("- 初始细胞数：%d", n_before),
  sprintf("- 去双胞后细胞数：%d（移除 %d，%.2f%%）", n_after_doublet, n_before - n_after_doublet, 100 * (n_before - n_after_doublet) / n_before),
  sprintf("- 去污染后细胞数：%d（再移除 %d，%.2f%%）", n_after_contam, removed_contam, 100 * removed_contam / n_after_doublet),
  "",
  "## 簇级非NK占比",
  paste0("```\n", paste(capture.output(print(df_cluster)), collapse = "\n"), "\n```"),
  "",
  "## 已保存图件",
  "- DoubletScore_Density.png",
  "- DoubletScore_Density_bySample.png（若按样本分组）",
  "- Cluster_nonNK_fraction.png",
  "- UMAP_filtered_by_SingleR.png",
  "- UMAP_filtered_clusters_by_Timepoint.png（若含timepoint）",
  "- NK_UCell_by_SingleR_label.png",
  "- Cluster_UCell_signature_means.png",
  "",
  "## R会话信息",
  paste0("```\n", paste(capture.output(utils::sessionInfo()), collapse = "\n"), "\n```")
)
writeLines(report_lines, con = report_path)
message("报告已生成：", report_path)

message("全部流程完成。")

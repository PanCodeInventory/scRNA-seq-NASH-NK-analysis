#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(optparse)
})

# ---------------- CLI ----------------
option_list <- list(
  make_option(c("--rds-in"), type = "character", help = "Input Seurat RDS path", default = "2_DataProcessing/RDS/nk.integrated.singleR_annotated.rds"),
  make_option(c("--out-rds"), type = "character", help = "Output Seurat RDS path", default = "2_DataProcessing/RDS/nk.integrated.singleR_annotated.noCluster6.rds"),
  make_option(c("--rm-clusters"), type = "character", help = "Comma-separated cluster labels to remove (compared against seurat_clusters/Idents as character)", default = "6"),
  make_option(c("--dims"), type = "integer", help = "Number of PCA dims for neighbors/UMAP", default = 10L),
  make_option(c("--resolution"), type = "double", help = "Louvain/Leiden clustering resolution", default = 0.3),
  make_option(c("--plots-dir"), type = "character", help = "Directory to write plots (UMAP)", default = "2_DataProcessing/3_Tuning/plots"),
  make_option(c("--logs-dir"), type = "character", help = "Directory to write logs/config", default = "2_DataProcessing/3_Tuning/logs"),
  make_option(c("--umap-width"), type = "double", help = "UMAP figure width (inches)", default = 9),
  make_option(c("--umap-height"), type = "double", help = "UMAP figure height (inches)", default = 7),
  make_option(c("--dpi"), type = "integer", help = "UMAP figure dpi", default = 300),
  make_option(c("--timepoint-order"), type = "character", help = "Comma-separated timepoint order for faceting", default = "0W_NCD,1W_MCD,2W_MCD,6W_MCD")
)
opt <- parse_args(OptionParser(option_list = option_list))

# ---------------- Helpers ----------------
split_commas <- function(x) {
  if (is.null(x) || is.na(x) || x == "") character(0) else trimws(strsplit(x, ",")[[1]])
}

ensure_dirs <- function(...) {
  dirs <- list(...)
  for (d in dirs) dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

safe_idents_vec <- function(obj) {
  if ("seurat_clusters" %in% colnames(obj@meta.data)) {
    as.character(obj@meta.data$seurat_clusters)
  } else {
    as.character(Idents(obj))
  }
}

detect_timepoint_col <- function(md) {
  cands <- c("timepoint", "Timepoint", "time", "Time", "orig.ident", "group", "condition", "sample", "dataset")
  cands <- cands[cands %in% colnames(md)]
  if (length(cands) > 0) cands[[1]] else NA_character_
}

detect_singleR_col <- function(md) {
  cands <- c("SingleR.pruned.labels", "SingleR.labels", "singleR.pruned.labels", "singleR_labels", "SingleR_label")
  cands <- cands[cands %in% colnames(md)]
  if (length(cands) > 0) cands[[1]] else NA_character_
}

assay_has_data <- function(obj, assay) {
  if (!(assay %in% Assays(obj))) return(FALSE)
  ok <- tryCatch({
    m <- GetAssayData(obj, assay = assay, slot = "data")
    is.matrix(m) || inherits(m, "dgCMatrix")
  }, error = function(e) FALSE)
  if (!ok) return(FALSE)
  dims_ok <- tryCatch({
    m <- GetAssayData(obj, assay = assay, slot = "data")
    ncol(m) > 0 && nrow(m) > 0
  }, error = function(e) FALSE)
  dims_ok
}

var_features_n <- function(obj, assay) {
  DefaultAssay(obj) <- assay
  tryCatch(length(VariableFeatures(obj)), error = function(e) 0L)
}

pick_assay <- function(obj) {
  # Priority: integrated (if has data & VFs) -> SCT (if VFs>0) -> RNA (fallback; will normalize/scale)
  if ("integrated" %in% Assays(obj)) {
    if (assay_has_data(obj, "integrated") && var_features_n(obj, "integrated") > 0) {
      return(list(assay = "integrated", mode = "ready"))
    }
  }
  if ("SCT" %in% Assays(obj)) {
    if (var_features_n(obj, "SCT") > 0) {
      # SCT often has VFs; data layer may be empty but Seurat uses SCT layers internally
      return(list(assay = "SCT", mode = "ready"))
    }
  }
  if ("RNA" %in% Assays(obj)) {
    return(list(assay = "RNA", mode = "prep")) # need Normalize/FindVariableFeatures/Scale
  }
  # Fallback to any assay present
  asy <- Assays(obj)
  if (length(asy) > 0) return(list(assay = asy[[1]], mode = "unknown"))
  stop("No assays found in object.")
}

safe_dims <- function(obj, requested) {
  if (!("pca" %in% names(obj@reductions))) return(1:requested)
  nd <- tryCatch(ncol(Embeddings(obj, "pca")), error = function(e) requested)
  1:min(requested, nd)
}

save_plot <- function(p, path_png, path_pdf, width, height, dpi) {
  ggsave(path_png, p, width = width, height = height, dpi = dpi, units = "in")
  ggsave(path_pdf, p, width = width, height = height, dpi = dpi, units = "in")
}

# ---------------- Main ----------------
main <- function(opt) {
  message(sprintf("[INFO] Input: %s", opt$`rds-in`))
  if (!file.exists(opt$`rds-in`)) stop(sprintf("Input RDS not found: %s", opt$`rds-in`))
  ensure_dirs(dirname(opt$`out-rds`), opt$`plots-dir`, opt$`logs-dir`)

  ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
  run_cfg <- file.path(opt$`logs-dir`, sprintf("run_config_remove_clusters_%s.txt", ts))
  sel_params <- file.path(opt$`logs-dir`, "selected_params.txt")
  sess_file <- file.path(opt$`logs-dir`, sprintf("sessionInfo_%s.txt", ts))

  obj <- readRDS(opt$`rds-in`)

  # Backup original clusters to seurat_clusters_orig
  md <- obj@meta.data
  if ("seurat_clusters" %in% colnames(md)) {
    obj$seurat_clusters_orig <- as.character(md$seurat_clusters)
  } else {
    obj$seurat_clusters_orig <- as.character(Idents(obj))
  }

  clusters_vec <- safe_idents_vec(obj)
  rm_clusters <- split_commas(opt$`rm-clusters`)
  if (length(rm_clusters) == 0) stop("No clusters specified to remove.")
  message(sprintf("[INFO] Removing clusters: %s", paste(rm_clusters, collapse = ", ")))
  to_keep <- ! (as.character(clusters_vec) %in% rm_clusters)
  cells_keep <- colnames(obj)[to_keep]
  n_removed <- sum(!to_keep)
  message(sprintf("[INFO] Cells before: %d; after: %d; removed: %d", ncol(obj), length(cells_keep), n_removed))

  if (length(cells_keep) == 0) stop("All cells removed; check rm-clusters specification.")
  obj2 <- subset(obj, cells = cells_keep)

  # Ensure seurat_clusters column exists (pre-recluster, reflects old labels minus removed)
  if (!("seurat_clusters" %in% colnames(obj2@meta.data))) {
    obj2$seurat_clusters <- as.character(Idents(obj2))
  } else {
    obj2$seurat_clusters <- as.character(obj2$seurat_clusters)
    obj2$seurat_clusters <- droplevels(factor(obj2$seurat_clusters)) %>% as.character()
  }

  # Choose assay and prepare if needed
  choice <- pick_assay(obj2)
  DefaultAssay(obj2) <- choice$assay
  message(sprintf("[INFO] Using assay: %s (mode=%s)", choice$assay, choice$mode))
  if (choice$assay == "RNA" && choice$mode %in% c("prep", "unknown")) {
    message("[INFO] Preparing RNA assay: Normalize -> FindVariableFeatures -> ScaleData")
    obj2 <- NormalizeData(obj2, verbose = FALSE)
    obj2 <- FindVariableFeatures(obj2, nfeatures = 3000, verbose = FALSE)
    obj2 <- ScaleData(obj2, verbose = FALSE)
  }

  # PCA -> Neighbors -> Clusters -> UMAP
  vf_n <- tryCatch(length(VariableFeatures(obj2)), error = function(e) 0L)
  if (vf_n == 0) {
    message("[WARN] VariableFeatures is empty; attempting to compute on RNA assay.")
    DefaultAssay(obj2) <- "RNA"
    obj2 <- NormalizeData(obj2, verbose = FALSE)
    obj2 <- FindVariableFeatures(obj2, nfeatures = 3000, verbose = FALSE)
    obj2 <- ScaleData(obj2, verbose = FALSE)
  }

  message(sprintf("[INFO] Running PCA/Neighbors/Clusters/UMAP with dims=1:%d, resolution=%.3f", opt$dims, opt$resolution))
  obj2 <- RunPCA(obj2, features = VariableFeatures(obj2), npcs = max(10L, opt$dims), verbose = FALSE)
  dims_use <- safe_dims(obj2, opt$dims)
  obj2 <- FindNeighbors(obj2, reduction = "pca", dims = dims_use, verbose = FALSE)
  obj2 <- FindClusters(obj2, resolution = opt$resolution, verbose = FALSE)
  obj2 <- RunUMAP(obj2, reduction = "pca", dims = dims_use, verbose = FALSE)

  # Persist new clusters into meta.data
  obj2$seurat_clusters <- as.character(Idents(obj2))

  # Save new object
  saveRDS(obj2, file = opt$`out-rds`)
  message(sprintf("[INFO] Saved new object: %s", opt$`out-rds`))

  # UMAP plots
  md2 <- obj2@meta.data
  tp_col <- detect_timepoint_col(md2)
  if (is.na(tp_col)) {
    message("[WARN] No timepoint-like column detected; UMAP by timepoint will be skipped.")
  }
  sr_col <- detect_singleR_col(md2)
  if (is.na(sr_col)) {
    message("[WARN] No SingleR label column detected; UMAP by SingleR will be skipped.")
  }

  tp_order <- split_commas(opt$`timepoint-order`)
  # 1) UMAP faceted by timepoint (colored by clusters)
  if (!is.na(tp_col)) {
    # Ensure factor order if provided
    md2[[tp_col]] <- as.character(md2[[tp_col]])
    if (length(tp_order) > 0) {
      md2[[tp_col]] <- factor(md2[[tp_col]], levels = tp_order)
    } else {
      md2[[tp_col]] <- factor(md2[[tp_col]])
    }
    obj2@meta.data <- md2

    p_tp <- tryCatch({
      DimPlot(obj2, reduction = "umap", group.by = "seurat_clusters", split.by = tp_col, label = TRUE, label.size = 3) +
        ggtitle(sprintf("UMAP by timepoint (split=%s)", tp_col))
    }, error = function(e) {
      message(sprintf("[WARN] Failed to create UMAP by timepoint: %s", e$message))
      NULL
    })
    if (!is.null(p_tp)) {
      out_png <- file.path(opt$`plots-dir`, "UMAP_noCluster6_byTimepoint.png")
      out_pdf <- file.path(opt$`plots-dir`, "UMAP_noCluster6_byTimepoint.pdf")
      save_plot(p_tp, out_png, out_pdf, width = opt$`umap-width`, height = opt$`umap-height`, dpi = opt$dpi)
      message(sprintf("[INFO] Saved: %s ; %s", out_png, out_pdf))
    }
  }

  # 2) UMAP colored by SingleR labels (no split)
  if (!is.na(sr_col)) {
    p_sr <- tryCatch({
      DimPlot(obj2, reduction = "umap", group.by = sr_col, label = FALSE) +
        ggtitle(sprintf("UMAP by SingleR (%s)", sr_col))
    }, error = function(e) {
      message(sprintf("[WARN] Failed to create UMAP by SingleR: %s", e$message))
      NULL
    })
    if (!is.null(p_sr)) {
      out_png <- file.path(opt$`plots-dir`, "UMAP_noCluster6_bySingleR.png")
      out_pdf <- file.path(opt$`plots-dir`, "UMAP_noCluster6_bySingleR.pdf")
      save_plot(p_sr, out_png, out_pdf, width = opt$`umap-width`, height = opt$`umap-height`, dpi = opt$dpi)
      message(sprintf("[INFO] Saved: %s ; %s", out_png, out_pdf))
    }
  }

  # Write run config and session info
  cfg <- tibble::tibble(
    key = c("rds_in", "out_rds", "rm_clusters", "dims", "resolution", "assay_used", "timepoint_col", "singleR_col"),
    value = c(opt$`rds-in`, opt$`out-rds`, paste(split_commas(opt$`rm-clusters`), collapse = ","), opt$dims, opt$resolution, DefaultAssay(obj2), ifelse(is.na(tp_col), "", tp_col), ifelse(is.na(sr_col), "", sr_col))
  )
  readr::write_tsv(cfg, run_cfg)
  readr::write_lines(sprintf("dims=%d\nresolution=%.3f", opt$dims, opt$resolution), sel_params)
  readr::write_lines(paste(capture.output(sessionInfo()), collapse = "\n"), sess_file)

  message(sprintf("[INFO] Done. Config: %s ; Session: %s", run_cfg, sess_file))
}

tryCatch(
  { main(opt) },
  error = function(e) {
    message("[ERROR] ", e$message)
    quit(status = 1)
  }
)

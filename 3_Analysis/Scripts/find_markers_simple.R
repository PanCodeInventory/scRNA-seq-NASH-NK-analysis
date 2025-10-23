#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
})

# Paths
rds_path <- "2_DataProcessing/RDS/nk.integrated.singleR_annotated.noCluster6.tuned.rds"
outdir   <- "3_Analysis/1.ClusterAnalysis"
data_dir <- file.path(outdir, "data")
logs_dir <- file.path(outdir, "logs")

dirs <- c(outdir, data_dir, logs_dir)
for (d in dirs) dir.create(d, showWarnings = FALSE, recursive = TRUE)

log_file <- file.path(logs_dir, sprintf("find_markers_simple_%s.txt", format(Sys.time(), "%Y%m%d_%H%M%S")))
con <- file(log_file, open = "wt")
wlog <- function(...) writeLines(paste0(...), con = con)

tryCatch({
  wlog(sprintf("[INFO] Loading RDS: %s", rds_path))
  obj <- readRDS(rds_path)
  wlog(sprintf("[INFO] Cells: %d; meta.data cols: %d", ncol(obj), ncol(obj@meta.data)))
  wlog(sprintf("[INFO] meta.data columns: %s", paste(colnames(obj@meta.data), collapse = ", ")))

  # Set Idents from seurat_clusters if available, else fallback to existing Idents
  if ("seurat_clusters" %in% colnames(obj@meta.data)) {
    cl <- as.character(obj@meta.data$seurat_clusters)
    if (length(cl) != ncol(obj)) stop("Cluster vector length mismatch with number of cells.")
    names(cl) <- colnames(obj)
    Idents(obj) <- factor(cl)
    wlog(sprintf("[INFO] Idents set from seurat_clusters; unique clusters: %d", length(unique(cl))))
  } else {
    wlog("[WARN] seurat_clusters not found; using existing Idents(obj).")
    Idents(obj) <- Idents(obj)
    wlog(sprintf("[INFO] unique clusters in Idents(obj): %d", length(unique(as.character(Idents(obj))))))
  }

  # DefaultAssay: prefer 'integrated' then 'SCT' then 'RNA'
  if ("integrated" %in% Assays(obj)) {
    DefaultAssay(obj) <- "integrated"
    wlog("[INFO] DefaultAssay set to 'integrated'")
  } else if ("SCT" %in% Assays(obj)) {
    DefaultAssay(obj) <- "SCT"
    wlog("[INFO] DefaultAssay set to 'SCT'")
  } else if ("RNA" %in% Assays(obj)) {
    DefaultAssay(obj) <- "RNA"
    wlog("[INFO] DefaultAssay set to 'RNA'")
  } else {
    stop("No suitable assay found (integrated/SCT/RNA).")
  }

  da <- DefaultAssay(obj)
  # Check data dims and variable features
  data_dims <- tryCatch({
    m <- GetAssayData(obj, assay = da, slot = "data")
    c(nrow = nrow(m), ncol = ncol(m))
  }, error = function(e) c(nrow = NA_integer_, ncol = NA_integer_))
  vf_len <- tryCatch(length(VariableFeatures(obj)), error = function(e) NA_integer_)
  wlog(sprintf("[INFO] Assay '%s' data dims: %sx%s; VariableFeatures: %s", da, data_dims["nrow"], data_dims["ncol"], vf_len))

  # If RNA chosen and data empty, prepare RNA
  if (da == "RNA" && (is.na(data_dims["nrow"]) || data_dims["nrow"] == 0 || is.na(vf_len) || vf_len == 0)) {
    wlog("[INFO] Preparing RNA assay: NormalizeData, FindVariableFeatures, ScaleData")
    obj <- NormalizeData(obj, assay = "RNA")
    obj <- FindVariableFeatures(obj, assay = "RNA")
    obj <- ScaleData(obj, assay = "RNA", features = VariableFeatures(obj))
    DefaultAssay(obj) <- "RNA"
    vf_len <- tryCatch(length(VariableFeatures(obj)), error = function(e) NA_integer_)
    wlog(sprintf("[INFO] RNA VariableFeatures after prepare: %s", vf_len))
  }

  # Run FindAllMarkers
  wlog("[INFO] Running FindAllMarkers (wilcox, min.pct=0.1, logfc.threshold=0.25, only.pos=TRUE)")
  markers <- FindAllMarkers(
    object = obj,
    only.pos = TRUE,
    min.pct = 0.1,
    logfc.threshold = 0.25,
    test.use = "wilcox"
  )
  if (nrow(markers) == 0) stop("FindAllMarkers returned 0 rows.")
  # Ensure avg_log2FC column
  if (!("avg_log2FC" %in% colnames(markers)) && ("avg_logFC" %in% colnames(markers))) {
    markers <- markers %>% mutate(avg_log2FC = avg_logFC)
  }

  # Save all markers
  all_csv <- file.path(data_dir, "markers_all_clusters.csv")
  write_csv(markers, all_csv)
  wlog(sprintf("[INFO] Wrote: %s", all_csv))

  # Top10 per cluster
  top_df <- markers %>%
    filter(!is.na(p_val_adj), p_val_adj <= 0.05) %>%
    group_by(cluster) %>%
    arrange(desc(avg_log2FC), p_val_adj, .by_group = TRUE) %>%
    slice_head(n = 10) %>%
    ungroup()
  top_csv <- file.path(data_dir, "markers_top10_per_cluster.csv")
  write_csv(top_df, top_csv)
  wlog(sprintf("[INFO] Wrote: %s", top_csv))

  wlog("[INFO] Done.")
}, error = function(e) {
  wlog(paste0("[ERROR] ", e$message))
  close(con)
  quit(status = 1)
})

close(con)

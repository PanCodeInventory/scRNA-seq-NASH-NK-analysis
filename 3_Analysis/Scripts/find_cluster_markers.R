#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
  library(optparse)
})

# ---------- Helpers ----------
ensure_idents <- function(obj, cluster_col = "") {
  if (cluster_col != "") {
    if (!cluster_col %in% colnames(obj@meta.data)) {
      stop(sprintf("Specified --cluster-col '%s' not found in meta.data.", cluster_col))
    }
    Idents(obj) <- obj@meta.data[[cluster_col]]
  } else {
    id_vec <- Idents(obj)
    if (length(unique(as.character(id_vec))) <= 1) {
      candidates <- c("seurat_clusters", "clusters", "cluster", "Cluster")
      found <- candidates[candidates %in% colnames(obj@meta.data)]
      if (length(found) > 0) Idents(obj) <- obj@meta.data[[found[[1]]]]
    }
  }
  obj
}

assay_exists <- function(obj, assay) {
  assay %in% Assays(obj)
}

has_assay_data <- function(obj, assay) {
  ok <- FALSE
  try({
    mat <- GetAssayData(obj, assay = assay, slot = "data")
    ok <- (!is.null(mat)) && (nrow(mat) > 0) && (ncol(mat) > 0)
  }, silent = TRUE)
  ok
}

has_variable_features <- function(obj, assay) {
  ok <- FALSE
  try({
    DefaultAssay(obj) <- assay
    vf <- VariableFeatures(obj)
    ok <- length(vf) > 0
  }, silent = TRUE)
  ok
}

prepare_rna_if_needed <- function(obj) {
  DefaultAssay(obj) <- "RNA"
  if (!has_assay_data(obj, "RNA") || !has_variable_features(obj, "RNA")) {
    message("[INFO] Preparing RNA assay: NormalizeData, FindVariableFeatures, ScaleData")
    obj <- NormalizeData(obj, assay = "RNA")
    obj <- FindVariableFeatures(obj, assay = "RNA")
    obj <- ScaleData(obj, assay = "RNA", features = VariableFeatures(obj))
  }
  obj
}

choose_assay <- function(obj, priorities = c("integrated", "SCT", "RNA")) {
  pri <- priorities[priorities != ""]
  available <- pri[pri %in% Assays(obj)]
  for (a in available) {
    if (a == "SCT") {
      if (has_assay_data(obj, "SCT") && has_variable_features(obj, "SCT")) {
        return(list(assay = "SCT", prepared = TRUE))
      } else {
        message("[WARN] SCT assay lacks data or variable features; will consider fallback.")
      }
    } else if (a == "RNA") {
      return(list(assay = "RNA", prepared = FALSE))
    } else {
      # other assay types (e.g. 'integrated'): require non-empty data layer
      if (has_assay_data(obj, a)) {
        return(list(assay = a, prepared = TRUE))
      }
    }
  }
  # Fallback to RNA if present
  if ("RNA" %in% Assays(obj)) {
    return(list(assay = "RNA", prepared = FALSE))
  }
  stop("No suitable assay found among priorities or available assays.")
}

write_run_logs <- function(logs_dir, config_lines) {
  ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
  cfg_file <- file.path(logs_dir, sprintf("run_config_%s.txt", ts))
  ses_file <- file.path(logs_dir, sprintf("sessionInfo_%s.txt", ts))
  writeLines(config_lines, con = cfg_file)
  writeLines(capture.output(sessionInfo()), con = ses_file)
  invisible(list(config = cfg_file, session = ses_file))
}

unify_logfc_column <- function(df) {
  if (!("avg_log2FC" %in% colnames(df)) && ("avg_logFC" %in% colnames(df))) {
    df <- df %>% mutate(avg_log2FC = avg_logFC)
  }
  df
}

# ---------- CLI ----------
option_list <- list(
  make_option(c("--rds"), type = "character", default = "1_Files/RDS/nk1.1_integrated.tuned.rds",
              help = "Input Seurat RDS path [default %default]"),
  make_option(c("--outdir"), type = "character", default = "3_Analysis/1.ClusterAnalysis",
              help = "Output root directory for data/plots/logs [default %default]"),
  make_option(c("--assay-priority"), type = "character", default = "integrated,SCT,RNA",
              help = "Comma-separated assay priorities (e.g. 'integrated,SCT,RNA')"),
  make_option(c("--cluster-col"), type = "character", default = "",
              help = "Cluster column to set as Idents; use current Idents if empty"),
  make_option(c("--only-pos"), type = "logical", default = TRUE,
              help = "Find positive markers only [default %default]"),
  make_option(c("--min-pct"), type = "double", default = 0.1,
              help = "Minimum fraction of cells expressing the feature [default %default]"),
  make_option(c("--logfc-threshold"), type = "double", default = 0.25,
              help = "Log2 fold-change threshold [default %default]"),
  make_option(c("--test-use"), type = "character", default = "wilcox",
              help = "Statistical test (e.g. 'wilcox', 'MAST') [default %default]"),
  make_option(c("--fdr"), type = "double", default = 0.05,
              help = "Adjusted p-value (FDR) threshold for TopN selection [default %default]"),
  make_option(c("--topn"), type = "integer", default = 10,
              help = "Top N markers per cluster to export [default %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

# ---------- Main ----------
main <- function(opt) {
  if (!file.exists(opt$rds)) {
    stop(sprintf("Input RDS not found: %s", opt$rds))
  }

  outroot <- opt$outdir
  data_dir  <- file.path(outroot, "data")
  plots_dir <- file.path(outroot, "plots")
  logs_dir  <- file.path(outroot, "logs")
  dirs <- c(outroot, data_dir, plots_dir, logs_dir)
  for (d in dirs) dir.create(d, showWarnings = FALSE, recursive = TRUE)

  message(sprintf("[INFO] Loading RDS: %s", opt$rds))
  obj <- readRDS(opt$rds)

  # Robust clusters (avoid ensure_idents early error)
  cluster_vec <- NULL
  if (opt$cluster_col != "" && opt$cluster_col %in% colnames(obj@meta.data)) {
    cluster_vec <- as.character(obj@meta.data[[opt$cluster_col]])
    message(sprintf("[DEBUG] Using provided cluster-col '%s'", opt$cluster_col))
  } else if ("seurat_clusters" %in% colnames(obj@meta.data)) {
    cluster_vec <- as.character(obj@meta.data$seurat_clusters)
    message("[DEBUG] Using meta.data$seurat_clusters")
  } else {
    cluster_vec <- as.character(Idents(obj))
    message("[DEBUG] Using Idents(obj)")
  }
  if (is.null(cluster_vec) || length(cluster_vec) == 0) {
    stop("Cluster vector is empty; please provide --cluster-col or check meta.data/Idents.")
  }
  message(sprintf("[DEBUG] cluster_vec length: %d; unique clusters: %d", length(cluster_vec), length(unique(cluster_vec))))
  if (length(unique(cluster_vec)) <= 1) {
    stop("Only one unique cluster detected; differential testing would be trivial.")
  }
  # Ensure length and naming align with object cells, then set Idents
  if (length(cluster_vec) != ncol(obj)) {
    stop(sprintf("Cluster vector length (%d) != number of cells (%d). Please check meta.data/cluster column.",
                 length(cluster_vec), ncol(obj)))
  }
  # Name the cluster vector to align with colnames(obj)
  names(cluster_vec) <- colnames(obj)
  Idents(obj) <- factor(cluster_vec)
  message("[DEBUG] Idents set from cluster_vec; levels: ", paste(levels(Idents(obj)), collapse = ", "))

  # Assay selection
  pri <- trimws(unlist(strsplit(opt$assay_priority, ",")))
  sel <- choose_assay(obj, priorities = pri)
  assay <- sel$assay
  message(sprintf("[INFO] Selected assay: %s", assay))
  if (assay == "RNA" && !sel$prepared) {
    obj <- prepare_rna_if_needed(obj)
  } else {
    # For non-RNA assays, ensure DefaultAssay set to selected
    DefaultAssay(obj) <- assay
  }

  # Optional: set parallel if future available — force sequential to avoid environment issues
  if (requireNamespace("future", quietly = TRUE)) {
    suppressPackageStartupMessages(library(future))
    plan(sequential)
    message("[INFO] Using future::sequential plan for consistency.")
  }

  # Debug assay state before DE
  message(sprintf("[DEBUG] DefaultAssay: %s", DefaultAssay(obj)))
  da <- DefaultAssay(obj)
  data_dims <- tryCatch({
    m <- GetAssayData(obj, assay = da, slot = "data")
    c(nrow = nrow(m), ncol = ncol(m))
  }, error = function(e) c(nrow = NA_integer_, ncol = NA_integer_))
  vf_len <- tryCatch(length(VariableFeatures(obj)), error = function(e) NA_integer_)
  message(sprintf("[DEBUG] Assay '%s' data dims: %sx%s; VariableFeatures: %s", da, data_dims["nrow"], data_dims["ncol"], vf_len))

  # Run FindAllMarkers with robust tryCatch and fallback
  message("[INFO] Running FindAllMarkers ...")
  markers <- NULL
  fm_err <- NULL
  tryCatch({
    markers <- FindAllMarkers(
      object = obj,
      only.pos = opt$only_pos,
      min.pct = opt$min_pct,
      logfc.threshold = opt$logfc_threshold,
      test.use = opt$test_use
    )
  }, error = function(e) {
    fm_err <<- e$message
    message(sprintf("[WARN] FindAllMarkers failed on assay '%s': %s", da, e$message))
  })

  # Fallback attempts if needed
  if (is.null(markers)) {
    # Try SCT if available and not current
    if ("SCT" %in% Assays(obj) && da != "SCT") {
      message("[INFO] Fallback: switching to SCT assay")
      DefaultAssay(obj) <- "SCT"
      vf_len <- tryCatch(length(VariableFeatures(obj)), error = function(e) NA_integer_)
      message(sprintf("[DEBUG] SCT VariableFeatures: %s", vf_len))
      tryCatch({
        markers <- FindAllMarkers(
          object = obj,
          only.pos = opt$only_pos,
          min.pct = opt$min_pct,
          logfc.threshold = opt$logfc_threshold,
          test.use = opt$test_use
        )
      }, error = function(e) {
        message(sprintf("[WARN] FindAllMarkers failed on SCT: %s", e$message))
      })
    }
  }
  if (is.null(markers)) {
    # Try RNA (prepare if needed)
    if ("RNA" %in% Assays(obj) && da != "RNA") {
      message("[INFO] Fallback: switching to RNA assay and preparing if needed")
      obj <- prepare_rna_if_needed(obj)
      DefaultAssay(obj) <- "RNA"
      vf_len <- tryCatch(length(VariableFeatures(obj)), error = function(e) NA_integer_)
      message(sprintf("[DEBUG] RNA VariableFeatures after prepare: %s", vf_len))
      tryCatch({
        markers <- FindAllMarkers(
          object = obj,
          only.pos = opt$only_pos,
          min.pct = opt$min_pct,
          logfc.threshold = opt$logfc_threshold,
          test.use = opt$test_use
        )
      }, error = function(e) {
        message(sprintf("[ERROR] FindAllMarkers failed on RNA: %s", e$message))
        fm_err <<- e$message
      })
    }
  }

  if (is.null(markers)) {
    stop(sprintf("FindAllMarkers failed across assays. Last error: %s", fm_err))
  }
  if (nrow(markers) == 0) {
    stop("FindAllMarkers returned 0 rows; check parameters and assay preparation.")
  }

  markers <- unify_logfc_column(markers)

  all_csv <- file.path(data_dir, "markers_all_clusters.csv")
  write_csv(markers, all_csv)
  message(sprintf("[INFO] Wrote: %s", all_csv))

  # TopN per cluster with FDR filter
  topn <- opt$topn
  fdr_thr <- opt$fdr
  top_df <- markers %>%
    filter(!is.na(p_val_adj), p_val_adj <= fdr_thr) %>%
    group_by(cluster) %>%
    arrange(desc(avg_log2FC), p_val_adj, .by_group = TRUE) %>%
    slice_head(n = topn) %>%
    ungroup()

  top_csv <- file.path(data_dir, "markers_top10_per_cluster.csv")
  write_csv(top_df, top_csv)
  message(sprintf("[INFO] Wrote: %s", top_csv))

  # Logs
  config_lines <- c(
    sprintf("rds: %s", opt$rds),
    sprintf("outdir: %s", opt$outdir),
    sprintf("assay_priority: %s", opt$assay_priority),
    sprintf("selected_assay: %s", assay),
    sprintf("only_pos: %s", opt$only_pos),
    sprintf("min_pct: %s", opt$min_pct),
    sprintf("logfc_threshold: %s", opt$logfc_threshold),
    sprintf("test_use: %s", opt$test_use),
    sprintf("fdr: %s", opt$fdr),
    sprintf("topn: %s", opt$topn),
    sprintf("n_cells: %s", ncol(obj)),
    sprintf("n_clusters: %s", length(unique(as.character(Idents(obj)))))
  )
  write_run_logs(logs_dir, config_lines)

  message("[INFO] Done.")
}

tryCatch(
  { main(opt) },
  error = function(e) {
    message("[ERROR] ", e$message)
    quit(status = 1)
  }
)

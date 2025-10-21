#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(readr)
  library(optparse)
})

# ---------- CLI ----------
option_list <- list(
  make_option(c("--rds"), type = "character", default = "1_Files/RDS/nk1.1_integrated.tuned.rds",
              help = "Input Seurat RDS path [default %default]"),
  make_option(c("--outdir"), type = "character", default = "3_Analysis/1.ClusterAnalysis",
              help = "Output root directory for logs/data [default %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

# ---------- Helpers ----------
safe_len <- function(x) if (is.null(x)) 0L else length(x)

detect_timepoint_candidates <- function(md) {
  candidates <- c("timepoint", "Timepoint", "time", "Time", "orig.ident", "group", "condition", "sample", "dataset")
  candidates[candidates %in% colnames(md)]
}

detect_timepoint_col <- function(md, preferred = "") {
  if (preferred != "" && preferred %in% colnames(md)) return(preferred)
  cand <- detect_timepoint_candidates(md)
  if (length(cand) == 0) return(NA_character_)
  cand[[1]]
}

safe_idents <- function(obj) {
  # Prefer seurat_clusters column if present; else use Idents(obj)
  if ("seurat_clusters" %in% colnames(obj@meta.data)) {
    as.character(obj@meta.data$seurat_clusters)
  } else {
    as.character(Idents(obj))
  }
}

assay_summary <- function(obj, assay) {
  s <- list(assay = assay)
  s$exists <- assay %in% Assays(obj)
  if (!s$exists) return(s)
  DefaultAssay(obj) <- assay
  # Try data/scale slots
  s$data_dims <- tryCatch({
    m <- GetAssayData(obj, assay = assay, slot = "data")
    c(nrow = nrow(m), ncol = ncol(m))
  }, error = function(e) c(nrow = NA_integer_, ncol = NA_integer_))
  s$counts_dims <- tryCatch({
    m <- GetAssayData(obj, assay = assay, slot = "counts")
    c(nrow = nrow(m), ncol = ncol(m))
  }, error = function(e) c(nrow = NA_integer_, ncol = NA_integer_))
  s$scale_dims <- tryCatch({
    m <- GetAssayData(obj, assay = assay, slot = "scale.data")
    c(nrow = nrow(m), ncol = ncol(m))
  }, error = function(e) c(nrow = NA_integer_, ncol = NA_integer_))
  s$var_features_n <- tryCatch(length(VariableFeatures(obj)), error = function(e) NA_integer_)
  s
}

reductions_summary <- function(obj) {
  reds <- names(obj@reductions)
  lapply(reds, function(r) {
    x <- obj@reductions[[r]]
    list(name = r,
         key = tryCatch(x@key, error = function(e) NA_character_),
         dims = tryCatch(ncol(Embeddings(x)), error = function(e) NA_integer_),
         cells = tryCatch(nrow(Embeddings(x)), error = function(e) NA_integer_))
  })
}

graphs_summary <- function(obj) {
  g <- names(obj@graphs)
  lapply(g, function(n) {
    m <- obj@graphs[[n]]
    list(name = n,
         dims = tryCatch(dim(m), error = function(e) c(NA_integer_, NA_integer_)))
  })
}

meta_overview <- function(md) {
  tibble(
    column = colnames(md),
    class  = sapply(md, function(x) paste(class(x), collapse = "|")),
    n_na   = sapply(md, function(x) sum(is.na(x))),
    n_unique = sapply(md, function(x) length(unique(as.character(x))))
  )
}

# ---------- Main ----------
main <- function(opt) {
  if (!file.exists(opt$rds)) {
    stop(sprintf("Input RDS not found: %s", opt$rds))
  }
  outroot <- opt$outdir
  data_dir <- file.path(outroot, "data")
  logs_dir <- file.path(outroot, "logs")
  dirs <- c(outroot, data_dir, logs_dir)
  for (d in dirs) dir.create(d, showWarnings = FALSE, recursive = TRUE)

  ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
  log_file <- file.path(logs_dir, sprintf("inspect_rds_%s.txt", ts))
  tp_csv   <- file.path(data_dir,  sprintf("inspect_timepoint_counts_%s.csv", ts))
  cl_csv   <- file.path(data_dir,  sprintf("inspect_cluster_counts_%s.csv", ts))
  meta_csv <- file.path(data_dir,  sprintf("inspect_meta_overview_%s.csv", ts))

  con <- file(log_file, open = "wt")

  wlog <- function(...) writeLines(paste0(...), con = con)

  wlog(sprintf("RDS: %s", opt$rds))
  obj <- readRDS(opt$rds)
  wlog(sprintf("Object class: %s", paste(class(obj), collapse = "|")))
  wlog(sprintf("Cells (ncol): %d; Features (nrow RNA default): %s", ncol(obj), tryCatch(nrow(obj), error = function(e) NA_integer_)))

  # Assays
  assays <- Assays(obj)
  wlog(sprintf("Assays: %s", paste(assays, collapse = ", ")))
  wlog(sprintf("DefaultAssay: %s", DefaultAssay(obj)))

  summaries <- lapply(assays, function(a) assay_summary(obj, a))
  for (s in summaries) {
    wlog(sprintf("Assay '%s' exists=%s; data_dims=%s; counts_dims=%s; scale_dims=%s; var_features_n=%s",
                 s$assay, s$exists,
                 paste(s$data_dims, collapse = "x"),
                 paste(s$counts_dims, collapse = "x"),
                 paste(s$scale_dims, collapse = "x"),
                 s$var_features_n))
  }

  # Reductions
  red_sum <- reductions_summary(obj)
  if (length(red_sum) > 0) {
    wlog("Reductions:")
    for (r in red_sum) {
      wlog(sprintf(" - %s: key=%s, dims=%s, cells=%s", r$name, r$key, r$dims, r$cells))
    }
  } else {
    wlog("Reductions: None")
  }

  # Graphs
  gr_sum <- graphs_summary(obj)
  if (length(gr_sum) > 0) {
    wlog("Graphs:")
    for (g in gr_sum) {
      wlog(sprintf(" - %s: dims=%s", g$name, paste(g$dims, collapse = "x")))
    }
  } else {
    wlog("Graphs: None")
  }

  # Idents and clusters
  id_len_before <- safe_len(Idents(obj))
  wlog(sprintf("Idents length before: %d", id_len_before))
  clusters_vec <- safe_idents(obj)
  wlog(sprintf("Cluster vector length: %d; unique clusters: %d", length(clusters_vec), length(unique(clusters_vec))))
  cl_tbl <- tibble(cluster = clusters_vec) %>% count(cluster, name = "count") %>% arrange(desc(count))
  write_csv(cl_tbl, cl_csv)
  wlog(sprintf("Cluster counts written: %s", cl_csv))

  # Meta overview
  md <- obj@meta.data
  wlog(sprintf("meta.data dims: %dx%d", nrow(md), ncol(md)))
  wlog(sprintf("meta.data columns: %s", paste(colnames(md), collapse = ", ")))
  meta_ov <- meta_overview(md)
  write_csv(meta_ov, meta_csv)
  wlog(sprintf("Meta overview written: %s", meta_csv))

  # Timepoint detection and counts
  tp_cands <- detect_timepoint_candidates(md)
  wlog(sprintf("Timepoint candidates present: %s", paste(tp_cands, collapse = ", ")))
  tp_col <- detect_timepoint_col(md, preferred = "")
  wlog(sprintf("Selected timepoint column: %s", tp_col))
  if (!is.na(tp_col)) {
    tp_vals <- md[[tp_col]]
    # coerce to character safely
    tp_chr <- as.character(tp_vals)
    tp_tbl <- tibble(timepoint = tp_chr) %>%
      mutate(is_na = is.na(timepoint) | timepoint == "") %>%
      count(timepoint, is_na, name = "count") %>%
      arrange(desc(count))
    write_csv(tp_tbl, tp_csv)
    wlog(sprintf("Timepoint counts written: %s", tp_csv))
    wlog(sprintf("Unique timepoints (raw): %s", paste(unique(tp_chr), collapse = ", ")))
    wlog(sprintf("NA/empty count: %d", sum(is.na(tp_chr) | tp_chr == "")))
  } else {
    wlog("No suitable timepoint column found.")
  }

  # Potential pitfalls summary
  wlog("Potential pitfalls summary:")
  if (!is.na(tp_col)) {
    na_tp <- sum(is.na(md[[tp_col]]) | as.character(md[[tp_col]]) == "")
    if (na_tp > 0) {
      wlog(sprintf(" - Timepoint column '%s' contains NA/empty values: %d", tp_col, na_tp))
    } else {
      wlog(sprintf(" - Timepoint column '%s' has no NA/empty", tp_col))
    }
  }
  if (length(unique(clusters_vec)) <= 1) {
    wlog(" - Only one unique cluster; plots may be trivial.")
  } else {
    wlog(" - Multiple clusters detected; OK.")
  }

  # Session info
  wlog("SessionInfo:")
  wlog(paste(capture.output(sessionInfo()), collapse = "\n"))

  close(con)

  message(sprintf("[INFO] Inspection completed. Logs: %s", log_file))
}

tryCatch(
  { main(opt) },
  error = function(e) {
    message("[ERROR] ", e$message)
    quit(status = 1)
  }
)

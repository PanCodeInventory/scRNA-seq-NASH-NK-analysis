#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(scales)
  library(optparse)
})

# ---------- Helpers ----------
stop_if_missing_cols <- function(md, cols) {
  miss <- cols[!cols %in% colnames(md)]
  if (length(miss) > 0) stop(sprintf("Missing columns: %s", paste(miss, collapse = ", ")))
}

detect_timepoint_col <- function(obj, preferred = "") {
  md <- obj@meta.data
  if (preferred != "" && preferred %in% colnames(md)) return(preferred)
  candidates <- c("timepoint", "Timepoint", "time", "Time", "orig.ident", "group", "condition", "sample", "dataset")
  found <- candidates[candidates %in% colnames(md)]
  if (length(found) == 0) {
    stop("No suitable timepoint column found. Please provide --timepoint-col explicitly.")
  }
  return(found[[1]])
}

derive_timepoint_order <- function(values, user_order = "") {
  v <- as.character(values)
  # drop NA and empty strings to avoid logical NA in indexing
  u <- unique(v[!is.na(v) & v != ""])
  if (length(u) == 0) return(character(0))
  if (user_order != "") {
    ord <- trimws(unlist(strsplit(user_order, ",")))
    ord <- ord[ord != ""]
    # keep only those present; append any missing at the end by their original order
    ord_fin <- unique(c(ord[ord %in% u], u[!(u %in% ord)]))
    return(ord_fin)
  }
  default_ord <- c("NCD", "MCD-1W", "MCD-2W", "MCD-6W")
  if (all(default_ord %in% u)) return(default_ord)
  # fallback: keep encounter order without NA
  return(u)
}

ensure_idents <- function(obj, cluster_col = "") {
  if (cluster_col != "") {
    if (!cluster_col %in% colnames(obj@meta.data)) {
      stop(sprintf("Specified --cluster-col '%s' not found in meta.data.", cluster_col))
    }
    Idents(obj) <- obj@meta.data[[cluster_col]]
  } else {
    # if Idents not set meaningfully, try common columns
    id_vec <- Idents(obj)
    if (length(unique(as.character(id_vec))) <= 1) {
      candidates <- c("seurat_clusters", "clusters", "cluster", "Cluster")
      found <- candidates[candidates %in% colnames(obj@meta.data)]
      if (length(found) > 0) Idents(obj) <- obj@meta.data[[found[[1]]]]
    }
  }
  obj
}

make_line_plot <- function(df_prop, topk = 12) {
  # choose top clusters by mean proportion
  avg <- df_prop %>%
    group_by(cluster) %>%
    summarise(mean_prop = mean(proportion, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(mean_prop))
  keep <- if (isTRUE(topk > 0)) head(avg$cluster, topk) else sort(unique(df_prop$cluster))
  dfp <- df_prop %>% filter(cluster %in% keep)

  ggplot(dfp, aes(x = timepoint, y = proportion, color = cluster, group = cluster)) +
    geom_line(linewidth = 0.7, alpha = 0.9) +
    geom_point(size = 1.5) +
    scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
    labs(x = "Timepoint", y = "Proportion", color = "Cluster", title = "Cluster proportion over time") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "right",
          panel.grid.minor = element_blank())
}

make_stacked_bar <- function(df_counts) {
  ggplot(df_counts, aes(x = timepoint, y = count, fill = cluster)) +
    geom_bar(stat = "identity", position = "fill") +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    labs(x = "Timepoint", y = "Composition (proportion)", fill = "Cluster",
         title = "Cluster composition per timepoint (stacked)") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "right",
          panel.grid.minor = element_blank())
}

write_run_logs <- function(logs_dir, config_lines) {
  ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
  cfg_file <- file.path(logs_dir, sprintf("run_config_%s.txt", ts))
  ses_file <- file.path(logs_dir, sprintf("sessionInfo_%s.txt", ts))
  writeLines(config_lines, con = cfg_file)
  writeLines(capture.output(sessionInfo()), con = ses_file)
  invisible(list(config = cfg_file, session = ses_file))
}

# ---------- CLI ----------
option_list <- list(
  make_option(c("--rds"), type = "character", default = "1_Files/RDS/nk1.1_integrated.tuned.rds",
              help = "Input Seurat RDS path [default %default]"),
  make_option(c("--outdir"), type = "character", default = "3_Analysis/1.ClusterAnalysis",
              help = "Output root directory for data/plots/logs [default %default]"),
  make_option(c("--timepoint-col"), type = "character", default = "",
              help = "Timepoint column name in meta.data; auto-detect if empty"),
  make_option(c("--cluster-col"), type = "character", default = "",
              help = "Cluster column to set as Idents; use current Idents if empty"),
  make_option(c("--timepoint-order"), type = "character", default = "",
              help = "Comma-separated order for timepoints (e.g. 'NCD,MCD-1W,MCD-2W,MCD-6W')"),
  make_option(c("--topk"), type = "integer", default = 12,
              help = "Top K clusters to show in line plot by mean proportion; 0=all [default %default]"),
  make_option(c("--formats"), type = "character", default = "png,pdf",
              help = "Comma-separated output formats for plots [default %default]"),
  make_option(c("--width"), type = "double", default = 9,
              help = "Plot width (inches) [default %default]"),
  make_option(c("--height"), type = "double", default = 6,
              help = "Plot height (inches) [default %default]"),
  make_option(c("--dpi"), type = "integer", default = 300,
              help = "Plot DPI [default %default]")
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
  message(sprintf("[DEBUG] meta.data rows: %d, cols: %d", nrow(obj@meta.data), ncol(obj@meta.data)))
  # Clusters vector: prefer seurat_clusters from meta.data; fallback to Idents(obj)
  cluster_vec <- if ("seurat_clusters" %in% colnames(obj@meta.data)) as.character(obj@meta.data$seurat_clusters) else as.character(Idents(obj))
  message(sprintf("[DEBUG] cluster_vec length: %d, unique clusters: %d", length(cluster_vec), length(unique(cluster_vec))))
  if (is.null(cluster_vec) || length(cluster_vec) == 0) {
    stop("Cluster vector is empty; please specify --cluster-col explicitly or check the object.")
  }
  if (length(unique(cluster_vec)) <= 1) {
    message("[WARN] Only one unique cluster detected; proceeding may produce trivial plots.")
  }

  # Timepoint detection and ordering (robust with tryCatch)
  timepoint_block_success <- TRUE
  tryCatch({
    tp_col <- detect_timepoint_col(obj, opt$timepoint_col)
    tp_vals <- obj@meta.data[[tp_col]]
    tp_chr <- as.character(tp_vals)
    present <- unique(tp_chr[!is.na(tp_chr) & tp_chr != ""])
    message(sprintf("[DEBUG] timepoint_col: %s; unique values (raw): %s", tp_col, paste(unique(tp_chr), collapse = ", ")))
    message(sprintf("[DEBUG] present timepoints (non-empty): %s", paste(present, collapse = ", ")))
    tp_order <- derive_timepoint_order(tp_vals, opt$timepoint_order)
    message(sprintf("[DEBUG] derived timepoint order (clean): %s", paste(tp_order, collapse = ", ")))
    # ensure final order only includes present values, then append any remaining present values by encounter order
    tp_ord_fin <- unique(c(tp_order[tp_order %in% present], present[!(present %in% tp_order)]))
    message(sprintf("[DEBUG] final timepoint order used: %s", paste(tp_ord_fin, collapse = ", ")))
    if (length(tp_ord_fin) == 0) {
      tp_ord_fin <- present
      message(sprintf("[WARN] Final timepoint order empty; fallback to present values: %s", paste(tp_ord_fin, collapse = ", ")))
    }

    cell_names <- rownames(obj@meta.data)
    if (is.null(cell_names) || length(cell_names) == 0) {
      cell_names <- colnames(obj)
    }
    # Build meta using raw character timepoint, then filter and factor
    keep_mask <- !is.na(tp_chr) & tp_chr != ""
    df_meta <- data.frame(
      cell = cell_names,
      cluster = cluster_vec,
      timepoint = tp_chr,
      stringsAsFactors = FALSE
    )
    df_meta <- df_meta[keep_mask, , drop = FALSE]
    tp_factor <- factor(df_meta$timepoint, levels = tp_ord_fin)
    df_meta$timepoint <- tp_factor
    message(sprintf("[DEBUG] df_meta dims: %dx%d; first timepoints: %s", nrow(df_meta), ncol(df_meta), paste(head(unique(df_meta$timepoint)), collapse = ", ")))
  }, error = function(e) {
    timepoint_block_success <<- FALSE
    message(sprintf("[ERROR] timepoint detection block failed: %s", e$message))
    # Fallback: select a reasonable timepoint column and build meta safely
    tp_col <<- if ("timepoint" %in% colnames(obj@meta.data)) "timepoint" else if ("orig.ident" %in% colnames(obj@meta.data)) "orig.ident" else colnames(obj@meta.data)[1]
    tp_vals <<- obj@meta.data[[tp_col]]
    tp_chr <<- as.character(tp_vals)
    cell_names <<- rownames(obj@meta.data)
    if (is.null(cell_names) || length(cell_names) == 0) {
      cell_names <<- colnames(obj)
    }
    present <<- unique(tp_chr[!is.na(tp_chr) & tp_chr != ""])
    if (length(present) == 0) present <<- unique(tp_chr)
    tp_ord_fin <<- present
    keep_mask <<- !is.na(tp_chr) & tp_chr != ""
    df_meta <<- data.frame(
      cell = cell_names,
      cluster = cluster_vec,
      timepoint = tp_chr,
      stringsAsFactors = FALSE
    )
    df_meta <<- df_meta[keep_mask, , drop = FALSE]
    df_meta$timepoint <<- factor(df_meta$timepoint, levels = tp_ord_fin)
    message(sprintf("[WARN] Fallback path used. timepoint_col=%s; levels=%s", tp_col, paste(tp_ord_fin, collapse = ", ")))
  })

  # Counts and proportions
  df_counts <- df_meta %>% count(timepoint, cluster, name = "count")
  message(sprintf("[DEBUG] df_counts rows: %d; unique timepoints: %d; unique clusters: %d",
                  nrow(df_counts), length(unique(df_counts$timepoint)), length(unique(df_counts$cluster))))
  df_totals <- df_counts %>% group_by(timepoint) %>%
    summarise(total = sum(count), .groups = "drop")
  df_prop <- df_counts %>%
    left_join(df_totals, by = "timepoint") %>%
    mutate(proportion = ifelse(total > 0, count / total, NA_real_))
  message(sprintf("[DEBUG] df_prop rows: %d; any NA proportion: %s",
                  nrow(df_prop), any(is.na(df_prop$proportion))))

  # Write CSVs
  counts_csv <- file.path(data_dir, "cluster_counts_by_timepoint.csv")
  props_csv  <- file.path(data_dir, "cluster_proportions_by_timepoint.csv")
  write_csv(df_counts %>% arrange(timepoint, cluster), counts_csv)
  write_csv(df_prop %>% arrange(timepoint, cluster),  props_csv)
  message(sprintf("[INFO] Wrote: %s", counts_csv))
  message(sprintf("[INFO] Wrote: %s", props_csv))

  # Plots
  formats <- trimws(unlist(strsplit(opt$formats, ",")))
  formats <- formats[formats != ""]

  p_line <- make_line_plot(df_prop, topk = opt$topk)
  for (fmt in formats) {
    out <- file.path(plots_dir, sprintf("cluster_proportion_lineplot.%s", fmt))
    ggsave(out, p_line, width = opt$width, height = opt$height, dpi = opt$dpi, units = "in", bg = "white")
    message(sprintf("[INFO] Saved: %s", out))
  }

  p_bar <- make_stacked_bar(df_counts)
  for (fmt in formats) {
    out <- file.path(plots_dir, sprintf("cluster_composition_stackedbar.%s", fmt))
    ggsave(out, p_bar, width = opt$width, height = opt$height, dpi = opt$dpi, units = "in", bg = "white")
    message(sprintf("[INFO] Saved: %s", out))
  }

  # Logs (robust)
  timepoint_order_str <- ""
  if (exists("tp_ord_fin")) {
    timepoint_order_str <- paste(tp_ord_fin, collapse = ",")
  } else if (exists("tp_order")) {
    timepoint_order_str <- paste(tp_order, collapse = ",")
  } else {
    # last resort: use levels from df_meta if present
    if (exists("df_meta") && "timepoint" %in% colnames(df_meta)) {
      timepoint_order_str <- paste(levels(df_meta$timepoint), collapse = ",")
    }
  }
  config_lines <- c(
    sprintf("rds: %s", opt$rds),
    sprintf("outdir: %s", opt$outdir),
    sprintf("timepoint_col: %s", if (exists("tp_col")) tp_col else ""),
    sprintf("cluster_source: %s", if (opt$cluster_col == "") "Idents(obj)" else opt$cluster_col),
    sprintf("timepoint_order_used: %s", timepoint_order_str),
    sprintf("topk: %s", opt$topk),
    sprintf("formats: %s", paste(formats, collapse = ",")),
    sprintf("n_cells: %s", ncol(obj)),
    sprintf("n_clusters: %s", length(unique(cluster_vec)))
  )
  tryCatch(
    { write_run_logs(logs_dir, config_lines) },
    error = function(e) {
      message(sprintf("[WARN] write_run_logs failed: %s", e$message))
    }
  )

  message("[INFO] Done.")
}

tryCatch(
  { main(opt) },
  error = function(e) {
    message("[ERROR] ", e$message)
    quit(status = 1)
  }
)

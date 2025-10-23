#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(scales)
  library(optparse)
  library(viridis)
  library(RColorBrewer)
})

# ---------- 功能函数 ----------
create_stacked_bar_with_labels <- function(df_props,
                                          min_label_threshold = 0.02,
                                          color_palette = "viridis",
                                          title = "Cluster Composition by Timepoint",
                                          x_label = "Timepoint",
                                          y_label = "Proportion",
                                          label_size = 3,
                                          legend_position = "right",
                                          top_n = 5) {
  
  # 计算累积比例用于标签定位，按cluster顺序排列（0,1,2,3,4,5,6）
  df_labeled <- df_props %>%
    group_by(timepoint) %>%
    arrange(cluster) %>%
    mutate(
      # 重新计算累积比例，确保标签位置正确
      cumsum_low = lag(cumsum(proportion), default = 0),
      cumsum_high = cumsum(proportion),
      label_y = (cumsum_low + cumsum_high) / 2,  # 标签位置在每个堆叠的中间
      cluster = as.character(cluster)
    ) %>%
    ungroup()
  
  # 只为前5个cluster添加标签
  df_labeled <- df_labeled %>%
    group_by(timepoint) %>%
    mutate(
      rank = rank(-proportion, ties.method = "first"),
      label_text = ifelse(rank <= top_n & proportion >= min_label_threshold,
                         percent(proportion, accuracy = 0.1),
                         "")
    ) %>%
    ungroup() %>%
    select(-rank)
  
  # 设置时间点顺序
  timepoint_order <- c("0W_NCD", "1W_MCD", "2W_MCD", "6W_MCD")
  df_labeled$timepoint <- factor(df_labeled$timepoint, levels = timepoint_order)
  
  # 创建颜色映射
  if (color_palette == "viridis") {
    colors <- viridis(length(unique(df_labeled$cluster)))
  } else if (color_palette == "Set3") {
    colors <- brewer.pal(max(3, length(unique(df_labeled$cluster))), "Set3")
  } else if (color_palette == "rainbow") {
    colors <- rainbow(length(unique(df_labeled$cluster)))
  } else {
    colors <- color_palette
  }
  
  # 创建堆叠柱状图
  p <- ggplot(df_labeled, aes(x = timepoint, y = proportion, fill = cluster)) +
    geom_bar(stat = "identity", position = "fill", color = "white", linewidth = 0.5) +
    geom_text(aes(y = label_y, label = label_text),
              size = label_size + 1,
              color = "white",
              fontface = "bold",
              position = position_identity(),
              show.legend = FALSE,
              check_overlap = FALSE) +
    scale_y_continuous(labels = percent_format(accuracy = 1), 
                       limits = c(0, 1),
                       expand = c(0, 0)) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_fill_manual(values = colors) +
    labs(title = title,
         x = x_label,
         y = y_label,
         fill = "Cluster") +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = legend_position,
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      panel.border = element_rect(fill = NA, color = "gray80", linewidth = 0.5)
    )
  
  return(p)
}

# ---------- 命令行参数 ----------
option_list <- list(
  make_option(c("--input"), type = "character", 
              default = "3_Analysis/1.ClusterAnalysis/data/cluster_proportions_by_timepoint.csv",
              help = "输入CSV文件路径 [default %default]"),
  make_option(c("--output"), type = "character", 
              default = "3_Analysis/1.ClusterAnalysis/plots",
              help = "输出目录 [default %default]"),
  make_option(c("--formats"), type = "character", default = "png,pdf",
              help = "输出格式，逗号分隔 [default %default]"),
  make_option(c("--width"), type = "double", default = 10,
              help = "图表宽度(英寸) [default %default]"),
  make_option(c("--height"), type = "double", default = 6,
              help = "图表高度(英寸) [default %default]"),
  make_option(c("--dpi"), type = "integer", default = 300,
              help = "图表DPI [default %default]"),
  make_option(c("--min_threshold"), type = "double", default = 0.02,
              help = "显示标签的最小比例阈值 [default %default]"),
  make_option(c("--color_palette"), type = "character", default = "viridis",
              help = "颜色方案 (viridis, Set3, rainbow) [default %default]"),
  make_option(c("--title"), type = "character", 
              default = "Cluster Composition by Timepoint",
              help = "图表标题 [default %default]"),
  make_option(c("--label_size"), type = "integer", default = 3,
              help = "标签字体大小 [default %default]"),
  make_option(c("--top_n"), type = "integer", default = 5,
              help = "每个时间点显示的cluster数量 [default %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

# ---------- 主函数 ----------
main <- function(opt) {
  # 检查输入文件
  if (!file.exists(opt$input)) {
    stop(sprintf("输入文件不存在: %s", opt$input))
  }
  
  # 创建输出目录
  dir.create(opt$output, showWarnings = FALSE, recursive = TRUE)
  
  # 读取数据
  message(sprintf("[INFO] 读取数据: %s", opt$input))
  df_props <- read_csv(opt$input, show_col_types = FALSE)
  
  # 检查数据结构
  required_cols <- c("timepoint", "cluster", "proportion")
  missing_cols <- setdiff(required_cols, colnames(df_props))
  if (length(missing_cols) > 0) {
    stop(sprintf("缺少必需的列: %s", paste(missing_cols, collapse = ", ")))
  }
  
  message(sprintf("[INFO] 数据包含 %d 个时间点, %d 个cluster", 
                  length(unique(df_props$timepoint)), 
                  length(unique(df_props$cluster))))
  
  # 创建图表
  message("[INFO] 创建带百分比标签的堆叠柱状图")
  p <- create_stacked_bar_with_labels(
    df_props = df_props,
    min_label_threshold = opt$min_threshold,
    color_palette = opt$color_palette,
    title = opt$title,
    label_size = opt$label_size,
    top_n = opt$top_n
  )
  
  # 保存图表
  formats <- trimws(unlist(strsplit(opt$formats, ",")))
  formats <- formats[formats != ""]
  
  for (fmt in formats) {
    output_file <- file.path(opt$output, sprintf("cluster_composition_stackedbar_with_labels.%s", fmt))
    ggsave(output_file, p, 
           width = opt$width, 
           height = opt$height, 
           dpi = opt$dpi, 
           units = "in", 
           bg = "white")
    message(sprintf("[INFO] 保存图表: %s", output_file))
  }
  
  # 显示数据摘要
  message("[INFO] 数据摘要:")
  summary_table <- df_props %>%
    group_by(timepoint) %>%
    summarise(
      total_cells = sum(count),
      n_clusters = n(),
      .groups = "drop"
    )
  print(summary_table)
  
  message("[INFO] 完成!")
}

# 执行主函数
tryCatch(
  { main(opt) },
  error = function(e) {
    message(sprintf("[ERROR] %s", e$message))
    quit(status = 1)
  }
)
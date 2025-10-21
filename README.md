# scRNA-seq NASH NK细胞分析项目

## 项目概述

本项目专注于NASH（非酒精性脂肪性肝炎）疾病模型中NK细胞的单细胞RNA测序分析。通过分析不同时间点（0周、1周、2周、6周）的NK细胞样本，研究NASH疾病进程中NK细胞的变化规律。

## 项目结构

```
scRNA-seq/
├── README.md                 # 项目说明文档
├── .gitignore               # Git忽略文件配置
├── Files/                   # 项目数据和分析文件
│   ├── Filter Files/        # 过滤后的10x Genomics数据
│   │   ├── NCD_NK1.1/       # 正常对照组NK1.1样本
│   │   ├── MCD-1W_NK1.1/    # 1周MCD组NK1.1样本
│   │   ├── MCD-2W_NK1.1/    # 2周MCD组NK1.1样本
│   │   └── MCD-6W_NK1.1/    # 6周MCD组NK1.1样本
│   ├── UMAP/                # UMAP分析相关文件
│   │   ├── scripts/         # R分析脚本
│   │   ├── Results/         # 分析结果
│   │   │   ├── plots/       # 可视化图表
│   │   │   ├── data/        # 结果数据
│   │   │   └── rds/         # R对象文件
│   │   ├── tuning/          # 参数调优相关文件
│   │   └── *.md             # 分析文档
│   └── Annotation/          # 注释文件
└── Demos/                   # 演示文件
```

## 样本信息

### 实验分组
- **NCD (正常对照组)**: 0周时间点
- **MCD (模型组)**: 1周、2周、6周时间点

### 细胞类型
- **NK1.1细胞**: 自然杀伤细胞，标记为NK1.1
- **CD45.2细胞**: 白细胞共同抗原标记细胞

## 分析流程

### 1. 数据预处理
- 10x Genomics原始数据质控
- 细胞过滤和质量控制
- 标准化处理

### 2. 数据整合
- 使用Seurat包进行SCTransform标准化
- 锚点整合多个样本
- 创建整合后的Seurat对象

### 3. 降维与聚类
- 主成分分析（PCA）
- UMAP降维
- 细胞聚类分析

### 4. 可视化分析
- 生成各时间点的UMAP图
- 聚类结果可视化
- 趋势分析图表

## 主要脚本

- `Files/UMAP/scripts/generate_umap_nk.R` - 主要UMAP分析脚本
- `Files/UMAP/scripts/generate_umap_nk_post.R` - 后处理脚本
- `Files/UMAP/scripts/tune_nk_dims_resolution.R` - 参数调优脚本
- `Files/UMAP/scripts/generate_umap_nk_equalized.R` - 均衡化分析脚本

## 主要结果

### 图表文件
- UMAP聚类图（各时间点）
- 聚类组成趋势图
- 肘部图（PCA维度选择）
- 调优热图

### 数据文件
- 聚类组成趋势数据
- 调优指标数据
- 整合后的R对象

## 技术栈

- **R语言**: 主要分析语言
- **Seurat**: 单细胞分析包
- **ggplot2**: 数据可视化
- **patchwork**: 图形组合

## 使用说明

1. 确保安装了必要的R包：
   ```r
   install.packages(c("Seurat", "dplyr", "ggplot2", "patchwork"))
   ```

2. 运行主要分析脚本：
   ```r
   source("Files/UMAP/scripts/generate_umap_nk.R")
   ```

3. 查看分析文档：
   - `Files/UMAP/guidedoc.md` - 详细分析指南
   - `Files/UMAP/guidedoc_resolution.md` - 分辨率调优指南
   - `Files/UMAP/guidedoc_UMAP_Normalization.md` - 标准化指南

## 注意事项

- 原始数据文件较大，已配置.gitignore排除
- 结果文件（图片、R对象等）也排除在版本控制外
- 仅保留脚本和文档文件进行版本控制
- 建议在运行分析前确保有足够的存储空间

## 更新

- 2025-10-20
  - 新增脚本：`Files/UMAP/scripts/remove_doublets_and_contaminants.R`（集成 scDblFinder 去双胞、SingleR 自动注释、UCell 签名评分、去污染规则与重分析的主流程）
  - 修复 Seurat v5 多层 assay 转换为 SCE 的问题，增强元数据行名对齐与日志/报告目录创建的健壮性
  - 前台运行并迭代修复，产出早期日志与统计

- 2025-10-20 深夜
  - 新增特异性 SingleR 修复脚本：`Files/UMAP/scripts/singleR_annotation_fix.R`（从 counts 构建 SCE，scater::logNormCounts 生成 logcounts，显式以 logcounts 作为 SingleR 输入，规避 data 层为空告警）
  - 产出对象与文档：
    - `Files/Doublet_Removed/RDS/nk.integrated.singleR_annotated.rds`
    - `Files/Doublet_Removed/reports/singleR_fix_report.md`
    - `Files/Doublet_Removed/plots/SingleR_label_barplot.png`

- 2025-10-21
  - 将 SingleR 修复策略集成至主流程并完成全流程清理与重分析（保留 NK/ILC）
  - 去双胞结果：移除 312 个细胞（约 1.61%）
  - 生成清理后对象与报告：
    - `Files/Doublet_Removed/RDS/nk.integrated.filtered.rds`
    - `Files/Doublet_Removed/RDS/nk.integrated.doublet_scored.rds`
    - `Files/Doublet_Removed/reports/cleaning_report.md`
  - 生成图件：
    - `Files/Doublet_Removed/plots/DoubletScore_Density.png`
    - `Files/Doublet_Removed/plots/DoubletScore_Density_bySample.png`
    - `Files/Doublet_Removed/plots/Cluster_nonNK_fraction.png`
    - `Files/Doublet_Removed/plots/UMAP_filtered_by_SingleR.png`
    - `Files/Doublet_Removed/plots/UMAP_filtered_clusters_by_Timepoint.png`（若含 timepoint）
    - `Files/Doublet_Removed/plots/NK_UCell_by_SingleR_label.png`
    - `Files/Doublet_Removed/plots/Cluster_UCell_signature_means.png`

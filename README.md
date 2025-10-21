# scRNA-seq NASH NK细胞分析项目

## 项目概述

本项目专注于 NASH（非酒精性脂肪性肝炎）疾病模型中 NK 细胞的单细胞 RNA 测序分析。通过分析不同时间点（0 周、1 周、2 周、6 周）的 NK 细胞样本，研究 NASH 疾病进程中 NK 细胞的变化规律。

## 项目结构（当前主路径：2_DataProcessing/*）

```
scRNA-seq/
├── README.md
├── .gitignore
├── 1_Files/                          # 原始/预处理数据（按分组）
│   ├── NK1.1/
│   └── CD45.2/
├── 2_DataProcessing/                 # 主数据处理管线（脚本与产物）
│   ├── 1_Samples_Merging/
│   │   ├── Scripts & guidedoc*.md
│   │   └── Results/{data,plots,rds}
│   ├── 2_Doublet_Removed/            # 去双胞/注释/清理产出（报告与图件）
│   │   ├── RDS/                      # R 对象（清理后 / 单细胞注释 / noNKT / tuned）
│   │   ├── plots/                    # UMAP 等图件
│   │   └── reports/                  # 报告（cleaning / singleR 修复 / noNKT）
│   ├── 3_UMAP-Tuning/                # UMAP 调参与选择产出（metrics/plots/logs）
│   │   ├── data/                     # 调参指标与候选 CSV
│   │   ├── plots/                    # 调参热力图与候选 UMAP
│   │   └── logs/                     # 运行配置与会话信息
│   └── Scripts/                      # 脚本（生成、清理、调参）
│       ├── generate_umap_nk.R
│       ├── generate_umap_nk_post.R
│       ├── remove_doublets_and_contaminants.R
│       ├── singleR_annotation_fix.R
│       ├── remove_NKT_cells.R                  # 新增：在已注释对象上剔除 NKT
│       └── tune_noNKT_dims_resolution.R        # 新增：基于 noNKT 对象进行 dims × resolution 调参并重跑
├── 2_Filter/                         # 可选镜像产出目录（按你的偏好保留）
│   └── 2_Doublet_Removed/{RDS,plots,reports}
└── 3_Analysis/                       # 下游分析（待扩展）
```

历史路径兼容（Files/*）说明：
- 早期版本产物位于 `Files/UMAP/*` 与 `Files/Doublet_Removed/*`。当前主路径已迁移至 `2_DataProcessing/*`，新产物与脚本请以该路径为准。

## 样本与分组信息
- 分组：NCD（0W）与 MCD（1W/2W/6W）
- 细胞类型：NK1.1（自然杀伤细胞）与 CD45.2（白细胞共同抗原）

## 分析流程概览

1) 样本合并与整合（SCTransform + Anchors）
2) 去双胞（scDblFinder，按样本/时间点分组）
3) 自动注释（SingleR，logcounts 修复策略）
4) UCell 基因签名评分（NK/T/B/Myeloid/DC/Plasma/Endothelium/Fibroblast/Hepatocyte）
5) 去污染（细胞级阈值 + 簇级非 NK 占比阈值）
6) 重跑降维/聚类/UMAP（兼容 SCT/RNA 多模型，必要时回退）
7) NKT 剔除（基于 SingleR 标签严格规则）
8) dims × resolution 调参（UMAP/聚类）与最终参数选择
9) 按最终参数生成分面 UMAP（timepoint）与标签 UMAP（SingleR）

## 关键技术与兼容策略

- 平台/框架：R、Seurat 5.x、SingleCellExperiment、scDblFinder、SingleR、celldex、scater、UCell、ggplot2、patchwork
- Seurat v5 多层 assay 差异：
  - 优先使用 layer 接口获取数据，必要时回退 slot 接口；兼容 RNA/SCT/integrated
  - 避免直接 `as.SingleCellExperiment(seu)`，显式构建 SCE 并确保 counts/logcounts 与 colData 对齐
- SingleR 修复策略：
  - 从 Seurat 提取 counts → SCE → `scater::logNormCounts` 生成 logcounts → SingleR 显式 `assay.type="logcounts"`
- 去污染判定：
  - 细胞级：NK_UCell ≥ P60 且 Δ(NK − max(others)) ≥ 0.05，或 SingleR 注释命中 NK/ILC
  - 簇级：簇内非 NK 注释占比 ≥ 0.7 则整体剔除
- 降维稳健性：
  - RNA：NormalizeData → FindVariableFeatures → ScaleData → RunPCA
  - SCT：直接使用 VariableFeatures；若 SCT 无 VF 则回退 RNA 并自动计算 HVG

## 使用说明（一键运行关键步骤）

1) 在已完成 SingleR 注释的对象上剔除 NKT
- 输入：`2_DataProcessing/2_Doublet_Removed/RDS/nk.integrated.singleR_annotated.rds`
- 运行：
  ```bash
  Rscript 2_DataProcessing/Scripts/remove_NKT_cells.R
  ```
- 产出：
  - `2_DataProcessing/2_Doublet_Removed/RDS/nk.integrated.noNKT.rds`
  - `2_DataProcessing/2_Doublet_Removed/reports/noNKT_removal_report.md`
  - `2_DataProcessing/2_Doublet_Removed/reports/removed_NKT_cell_ids.csv`
  - `2_DataProcessing/2_Doublet_Removed/plots/NKT_removal_label_counts_before_after.png`

2) 基于 noNKT 对象进行 dims × resolution 调参并重跑
- 输入：`2_DataProcessing/2_Doublet_Removed/RDS/nk.integrated.noNKT.rds`
- 运行：
  ```bash
  Rscript 2_DataProcessing/Scripts/tune_noNKT_dims_resolution.R
  ```
- 产出：
  - 调参指标：`2_DataProcessing/3_UMAP-Tuning/data/nk_noNKT_tuning_metrics.csv`
  - 候选组合：`2_DataProcessing/3_UMAP-Tuning/data/nk_noNKT_tuning_best_per_dims.csv`
               `2_DataProcessing/3_UMAP-Tuning/data/nk_noNKT_tuning_top_candidates.csv`
  - 热力图/分面：`2_DataProcessing/3_UMAP-Tuning/plots/heatmap_*_noNKT.(png|pdf)`
                 `2_DataProcessing/3_UMAP-Tuning/plots/UMAP_noNKT_tuning_dims*_res*_byTimepoint.(png|pdf)`
  - 最终对象：`2_DataProcessing/2_Doublet_Removed/RDS/nk.integrated.noNKT.tuned.rds`
  - 最终图件：`2_DataProcessing/3_UMAP-Tuning/plots/UMAP_noNKT_final_byTimepoint.png`
               `2_DataProcessing/3_UMAP-Tuning/plots/UMAP_noNKT_final_bySingleR.png`
  - 参数与日志：`2_DataProcessing/3_UMAP-Tuning/logs/selected_params.txt`
                `2_DataProcessing/3_UMAP-Tuning/logs/run_config_*.txt`
                `2_DataProcessing/3_UMAP-Tuning/logs/sessionInfo_*.txt`

3) 参数建议与可选调整
- 若希望减少碎簇并增强连通：
  - 将 UMAP 参数提高平滑度（例如 n.neighbors=50、min.dist=0.5）
  - 降低聚类分辨率（例如 res≤0.4）
- 若需更严格的 NKT 排除：
  - 在 NKT 剔除时引入 NK_UCell − T_UCell ≥ 0.10 的差值约束
  - 扩充 SingleR 标签排除同义词（如 “NK T cell”, “natural killer T” 等）

## 主要脚本（当前有效）
- `2_DataProcessing/Scripts/remove_doublets_and_contaminants.R`：去双胞 + 注释 + UCell + 去污染 + 重分析主流程
- `2_DataProcessing/Scripts/singleR_annotation_fix.R`：SingleR 空 data 层修复（counts→logNormCounts→SingleR）
- `2_DataProcessing/Scripts/remove_NKT_cells.R`：在已注释对象上剔除 NKT 并生成报告与图件
- `2_DataProcessing/Scripts/tune_noNKT_dims_resolution.R`：基于 noNKT 对象进行 dims × resolution 调参、选择并重跑生成最终产物
- 历史脚本（仍可参考）：`Files/UMAP/scripts/*`

## 主要结果（样例）
- 去双胞/清理整体：
  - `2_DataProcessing/2_Doublet_Removed/reports/cleaning_report.md`
  - `2_DataProcessing/2_Doublet_Removed/plots/UMAP_filtered_by_SingleR.png`
  - `2_DataProcessing/2_Doublet_Removed/plots/UMAP_filtered_clusters_by_Timepoint.png`
- NKT 剔除：
  - `2_DataProcessing/2_Doublet_Removed/reports/noNKT_removal_report.md`
  - `2_DataProcessing/2_Doublet_Removed/plots/NKT_removal_label_counts_before_after.png`
- 调参与最终：
  - `2_DataProcessing/3_UMAP-Tuning/data/nk_noNKT_tuning_metrics.csv`
  - `2_DataProcessing/3_UMAP-Tuning/plots/UMAP_noNKT_final_byTimepoint.png`
  - `2_DataProcessing/3_UMAP-Tuning/plots/UMAP_noNKT_final_bySingleR.png`

## 技术栈
- R、Seurat、SingleCellExperiment、scDblFinder、SingleR、celldex、scater、UCell、ggplot2、patchwork

## 注意事项
- `.gitignore` 默认排除大体量数据与图件/RDS 等产物；脚本与文档纳入版本控制
- Seurat 版本差异可能影响 `FindClusters` 图名称；已在脚本内自动选择可用 `graph.name`

## 更新

- 2025-10-20
  - 新增脚本：`Files/UMAP/scripts/remove_doublets_and_contaminants.R`（集成 scDblFinder 去双胞、SingleR 自动注释、UCell 签名评分、去污染规则与重分析的主流程）
  - 修复 Seurat v5 多层 assay 转换为 SCE 的问题，增强元数据行名对齐与日志/报告目录创建的健壮性

- 2025-10-20 深夜
  - 新增特异性 SingleR 修复脚本：`Files/UMAP/scripts/singleR_annotation_fix.R`（从 counts 构建 SCE，scater::logNormCounts 生成 logcounts，显式以 logcounts 作为 SingleR 输入，规避 data 层为空告警）
  - 产出对象与文档：
    - `Files/Doublet_Removed/RDS/nk.integrated.singleR_annotated.rds`
    - `Files/Doublet_Removed/reports/singleR_fix_report.md`
    - `Files/Doublet_Removed/plots/SingleR_label_barplot.png`

- 2025-10-21（清理与重分析）
  - 将 SingleR 修复策略集成至主流程并完成全流程清理与重分析（保留 NK/ILC）
  - 去双胞结果：移除 312 个细胞（约 1.61%）
  - 生成清理后对象与报告（历史路径 Files/*）：
    - `Files/Doublet_Removed/RDS/nk.integrated.filtered.rds`
    - `Files/Doublet_Removed/RDS/nk.integrated.doublet_scored.rds`
    - `Files/Doublet_Removed/reports/cleaning_report.md`
  - 图件（历史路径 Files/*）：DoubletScore、Cluster_nonNK_fraction、UMAP 等

- 2025-10-21（NKT 剔除 + UMAP 调参与最终）
  - 新增脚本：`2_DataProcessing/Scripts/remove_NKT_cells.R`（严格规则排除 NKT）
  - 新增脚本：`2_DataProcessing/Scripts/tune_noNKT_dims_resolution.R`（调参与最终降维/聚类/UMAP）
  - 产出：
    - NKT 剔除：19126 → 19007（移除 119，0.62%）
    - 调参指标与候选：`2_DataProcessing/3_UMAP-Tuning/data/*`
    - 最终参数：dims=10、res=0.3（见 `selected_params.txt`）
    - 最终对象与图：`2_DataProcessing/2_Doublet_Removed/RDS/nk.integrated.noNKT.tuned.rds`；`3_UMAP-Tuning/plots/*`

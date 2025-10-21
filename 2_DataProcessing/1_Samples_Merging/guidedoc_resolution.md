# NK UMAP 分辨率（resolution）与主成分（dims）调参与脚本修改指南

目的
- 通过独立调参脚本系统探索合适的 PCA 主成分个数（dims）与聚类分辨率（resolution），在不影响正式产出脚本的前提下形成明确的参数选择依据。
- 在确定最佳组合后，对正式脚本 `Files/UMAP/scripts/generate_umap_nk.R` 进行最小侵入式修改，仅替换相关参数并保持现有输出结构与命名。

总体思路
- 将“参数探索”和“正式产出”分离：
  1) 新建调参脚本 `Files/UMAP/scripts/tune_nk_dims_resolution.R`，读取现有整合对象 RDS，跑一组 dims×resolution 网格，统计/绘制评价指标并输出候选组合；
  2) 根据候选结果选定最终的 dims* 与 res*；
  3) 在正式脚本中仅替换 FindNeighbors/RunUMAP 的 dims 与 FindClusters 的 resolution，其余逻辑与输出文件保持不变。

一、调参脚本设计与使用
- 文件：`Files/UMAP/scripts/tune_nk_dims_resolution.R`
- 输入：整合对象 RDS（默认路径：`Files/Filter Files/RDS/nk.integrated.rds`）
- 主要流程（效率优先）：
  - 对每个 dims 值：仅运行一次 FindNeighbors 和 RunUMAP；
  - 对该 dims 下的每个 res 值：运行一次 FindClusters，并记录指标；
  - 指标包括：簇数（n_clusters）、最小簇占比（min_cluster_frac）、簇占比分位数（q95_cluster_frac）、中位轮廓系数（median_silhouette，可选）、负轮廓占比（pct_silhouette_negative，可选）；
  - 输出热力图（median silhouette、n_clusters）与自动筛选的少量候选组合的 UMAP 分面图（按 timepoint 分面）。
- 输出目录（与正式产出隔离）：
  - `Files/UMAP/Results/tuning/plots`：调参热力图与少量候选组合的分面 UMAP
  - `Files/UMAP/Results/tuning/data`：指标表 `nk_tuning_metrics.csv`、候选表 `nk_tuning_best_per_dims.csv` 与 `nk_tuning_top_candidates.csv`
  - `Files/UMAP/Results/tuning/logs`：运行配置 `run_config_*.txt` 与会话信息 `sessionInfo_*.txt`
- 推荐选择规则（可调整）：
  - 在满足 `6 ≤ n_clusters ≤ 20`、`min_cluster_frac ≥ 0.5%` 的组合中；
  - 优先选择 `median_silhouette` 较高者，如并列则选择 `pct_silhouette_negative` 更低者；
  - 结合实际生物学解释与下游可视化效果作最终决策。

二、确定参数后的正式脚本修改指南
- 原脚本：`Files/UMAP/scripts/generate_umap_nk.R`
- 修改目标（仅两处核心参数）：
  - 将 PCA/UMAP/Neighbors 所用维度统一为选定的 `dims_to_use <- D`（如 20/25/30）；
  - 将 FindClusters 的分辨率改为选定的 `resolution_to_use <- R`（如 0.6/0.8）。
- 建议在脚本顶部引入参数块：
  ```r
  # 固定调参结果
  dims_to_use <- 25           # 例如：来自调参脚本推荐
  resolution_to_use <- 0.8    # 例如：来自调参脚本推荐
  ```
- 将以下调用替换为使用上述变量：
  ```r
  nk.integrated <- FindNeighbors(nk.integrated, dims = 1:dims_to_use)
  nk.integrated <- RunUMAP(nk.integrated, dims = 1:dims_to_use)
  nk.integrated <- FindClusters(nk.integrated, resolution = resolution_to_use)
  ```
- group.by 说明：
  - 当前脚本仅运行一次 FindClusters，使用 `group.by = "seurat_clusters"` 可以保持不变；
  - 若未来扩展为多分辨率循环，应改为使用 `group.by = paste0("integrated_snn_res.", resolution_to_use)` 避免覆盖。
- timepoint 因子顺序建议固定：
  ```r
  nk.integrated$timepoint <- factor(
    nk.integrated$timepoint,
    levels = c("0W_NCD","1W_MCD","2W_MCD","6W_MCD")
  )
  ```
- 其他保持不变：
  - 输出目录与命名维持现状（因为正式产出仅对应单一最佳参数组合，不需要在文件名中加入 res/dims）；
  - RDS 的保存与加载逻辑保持现状（或将最终对象另存为快照以便复查）。

三、推荐执行流程
1) 执行调参脚本：
   - 打开 `Files/UMAP/scripts/tune_nk_dims_resolution.R`，根据机器资源调整 `dims_grid` 与 `res_grid` 的大小；
   - 运行脚本，生成指标表与热力图，并输出少量候选组合的 UMAP 分面图；
   - 结合统计指标与图件，选定最终 `dims_to_use` 与 `resolution_to_use`。
2) 修改正式脚本：
   - 在 `generate_umap_nk.R` 顶部加入参数块（dims_to_use/resolution_to_use）；
   - 将 RunUMAP/FindNeighbors/FindClusters 的 dims/res 替换为上述变量；
   - 保持其他绘图与统计逻辑不变。
3) 运行正式脚本：
   - 正式脚本按选定参数生成出版级别图件（PNG/PDF/SVG），与既有产出目录结构保持一致；
   - 如需记录参数决策，可在 `Files/UMAP/Results/scripts` 记录一份简短的 `run_config_*.txt`。

四、注意事项与常见问题
- 计算复用：
  - 调参脚本中每个 dims 仅运行一次 Neighbors/UMAP，多个 res 仅运行 FindClusters，避免重复计算；
  - 正式脚本仅运行一次，参数固定。
- 轮廓系数计算：
  - 需要 R 包 `cluster`，未安装则自动跳过 silhouette 指标（仍保留基础统计）。
- 颜色与标签：
  - 分辨率升高可能导致簇数增多，标签可能拥挤，可适当调整 `label.size` 或关闭自动标签；
  - 统一调色板可选（非必须）。
- 命名与小数点：
  - 调参脚本的文件名中分辨率默认以 `res0p6` 风格（将小数点替换为 `p`）。
  - 正式脚本无需在文件名中标注 res/dims（因为最终仅对应一个固定组合）。
- 版本兼容：
  - 若 Seurat 版本变更导致 meta 列命名差异（如 `integrated_snn_res.R`），绘图时需留意 `group.by`。
- 复现性：
  - 固定 `set.seed(1234)` 与 `future::plan("sequential")`（如需并行，再评估内存与对象导出限制）。

五、示例修改片段（从默认 20/0.5 改为 25/0.8）
```r
# 新增：在脚本顶部定义参数
dims_to_use <- 25
resolution_to_use <- 0.8

# 替换原有固定参数
nk.integrated <- FindNeighbors(nk.integrated, dims = 1:dims_to_use)
nk.integrated <- FindClusters(nk.integrated, resolution = resolution_to_use)
nk.integrated <- RunUMAP(nk.integrated, dims = 1:dims_to_use)

# 出图部分无需修改（单一分辨率），仍可使用 group.by = "seurat_clusters"
p_split <- DimPlot(
  nk.integrated,
  reduction = "umap",
  group.by = "seurat_clusters",
  split.by = "timepoint",
  label = TRUE,
  repel = TRUE,
  ncol = 2
) + theme(strip.text.x = element_text(size = 12))
```

六、结论
- 通过独立调参脚本先探索，再在正式脚本中仅替换核心参数，能同时满足性能、可追溯与产出一致性要求；
- 正式脚本保持原有目录与命名，避免与调参产物混淆；
- 若未来需要比较多个分辨率，建议新增“多分辨率循环版”的探索脚本，不直接修改正式产出脚本结构。

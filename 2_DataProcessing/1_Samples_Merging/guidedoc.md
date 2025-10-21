# 项目名称: 为NASH各时间点NK细胞生成独立的UMAP图

# 编程环境: R with Seurat and plotting packages

# 输入文件:
/home/harry/NASH/scRNA-seq/Files/Filter Files
NCD_NK1.1
MCD-1W_NK1.1
MCD-2W_NK1.1
MCD-6W_NK1.1

# 输出文件
结果图片保存在/home/harry/NASH/scRNA-seq/Files/UMAP/Results/plots下
代码文件保存在 scRNA-seq/Files/UMAP/Results/scripts 下

指令1：环境设置、数据加载与整合
"请编写R代码，使用Seurat包完成以下任务：

加载 Seurat, dplyr, ggplot2, patchwork 包。
读取四个已质控的NK1.1样本的10x Genomics数据，创建Seurat对象。在创建对象时，直接添加 timepoint 元数据，值分别为: '0W_NCD', '1W_MCD', '2W_MCD', '6W_MCD'。将这四个对象存入一个列表。
跳过质控步骤。
使用SCTransform和锚点整合工作流，将列表中的四个对象整合为一个名为 nk.integrated 的Seurat对象。
将 nk.integrated 对象的默认分析流程（DefaultAssay）设置为 'integrated'。"
指令2：降维、聚类与UMAP计算
"接上一步，我们对整合后的 nk.integrated 对象进行细胞聚类和UMAP降维：

运行主成分分析 (PCA)：RunPCA(nk.integrated, verbose = FALSE)。
使用ElbowPlot辅助判断用于下游分析的PC维度。
构建细胞近邻图：FindNeighbors(nk.integrated, dims = 1:20) (请使用前20个PCs，或根据ElbowPlot的结果调整)。
识别细胞簇：FindClusters(nk.integrated, resolution = 0.5) (resolution=0.5是常用起始值)。
计算统一的UMAP坐标：RunUMAP(nk.integrated, dims = 1:20)。"
指令3：生成并保存出版级UMAP图
"现在，我们将生成最终的UMAP图，并以三种格式保存。请编写代码完成以下操作：

首先，生成一张总览图 (Overview Plot)，展示所有细胞并按seurat_clusters着色，以便查看整体聚类效果。
<R>
p_overview <- DimPlot(nk.integrated, reduction = "umap", group.by = "seurat_clusters", label = TRUE, label.size = 5, repel = TRUE) + 
              ggtitle("All NK Cells Clustered")
print(p_overview)
接下来，生成一张按时间点分割的关键图 (Split Plot)。这张图是核心成果，图中细胞按其所属的cluster着色。
<R>
p_split <- DimPlot(nk.integrated, reduction = "umap", group.by = "seurat_clusters", split.by = "timepoint", label = TRUE, repel = TRUE, ncol = 2) + 
           theme(strip.text.x = element_text(size = 12)) # 调整分面标题大小
print(p_split)
以三种格式保存这张关键的分割图 (p_split)，确保高分辨率和矢量格式，以满足出版要求。
请使用以下代码模板进行保存:

<R>
# 定义文件名和尺寸 (inches)
output_filename <- "UMAP_NK_Clusters_by_Timepoint"
plot_width <- 12
plot_height <- 10
# 1. 保存为 PNG 格式 (高分辨率栅格图)
# dpi: dots per inch，300是印刷标准，500更佳
ggsave(
  filename = paste0(output_filename, ".png"),
  plot = p_split,
  width = plot_width,
  height = plot_height,
  dpi = 500
)
# 2. 保存为 PDF 格式 (矢量图，无限缩放，最适合出版)
# PDF格式会保留所有文字和图形元素为矢量对象
ggsave(
  filename = paste0(output_filename, ".pdf"),
  plot = p_split,
  width = plot_width,
  height = plot_height,
  device = 'pdf'
)
# 3. 保存为 SVG 格式 (矢量图，适合网页和进一步编辑)
# SVG (Scalable Vector Graphics) 格式同样是矢量，可用 Adobe Illustrator 等软件编辑
ggsave(
  filename = paste0(output_filename, ".svg"),
  plot = p_split,
  width = plot_width,
  height = plot_height,
  device = 'svg'
)
print("All plots have been saved in PNG, PDF, and SVG formats.")
"


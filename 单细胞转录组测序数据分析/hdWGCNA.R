library(Seurat)
library(tidyverse)
library(cowplot)
library(patchwork)
library(WGCNA)
library(hdWGCNA)
library(igraph)
load("Data/pbmc.Rda")

# 初始化参数
{
  theme_set(theme_cowplot())
  set.seed(42)
  enableWGCNAThreads(nThreads = 7)
}

DimPlot(pbmc, group.by='celltype', label=TRUE) 

{
  # 创建WGCNA对象
  pbmc <- SetupForWGCNA(
    pbmc,
    gene_select = "fraction", 
    fraction = 0.05, # 选择至少在5%的细胞中表达的基因
    wgcna_name = "tutorial" 
  )
}

# 构造metacells
pbmc <- MetacellsByGroups(
  pbmc = pbmc,
  group.by = c("celltype", "seurat_clusters"), 
  reduction = 'harmony', 
  k = 25,  # KNN参数
  max_shared = 10, # 两个metacell可共享的最大细胞数
  ident.group = 'celltype' 
)

{
  # 标准化metacell表达矩阵
  pbmc <- NormalizeMetacells(pbmc)
}

# 创建表达矩阵用于WGCNA
pbmc <- SetDatExpr(
  pbmc,
  group_name = "T_cells", # 挑选感兴趣的细胞类型
  group.by='celltype', 
  assay = 'SCT',
  slot = 'data' 
)

#选取软阈值
pbmc <- TestSoftPowers(
  pbmc,
  networkType = 'signed' # 此外，还可选择“unsigned”或“signed hybrid”参数
)

# 1. unsigned：无符号网络，即不考虑基因表达的正负号，只考虑它们之间的关联性。
# 2. signed：带符号网络，即同时考虑基因表达的正负号和大小，可以反映基因之间的正负调控
# 关系。
# 3. hybrid：混合网络，即将无符号网络和带符号网络结合起来，利用它们之间的优势来提高网
# 络分析的准确性和可靠性。


plot_list <- PlotSoftPowers(pbmc)

P1=wrap_plots(plot_list, ncol=2)
P1

{
  power_table <- GetPowerTable(pbmc)
  head(power_table)
}
# 构造共表达网络

{
  pbmc <- ConstructNetwork(
    pbmc, soft_power=7, # 上图中黑色圆圈中的最佳软阈值
    setDatExpr=FALSE,
    tom_name = 'Vascular Myocyte' 
  )
}
P2=PlotDendrogram(pbmc, main='Vascular Myocyte hdWGCNA Dendrogram')
P2
# 灰色模块表示无法组成任何共表达模块的基因，这些基因不纳入后续的分析

# 计算模块特征值ME: Module Eigengenes，模块特征基因（值）
# 计算所有的MEs，比较耗时
pbmc <- ModuleEigengenes(
  pbmc,
  group.by.vars="seurat_clusters"# 根据其去批次
) #时间较长


# 计算模块连接性
pbmc <- ModuleConnectivity(
  pbmc,
  group.by = 'celltype', 
  group_name = 'Vascular Myocyte' # 感兴趣的细胞类型的kME
)

# 重命名module，
pbmc <- ResetModuleNames(
  pbmc,
  new_name = "Vascular Myocyte"
)

# 可视化每个模块中的基因的kME
P3=PlotKMEs(pbmc, ncol=5, n_hubs = 20)
P3

{
  # 获取模块分配表
  modules <- GetModules(pbmc)%>%subset(module!="grey")
}
# 获取那些与每个模块高度连接的基因（hub genes）
hub_df <- GetHubGenes(pbmc, n_hubs = 100)

write.csv(hub_df,"Data/hub_df_gene.csv")

# 保存一下hdWGCNA结果
save(pbmc, file = 'hdWGCNA.rda')


# 两种评分方式
# 一种是seurat自带的AddModuleScore评分，另一种是UCell评分。
# 利用hub genes进行打分
pbmc <- ModuleExprScore(
  pbmc,
  n_genes = 25,
  method='Seurat'#AddModuleScore
)


# 可视化
plot_list <- ModuleFeaturePlot(
  pbmc,
  features='hMEs', # 可选择MEs、hMEs、scores、average
  order=TRUE,
)

p4=wrap_plots(plot_list, ncol=3)
p4

# 模块相关性
{
  ModuleCorrelogram(pbmc)
}
# 下面两种图很重要，直接把模块与亚群，也就是表型相联系，
{
  # get hMEs from seurat object
  MEs <- GetMEs(pbmc, harmonized=TRUE)
  mods <- colnames(MEs); mods <- mods[mods != 'grey']
  
  # add hMEs to Seurat meta-data:
  pbmc@meta.data <- cbind(pbmc@meta.data, MEs)
}
# plot with Seurat's DotPlot function
p <- DotPlot(pbmc, features=mods, group.by = 'celltype')

# flip the x/y axes, rotate the axis labels, and change color scheme:
p6 <- p +
  coord_flip() +
  RotatedAxis() +
  scale_color_gradient2(high='red', mid='grey95', low='blue')

# plot output
p6

p <- VlnPlot(
  pbmc,
  features = 'Vascular Myocyte3', # 换成自己的模块
  group.by = 'celltype',
  pt.size = 0 # don't show actual data points
)

# add box-and-whisker plots on top:
p <- p + geom_boxplot(width=.25, fill='white')

# change axis labels and remove legend:
p7 <- p + xlab('') + ylab('hME') + NoLegend()
p7

# 核心基因的交互网络
ModuleNetworkPlot(pbmc) # 输出结果在ModuleNetworks文件夹

# 可视化多个网络模块
HubGeneNetworkPlot(
  pbmc,
  vertex.label.cex = 0.8,   # 控制模块标签字体大小
  hub.vertex.size = 8, # hub基因节点的大小
  other.vertex.size = 3, # 其他基因节点的大小
  n_hubs = 2, # 用于可视化的
  n_other = 5, # 随机选取的gene 
  edge_prop = 0.55, # 采样的边数
  mods = "all"
)
dev.off()





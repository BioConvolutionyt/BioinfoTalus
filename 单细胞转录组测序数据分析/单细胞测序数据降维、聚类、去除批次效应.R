library(Seurat)
library(data.table)
library(stringr)
library(tibble)
library(ggplot2)
library(patchwork)
library(clustree)
library(harmony)
load("Data/pbmc.Rda")

# 若没跑SCT则先运行SCTransform
pbmc <- SCTransform(pbmc)

DefaultAssay(pbmc) = "SCT"
Idents(pbmc) = "orig.ident"


# 进行降维、聚类
#PCA
pbmc <- RunPCA(pbmc, verbose = F)


#选取合适的PC
#主成分累积贡献大于90%,选择拐点
ElbowPlot(pbmc, ndims = 50)

#确定与每个 PC 的百分比   
pct <- pbmc [["pca"]]@stdev / sum( pbmc [["pca"]]@stdev) * 100

#计算每个 PC 的累计百分比
cumu <- cumsum(pct)
cumu

#设置PC
pcs = 1:41

# 去除批次效应
# 去除批次效应的算法很多，最常用的是Harmony
# 每种算法去批次的力度不同，Harmony较弱，CCA较强，批次效应明显时可尝试用CCA，但可能导致矫枉过正

# 去除批次效应前可先看看批次效应是否明显
DimPlot(pbmc, reduction = "umap")
dev.off()

#分割
DefaultAssay(pbmc) = "RNA"
pbmc[["RNA"]] <- split(pbmc[["RNA"]], f = pbmc$orig.ident)

# Harmony
pbmc <- IntegrateLayers(
  object = pbmc, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE)

# 其他去除批次效应的算法
# # CCA
# pbmc <- IntegrateLayers(
#   object = pbmc, method = CCAIntegration,
#   orig.reduction = "pca", new.reduction = "cca",
#   verbose = TRUE)
# 
# # RPCA
# pbmc <- IntegrateLayers(
#   object = pbmc, method = RPCAIntegration,
#   orig.reduction = "pca", new.reduction = "rpca",
#   verbose = FALSE)
# 
# # JointPCA
# pbmc <- IntegrateLayers(
#   object = pbmc, method = JointPCAIntegration,
#   orig.reduction = "pca", new.reduction = "JointPCA",
#   verbose = FALSE)

#将Layers融合
pbmc = JoinLayers(pbmc)

table(pbmc@meta.data$orig.ident)

#选取合适的分辨率
#从0.1-2的resolution结果均运行一遍
DefaultAssay(pbmc) = "SCT"
seq = seq(0.1,2,by=0.1)
pbmc <- FindNeighbors(pbmc, reduction = "harmony", dims = pcs) # reduction参数设置取决于先前去除批次效应的算法
for (res in seq){
  pbmc = FindClusters(pbmc, resolution = res, algorithm = 4) # algorithm = 4表示使用Leiden算法，效果比默认算法好
}

#画图
p1 = clustree(pbmc, prefix = "SCT_snn_res.")+coord_flip()
p = p1+plot_layout(widths = c(3,1))
p


# 基于聚类结果选择合适的分辨率
# 分辨率的选择不宜太大或太小，一般选择随着分辨率的增加，聚类情况趋于稳定时的分辨率，后续可能需要根据细胞注释的情况再做调整
# 降维、聚类的reduction参数设置取决于先前去除批次效应的算法
# 去除批次效应的算法可能会导致重复的点，使得tSNE过程报错
pbmc <- FindNeighbors(pbmc, reduction = "harmony",  dims = pcs) %>% FindClusters(resolution = 0.6, algorithm = 4) 
pbmc <- RunUMAP(pbmc, reduction = "harmony", dims = pcs, reduction.name = "umap.harmony")
pbmc <- RunTSNE(pbmc, reduction = "harmony", dims = pcs, reduction.name = "tsne.harmony")


save(pbmc, file="./Data/pbmc.Rda")

#画图
DimPlot(pbmc, reduction = "umap.harmony", label = T)
dev.off()

DimPlot(pbmc,reduction = "umap.harmony",label = F,group.by = "orig.ident")
dev.off()



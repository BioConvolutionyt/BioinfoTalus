library(Seurat)
library(clustree)
library(ggplot2)
library(patchwork)
library(tidyverse)


####空转数据读取####
# 空间转录组数据的读取方式很多，这里仅展示最简单的一种
# 确保目录结构中包含：filtered_feature_bc_matrix.h5
# spatial/ 子目录（含 tissue_positions.csv、scalefactors_json.json、tissue_lowres_image.png）
# 需要SeuratObject 5.3.0
stRNA <- Load10X_Spatial(data.dir = "Data/10x Visium/GSM8594561",
                         filename = "filtered_feature_bc_matrix.h5", # 过滤后的文件
                         assay = "Spatial", #存储数据的 assay 名称
                         slice = "GSM8594561" #空间切片命名
)


SpatialFeaturePlot(stRNA, 
                   features = "nCount_Spatial", 
                   pt.size = 5, # 散点大小
                   alpha = 0.6 # 透明度
                   )

# 空转数据一般可不进行质量控制
stRNA[["percent.mt"]] <- PercentageFeatureSet(stRNA, pattern = "^MT-")
stRNA <- SCTransform(stRNA, assay = "Spatial")

save(stRNA, file="Data/stRNA.Rda")


####降维，聚类####
DefaultAssay(stRNA) = "SCT"

#PCA
stRNA <- RunPCA(stRNA, verbose = F)
#选取合适的PC
#主成分累积贡献大于90%,选择拐点
ElbowPlot(stRNA, ndims = 50)

#确定与每个 PC 的百分比   
pct <- stRNA [["pca"]]@stdev / sum( stRNA [["pca"]]@stdev) * 100

#计算每个 PC 的累计百分比
cumu <- cumsum(pct)
cumu

#设置PC
pcs = 1:38

#选取合适的分辨率
#从0.1-2的resolution结果均运行一遍
seq = seq(0.1,2,by=0.1)
stRNA <- FindNeighbors(stRNA,  dims = pcs) 
for (res in seq){
  stRNA = FindClusters(stRNA, resolution = res, algorithm = 4)
}

#画图
p1 = clustree(stRNA, prefix = "SCT_snn_res.")+coord_flip()
p = p1+plot_layout(widths = c(3,1))
p

#降维聚类，reduction参数设置取决于先前去除批次效应的算法
stRNA <- FindNeighbors(stRNA, reduction = "pca",  dims = pcs) %>% FindClusters(resolution = 0.8, algorithm = 4) #需要根据聚类图选择合适的分别率
stRNA <- RunUMAP(stRNA, reduction = "pca", dims = pcs, reduction.name = "umap")
stRNA <- RunTSNE(stRNA, reduction = "pca", dims = pcs, reduction.name = "tsne")

#画图
DimPlot(stRNA, reduction = "umap", label = T)
dev.off()

save(stRNA, file="Data/stRNA.Rda")


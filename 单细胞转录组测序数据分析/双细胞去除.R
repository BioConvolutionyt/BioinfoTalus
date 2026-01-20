library(Seurat)
library(data.table)
library(stringr)
library(tibble)
library(ggplot2)
library(patchwork)
library(DoubletFinder)
library(clustree)
load("./Data/pbmc.Rda")

#数据标准化、寻找可变特征、数据缩放
pbmc = pbmc %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData()

# 后续进行的降维和聚类仅因为双细胞去除的算法需要相应的结果
# 去除双细胞后还会重新进行正式的降维、聚类

#PCA
pbmc <- RunPCA(pbmc, verbose = F)

#主成分分析图形
#DimPlot(object = pbmc, reduction = "pca")
#dev.off()

#选取合适的PC
#主成分累积贡献大于90%,选择拐点
ElbowPlot(pbmc, ndims = 50)

#确定与每个 PC 的百分比   
pct <- pbmc [["pca"]]@stdev / sum( pbmc [["pca"]]@stdev) * 100

#计算每个 PC 的累计百分比
cumu <- cumsum(pct)
cumu 

#设置PC（主成分数量），一般来说只要大于30，不同主成分数量对后续结果影响不大
pcs = 1:41  # 根据cumu选择累计贡献率大于90%时的主成分数量

#选取合适的分辨率
#从0.1-2的resolution结果均运行一遍
seq = seq(0.1,2,by=0.1)
pbmc <- FindNeighbors(pbmc,  dims = pcs) 
for (res in seq){
  pbmc = FindClusters(pbmc, resolution = res, algorithm = 4) # algorithm = 4表示使用Leiden算法，效果比默认算法好
}

#画图
p1 = clustree(pbmc, prefix = "RNA_snn_res.")+coord_flip()
p = p1+plot_layout(widths = c(3,1))
p

gc() # 及时清理内存

# 基于聚类结果选择合适的分辨率
# 分辨率的选择不宜太大或太小，一般选择随着分辨率的增加，聚类情况趋于稳定时的分辨率
# PCA、UMAP、寻找邻居、聚类
pbmc = pbmc %>% 
  RunPCA() %>% 
  RunUMAP(dims = pcs) %>%  
  RunTSNE(dims = pcs) %>% 
  FindNeighbors(dims = pcs) %>% 
  FindClusters(resolution = 0.6, algorithm = 4) # 这里选择分辨率为0.6

#首先获得最佳的pK值
#pK表示领域大小
sweep.res.list <- paramSweep(pbmc, PCs = pcs, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
pk_best = bcmvn %>% 
  dplyr::arrange(desc(BCmetric)) %>% 
  dplyr::pull(pK) %>% 
  .[1] %>% as.character() %>% as.numeric()

#然后估算出双细胞群中，homotypic doublets的比例
annotations <- pbmc$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
print(homotypic.prop)

#双细胞占比为5%左右
nExp_poi <- round(0.05*nrow(pbmc@meta.data))        
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop)) 

# 模拟出的artificial doublet数量。不同取值对识别结果影响不大，默认为0.25
pbmc <- doubletFinder(pbmc, PCs = pcs, 
                      pN = 0.25, pK = pk_best, nExp = nExp_poi.adj, 
                      reuse.pANN = NULL, sct = FALSE)
# 若出现：错误与xtfrm.data.frame(x)：无法xtfrm数据帧
# 将上述的reuse.pANN＝FALSE改成NULL

# 将列名改为"Double_score"和"Is_Double"
colnames(pbmc@meta.data)[length(colnames(pbmc@meta.data))-1] <- "Double_score"
colnames(pbmc@meta.data)[length(colnames(pbmc@meta.data))] <- "Is_Double"

# 绘制DoubletFinder分类的tsne图
# DimPlot(pbmc, reduction = "tsne", group.by = "Is_Double")
# 绘制双细胞分类的小提琴图
# VlnPlot(pbmc, group.by = "Is_Double", features = c("nCount_RNA", "nFeature_RNA"),pt.size = 0, ncol = 2)

# 过滤非单细胞数据
pbmc <- subset(pbmc, Is_Double == "Singlet")
# 去除无用的中间结果
pbmc@meta.data <- pbmc@meta.data[, !grepl("^RNA_", names(pbmc@meta.data))]
pbmc@meta.data[, "Is_Double"] <- NULL

save(pbmc, file="Data/pbmc.Rda")



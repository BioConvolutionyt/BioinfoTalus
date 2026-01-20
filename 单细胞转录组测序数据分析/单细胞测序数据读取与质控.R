library(Seurat)
library(data.table)
library(stringr)
library(tibble)


####数据读取####
# 10X数据的读取
# 每个样本对应的文件夹下需要包含barcodes，genes (features)，matrix三个文件，需要预先整理好
# 文件的命名必须为barcodes，genes/features，matrix
# 文件类型可以是.gz也可以是解压过的
samples <- list.files("Data/10x Single Cell")

# 创建一个空的列表
seurat_list <- list()

#读取数据并创建Seurat对象
for (sample in samples) {
  #文件路径
  data.path <- paste0("Data/10x Single Cell/", sample)
  #读取10x数据
  seurat_data <- Read10X(data.dir = data.path)
  #创建Seurat对象，删除，小于200个基因表达的细胞，小于3个细胞表达的基因
  seurat_obj <- CreateSeuratObject(counts=seurat_data, project=sample, min.features=200, min.cells=3)
  #添加到列表中
  seurat_list <- append(seurat_list, seurat_obj)
}

#合并
seurat_combined <- merge(seurat_list[[1]], 
                         y = seurat_list[-1],
                         add.cell.ids = samples)

#将Layers融合
pbmc = JoinLayers(seurat_combined)


####质量控制####
# 计算线粒体基因的百分比
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
# 计算核糖体基因的百分比
pbmc[["percent.rb"]] <- PercentageFeatureSet(pbmc, pattern = "^RP")
# 计算红细胞基因的百分比
# HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
# HB.genes <- CaseMatch(HB.genes, rownames(pbmc))
# pbmc[["percent.hb"]]<-PercentageFeatureSet(pbmc, features=HB.genes) 
# 绘制小提琴图，展示每个细胞中的基因数（nFeature_RNA）、检测到的分子数（nCount_RNA）以及线粒体基因的百分比（percent.mt）
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA","percent.mt", "percent.rb"), 
        ncol = 2) #,pt.size = 0


# nCount表示每个细胞中检测到的读段数
# nFeature表示每个细胞中检测到的基因数
# 上述两者，若过小则可能为异常的测序结果，若过大则可能为双细胞或多细胞
# 绘制散点图，比较不同特征之间的关系，离群散点可能为异常细胞
FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.rb")+ RotatedAxis()

FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+ RotatedAxis()


# 查看nCount与nFeature的分布
quantile(pbmc$nFeature_RNA, seq(0.01, 0.1, 0.01))
quantile(pbmc$nFeature_RNA, seq(0.9, 1, 0.01))
quantile(pbmc$nCount_RNA, seq(0.01, 0.1, 0.01))
quantile(pbmc$nCount_RNA, seq(0.9, 1, 0.01))
quantile(pbmc$percent.mt, seq(0.9, 1, 0.01))
quantile(pbmc$percent.rb, seq(0.9, 1, 0.01))


# 过滤，过滤条件因数据而异（建议查阅相关文献来确定特定数据的质控方式，以下是相对通用的模式）
sub1 <- pbmc$nCount_RNA >= 1000 # 过滤掉低于1000的细胞
sub2 <- pbmc$nFeature_RNA >= 200 & pbmc$nFeature_RNA <= 8000 # 过滤掉基因数小于200或大于8000的细胞
sub3 <- pbmc$percent.mt <= 20 # 过滤掉线粒体基因百分比大于20%的细胞
# sub4 <- pbmc$percent.rb <= 20 # 过滤掉核糖体基因百分比大于20%的细胞
# sub5 <- pbmc$percent.hb <= 1
sub <- sub1 & sub2 & sub3 # & sub4 & sub5 # 合并过滤条件
pbmc <- pbmc[, sub] # 应用过滤条件，保留符合条件的细胞

# 重新绘制
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA","percent.mt", "percent.rb"), 
        ncol = 2)

save(pbmc, file="Data/pbmc.Rda")


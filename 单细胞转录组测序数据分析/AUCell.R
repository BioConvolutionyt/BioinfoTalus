library(Seurat)
library(GSVA)
library(GSEABase)
library(AUCell)
library(doParallel)
library(future)
load("./Data/pbmc.Rda")


DefaultAssay(pbmc) <- "RNA"
#标准化
pbmc <- NormalizeData(pbmc)
dim(pbmc)

#获取表达矩阵
expr = as.matrix(GetAssayData(object = pbmc@assays$RNA, layer = "data"))
# 获取平均表达量
# expr = as.matrix(AverageExpression(pbmc, assays = "RNA", layer = "data", group.by = "celltype")[[1]])
dim(expr)

#筛选
expr <- expr[rowSums(expr)>0,]
dim(expr)

# 以Hallmark基于集为例
# 基因集需要为gmt文件
geneSets = getGmt("Data/msigdb_v7.0_GMTs/h.all.v7.0.symbols.gmt", geneIdType=SymbolIdentifier())

# AUCell，基于排序计算评分，可比性强
cells_rankings <- AUCell_buildRankings(expr) # 会绘制细胞的基因表达情况，横轴是表达的基因的数量，纵轴是细胞数量
plan(multisession, workers = 10)  # 使用 doParallel 兼容的后端
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)
Result <- getAUC(cells_AUC)
dev.off()

# 得分进行标准化
normalize=function(x){return((x-min(x))/(max(x)-min(x)))}
Result=normalize(Result)

#保存结果
Output=rbind(id=colnames(Result), Result)
write.table(Output, file="Data/ssgseaOut_AUCell.txt", sep="\t", quote=F, col.names=F)

#添加到meta.data中
pbmc <- AddMetaData(pbmc, metadata = t(Result))


#画图
#绘制单个个基于集得分
FeaturePlot(pbmc, features=rownames(Result)[3], reduction="umap", ncol=1) #features=rownames(Result)
DimPlot(pbmc, label=T, reduction="umap")
dev.off()


#保存结果
save(pbmc, "Data/pbmc.Rda")
write.csv(hmExp, "Data/geneset_score.csv")



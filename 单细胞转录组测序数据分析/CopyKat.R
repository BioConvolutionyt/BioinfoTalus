library(Seurat)
library(copykat)
load("Data/pbmc.Rda")


#只能单样本运行
table(pbmc@meta.data$orig.ident)
pbmc <-  subset(pbmc, orig.ident %in% "GSM9113377")
table(pbmc@meta.data$celltype)

#提取肿瘤来源细胞和部分正常参考细胞
pbmc_turmor <- subset(pbmc, celltype %in% 'Epithelial_cells') # 提取肿瘤候选细胞
pbmc_normal <- subset(pbmc, celltype %in% c("T_cells ","Plasma_cells", "B_cells")) # 取正常参考细胞，一般选免疫细胞

#合并
pbmc <- merge(pbmc_turmor, pbmc_normal)

#融合
DefaultAssay(pbmc) = "RNA"
pbmc = JoinLayers(pbmc)

#提取用于分析的表达矩阵
counts <- as.matrix(GetAssayData(object = pbmc@assays$RNA, layer = "counts"))

#设置正常参考细胞
ref <- colnames(pbmc_normal)

#运行
res <- copykat(rawmat=counts,
               ngene.chr=5,
               norm.cell.names=ref,
               sam.name="all",
               n.cores=1)

#无参考细胞时仅需运行这句
#res <- copykat(rawmat=counts,ngene.chr=5,sam.name="all",n.cores=10)

#保存结果
save(res, "Data/copykat.res.Rda")

#读入
malignant <- read.delim("all_copykat_prediction.txt")

#转化为数据框
malignant <- data.frame(copykat.pred = malignant$copykat.pred, row.names = malignant$cell.names)

#将结果添加到meta.data
pbmc <- AddMetaData(pbmc, metadata = malignant)
table(pbmc$copykat.pred)

#修改
pbmc$copykat <- recode(pbmc@meta.data$copykat.pred,
                       "not.defined" = "diploid")
table(pbmc$copykat)
#画图
DimPlot(pbmc, group.by = "copykat", cols = c("red", "blue", "gray50"))

save(pbmc, file="Data/pbmc.Rda")



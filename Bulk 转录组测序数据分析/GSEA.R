library(org.Hs.eg.db) #org.Hs.eg.db包主要注释基因:用于不同数据库ID间的转化
library(clusterProfiler)
library(tidyverse)
library(enrichplot)
library(RColorBrewer)
library(GseaVis)


load("Data/DEG_results.Rda")
DEG <- arrange(DEG, padj)

DEG <- rownames_to_column(DEG, "Gene")

geneList = DEG[,"log2FoldChange"]
names(geneList) = as.character(DEG[,'Gene'])
head(geneList)
geneList = sort(geneList, decreasing = TRUE)
head(geneList)

#GSEA基因集：https://zhuanlan.zhihu.com/p/504101161

msigdb_GMTs <- "Data/msigdb_v7.0_GMTs"
msigdb <- "h.all.v7.0.symbols.gmt"    
#读取上面指定的gmt文件
kegmt <- read.gmt(file.path(msigdb_GMTs,msigdb))

set.seed(42) #设置种子
gsea <-GSEA(geneList, TERM2GENE = kegmt) #GSEA分析
#转换成数据框
gsea_result_df <- as.data.frame(gsea)
save(gsea,gsea_result_df,file = "Data/GSEA_gene_h.all.rda")

#绘图
# 绘制上调的通路
gsea_up <- gsea_result_df[gsea_result_df$NES>0, ]
# pdf("h.gsea_up.pdf",11, 7)
gseaNb(gsea,
       geneSetID = as.vector(gsea_up$ID)[1:min(5, length(as.vector(gsea_up$ID)))],
       curveCol = rev(brewer.pal(n = 5, name = "Set1")),
       lineSize = 1.5,
       # htCol = c("#3287BC", "#EB6046"),
       htAlpha = 1,
       htHeight = 0.5,
       rmPrefix= F,
       base_size = 14,
       legend.position = c(0.85,0.75),
       addPval=F,
       pvalX = 0.9,
       pvalY = 0.9)
dev.off()

# 绘制下调的通路
gsea_down <- gsea_result_df[gsea_result_df$NES<0, ]
# pdf("h.gsea_down.pdf",11, 7)
gseaNb(gsea,
       geneSetID = as.vector(gsea_down$ID)[1:min(5, length(as.vector(gsea_down$ID)))],
       curveCol = rev(brewer.pal(n = 5, name = "Set1")),
       lineSize = 1.5,
       # htCol = c("#3287BC", "#EB6046"),
       htAlpha = 1,
       htHeight = 0.5,
       rmPrefix= F,
       base_size = 13,
       legend.position = c(0.85,0.7),
       addPval=F,
       pvalX = 0.9,
       pvalY = 0.9)
dev.off()

#换C7跑
msigdb_GMTs <- "Data/msigdb_v7.0_GMTs"
msigdb <- "c7.all.v7.0.symbols.gmt"   
#读取上面指定的gmt文件
kegmt <- read.gmt(file.path(msigdb_GMTs,msigdb))

set.seed(42) #设置种子
gsea <-GSEA(geneList,TERM2GENE = kegmt) #GSEA分析
#转换成数据框
gsea_result_df <- as.data.frame(gsea)
save(gsea,gsea_result_df,file = "Data/GSEA_gene_c7.rda")

#绘图
# 绘制上调的通路
gsea_up <- gsea_result_df[gsea_result_df$NES>0, ]
# pdf("c7.gsea_up.pdf",11, 7)
gseaNb(gsea,
       geneSetID = as.vector(gsea_up$ID)[1:min(5, length(as.vector(gsea_up$ID)))],
       curveCol = rev(brewer.pal(n = 5, name = "Set1")),
       lineSize = 1.5,
       # htCol = c("#3287BC", "#EB6046"),
       htAlpha = 1,
       htHeight = 0.5,
       rmPrefix= F,
       base_size = 14,
       legend.position = c(0.8,0.75),
       addPval=F,
       pvalX = 0.9,
       pvalY = 0.9)
dev.off()
# 绘制下调的通路
gsea_down <- gsea_result_df[gsea_result_df$NES<0, ]
# pdf("c7.gsea_down.pdf",11, 7)
gseaNb(gsea,
       geneSetID = as.vector(gsea_down$ID)[1:min(5, length(as.vector(gsea_down$ID)))],
       curveCol = rev(brewer.pal(n = 5, name = "Set1")),
       lineSize = 1.5,
       # htCol = c("#3287BC", "#EB6046"),
       htAlpha = 1,
       htHeight = 0.5,
       rmPrefix= F,
       base_size = 13,
       legend.position = c(0.8,0.7),
       addPval=F,
       pvalX = 0.9,
       pvalY = 0.9)
dev.off()





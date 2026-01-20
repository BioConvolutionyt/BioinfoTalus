library(tidyverse)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ggnewscale)
library(ggplot2)
library(pathview)

# 以差异基因为例进行富集分析
load("Data/DEG_results.Rda")
# 筛选差异显著的基因
DEG <- DEG[DEG$padj < 0.05 & DEG$log2FoldChange > 2,]

# 换成别的数据时，数据中需要包含名以SYMBOL为列名记录的基因名称
# 将基因SYMBOL转换为富集分析所需的ENTREZID
DEG <- rownames_to_column(DEG,"SYMBOL")
genelist <- bitr(DEG$SYMBOL, fromType="SYMBOL",
                 toType="ENTREZID", OrgDb='org.Hs.eg.db')
DEG <- inner_join(DEG,genelist,by="SYMBOL")


# GO
ego <- enrichGO(gene = DEG$ENTREZID,
                OrgDb = org.Hs.eg.db, 
                ont = "all",
                pAdjustMethod = "BH",
                minGSSize = 1,
                pvalueCutoff =0.05, 
                qvalueCutoff =0.05,
                readable = TRUE)

ego_res <- ego@result
save(ego,ego_res,file = "Data/GO_results.Rda")


# KEGG
# 需要联网
options(timeout = 300)
kk <- enrichKEGG(gene = DEG$ENTREZID,
                 organism = 'hsa', # 物种
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.1)
kk_res <- kk@result
save(kk,kk_res,file = "Data/KEGG_results.Rda")

#网络图
# 生成绘图所需的foldChange List
# 没有foldChange数据就不执行如下代码，且绘图时不传入foldChange参数
List = DEG$log2FoldChange
names(List)= DEG$ENTREZID
head(List)
List = sort(List,decreasing = T)

#GO
cnetplot(ego, foldChange = List, circular = TRUE, colorEdge = TRUE)
#KEGG
cnetplot(kk, foldChange = List, circular = TRUE, colorEdge = TRUE)

# 气泡图和条形图
# GO
ego_df <- as.data.frame(ego)
dotplot(ego, split="ONTOLOGY",showCategory=5)+facet_grid(ONTOLOGY~., scale="free") #用oncology来分框，每个框里而显示5个条口
barplot(ego, split="ONTOLOGY",showCategory=5)+facet_grid(ONTOLOGY~., scale="free")
# KEGG
kk_df <- as.data.frame(kk)
dotplot(kk, showCategory=15, orderBy="GeneRatio", label_format=70)
barplot(kk, showCategory=15, drop=TRUE, label_format=70)

# KEGG通路可视化
pathway <- pathview(gene.data=DEG$SYMBOL, # 基因列表
                  pathway.id="hsa04820", # 需要可视化的通路名
                  species="hsa", # 物种
                  gene.idtype="SYMBOL", # 基因ID类型 
                  kegg.native=TRUE)

# KEGG通路+差异表达可视化
deg_kegg <- DEG[, c("SYMBOL", "log2FoldChange")]
deg_kegg <- column_to_rownames(deg_kegg, "SYMBOL")
pathway<-pathview(gene.data=deg_kegg, # 基因列表
                  pathway.id="hsa04820", # 通路名
                  species="hsa", # 物种
                  gene.idtype="SYMBOL", # 基因ID类型 
                  kegg.native=TRUE)


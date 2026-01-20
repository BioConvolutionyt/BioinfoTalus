library(tidyverse)
library(DESeq2) # 做差异分析的包，也可用edgeR，limma做差异分析
library(ggpubr)
library(ggthemes)

# 注意：差异分析必须用counts数据
load("Data/TCGA_counts_matrix.Rda")

# 定义分组条件，这里以是否是肿瘤样本为分组条件（01A后缀的样本为肿瘤，11A为正常）
conditions = data.frame(group = factor(ifelse(endsWith(colnames(counts_matrix), "01A"), 1, 0), levels=c(1, 0)))

# 对所有数值列进行四舍五入取整
for (col in names(counts_matrix)) {
  if (is.numeric(counts_matrix[[col]])) {
    counts_matrix[[col]] <- round(counts_matrix[[col]])
  }
}

# 差异分析准备工作
dds <- DESeqDataSetFromMatrix(
  countData = counts_matrix,
  colData = conditions,
  design = ~ group)

# 开始差异分析
dds <- DESeq(dds)
resultsNames(dds)
# 提取结果
DEG <- results(dds)
DEG <- as.data.frame(DEG)
DEG <- DEG[rowSums(is.na(DEG)) == 0, ] # 去除包含空值NA的行
# 筛选差异显著的基因
DEG <- DEG[DEG$padj < 0.05, ] # 这里仅针对调整后的p值做筛选，后续可按需对log2FC做筛选
save(DEG, file="Data/DEG_results.Rda")


# 火山图
#添加上下调信息
logFC_cutoff <- 1
type1 = (DEG$padj < 0.05)&(DEG$log2FoldChange < -logFC_cutoff)
type2 = (DEG$padj < 0.05)&(DEG$log2FoldChange > logFC_cutoff)
DEG$change = ifelse(type1,"DOWN",ifelse(type2,"UP","NOT"))
DEG$logP <- -log10(DEG$padj)

# 添加基因标签信息
DEG$Label = ""  # 新加一列label
# 对差异基因的p值进行从小到大的排序
DEG <- DEG[order(DEG$padj), ]
# 将行名赋给DEG$Gene
DEG$Gene <- rownames(DEG)
# 高表达的基因中，选择fdr值最小的5个
up.genes <- head(DEG$Gene[which(DEG$change == "UP")], 5)
# 低表达的基因中，选择fdr值最小的5个
down.genes <- head(DEG$Gene[which(DEG$change == "DOWN")], 5)
# 将up.genes和down.genes合并，并加入到Label中
DEG.top5.genes <- c(as.character(up.genes), as.character(down.genes))
# 在DEG$Label中为DEG.top5.genes对应的基因添加标签
DEG$Label[match(DEG.top5.genes, DEG$Gene)] <- DEG.top5.genes

ggscatter(DEG, x = "log2FoldChange", y = "logP",
          color = "change",
          palette = c("blue", "black", "red"),
          size = 1,
          label = DEG$Label,
          font.label = 8,
          repel = TRUE,
          xlab = "log2FoldChange",
          ylab = "-log10(Adjust P-value)") +
  theme_base() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-logFC_cutoff, logFC_cutoff), linetype = "dashed")




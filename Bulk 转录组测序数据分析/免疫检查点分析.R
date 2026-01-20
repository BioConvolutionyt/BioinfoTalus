library(limma)
library(reshape2)
library(ggpubr)
library(ggplot2)
library(tidyverse)

# 基因表达数据
load("Data/TCGA_tpms_matrix.Rda")
# 选出肿瘤样本
tpms_matrix <- tpms_matrix[, endsWith(colnames(tpms_matrix), "01A")]
colnames(tpms_matrix) <- substr(colnames(tpms_matrix), 1, 12)

rt <- as.data.frame(t(tpms_matrix))
# 免疫检查点相关基因
imm_genes <- read.table("Data/ICGs.txt", sep="\t", header = F)
# 分组数据
# 以免疫评分为例对样本进行高低风险组分组
ESTIMATE_result <- read.table("Data/TCGA_ESTIMATE_result.txt",sep = "\t", row.names = 1, header = TRUE)
risk <- data.frame(ifelse(ESTIMATE_result$ImmuneScore > median(ESTIMATE_result$ImmuneScore), "Low", "High"))

# 将免疫检查点相关基因与数据中包含的基因取交集
inter_gene <- intersect(imm_genes$V1, colnames(rt))
rt <- rt[, inter_gene]
rt <- rt[rownames(ESTIMATE_result), ] # 合并分组数据前先对齐索引
rt$risk <- risk[, 1]


# 保留符合P值筛选条件的基因
pvalue.sig = 0.05
sigGene = c()

for(i in colnames(rt)[1:(ncol(rt)-1)]) {
  if(sd(rt[,i]) < 0.001) { next }
  wilcoxTest = wilcox.test(rt[,i] ~ rt[,"risk"])
  pvalue = wilcoxTest$p.value
  if(wilcoxTest$p.value < pvalue.sig) {
    sigGene = c(sigGene, i)
  }
}
sigGene = c(sigGene, "risk")
rt = rt[, sigGene]

rt = melt(rt, id.var = c("risk"))
colnames(rt) = c("risk", "gene", "expression")
#最终的 rt 是一个包含三列（risk、gene、expression）的长格式数据框
#每行表示某个样本在某个经 Wilcoxon 检验显著的免疫检查点基因上的表达量及其分组（High/Low）

save(rt, file="Data/immcheckpoint.Rda")


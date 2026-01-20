library(tidyverse)
library(glmnet)

#读取文件
load("Data/DEG_results.Rda")
load("Data/TCGA_tpms_matrix.Rda")
surv <- read.csv("Data/clinical.csv", row.names = 1)

# 选出肿瘤样本
tpms_matrix <- tpms_matrix[, endsWith(colnames(tpms_matrix), "01A")]
colnames(tpms_matrix) <- substr(colnames(tpms_matrix), 1, 12)
tpms_matrix <- as.data.frame(t(tpms_matrix))
# 去除生存时间缺失的样本
surv <- surv[surv$OS.time != "'--", ]
# 筛选显著差异基因
DEG <- DEG[abs(DEG$log2FoldChange) > 2, ]

# 合并生存数据与表达矩阵
inter_sample <- intersect(rownames(surv), rownames(tpms_matrix))
surv <- surv[inter_sample, ]
tpms_matrix <- tpms_matrix[inter_sample, ]
surv.expr <- cbind(surv[,c("OS.time", "OS")], tpms_matrix[,rownames(DEG)])
# 将生存数据类型转换成数值型
surv.expr[,c("OS.time", "OS")] <- lapply(surv.expr[,c("OS.time", "OS")], as.numeric)

# 使用交叉验证的 LASSO（alpha=1 表示 LASSO）来选择特征
cv_fit <-cv.glmnet(x=as.matrix(surv.expr[,3:ncol(surv.expr)]), y=surv.expr$OS.time, alpha=1)
# 绘制交叉验证结果（随 lambda 的变化情况）
plot(cv_fit)

fit <-glmnet(x=as.matrix(surv.expr[,3:ncol(surv.expr)]), y=surv.expr$OS.time, alpha=1)
# 绘制系数随 lambda 变化的图
plot(fit,xvar ="lambda")

# 从交叉验证结果中获取达到最优性能的 lambda 对应的系数向量
coefficient <- coef(cv_fit,s="lambda.min")

# 提取关键基因
Active.index<- which(as.numeric(coefficient)!=0)
Active.coefficient<- as.numeric(coefficient)[Active.index]
key_gene<- rownames(coefficient)[Active.index][-1]
key_gene <- as.data.frame(key_gene)
colnames(key_gene) <- "SYMBOL"

save(key_gene, file = "Data/key_gene.Rda")



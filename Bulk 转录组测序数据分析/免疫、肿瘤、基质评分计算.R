library(tidyverse)
library(utils)
library(estimate)

load("Data/TCGA_tpms_matrix.Rda")

# 选出肿瘤样本
tpms_matrix <- tpms_matrix[, endsWith(colnames(tpms_matrix), "01A")]
colnames(tpms_matrix) <- substr(colnames(tpms_matrix), 1, 12)
write.table(tpms_matrix, "Data/TCGA_turmor_exp.txt", sep="\t", quote=F)


#计算评分
filterCommonGenes(input.f = "Data/TCGA_turmor_exp.txt",   
                  output.f = "Data/TCGA_turmor_exp.gct",   #输出文件名
                  id = "GeneSymbol")   #行名为gene symbol
estimateScore("Data/TCGA_turmor_exp.gct",   #刚才的输出文件名
              "Data/TCGA_estimate_score.txt",   #新的输出文件名（即估计的结果文件）
              platform="affymetrix")   #默认平台

#提取结果并整理
ESTIMATE_result <- read.table("Data/TCGA_estimate_score.txt", sep = "\t", row.names = 1,check.names = F,stringsAsFactors = F,header = T)
ESTIMATE_result <- ESTIMATE_result[,-1]  
colnames(ESTIMATE_result) <- ESTIMATE_result[1,]   
ESTIMATE_result <- as.data.frame(t(ESTIMATE_result[-1,]))
rownames(ESTIMATE_result) <- colnames(tpms_matrix)

#保存结果
write.table(ESTIMATE_result, file = "Data/TCGA_ESTIMATE_result.txt",sep = "\t",row.names = T,col.names = NA,quote = F)



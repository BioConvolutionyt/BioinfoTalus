library(GSVA)
library(BiocParallel)
library(tidyverse)

# 读入表达矩阵
load("Data/TCGA_tpms_matrix.Rda")
# 读入待打分的基因集（可以自定义，本质上为列表，列表的names属性为条目名称，对应元素为基因名向量）
load("Data/ssGSEA data packages/TISIDB肿瘤浸润淋巴细胞基因集.rdata")

exp_matrix <- as.matrix(tpms_matrix)

# 设定参数
param <- ssgseaParam(
  exprData = exp_matrix, # 表达矩阵
  geneSets = tisidb_cell # 基因集
  )

# 如果要用GSVA则运行下面的代码
# param <- gsvaParam(exp_matrix, 
#                      tisidb_cell, 
#                      kcdf='Gaussian',
#                      absRanking=TRUE)


# 执行ssGSEA
ssGSEA_matrix <- gsva(
  param,
  verbose = TRUE)

ssGSEA_matrix <- as.data.frame(ssGSEA_matrix)
save(ssGSEA_matrix, file="Data/ssGSEA_matrix.Rda")



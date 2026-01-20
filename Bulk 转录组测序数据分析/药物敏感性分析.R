library(limma)       # 用于基因表达数据的分析
library(oncoPredict) # oncoPredict包，用于预测药物敏感性
library(parallel)    # 并行计算，提高计算效率
set.seed(42)      # 设置随机种子，确保结果的可重复性


load("Data/TCGA_tpms_matrix.Rda")
# 选出肿瘤样本
tpms_matrix <- tpms_matrix[, endsWith(colnames(tpms_matrix), "01A")]
colnames(tpms_matrix) <- substr(colnames(tpms_matrix), 1, 12)

# 将数据框转换为矩阵
rt <- as.matrix(tpms_matrix) #行为基因列为样本

# 读取药物敏感性数据
GDSC2_Expr <- readRDS(file = 'Data/GDSC2_Expr.rds')  # 读取GDSC2表达数据
GDSC2_Res <- readRDS(file = 'Data/GDSC2_Res.rds')    # 读取GDSC2药物反应数据
# 对药物反应数据进行指数变换
GDSC2_Res <- exp(GDSC2_Res)

# 预测药物敏感性
calcPhenotype(
  trainingExprData = GDSC2_Expr,    # 训练组的基因表达数据
  trainingPtype = GDSC2_Res,        # 训练组的药物敏感性数据
  testExprData = rt,                # 测试组的基因表达数据
  batchCorrect = 'eb',              # 使用经验贝叶斯方法进行批次效应校正
  powerTransformPhenotype = TRUE,   # 对药物敏感性数据进行幂变换
  removeLowVaryingGenes = 0.95,      # 去除表达波动小于0.3的基因
  minNumSamples = 10,               # 样本数目少于10的基因将被去除
  printOutput = TRUE,               # 输出结果
  removeLowVaringGenesFrom = 'rawData' # 从原始数据中去除波动小的基因
)

# 加载绘图所需的包
library(ggplot2)
library(ggpubr)

# 读取预测的药物敏感性数据和风险数据
sen <- read.csv("calcPhenotype_Output/DrugPredictions.csv", row.names = 1)  # 药物敏感性预测结果

# 以免疫评分为例对样本进行高低风险组分组
ESTIMATE_result <- read.table("Data/TCGA_ESTIMATE_result.txt",sep = "\t", row.names = 1, header = TRUE)
risk <- data.frame(ifelse(ESTIMATE_result$ImmuneScore > median(ESTIMATE_result$ImmuneScore), "Low", "High"))
colnames(risk) <- "risk"

# 去除药物敏感性数据列名中的数字后缀（可以不运行）
colnames(sen) = gsub("(.*)\\_(\\d+)", "\\1", colnames(sen))

# 将风险数据和药物敏感性数据合并
rt <- cbind(risk, sen)

# 设置p值过滤阈值
pFilter <- 0.05

# 将风险分组转换为因子类型，并指定分组顺序
rt$risk <- factor(rt$risk, levels = c("Low", "High"))
type <- levels(factor(rt[,"risk"]))     # 获取风险分组的不同水平
comp <- combn(type, 2)                  # 生成两两组合比较
my_comparisons <- list()                # 初始化比较列表
for(i in 1:ncol(comp)) {
  my_comparisons[[i]] <- comp[, i]      # 将组合添加到比较列表中
}

# 对每种药物进行药物敏感性差异分析和绘图
for(drug in colnames(rt)[2:ncol(rt)]) {
  rt1 <- rt[, c(drug, "risk")]          # 提取药物敏感性和风险数据
  colnames(rt1) <- c("Drug", "Risk")    # 重命名列名
  rt1 <- na.omit(rt1)                   # 去除缺失值
  rt1$Drug <- log2(rt1$Drug + 1)        # 对药物敏感性数据进行log2变换
  
  # 使用Wilcoxon秩和检验比较不同风险组间的药物敏感性差异
  test <- wilcox.test(Drug ~ Risk, data = rt1)
  diffPvalue <- test$p.value            # 获取p值
  
  # 如果p值小于设定的阈值，绘制箱线图
  if(diffPvalue < pFilter) {
    boxplot <- ggboxplot(rt1, x = "Risk", y = "Drug", fill = "Risk",
                         xlab = "Risk", ylab = paste0(drug, " sensitivity"),
                         legend.title = "Risk",
                         palette = c("#cf203a", "#209bcd")) + 
      stat_compare_means(comparisons = my_comparisons) # 添加比较统计
    
    # 保存箱线图为PDF文件
    pdf(file = paste0("drugSensitivity.", drug, ".pdf"), width = 5, height = 4.5)
    print(boxplot)
    dev.off()  # 关闭PDF设备
  }
}




library(survival)
library(forestplot) # 用于绘制森林图
library(survminer)
library(tidyverse)

# 读入临床数据
surv <- read.csv("Data/clinical.csv", row.names = 1)
# 将生存数据类型转换成数值型
surv[,c("OS.time", "OS", "Age")] <- lapply(surv[,c("OS.time", "OS", "Age")], as.numeric)

# 以免疫评分为例进行多因素COX回归
# 合并免疫评分与生存数据
ESTIMATE_result <- read.table("Data/TCGA_ESTIMATE_result.txt",sep = "\t", row.names = 1, header = TRUE)
inter_sample <- intersect(rownames(ESTIMATE_result), rownames(surv))
surv <- surv[inter_sample, ]
ESTIMATE_result <- ESTIMATE_result[inter_sample, ]
surv$ImmuneScore <- ESTIMATE_result$ImmuneScore

rt <- filter(surv, if_all(everything(), ~ .x != "'--")) # 去除缺失值

# 去除细分亚型，因为亚型样本少，可能导致模型不稳定或无法收敛
rt$`T` <- gsub("a","", rt$`T`)
rt$`T` <- gsub("b","", rt$`T`)
rt$N <- gsub("a","", rt$N)
rt$N <- gsub("b","", rt$N)
rt$N <- gsub("c","", rt$N)
rt$stage <- gsub("A","",rt$stage)
rt$stage <- gsub("B","",rt$stage)
rt$stage <- gsub("C","",rt$stage)
# 去除异常变量（根据summary(multiCox)结果判断哪个变量异常，如果有的话则去除）
# 因为NX在数据中数量较少，会导致异常，故去除
rt <- rt[rt$N != "NX", ]


# 进行多因素Cox回归分析（所有变量同时分析）
multiCox=coxph(Surv(OS.time, OS) ~ ., data = rt)

multiCoxSum=summary(multiCox)# 获取回归结果的摘要

# 绘制多因素森林图
# 先看是否存在异常值再绘图，否则会报错。若存在，则需在读入数据时去掉异常的分析对象重新分析
summary(multiCox)

# 绘制多因素分析的森林图
ggforest(multiCox,
         main="HR",                      # 图的标题
         cpositions=c(0.02,0.22,0.4),    # 自定义元素的位置
         fontsize=0.7,                   # 字体大小
         refLabel="reference",           # 参考值的标签
         noDigits=2)                     # 显示数值的小数位数
dev.off()

# 准备多因素分析的结果以保存到文件
result = data.frame()
result = cbind(
  HR=multiCoxSum$conf.int[,"exp(coef)"],  # 提取HR值
  HR.95L=multiCoxSum$conf.int[,"lower .95"],  # 提取HR的95%下限
  HR.95H=multiCoxSum$conf.int[,"upper .95"],  # 提取HR的95%上限
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])  # 提取p值
result = cbind(id=row.names(result), result)  # 将变量名作为一列

# 将多因素分析的结果保存
write.csv(result,"Data/Multivariate_result.csv")


# 另一种绘图风格
HT <- result

data = as.data.frame(HT)  # 将数据转换为矩阵
data[, "HR"] <- as.numeric(as.character(data[, "HR"]))
data[, "HR.95L"] <- as.numeric(as.character(data[, "HR.95L"]))
data[, "HR.95H"] <- as.numeric(as.character(data[, "HR.95H"]))
data[, "pvalue"] <- as.numeric(as.character(data[, "pvalue"]))

HR = data[,2:4]  # 提取HR值及其上下限


{
  hr = sprintf("%.3f", HR[,"HR"])  # 格式化HR值
  hrLow = sprintf("%.3f", HR[,"HR.95L"])  # 格式化HR的95%下限
  hrHigh = sprintf("%.3f", HR[,"HR.95H"])  # 格式化HR的95%上限
  pVal = data[,"pvalue"]  # 提取p值
  pVal = ifelse(pVal < 0.001, "<0.001", sprintf("%.3f", pVal))  # 格式化p值，小于0.001的显示为<0.001
  
  # 创建表格内容，包括变量名、p值和HR及其置信区间
  tabletext <- list(c(NA, rownames(HR)),
                    append("pvalue", pVal),
                    append("Hazard ratio", paste0(hr, "(", hrLow, "-", hrHigh, ")")))
}

last_row <- nrow(data)+2 # 根据需要调整最后一行号

# 初始化 hrzl_lines 列表
hrzl_lines <- list("1" = gpar(lwd=2, col="black"), 
                   "2" = gpar(lwd=1.5, col="black"))

# 使用[[ ]] 来动态添加元素
hrzl_lines[[as.character(last_row)]] <- gpar(lwd=2, col="black")


forestplot(labeltext=tabletext,
           mean=c(NA,as.numeric(data$HR)),
           lower=c(NA,as.numeric(data$HR.95L)), 
           upper=c(NA,as.numeric(data$HR.95H)),
           fn.ci_norm="fpDrawDiamondCI",# box类型选择钻石
           col=fpColors(box="#00A896", lines="#02C39A", zero = "black"),# box颜色
           boxsize=0.4,# box大小固定
           lwd.ci=1,
           ci.vertices.height = 0.1,ci.vertices=T,# 显示区间
           zero=1,# zero线横坐标
           lwd.zero=1.5,# zero线宽
           lwd.xaxis=2,
           xlog=TRUE,
           xlab="Hazard ratios",
           txt_gp=fpTxtGp(label=gpar(cex=1.2),# 各种字体大小设置
                          ticks=gpar(cex=0.85),
                          xlab=gpar(cex=1),
                          title=gpar(cex=1.5)),
           hrzl_lines=hrzl_lines,
           new_page = F
)
dev.off()








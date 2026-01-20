library(rms)# 用于构建和分析回归模型
library(survival)
library(tidyverse)

# 读入临床数据
surv <- read.csv("Data/clinical.csv", row.names = 1)
# 将生存数据类型转换成数值型
surv[,c("OS.time", "OS")] <- lapply(surv[,c("OS.time", "OS")], as.numeric)

# 以免疫评分为例构建Nomogram
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


# 初始化数据分布信息，方便后续模型拟合和绘图
ddist=datadist(rt)
options(datadist='ddist')
# 构建Cox比例风险回归模型
#分析生存时间和生存状态与其他变量之间的关系
colnames(rt)
res.cox = cph(Surv(OS.time, OS) ~ demographic.gender + demographic.race +
               ImmuneScore + `T`+ M + N, 
            x = T, y = T, surv = T, data=rt)
# 创建一个用于计算生存概率的函数
surv_cal=Survival(res.cox)
# 构建列线图（Nomogram），预测1年、2年和3年的生存概率
nom <- nomogram(res.cox, fun=list(function(x) surv_cal(1, x), function(x) surv_cal(2, x), 
                                  function(x) surv_cal(3, x)), lp=F, 
                funlabel=c("1-year survival", "2-year survival", "3-year survival"), 
                fun.at=seq(0.90, 0.1,by=-0.1))  
# 绘制列线图
plot(nom, 
     xfrac = 0.4, # 调整比例
     lwd = 2, # 线条宽度
     col.grid = gray(c(0.8, 0.95)), # 网格颜色
     labelsize = 1.2, # 标签文字大小
     cex.axis = 1.1, # 轴的文字大小
)
dev.off()

#输出列线图的风险得分
rt$Nom <- predict(res.cox, type = "lp")

# 输出风险得分到.csv文件
write.csv(rt["Nom"], "Data/Nomo_score.csv", row.names = TRUE)






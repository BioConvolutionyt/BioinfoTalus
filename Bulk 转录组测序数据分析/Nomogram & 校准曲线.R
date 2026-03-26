library(rms) # 用于构建和分析回归模型
library(survival)
library(tidyverse)
library(regplot)

# 读入临床数据
surv <- read.csv("Data/clinical.csv", row.names = 1, stringsAsFactors = FALSE)
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

# 将字符型变量转换成因子
rt[sapply(rt, is.character)] <- lapply(rt[sapply(rt, is.character)], as.factor)
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
     col.grid = gray(c(0.5, 0.85)), # 网格颜色
     labelsize = 1.2, # 标签文字大小
     cex.axis = 1, # 轴的文字大小
)
dev.off()

# 另一种绘图风格
regplot(res.cox,
        observation = TRUE,
        rank = "sd",
        failtime = c(365, 1095, 1825), # 对应1、3、5年
        points = TRUE)
dev.off()

#输出列线图的风险得分
rt$Nom <- predict(res.cox, type = "lp")

# 输出风险得分到.csv文件
write.csv(rt["Nom"], "Data/Nomo_score.csv", row.names = TRUE)


# 校准曲线（1，2，3年）
# 由于示例数据问题，此处绘图效果不佳
col <- c("#4DBBD5FF", "#E64B35FF", "#00A087FF") # 定义颜色
# 1年
f1 <- cph(Surv(OS.time, OS) ~ demographic.gender + demographic.race +
            ImmuneScore + `T`+ M + N, x=T, y=T, surv=T, data=rt, time.inc=1)
cal1 <- calibrate(f1, cmethod="KM", method="boot", u=1, m=(nrow(rt)/3), B=1000)
plot(cal1, xlim = c(0, 1), ylim = c(0, 1),
     xlab = "Nomogram-predicted OS (%)", 
     ylab = "Observed OS (%)", 
     lwd = 2, col = col[1], sub = FALSE)
points(cal1[, "mean.predicted"], 
       cal1[, "KM"], 
       pch = 16, cex = 1.3, col = col[1])
# 2年
f2 <- cph(Surv(OS.time, OS) ~ demographic.gender + demographic.race +
            ImmuneScore + `T`+ M + N, x=T, y=T, surv=T, data=rt, time.inc=2)
cal2 <- calibrate(f2, cmethod="KM", method="boot", u=2, m=(nrow(rt)/3), B=1000)
plot(cal2, xlim = c(0, 1), ylim = c(0, 1),
     xlab = "", 
     ylab = "", 
     lwd = 2, col = col[2], sub = FALSE, add=T) # 如果要分开绘制，add设置为FALSE
points(cal2[, "mean.predicted"], 
       cal2[, "KM"], 
       pch = 16, cex = 1.3, col = col[2])
# 3年
f3 <- cph(Surv(OS.time, OS) ~ demographic.gender + demographic.race +
            ImmuneScore + `T`+ M + N, x=T, y=T, surv=T, data=rt, time.inc=3)
cal3 <- calibrate(f3, cmethod="KM", method="boot", u=3, m=(nrow(rt)/3), B=1000)
plot(cal3, xlim = c(0, 1), ylim = c(0, 1),
     xlab = "", 
     ylab = "", 
     lwd = 2, col = col[3], sub = FALSE, add=T) # 如果要分开绘制，add设置为FALSE
points(cal3[, "mean.predicted"], 
       cal3[, "KM"], 
       pch = 16, cex = 1.3, col = col[3])

# 添加图例，标示不同时间点的校准曲线
legend("bottomright", 
       legend = c("1-Year OS", "2-Year OS", "3-Year OS"),
       col = col,
       lwd = 2, 
       pch = 16,        
       bty = "n",       
       cex = 1.1       # 调整图例字体大小
)
dev.off()



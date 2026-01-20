library(tidyverse)
library(utils)
library(estimate)
library(survminer)
library(survival)

surv <- read.csv("Data/clinical.csv", row.names = 1)
surv$OS.time <- as.numeric(surv$OS.time)/365 # 时间转换成以年为单位
surv <- surv[!is.na(surv$OS.time),]

# 以免疫评分为例进行高低风险组划分
# 合并免疫评分与生存数据
ESTIMATE_result <- read.table("Data/TCGA_ESTIMATE_result.txt",sep = "\t", row.names = 1, header = TRUE)
inter_sample <- intersect(rownames(ESTIMATE_result), rownames(surv))
surv <- surv[inter_sample, ]
ESTIMATE_result <- ESTIMATE_result[inter_sample, ]
surv$ImmuneScore <- ESTIMATE_result$ImmuneScore


# 基于survival寻找最优 cutpoint
res.cut <- surv_cutpoint(
  data = surv,
  time = "OS.time",        # 生存时间变量名
  event = "OS",     # 事件状态变量名（1=死亡，0=删失）
  variables = "ImmuneScore"     # 要寻找 cutpoint 的连续变量
)
cutpoint <- res.cut$cutpoint[[1]]

# 或者直接使用中位数作为cutpoint
# cutpoint <- median(surv$ImmuneScore)

surv$group <- ifelse(surv$ImmuneScore > cutpoint, "Low","High") # 注意这里是免疫评分较高的为低风险组
surv$group <- factor(surv$group, levels = c("Low","High"))

fitd <- survdiff(Surv(OS.time, OS) ~ group,
                 data      = surv,
                 na.action = na.exclude)
pValue <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)

# 拟合生存曲线
fit <- survfit(Surv(OS.time, OS)~ group, data = surv)
summary(fit)
p.lab <- paste0("P", ifelse(pValue < 0.001, " < 0.001", paste0(" = ",round(pValue, 3))))

survplot <- ggsurvplot(fit,
                       data = surv,
                       pval = p.lab,
                       conf.int = FALSE, # 显示置信区间
                       risk.table = TRUE, # 显示风险表
                       risk.table.col = "strata",
                       palette = "aaas", # 配色
                       legend.labs = c("Low", "High"), # 图例
                       size = 1,
                       break.time.by = 3, # x轴步长
                       surv.median.line = "none", # 限制垂直和水平的中位生存
                       ylab = "Survival probability (%)", # 修改y轴标签
                       xlab = "Time (Years)", # 修改x轴标签
                       ncensor.plot = FALSE, # 显示删失图块
                       risk.table.y.text = TRUE)
survplot$plot <- survplot$plot + 
  theme(panel.grid.major = element_line(colour = "gray", linetype = "dashed"),
        panel.grid.minor = element_line(colour = "gray", linetype = "dotted"))
print(survplot)
dev.off()
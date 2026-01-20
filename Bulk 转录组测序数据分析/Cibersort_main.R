library(e1071)
library(parallel)
library(preprocessCore)
library(tidyverse)
source("Bulk 转录组测序数据分析/CIBERSORT.R")

load("Data/TCGA_tpms_matrix.Rda")
sig_matrix <- "Data/LM22.txt"


# 选出肿瘤样本
tpms_matrix <- tpms_matrix[, endsWith(colnames(tpms_matrix), "01A")]
colnames(tpms_matrix) <- substr(colnames(tpms_matrix), 1, 12)
# 必须以这种方式写出的文件才能被CIBERSORT正确地处理
write.table(tpms_matrix, "Data/TCGA_tpms_matrix_cleaned.txt", sep = "\t", row.names = T, col.names = NA, quote = F)


exp = "Data/TCGA_tpms_matrix_cleaned.txt"   #肿瘤患者表达谱
res_cibersort <- CIBERSORT(sig_matrix, exp, perm=100, QN=TRUE)
res_cibersort <- res_cibersort[,1:22]   
ciber.res <- res_cibersort[,colSums(res_cibersort) > 0]   #去除丰度全为0的细胞
ciber.res <- as.data.frame(ciber.res)
save(ciber.res, file="Data/ciber.res.Rda")

# Cibersort彩虹图
mycol <- ggplot2::alpha(rainbow(ncol(ciber.res)), 0.7) #创建彩虹色板（带70%透明度）
par(bty="o", mgp = c(2.5,0.3,0), mar = c(2.1,4.1,2.1,10.1),tcl=-.25,las = 1,xpd = F)
barplot(as.matrix(t(ciber.res)),
        border = NA, # 柱子无边框
        names.arg = rep("",nrow(ciber.res)), # 无横坐标样本名
        yaxt = "n", # 先不绘制y轴
        ylab = "Relative percentage", # 修改y轴名称
        col = mycol) # 采用彩虹色板
axis(side = 2, at = c(0,0.2,0.4,0.6,0.8,1), # 补齐y轴添加百分号
     labels = c("0%","20%","40%","60%","80%","100%"))
legend(par("usr")[2]-60,
       par("usr")[4],
       legend = colnames(ciber.res),
       xpd = T,
       fill = mycol,
       cex = 0.65,
       border = NA,
       y.intersp = 1,
       x.intersp = 0.2,
       bty = "n")
dev.off()   #关闭画板


# 后续可对结果进行分组比较



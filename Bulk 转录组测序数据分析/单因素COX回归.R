library(tidyverse)
library(survival)
library(forestplot)


#读取文件
load("Data/key_gene.Rda")
load("Data/TCGA_tpms_matrix.Rda")
surv <- read.csv("Data/clinical.csv", row.names = 1)

# 选出肿瘤样本
tpms_matrix <- tpms_matrix[, endsWith(colnames(tpms_matrix), "01A")]
colnames(tpms_matrix) <- substr(colnames(tpms_matrix), 1, 12)
tpms_matrix <- as.data.frame(t(tpms_matrix))
# 去除生存时间缺失的样本
surv <- surv[surv$OS.time != "'--", ]

# 合并生存数据与表达矩阵
inter_sample <- intersect(rownames(surv), rownames(tpms_matrix))
surv <- surv[inter_sample, ]
tpms_matrix <- tpms_matrix[inter_sample, ]
surv.expr <- cbind(surv[,c("OS.time", "OS")], tpms_matrix[,key_gene$SYMBOL])
# 将所有数据类型成数值型
surv.expr[] <- lapply(surv.expr, as.numeric)


Coxoutput <- NULL 

for(i in 3:ncol(surv.expr)){
  g <- colnames(surv.expr)[i]
  cox <- coxph(Surv(OS.time,OS) ~ surv.expr[,i], data = surv.expr) # 单变量cox模型
  coxSummary = summary(cox)
  
  Coxoutput <- rbind.data.frame(Coxoutput,
                                data.frame(gene = g,
                                           HR = as.numeric(coxSummary$coefficients[,"exp(coef)"])[1],
                                           z = as.numeric(coxSummary$coefficients[,"z"])[1],
                                           pvalue = as.numeric(coxSummary$coefficients[,"Pr(>|z|)"])[1],
                                           lower = as.numeric(coxSummary$conf.int[,3][1]),
                                           upper = as.numeric(coxSummary$conf.int[,4][1]),
                                           stringsAsFactors = F),
                                stringsAsFactors = F)
}

Coxoutput <- arrange(Coxoutput,pvalue)

#筛选top基因
gene_sig <- Coxoutput[Coxoutput$pvalue < 0.005,] # 取出p值小于0.005的基因
save(gene_sig, file = "Data/gene_cox_sig.Rda")


# 绘制森林图
# 输入表格的制作
tabletext <- cbind(c("Gene",gene_sig$gene),
                   c("HR",format(round(as.numeric(gene_sig$HR),3),nsmall = 3)),
                   c("lower 95%CI",format(round(as.numeric(gene_sig$lower),3),nsmall = 3)),
                   c("upper 95%CI",format(round(as.numeric(gene_sig$upper),3),nsmall = 3)),
                   c("pvalue",format(round(as.numeric(gene_sig$p),3),nsmall = 3)))


last_row <- nrow(gene_sig)+2 # 根据需要调整最后一行号

# 初始化 hrzl_lines 列表
hrzl_lines <- list("1" = gpar(lwd=2, col="black"), 
                   "2" = gpar(lwd=1.5, col="black"))

# 使用[[ ]] 来动态添加元素
hrzl_lines[[as.character(last_row)]] <- gpar(lwd=2, col="black")


forestplot(labeltext=tabletext,
           mean=c(NA,as.numeric(gene_sig$HR)),
           lower=c(NA,as.numeric(gene_sig$lower)), 
           upper=c(NA,as.numeric(gene_sig$upper)),
           fn.ci_norm="fpDrawDiamondCI",# box类型选择钻石
           col=fpColors(box="#00A896", lines="#02C39A", zero = "black"),# box颜色
           boxsize=0.4,# box大小固定
           lwd.ci=1,
           ci.vertices.height = 0.1,ci.vertices=T,# 显示区间
           zero=1,# zero线横坐标
           lwd.zero=1.5,# zero线宽
           # xticks = c(0.5,1,1.5), # 横坐标刻度根据需要可随意设置
           lwd.xaxis=2,
           xlab="Hazard ratios",
           txt_gp=fpTxtGp(label=gpar(cex=1.2),# 各种字体大小设置
                          ticks=gpar(cex=0.85),
                          xlab=gpar(cex=1),
                          title=gpar(cex=1.5)),
           hrzl_lines=hrzl_lines,
           new_page = F
)
dev.off()

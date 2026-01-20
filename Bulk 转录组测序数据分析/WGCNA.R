library(WGCNA)
library(tidyverse)

####数据整理####
load("Data/TCGA_tpms_matrix.Rda")
datExpr0 = as.data.frame(t(tpms_matrix))

#主要看缺失值
gsg = goodSamplesGenes(datExpr0, verbose = 3)

gsg$allOK 
if (!gsg$allOK){
  # 把含有缺失值的基因或样本打印出来
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # 去掉那些缺失值
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

#样本过滤
sampleTree = hclust(dist(datExpr0), method = "average")

par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
dev.off()


#如果有异常值就需要去掉，根据聚类图自己设置cutHeight 参数的值
#单独的图
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2) +
  #想用哪里切，就把“h = 320”和“cutHeight = 320”中换成你的cutoff
  abline(h = 320, col = "red")

clust = cutreeStatic(sampleTree, cutHeight = 320, minSize = 10)
keepSamples = (clust==1)    #1是保留的
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
dim(datExpr)

#没有异常样本就不需要去除
datExpr = datExpr0
rm(datExpr0)

####表型信息整理####
# 示例数据中表型信息仅有Tumor/Normal
clinical <- data.frame(sampleID = rownames(datExpr),
                       Normal = ifelse(endsWith(rownames(datExpr), "01A"), 1, 0),
                       Tumor = ifelse(endsWith(rownames(datExpr), "01A"), 0, 1))
clinical <- column_to_rownames(clinical, "sampleID")

datTraits = clinical[c("Normal","Tumor")]
sampleTree2 = hclust(dist(datExpr), method = "average")
# 各个样本的表现: 白色表示低，红色为高，灰色为缺失
traitColors = numbers2colors(datTraits, signed = FALSE)
# 把样本聚类和表型绘制在一起
# pdf("Phenotypic clustering.pdf", 13, 7)
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")
dev.off()



####正式分析####
#软阈值的筛选
#设置一系列软阈值，范围是1-30之间
powers = c(1:10, seq(from = 12, to=30, by=2))

sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

#这个结果就是推荐的软阈值
sft$powerEstimat


cex1 = 0.9 #一般是0.9.不能低于0.85
# pdf("Soft threshold selection.pdf", 11, 6)
par(mfrow = c(1,2))
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
abline(h=cex1,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()


#进一步构建网络
power = sft$powerEstimate
power#如果前面没有得到推荐的软阈值，就要根据上面的图，自己选择

#需要时间
net = blockwiseModules(datExpr, power = power,
                       TOMType = "unsigned", 
                       minModuleSize = 50, # minModuleSize 默认30，参数设置最小模块的基因数，值越小，小的模块就会被保留下来
                       reassignThreshold = 0, 
                       mergeCutHeight = 0.6,# mergeCutHeight 默认0.25，设置合并相似性模块的距离，值越小，就越不容易被合并，保留下来的模块就越多。
                       deepSplit = 1 ,# deepSplit 默认2，调整划分模块的敏感度，值越大，越敏感，得到的模块就越多
                       numericLabels = TRUE,
                       pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "testTOM",
                       verbose = 3)


#此处展示得到了多少模块，每个模块里面有多少基因。
table(net$colors)
#如果结果不合理，可以适当调整参数进行修改

#具体细节可参考：https://zhuanlan.zhihu.com/p/34697561

mergedColors = labels2colors(net$colors)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
#每种颜色都代表着一种基因模块，即里面的基因功能是相似的
#灰色没意义

#保存每个模块对应的基因
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs
geneTree = net$dendrograms[[1]]
gm = data.frame(net$colors)
gm$color = moduleColors
head(gm)

genes = split(rownames(gm),gm$color)
save(genes,file = "Data/genes.Rdata")


#模块与表型的相关性
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)


#热图
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)

# pdf("Module-trait relationships.pdf", 8, 8)
par(mar = c(6, 9, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed (100),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.8,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()


#相关系数最好大于0.8，实在没有，0.6或0.7也行


#把gene module输出到文件
text <- unique(moduleColors)
for (i  in 1:length(text)) {
  y=t(assign(paste(text[i],"expr",sep = "."),datExpr[moduleColors==text[i]]))
  write.csv(y,paste(paste("./", text[i], sep=""),"csv",sep = "."),quote = F)
}


#GS与MM
#GS代表模块里的每个基因与形状的相关性
#MM代表每个基因和所在模块之间的相关性，表示是否与模块的趋势一致
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")


#第几列的表型是最关心的，下面的i就设置为几。
#与关心的表型相关性最高的模块赋值给下面的module。
traitData = datTraits 
i = 2 #替换自己想要的表型列
module = "turquoise"##替换自己的样本颜色
assign(colnames(traitData)[i],traitData[i])
instrait = eval(parse(text = colnames(traitData)[i]))
geneTraitSignificance = as.data.frame(cor(datExpr, instrait, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(instrait), sep="")
names(GSPvalue) = paste("p.GS.", names(instrait), sep="")
column = match(module, modNames) #找到目标模块所在列
moduleGenes = moduleColors==module #找到模块基因所在行
# par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()

# ggplot2风格
# 示例数据 - 替换为你的实际数据
mm_values <- abs(geneModuleMembership[moduleGenes, column])
gs_values <- abs(geneTraitSignificance[moduleGenes, 1])

# 创建一个数据框
df <- data.frame(MM = mm_values, GS = gs_values)

# 计算 Pearson 相关系数和 p 值
cor_test <- cor.test(df$MM, df$GS, method = "pearson")
r_val  <- round(cor_test$estimate, 3)   # 相关系数 r
p_val  <- cor_test$p.value              # p 值

# 格式化 p 值（避免科学计数法太小）
p_str <- ifelse(p_val < 0.001, "p < 0.001", paste("p =", round(p_val, 3)))
# 构建标签文本
label_text <- paste("r =", r_val, "\n", p_str)
# pdf("Module Membership vs. Gene Significance.pdf", 7, 6.5)
p <- ggplot(df, aes(x = MM, y = GS)) +
  geom_point(shape = 1, color = "turquoise", size = 2) + # 绘制散点图
  geom_smooth(method = "lm", se = TRUE, color = "black", size = 0.8) + # 添加线性拟合曲线及置信区间
  annotate("text", x = 0.1, y = Inf, label = label_text,
           hjust = -0.1, vjust = 1.5, size = 5)+
  labs(
    title = "Module Membership vs. Gene Significance",
    x = paste("Module Membership in", module, "module"),
    y = "Gene significance"
  ) +
  theme_classic(base_size = 14) # 设置主题和字体大小
p + theme(plot.title = element_text(hjust = 0.5))

dev.off()

# 注意：如果分析的数据为肿瘤数据则不建议去除细胞周期的影响！！！
# 由于肿瘤的发生与细胞周期相关的基因有一定联系，去除细胞周期可能导致错误的分析结果
# 一般不建议对肿瘤单细胞测序数据进行细胞周期去除
library(Seurat)
library(data.table)
library(stringr)
library(tibble)
load("./Data/pbmc.Rda")


#获取G2M期相关基因
g2m_genes <- cc.genes$g2m.genes
g2m_genes <- CaseMatch(search=g2m_genes, match=rownames(pbmc))

#获取S期相关基因
s_genes <- cc.genes$s.genes    
s_genes <- CaseMatch(search=s_genes, match=rownames(pbmc))

#细胞周期阶段评分
pbmc <- CellCycleScoring(pbmc, g2m.features=g2m_genes, s.features=s_genes)

# 查看细胞周期的影响情况
DimPlot(pbmc, group.by = "Phase")

# SCT
pbmc <- SCTransform(pbmc, vars.to.regress = c("S.Score", "G2M.Score")) # 将细胞周期的影响回归掉

save(pbmc, file="Data/pbmc.Rda")

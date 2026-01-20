library(dplyr)
library(Seurat)
library(tidyverse)
library(patchwork)
library(infercnv)
load("Data/pbmc.Rda")

table(pbmc$celltype)

exprMatrix <- as.matrix(GetAssayData(pbmc, layer ='counts'))
cellAnnota <- subset(pbmc@meta.data, select='celltype')
dim(exprMatrix)
write.table(cellAnnota,file ="./Data/groupFiles.txt",sep = '\t',col.names = F)

table(cellAnnota$celltype)


infercnv_obj = CreateInfercnvObject(raw_counts_matrix=exprMatrix,
                                    annotations_file="./Data/groupFiles.txt",
                                    delim="\t",
                                    gene_order_file= "geneLocate.txt",
                                    ref_group_names="T_cells") #infercnv分析拷贝数变异时参照的正常细胞

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, 
                             out_dir=  './Data/infercnv/' ,  
                             cluster_by_groups=T,   
                             hclust_method="ward.D2", 
                             plot_steps=F,
                             denoise=TRUE,#去噪音
                             HMM=F,
                             output_format="pdf",
                             num_threads=32,
                             write_expr_matrix = T
)

save(infercnv_obj, file="Data/infercnv_obj.Rda")

grp=read.table("./Data/infercnv/infercnv.observation_groupings.txt",sep = "",header = T)
obs=read.table("./Data/infercnv/infercnv.observations.txt", header = T,check.names = F)
max(obs)
min(obs)
obs[obs>0.8 & obs<0.93]=2
obs[obs>=0.93 & obs<0.95]=1
obs[obs>=0.95 & obs<1.05]=0
obs[obs>=1.05 & obs<1.07]=1
obs[obs>=1.07 & obs<1.2]=2

scores=as.data.frame(colSums(obs))


scores$celltype=grp$Dendrogram.Group
colnames(scores)=c("score","celltype")
save(scores, file="Data/inferCNV_scores.Rda")


library(RColorBrewer)
infercnv::plot_cnv(infercnv_obj, #上两步得到的infercnv对象
                   plot_chr_scale = T, #画染色体全长，默认只画出（分析用到的）基因
                   output_filename = "Data/Custom_plot",output_format = "pdf", #保存为pdf文件
                   custom_color_pal =  color.palette(c("#8DD3C7","white","#BC80BD"), c(2, 2))) #改颜色


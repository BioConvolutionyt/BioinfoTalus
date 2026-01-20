library(Seurat)
library(ggplot2)
library(patchwork)
library(ROGUE)
library(SingleR)
library(scRNAtoolVis)
library(SCpubr)
library(ggsci)
library(randomcoloR)
load("Data/pbmc.Rda")

####自动注释####
# 使用SingleR进行自动注释的准确性和认可度没有人工注释高，一般都进行人工注释
# 但自动注释结果可作为人工注释有效的辅助验证

#加载参考数据库
load("Data/ref_Human_all.RData")

#获取表达矩阵
testdata = GetAssayData(object = pbmc@assays$RNA, layer = "counts")

#获取clusters
clusters <- pbmc@meta.data$seurat_clusters

#运行singleR
cellpred <- SingleR(test = testdata, ref = ref_Human_all, clusters = clusters, assay.type.test = "logcounts", 
                    labels = ref_Human_all@colData@listData[["label.main"]], assay.type.ref = "logcounts")

#获取每个cluster的细胞类型
celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = F)

#添加到seurat.metadata对象中
pbmc@meta.data$SingleR = "NA"
for(i in 1:nrow(celltype)){
  pbmc@meta.data[which(pbmc$seurat_clusters == celltype$ClusterID[i]),'SingleR'] <- celltype$celltype[i]
}

# 画图
# 注释的细胞类型在每个参考细胞类型上的相似度热图，可用于评估注释效果
p = plotScoreHeatmap(cellpred)

p1 <- DimPlot(pbmc, group.by = "SingleR", label = T,reduction = "tsne")
p2 <- DimPlot(pbmc, group.by = "SingleR", label = T,reduction = "umap")
p <- p1 | p2
p

save(pbmc, file="Data/pbmc.Rda")

####人工注释####
#寻找marker
DefaultAssay(pbmc) = "RNA"
pbmc.markers.RNA <- FindAllMarkers(pbmc, only.pos = TRUE, logfc.threshold = 1)
save(pbmc.markers.RNA, file="Data/markers.RNA.Rda")

# 创建marker集合
# marker的选定可基于相关的文献
# 注意：不同组织，甚至相同组织不同疾病的细胞marker都不一定是通用的
# 这里marker的质量较高，每个cluster属于哪种细胞类型比较明确
# 但大部分情况不会这么顺利，需要不断换marker
marker <- list("Mast_cells" = c("TPSAB1", "GATA2", "MS4A2"),
               "Smooth_muscle_cells" = c("TAGLN", "ACTA2", "MYH11"),
               "Myeloid_cells" = c("LYZ", "CD14", "CD68", "CSF1R"),
               "Endothelia_cells" = c("PECAM1", "CLDN5", "VWF"),
               "Epithelial_cells"=c("KRT8", "CLDN4","TSPAN8"),
               "Fibroblasts"=c("COL1A1", "DCN", "LUM"),
               "B_cells"=c("HLA-DRA", "MS4A1"),
               "Plasma_cells"=c("MZB1", "JCHAIN"),
               "T_cells"=c("CD3D", "CD2", "TRBC2"),
               "NE" = c("FCGR3B", "CSF3R", "CXCR2", "G0S2")
)


# 气泡图
# 如果出现报错The facets argument of facet_grid() was deprecated .....，请更新Seurat版本至5.4.0
do_DotPlot(sample = pbmc, features = marker, dot.scale = 10,
           legend.length = 50, legend.framewidth = 2, font.size = 10)

#另一种颜色气泡图
DotPlot(pbmc, features = unique(as.vector(unlist(marker))), scale = FALSE, cols = "RdYlBu")+RotatedAxis()

# 如果不好注释，可调整聚类分辨率重新聚类（仅需要运行下面的代码即可）
# 例如某个簇同时表达两种细胞类型的marker，则可以尝试提高聚类分辨率，看能否将其聚为不同的两类再进行注释
DefaultAssay(pbmc) = "SCT"
pbmc <- FindClusters(pbmc, resolution = 0.8, algorithm = 4)
# 除了看marker的显著程度，还可以观察两个簇在umap降维散点图中的位置
# 如果两个簇相距较近则细胞类型可能比较相似，相距较远则细胞类型可能相差较大
DimPlot(pbmc, reduction = "umap.harmony", label = T)


# 输入注释结果（每个簇对应的细胞类型）
pbmc$celltype <- recode(pbmc@meta.data$seurat_clusters,
                        "1" = "Myeloid_cells",
                        "2" = "Epithelial_cells",
                        "3" = "NE",
                        "4" = "Myeloid_cells",
                        "5" = "T_cells",
                        "6" = "Epithelial_cells",
                        "7" = "T_cells",
                        "8" = "Plasma_cells",
                        "9" = "Epithelial_cells",
                        "10" = "Plasma_cells",
                        "11" = "B_cells",
                        "12" = "Endothelia_cells",
                        "13" = "Plasma_cells",
                        "14" = "Fibroblasts", 
                        "15" = "Mast_cells", 
                        "16" = "Epithelial_cells", 
                        "17" = "Smooth_muscle_cells",
)

Idents(pbmc) = "seurat_clusters"
Idents(pbmc) = "celltype"
# 展示注释结果
DimPlot(pbmc, reduction = "umap", label = T, label.size = 3.5)+theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),
                                                                                     legend.position = "right")
DimPlot(pbmc, reduction = "tsne", label = T, label.size = 3.5)+theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right")


cors <- pal_igv()(30) # 也可以自定义颜色
DimPlot(pbmc,reduction = "tsne",label = T,group.by="seurat_clusters",cols = cors)
DimPlot(pbmc,reduction = "umap",label = T,group.by="celltype",cols = cors)


# 后续为绘图代码，可选做
# 寻找每个细胞类型的markers
pbmc.markers.celltype <- FindAllMarkers(pbmc, only.pos = TRUE,logfc.threshold = 1,min.pct = 0.3)
save(pbmc.markers.celltype, file="Data/markers.celltype.Rda")

# 提取前5的marker
top5 <- pbmc.markers.celltype %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)

# 对这些基因进行scale data
markers = as.data.frame(top5[,"gene"])
pbmc <- ScaleData(pbmc, features = as.character(unique(markers$gene)))

#绘制热图
DoHeatmap(pbmc,
          features = as.character(unique(markers$gene)),
          group.by = "celltype")+
  ggsci::scale_colour_npg() +
  scale_fill_gradient2(low = '#0099CC', mid = 'white', high = '#CC0033', 
                       name = 'Z-score')

#定义颜色
colaa=distinctColorPalette(100)

#单基因小提琴图
VlnPlot(pbmc, features = "CD3D",group.by="celltype")+NoLegend()

# umap单个基因图
FeaturePlot(pbmc, features = "CD3D")

#多基因分布图
symbol=c("CD3D", "CD2", "TRBC2")

VlnPlot(pbmc, features = symbol,group.by = "celltype", stack=TRUE,cols = colaa)+ NoLegend()   
#分组分半小提琴图
p<-VlnPlot(pbmc, features = symbol,stack=T,pt.size=0,flip = T,split.by = 'celltype',
           group.by = "celltype",
           cols = c("#78C2C4","#C73E3A"),
           split.plot = T)+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(colour = 'black',size = 10,angle = 90),
        legend.position = 'none')
p

#堆叠柱状图
cell.prop<-as.data.frame(prop.table(table(pbmc@meta.data$celltype, pbmc@meta.data$orig.ident)))
colnames(cell.prop)<-c("cluster","group","proportion")

ggplot(cell.prop,aes(group,proportion,fill=cluster))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+
  theme_bw()+
  theme(axis.ticks.length=unit(0.5,'cm'))+
  guides(fill=guide_legend(title=NULL))

save(pbmc, file="Data/pbmc.Rda")




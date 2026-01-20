library(Seurat)
library(ggplot2)
library(slingshot)
library(RColorBrewer)
library(fields)
library(scop) # 一个端到端的单细胞组学分析 + 高级绘图包
library(zellkonverter)

load("Data/subgroups.Rda")

# 赋值给临时变量，以便于后续无需修改变量名
pbmc <- subgroups

# 换个颜色
nature_colors <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
                   "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85",
                   "#FFB000", "#A20056")

CellDimPlot(
  srt = pbmc,
  group.by = c("celltype"),
  reduction = "UMAP",
  palcolor = nature_colors,
  label = TRUE,
  label_insitu = TRUE,
  theme_use = "theme_blank"
)


# 转换为SingleCellExperiment对象并进行slingshot分析
# 使用PCA降维后的结果进行slingshot一般比UMPA更加稳定
reduced <- "pca"
sce <- as.SingleCellExperiment(pbmc)
reducedDim(sce, toupper(reduced)) <- Embeddings(pbmc, reduced)
sce <- slingshot(data = sce,
                 clusterLabels = 'celltype',
                 reducedDim = toupper(reduced),
                 # 分化起始点需要人为设定
                 start.clus = "Naïve_T_cells",
                 end.clus = "Effector_T_cells")

# 绘制轨迹图
sds <- SlingshotDataSet(sce)
save(sds, file = "Data/slingshot_results.Rda")
sds

# Seurat风格
# 提取轨迹信息
curves_list <- sds@curves
lineage_names <- names(curves_list)

get_curve_df <- function(curve, lineage, cut_frac = 0.1, dims = c(1, 2)) {
  s <- curve[["s"]]
  if (is.null(s) || nrow(s) < 2) return(NULL)
  
  df <- as.data.frame(s[, dims, drop = FALSE])
  colnames(df) <- c("Dim1", "Dim2")
  
  # 两端裁剪（防止端点过长/不稳定）
  cut_n <- floor(nrow(df) * cut_frac)
  if (cut_n > 0 && nrow(df) > 2 * cut_n) {
    df <- df[(cut_n + 1):(nrow(df) - cut_n), , drop = FALSE]
  }
  
  df$lineage <- lineage
  df$pt_order <- seq_len(nrow(df))   # 保证每条曲线内部按顺序连线
  df
}

curves_df_all <- do.call(
  rbind,
  Map(get_curve_df, curves_list, lineage_names)
)


p <- CellDimPlot(
  srt = pbmc,
  group.by = "celltype",
  palcolor = nature_colors,
  reduction = toupper(reduced),
  pt.size = 1,
  theme_use = "theme_blank"
)

p +
  ggnewscale::new_scale_color() +
  geom_path(
    data = curves_df_all,
    aes(Dim1, Dim2, group = lineage),
    colour = "black",
    linewidth = 2,
    arrow = arrow(length = unit(0.3, "cm"))
  ) +
  geom_path(
    data = curves_df_all,
    aes(Dim1, Dim2, group = lineage, color = lineage),
    linewidth = 1.0,
    arrow = arrow(length = unit(0.3, "cm"))
  ) +
  scale_color_brewer(palette = "Set2", name = "lineage")




# 伪时间可视化
pbmc$slingPseudotime <- sce$slingPseudotime_1

FeatureDimPlot(
  srt = pbmc,
  features = c("slingPseudotime"),
  reduction =toupper(reduced),
  palcolor = nature_colors,
  pt.size = 1,
  theme_use ="theme_blank",
  theme_args = list(
    legend.key.size = grid::unit(0.7, 'cm'), # 调整colorbar的宽度
    plot.title = element_blank()
  )
)+
  geom_path(data = curves_df, 
            aes(x = Dim1, y = Dim2), 
            group = 1,
            linewidth = 1,
            arrow = arrow(length = unit(0.3, "cm")))

subgroups <- pbmc
save(subgroups, file="./Data/subgroups.Rda")

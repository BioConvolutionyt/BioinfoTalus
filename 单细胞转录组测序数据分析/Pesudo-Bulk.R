library(Seurat)   
library(DESeq2)
library(limma) 
library(Matrix)
load("Data/pbmc.Rda")


# 主函数
# 函数接收两个Seurat对象
# 按照分组变量group.by构造伪变量进行两个Seurat对象之间的pseudo_bulk差异分析
# 若两个Seurat对象构造的伪变量不一致，intersect_groups参数决定是否取交集再进行后续分析
# assay用于指定参与差异分析的表达矩阵
# method用于指定差异分析方法，可选 edgeR（基于负二项分布），limma（基于线性模型）
pesudoBulkDE <- function(obj1, obj2,
                             group.by,
                             intersect_groups = TRUE,
                             method = c("edgeR", "limma")) {
  
  method <- match.arg(method)
  
  # 提取分组信息
  meta1 <- obj1@meta.data
  meta2 <- obj2@meta.data
  
  grps1 <- unique(as.character(meta1[[group.by]]))
  grps2 <- unique(as.character(meta2[[group.by]]))
  
  common_groups <- intersect(grps1, grps2)
  use_groups <- if (intersect_groups) common_groups else union(grps1, grps2)
  
  # 为单个对象构建伪样本计数和样本信息
  build_pseudo <- function(obj, obj_label) {
    counts <- GetAssayData(obj, layer = "counts", assay = DefaultAssay(obj))
    meta <- obj@meta.data
    genes <- rownames(counts)
    
    # 只对在 use_groups 中且对象中存在的分组进行伪样本构建
    groups_available <- intersect(use_groups, unique(as.character(meta[[group.by]])))
    if (length(groups_available) == 0) {
      return(list(counts = NULL, sample_info = data.frame(), genes = genes))
    }
    
    pseudo_counts <- NULL
    sample_info <- data.frame(sample_id = character(0),
                              object = character(0),
                              group = character(0),
                              stringsAsFactors = FALSE)
    
    for (g in groups_available) {
      cells_in_group <- rownames(meta)[which(meta[[group.by]] == g)]
      if (length(cells_in_group) == 0) next
      
      sub_counts <- counts[, cells_in_group, drop = FALSE]
      summed <- Matrix::rowSums(sub_counts)
      
      colname <- paste0(obj_label, "_", g)
      if (is.null(pseudo_counts)) {
        pseudo_counts <- matrix(summed, ncol = 1)
        rownames(pseudo_counts) <- genes
        colnames(pseudo_counts) <- colname
      } else {
        pseudo_counts <- cbind(pseudo_counts, summed)
        colnames(pseudo_counts)[ncol(pseudo_counts)] <- colname
      }
      
      sample_info <- rbind(sample_info,
                           data.frame(sample_id = colname,
                                      object = obj_label,
                                      group = as.character(g),
                                      stringsAsFactors = FALSE))
    }
    
    list(counts = pseudo_counts, sample_info = sample_info, genes = genes)
  }
  
  # 构建两个对象的伪样本
  pc1 <- build_pseudo(obj1, "Obj1")
  pc2 <- build_pseudo(obj2, "Obj2")
  
  # 取公共基因，确保两边矩阵同维
  common_genes <- intersect(rownames(pc1$counts), rownames(pc2$counts))
  
  counts1 <- pc1$counts[common_genes, , drop = FALSE]
  counts2 <- pc2$counts[common_genes, , drop = FALSE]
  counts_combined <- cbind(counts1, counts2)
  
  sample_info_combined <- rbind(pc1$sample_info, pc2$sample_info)
  
  # 调整列名以匹配 counts_combined 的列顺序
  colnames(counts_combined) <- paste0(sample_info_combined$object,
                                      "_",
                                      sample_info_combined$group)
  
  # 设计矩阵：object 与 group 两个因素的加性效应
  object_factor <- factor(sample_info_combined$object)
  group_factor <- factor(sample_info_combined$group)
  design <- model.matrix(~ object_factor + group_factor)
  
  coef_names <- colnames(design)
  object_coefs <- coef_names[grep("object_factor", coef_names)]
  coef_index <- if (length(object_coefs) > 0) {
    which(coef_names == object_coefs[length(object_coefs)])
  } else {
    2  # 回退值
  }
  
  # 差异分析
  if (method == "edgeR") {
    # edgeR 路径
    dge <- DGEList(counts = counts_combined)
    dge$samples$group <- object_factor
    dge <- estimateDisp(dge, design)
    fit <- glmFit(dge, design)
    lrt <- glmLRT(fit, coef = coef_index)
    res <- as.data.frame(lrt$table)
    res$adj.P.Val <- p.adjust(res$PValue, method = "BH")
    # 组织输出：gene 列 + 统计量
    res$gene <- rownames(lrt$table)
  } else {
    # limma 路径
    voom_res <- voom(counts_combined, design, plot = FALSE)
    fit <- lmFit(voom_res, design)
    fit <- eBayes(fit)
    res <- topTable(fit, coef = coef_index, number = Inf, sort.by = "P")
    res <- as.data.frame(res)
    res$gene <- rownames(res)
  }
  
  # 结果整理
  res <- res[, c("gene", setdiff(colnames(res), "gene"))]
  rownames(res) <- NULL
  
  return(res)
}


# 这里以Myeloid_cells与Epithelial_cells之间的差异基因分析为例
seurat_obj1 <- subset(pbmc, celltype=="Myeloid_cells")
seurat_obj2 <- subset(pbmc, celltype=="Epithelial_cells")

# 注意：分组变量group.by的选择非常关键！选择不当可能导致报错或无意义的结果
# 两个对象都需要在元数据中存在同名的group.by列
# 且该列的分组在两对象中应具有相同的生物学意义
# 通常会用能代表批次、样本或细胞类型的字段
# 若分析同类样本，不同细胞间的差异基因，如肿瘤样本中T细胞与B细胞间的差异基因，可选group.by = "orig.ident"
# 若分析非同类样本之间的差异基因，如肿瘤与正常组织间的差异基因，可选group.by = "celltype"
# 一般不用seurat_clusters做分组变量
deg_results <- pesudoBulkDE(seurat_obj1, 
                           seurat_obj2, 
                           group.by = "orig.ident", 
                           method = "edgeR")

# 筛选
# deg_result <- deg_result[deg_result$adj.P.Val < 0.05, ]

save(deg_results, file="Data/pesudoBulkDE.Rda")




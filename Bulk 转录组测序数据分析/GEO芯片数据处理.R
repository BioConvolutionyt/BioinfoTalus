# 不同平台的芯片数据注释方式可能存在差异，这里以GPL17586平台为例
# 该脚本仅适用于人类数据集
library(data.table)
library(stringr)
library(biomaRt)

# 读入平台注释文件
GPL <- read.delim("Data/GEO/GPL17586-45144.txt", header = TRUE, check.names = FALSE, row.names = 1, quote = "", comment.char = "")
# 读入表达矩阵
GEO_matrix <- read.table("Data/GEO/GSE74777_series_matrix.txt", sep = "\t", header = T, row.names = 1, check.names = F) 


# 筛选
GPL <- GPL[, c("gene_assignment", "mrna_assignment")]
GPL$gene_symbol <- ""

# 定义函数
# 判断基因名称是否合法
is_legal_gene <- function(x) {
  !is.na(x) &
    !grepl("^(OTTHUMG|LOC)\\d+$", x) &
    x != "---" &
    !grepl(" ", x)
}

# 提取从字符串中基因信息
extract_gene_info <- function(gene_assignment, mrna_assignment){
  
  ga <- strsplit(gene_assignment, "//", fixed = TRUE)[[1]]
  ga <- trimws(ga)
  
  if (length(ga) >= 2 && is_legal_gene(ga[2])) return(ga[2])
  if (length(ga) >= 7 && is_legal_gene(ga[7])) return(ga[7])
  
  ma <- strsplit(mrna_assignment, "//", fixed = TRUE)[[1]]
  ma <- trimws(ma)
  
  if (length(ma) >= 3 && ma[2] == "RefSeq") {
    gene <- sub(".*\\(([^)]*)\\).*", "\\1", ma[3])
    return(ifelse(is_legal_gene(gene), gene, ""))
  }
  
  if (length(ma) >= 3 && ma[2] == "ENSEMBL") {
    gene <- regmatches(ma[3], regexpr("ENSG\\d{11}(?:\\.\\d+)?", ma[3]))
    return(ifelse(is_legal_gene(gene), gene, ""))
  }
  
  ""
}

# 执行注释
GPL$gene_symbol <- vapply(
  seq_len(nrow(GPL)),
  function(i) {
    extract_gene_info(GPL$gene_assignment[i], GPL$mrna_assignment[i])
  },
  character(1)
)

# 筛选
GPL <- GPL[GPL$gene_symbol != "", ]
GPL_ensg <- GPL[startsWith(GPL$gene_symbol, "ENSG"), ]

# 对ENSEMBL ID进一步进行注释
mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")  # 人类
res <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters    = "ensembl_gene_id",
  values     = GPL_ensg$gene_symbol,
  mart       = mart
)
res <- res[res$hgnc_symbol != "", ]

idx <- match(GPL_ensg$gene_symbol, res$ensembl_gene_id)
valid <- !is.na(idx)

GPL_ensg$gene_symbol[valid] <- res$hgnc_symbol[idx[valid]]

# 对在线注释失败的ENSEMBL ID做进一步处理
GPL_ensg_blank <- GPL_ensg[startsWith(GPL_ensg$gene_symbol, "ENSG"), ]
annotation <- read.csv("Data/GEO/annotation.csv", header = T) # 来自TCGA，人工整理的注释文件

inter_ID <- intersect(annotation$gene_id, GPL_ensg_blank$gene_symbol)
annotation <- annotation[annotation$gene_id %in% inter_ID, ]
annotation <- unique(annotation)

ann_map <- setNames(annotation$gene_name, annotation$gene_id)
idx <- match(GPL_ensg_blank$gene_symbol, names(ann_map))
valid <- !is.na(idx)

GPL_ensg_blank$gene_symbol[valid] <- ann_map[idx[valid]]

# 将本地注释和在线注释结果合并
GPL_ensg[rownames(GPL_ensg_blank), "gene_symbol"] <- 
  GPL_ensg_blank$gene_symbol
# 将ENSEMBL ID注释结果与原始注释结果合并
GPL[rownames(GPL_ensg), "gene_symbol"] <- GPL_ensg$gene_symbol
GPL <- GPL[!startsWith(GPL$gene_symbol, "ENSG"), ]

# 将平台注释文件探针ID与芯片数据探针ID取交集
inter_probe <- intersect(rownames(GPL), rownames(GEO_matrix))
# 合并表达矩阵与注释结果
GPL <- GPL[inter_probe, ]
GEO_matrix <- GEO_matrix[inter_probe, ]
GEO_matrix <- cbind(gene_symbol = GPL$gene_symbol, GEO_matrix)

# 合并重复探针
dt <- as.data.table(GEO_matrix)
GEO_matrix <- dt[
  ,
  lapply(.SD, mean, na.rm = TRUE),
  by = gene_symbol
]

GEO_matrix <- column_to_rownames(GEO_matrix, "gene_symbol")

# 标准化（若数据上传者明确说明已做过标准化则可跳过这一步，这里直接跳过）
# boxplot(GEO_matrix, outline=FALSE, notch=T, las=2)
# GEO_matrix = normalizeBetweenArrays(exGEO_matrixp)

save(GEO_matrix, file="Data/GEO_matrix.Rda")


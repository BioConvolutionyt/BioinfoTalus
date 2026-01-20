library(rjson)
library(tidyverse)

#metadata文件名
json <- jsonlite::fromJSON("Data/TCGA/RNA seq/metadata.cart.2025-09-22.json")  
#View(json)

{
  sample_id <- sapply(json$associated_entities,function(x){x[,1]})
  file_sample <- data.frame(sample_id,file_name=json$file_name)  
}
#Counts文件夹名
count_file <- list.files("Data/TCGA/RNA seq/counts files",pattern = '[.]tsv$',recursive = TRUE)  
{
  count_file_name <- strsplit(count_file, split='/')
  count_file_name <- sapply(count_file_name,function(x){x[2]})
}

matrix = data.frame(matrix(nrow=60660,ncol=0))

#下面的修改样本例数
for (i in 1:nrow(json)){
  path = paste0('Data/TCGA/RNA seq/counts files/',count_file[i])  
  data<- read.delim(path,fill = TRUE,header = FALSE,row.names = 1)
  colnames(data)<-data[2,]
  data <-data[-c(1:6),]
  # 注意：需要counts则用data[3]，需要tpms则用data[6]
  data <- data[3]   #TCGA提供的数据类型包括 3：unstranded；4：stranded_first；5：stranded_second；6：tpm_unstranded；7：fpkm_unstranded；8：fpkm_uq_unstranded
  colnames(data) <- file_sample$sample_id[which(file_sample$file_name==count_file_name[i])]
  matrix <- cbind(matrix,data)
}

# 提取注释信息
tsv_files <- list.files(
  path = "Data/TCGA/RNA seq/counts files",
  pattern = "\\.tsv$",          
  recursive = TRUE,             
  full.names = TRUE             
)

names <- read.table(tsv_files[1], sep = "\t", header = T, row.names = 1)
names <- names[-c(1:4), c(1:2)]

# 获取矩阵和基因名的交集
same <- intersect(rownames(matrix), rownames(names))

# 筛选交集部分的数据
matrix <- matrix[same, ]
names <- names[same, ]
# 添加基因symbol到矩阵
matrix$symbol <- names[, 1]
matrix$type <- names[, 2]

matrix <- matrix[, c(ncol(matrix), 1:(ncol(matrix) - 1))]  # 调整列顺序
matrix <- matrix[, c(ncol(matrix), 1:(ncol(matrix) - 1))]  # 调整列顺序

# 提取蛋白质编码基因与长链非编码RNA基因
counts_matrix <- matrix[matrix$type %in% c("protein_coding"), ] # 需要提取其他类型的基因如lncRNA,可在这里修改
counts_matrix$type <- NULL

#基因symbol去重
counts_matrix <- counts_matrix %>%
  mutate(across(-symbol, ~ suppressWarnings(as.numeric(.x)))) %>%
  group_by(symbol) %>%
  summarize(across(where(is.numeric), max, na.rm = TRUE), .groups = "drop")
#设置行名
counts_matrix <- as.data.frame(counts_matrix)
rownames(counts_matrix) <- counts_matrix$symbol
counts_matrix$symbol <- NULL

#样本去重
colnames(counts_matrix) <- substr(colnames(counts_matrix), 1,16)
counts_matrix <- counts_matrix[,endsWith(colnames(counts_matrix), "A"), drop = FALSE]

# 如果使用TPM数据，建议进行对数化处理
# tpms_matrix <- log2(counts_matrix + 1)
# save(tpms_matrix, file="Data/TCGA_tpms_matrix.Rda")
save(counts_matrix, file="Data/TCGA_counts_matrix.Rda")



library(maftools)
library(dplyr)

# 数据整理
files <- list.files(path = "Data/TCGA/Mutation/gdc_download_20250923_074020.325351",
                    pattern='*.gz',recursive = TRUE)
all_mut <- data.frame()

for(file in files){
  mut<-read.delim(file, skip=7, header =T, fill=TRUE, sep = "\t")
  all_mut <- rbind(all_mut, mut)
}

#保留16位样本条形码用于排序（01A 优先判断）
all_mut$TSB16 <- substr(all_mut$Tumor_Sample_Barcode, 1, 16)

#保留前12位字符（case 级）
all_mut$Tumor_Sample_Barcode = substr(all_mut$Tumor_Sample_Barcode, 1, 12)

all_mut <- read.maf(all_mut)

mut_data <- all_mut@data %>%
  .[, c("Hugo_Symbol", "Variant_Classification", "Tumor_Sample_Barcode", "TSB16")] %>%  # 修改：加入 TSB16 以便排序
  as.data.frame() %>%
  mutate(
    Tumor_Sample_Barcode = substring(.$Tumor_Sample_Barcode, 1, 12),
    vial        = substr(.$TSB16, 16, 16)   #解析 vial（A/B/…）
  ) %>%
  arrange(
    Tumor_Sample_Barcode,
    vial,                 #优先按字母顺序
    TSB16                 #再按16位做稳定排序
  ) %>%
  select(Hugo_Symbol, Variant_Classification, Tumor_Sample_Barcode)  #去掉排序辅助列

# 提取基因
gene <- as.character(unique(mut_data$Hugo_Symbol))
# 提取样本
sample <- as.character(unique(mut_data$Tumor_Sample_Barcode))

# 创建 data.frame
mut_01 <- as.data.frame(
  matrix(0, length(gene), length(sample),
         dimnames = list(gene, sample))
)
# 记录突变类型数据
# mut_type <- as.data.frame(matrix("", length(gene), length(sample),
#                                  dimnames = list(gene, sample)))

# 将信息填入
for (i in 1:nrow(mut_data)){
  mut_01[as.character(mut_data[i, 1]), as.character(mut_data[i, 3])] <- 1
}
# for (i in 1:nrow(mut_data)) {
#   mut_type[as.character(mut_data[i, 1]), as.character(mut_data[i, 3])] <- as.character(mut_data[i, 2])
# }


write.csv(mut_data, "Data/mut_metadata.csv", row.names = FALSE)
write.csv(mut_01, "Data/TCGA_mut_01.csv", row.names = TRUE)
# write.csv(mut_type, "Data/TCGA_mut_type.csv", row.names = TRUE)
# write.mafSummary(all_mut, basename = "all_mut")
save(all_mut, file = "Data/all_mut.Rda")


# maftools瀑布图
oncoplot(
  maf = all_mut,
  top = 30,
  fontSize = 0.6,
  sortByAnnotation = TRUE
)
dev.off()









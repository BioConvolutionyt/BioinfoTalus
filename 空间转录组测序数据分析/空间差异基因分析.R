library(SPARK)
library(Seurat)
load("Data/stRNA.Rda")

# 提取原始计数矩阵和空间坐标
counts <- stRNA@assays[["SCT"]]@counts
coords <- GetTissueCoordinates(stRNA, image = "Spatial", scale = NULL)
# 对齐计数矩阵与坐标
counts <-counts[, rownames(coords)]
info <- data.frame(x = coords$x,
                   y = coords$y
)
rownames(info) <- rownames(coords)

# 创建SPARK对象（自动过滤低质量Spot）
spark <- CreateSPARKObject(
  counts = counts, 
  location = info[, 1:2],  # 只需x,y列
  percentage = 0.1,        # 相当于log(CPM/10k + 1)标准化
  min_total_counts = 10     # 过滤总UMI<10的Spot
)

# 计算并存储总UMI数
spark@lib_size <- apply(spark@counts, 2, sum)

# 估计参数（使用内置的lib_size）
spark <- spark.vc(
  object = spark, 
  covariates = NULL,       # 无额外协变量
  lib_size = spark@lib_size, # 使用对象内已计算的lib_size
  num_core = 5,            # 并行计算
  verbose = FALSE
)

# 计算p值
spark <- spark.test(
  object = spark,
  check_positive = TRUE,   # 检查空间正相关
  verbose = FALSE
)

# 查看结果
head(spark@res_mtest)
save(spark, file="Data/spark.Rda")



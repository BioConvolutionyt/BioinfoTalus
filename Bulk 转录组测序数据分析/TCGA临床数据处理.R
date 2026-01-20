library(readr)        
library(tidyverse)    
library(dplyr)        
library(data.table)   

# 读取临床数据
clin <- fread("Data/TCGA/RNA seq/clinical.cart.2025-09-22/clinical.tsv", data.table = FALSE)  # 使用data.table的fread函数

# 提取感兴趣的列，并删除重复项
# 常规分析中，下面列出的几列是最常用的临床信息，也可以根据研究需要添加其他信息，如吸烟史等。
clin_time <- clin %>%
  dplyr::select(  # 选择临床数据中的相关列
    cases.submitter_id,    # 案例提交者ID
    demographic.vital_status,         # 生存状态（Alive / Dead）
    demographic.days_to_death,        # 从诊断到死亡的天数
    diagnoses.days_to_last_follow_up, # 从诊断到最后一次随访的天数
    demographic.age_at_index,         # 诊断时的年龄
    demographic.gender,               # 性别
    demographic.race,
    diagnoses.ajcc_pathologic_t,    # 病理T分期
    diagnoses.ajcc_pathologic_m,    # 病理M分期
    diagnoses.ajcc_pathologic_n,    # 病理N分期
    diagnoses.ajcc_pathologic_stage # 病理分期
  ) %>%
  dplyr::filter(!duplicated(cases.submitter_id))  # 删除重复的病例ID，只保留每个ID的第一次记录

# 查看数据的前几行，确认数据的提取是否正确
head(clin_time)

# 创建一个新的变量“OS.time”和“OS”：
# - 若病人存活，使用“days_to_last_follow_up”列的数据表示生存时间
# - 若病人死亡，使用“days_to_death”列的数据表示生存时间
# - "OS"列用于表示生存状态，存活为0，死亡为1
clin_merge <- clin_time %>%
  dplyr::mutate(  # 新增两列，用于表示生存时间和生存状态
    OS.time = case_when(
      demographic.vital_status == "Alive" ~ diagnoses.days_to_last_follow_up,  # 存活使用最后一次随访天数
      demographic.vital_status == "Dead" ~ demographic.days_to_death             # 死亡使用死亡天数
    ),
    OS = case_when(
      demographic.vital_status == "Alive" ~ 0,   # 存活用0表示
      demographic.vital_status == "Dead" ~ 1      # 死亡用1表示
    )
  )

# 查看数据的前几行，确认“OS.time”和“OS”列的创建是否正确
head(clin_merge)

# 重命名某些列，使列名更加直观易懂
clin_merge <- clin_merge %>%
  dplyr::rename(
    Age = demographic.age_at_index,         # 将age_at_index重命名为Age
    T = diagnoses.ajcc_pathologic_t,      # 将pathologic_t重命名为T
    M = diagnoses.ajcc_pathologic_m,      # 将pathologic_m重命名为M
    N = diagnoses.ajcc_pathologic_n,      # 将pathologic_n重命名为N
    stage = diagnoses.ajcc_pathologic_stage # 将pathologic_stage重命名为stage
  )



# 删除第2、3、4列（demographic.vital_status, days_to_death, days_to_last_follow_up），并去除缺失值
clin <- clin_merge[,-c(2, 3, 4)]  
clin <- na.omit(clin)          # 删除任何含有缺失值的行

# 将最终处理后的数据保存为一个CSV文件
write.csv(clin, file = "Data/clinical.csv", row.names = FALSE, quote = FALSE)  # 将数据保存到“clinical.csv”文件中，不保留行名，并且不加引号



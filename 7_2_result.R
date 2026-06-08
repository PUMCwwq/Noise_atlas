##################################################################
############ -------- result_2_var_noise_gene -------#############
############ ------      disease vs health       -------##########
############ ------        Monod 结果可视化      -------########## 
##################################################################
library(dplyr)
library(knitr)
library(tidyr)
library(stats)
library(tibble)


################ ✅️✅️✅️✅️✅️✅️✅️✅️✅️✅️✅️
################ ✅️✅️✅️✅️✅️✅️✅️✅️✅️✅️✅️     Extrinsic
################ ✅️✅️✅️✅️✅️✅️✅️✅️✅️✅️✅️
##### （一）读取Extrinsic的结果，计算每组参数的均值
############ 这个数据集只有5个细胞类型--没有otherT !!!
# ['covid1', 'covid2', 'covid3', 'covid4', 'covid5', 'covid6', 'covid7', 'covid8', 'covid9', 'covid10', 'covid11'] 
rm(list=ls())
covid_files <- list.files(path = "/data/wuwq/noise/Monod_FLU/phys", pattern = "^B_covid\\d+_Extrinsic\\.csv$", full.names = TRUE) # ^代表以CD4开头
covid_data <- lapply(covid_files, read.csv, header = TRUE)
# 合并所有数据框为一个长格式数据框
combined_data <- do.call(rbind, covid_data)
# 计算每个基因的三个参数的均值
covid_summary <- combined_data %>%
  group_by(X) %>%
  summarise(
    mean_log10_b = mean(log10.b, na.rm = TRUE),
    se_log10_b = sd(log10.b, na.rm = TRUE) / sqrt(n()),
    mean_log10_beta = mean(log10.beta, na.rm = TRUE),
    se_log10_beta = sd(log10.beta, na.rm = TRUE) / sqrt(n()),
    mean_log10_gamma = mean(log10.gamma, na.rm = TRUE),
    se_log10_gamma = sd(log10.gamma, na.rm = TRUE) / sqrt(n())
  )
write.csv(covid_summary, "/data/wuwq/noise/Monod_FLU/phys/phys_Extrinsic_sum/B_covid.csv", row.names = FALSE)

# ['flu1', 'flu2', 'flu3', 'flu4', 'flu5']
rm(list=ls())
covid_files <- list.files(path = "/data/wuwq/noise/Monod_FLU/phys", pattern = "^B_flu\\d+_Extrinsic\\.csv$", full.names = TRUE) # ^代表以CD4开头
covid_data <- lapply(covid_files, read.csv, header = TRUE)
# 合并所有数据框为一个长格式数据框
combined_data <- do.call(rbind, covid_data)
# 计算每个基因的三个参数的均值
covid_summary <- combined_data %>%
  group_by(X) %>%
  summarise(
    mean_log10_b = mean(log10.b, na.rm = TRUE),
    se_log10_b = sd(log10.b, na.rm = TRUE) / sqrt(n()),
    mean_log10_beta = mean(log10.beta, na.rm = TRUE),
    se_log10_beta = sd(log10.beta, na.rm = TRUE) / sqrt(n()),
    mean_log10_gamma = mean(log10.gamma, na.rm = TRUE),
    se_log10_gamma = sd(log10.gamma, na.rm = TRUE) / sqrt(n())
  )
write.csv(covid_summary, "/data/wuwq/noise/Monod_FLU/phys/phys_Extrinsic_sum/B_flu.csv", row.names = FALSE)

# ['nor1', 'nor2', 'nor3', 'nor4']
rm(list=ls())
covid_files <- list.files(path = "/data/wuwq/noise/Monod_FLU/phys", pattern = "^B_nor\\d+_Extrinsic\\.csv$", full.names = TRUE) # ^代表以CD4开头
covid_data <- lapply(covid_files, read.csv, header = TRUE)
# 合并所有数据框为一个长格式数据框
combined_data <- do.call(rbind, covid_data)
# 计算每个基因的三个参数的均值
covid_summary <- combined_data %>%
  group_by(X) %>%
  summarise(
    mean_log10_b = mean(log10.b, na.rm = TRUE),
    se_log10_b = sd(log10.b, na.rm = TRUE) / sqrt(n()),
    mean_log10_beta = mean(log10.beta, na.rm = TRUE),
    se_log10_beta = sd(log10.beta, na.rm = TRUE) / sqrt(n()),
    mean_log10_gamma = mean(log10.gamma, na.rm = TRUE),
    se_log10_gamma = sd(log10.gamma, na.rm = TRUE) / sqrt(n())
  )
write.csv(covid_summary, "/data/wuwq/noise/Monod_FLU/phys/phys_Extrinsic_sum/B_nor.csv", row.names = FALSE)





##### （二）读取Extrinsic的结果，计算每组参数的均值  5种细胞类型--没有otherT
#########################################################################
# △log10b=(log10 b1)-(log10 b2)  → Diff_log10.b
# △log10 beta/k = (log10 beta1)-(log10 beta2)   Diff_log10.beta
# △log10 gamma/k = (log10 gamma1)-(log10 gamma2)   Diff_log10.gamma
# △log10 k = -1/2 * (△log10 beta/k + △log10 gamma/k)      Diff_log10.k   
# Diff_log10.k = -1/2 * (Diff_log10.beta+ Diff_log10.gamma)
##########################################################################

##### 1.✅️✅️✅️  covid组和 nor组    
### 1.1. log10.b     log10.beta     log10.gamma      先计算检验
rm(list=ls())
# 读取covid组数据
covid_files <- list.files(path = "/data/wuwq/noise/Monod_FLU/phys", pattern = "^otherT_covid\\d+_Extrinsic\\.csv$", full.names = TRUE)
covid_data <- lapply(covid_files, read.csv, header = TRUE)
covid_combined_data <- do.call(rbind, covid_data)
# 读取nor组数据
nor_files <- list.files(path = "/data/wuwq/noise/Monod_FLU/phys", pattern = "^otherT_nor\\d+_Extrinsic\\.csv$", full.names = TRUE)
nor_data <- lapply(nor_files, read.csv, header = TRUE)
nor_combined_data <- do.call(rbind, nor_data)
# 找到covid和nor组的交集基因
common_genes <- intersect(unique(covid_combined_data$X), unique(nor_combined_data$X))
# 准备数据框架存储结果
results <- data.frame(
  Gene = character(),
  Mean_Covid_log10.beta = numeric(),
  Mean_Nor_log10.beta = numeric(),
  Diff_log10.beta = numeric(),
  P_Value = numeric(),
  FDR = numeric(),
  stringsAsFactors = FALSE
)
# 进行秩和检验（Mann-Whitney U检验）
for (gene in common_genes) {
  # 提取covid和nor组对应基因的数据
  covid_gene_data <- covid_combined_data[covid_combined_data$X == gene, "log10.beta"]
  nor_gene_data <- nor_combined_data[nor_combined_data$X == gene, "log10.beta"]
  # 确保每组数据不为空
  if (length(covid_gene_data) > 0 && length(nor_gene_data) > 0) {
    # 计算均值
    mean_covid <- mean(covid_gene_data, na.rm = TRUE)
    mean_nor <- mean(nor_gene_data, na.rm = TRUE)
    diff <- mean_covid - mean_nor
    # 进行秩和检验（Mann-Whitney U检验）
    test_result <- wilcox.test(covid_gene_data, nor_gene_data)
    # 提取p值并进行FDR校正
    p_value <- test_result$p.value
    results <- rbind(results, data.frame(
      Gene = gene,
      Mean_Covid_log10.beta = mean_covid,
      Mean_Nor_log10.beta = mean_nor,
      Diff_log10.beta = diff,   # 此处对于组间计算时用减法，因为Diff_log10.b实际为log10(Covid/Nor)
      P_Value = p_value,
      FDR = p.adjust(p_value, method = "fdr")
    ))
  }
}
results <- results[order(results$FDR), ]
write.csv(results, "/data/wuwq/noise/Monod_FLU/phys/phys_Extrinsic_sum/covid_nor_wilcox/otherT_log10beta.csv", row.names = FALSE)

### 1.2. 汇总 △log10 b, log10 beta/k , △log10 gamma/k, △log10 k   →  即 Diff_log10.b  Diff_log10.beta  Diff_log10.gamma  Diff_log10.k
# △log10b=(log10 b1)-(log10 b2)  → Diff_log10.b
# △log10 beta/k = (log10 beta1)-(log10 beta2)   Diff_log10.beta
# △log10 gamma/k = (log10 gamma1)-(log10 gamma2)   Diff_log10.gamma
# △log10 k = -1/2 * (△log10 beta/k + △log10 gamma/k)      Diff_log10.k   
# Diff_log10.k = -1/2 * (Diff_log10.beta+ Diff_log10.gamma)

rm(list=ls())
# 读取上一步计算出的三个参数组间变化值
log10b <- read.csv("/data/wuwq/noise/Monod_FLU/phys/phys_Extrinsic_sum/covid_nor_wilcox/B_log10b.csv", header = TRUE)
log10beta <- read.csv("/data/wuwq/noise/Monod_FLU/phys/phys_Extrinsic_sum/covid_nor_wilcox/B_log10beta.csv", header = TRUE)
log10gamma <- read.csv("/data/wuwq/noise/Monod_FLU/phys/phys_Extrinsic_sum/covid_nor_wilcox/B_log10gamma.csv", header = TRUE)
# 合并log10b和log10beta表，基于Gene列
merged_data <- merge(
  log10b[, c("Gene", "Diff_log10.b")], 
  log10beta[, c("Gene", "Diff_log10.beta")], 
  by = "Gene"
)
# 再将log10gamma与前面合并结果结合，基于Gene列
covid_nor_phys <- merge(
  merged_data,
  log10gamma[, c("Gene", "Diff_log10.gamma")],
  by = "Gene"
)
# 新增一列Diff_log10.k
covid_nor_phys$Diff_log10.k <- -1/2 * (covid_nor_phys$Diff_log10.beta + covid_nor_phys$Diff_log10.gamma)
write.csv(covid_nor_phys, "/data/wuwq/noise/Monod_FLU/phys/phys_Extrinsic_sum/covid_nor_wilcox/allphys_B.csv", row.names = FALSE)



##### 2.✅️✅️✅         flu组 和 nor组
### 2.1. log10.b     log10.beta     log10.gamma      先计算检验
rm(list=ls())
# 读取flu组数据
flu_files <- list.files(path = "/data/wuwq/noise/Monod_FLU/phys", pattern = "^B_flu\\d+_Extrinsic\\.csv$", full.names = TRUE)
flu_data <- lapply(flu_files, read.csv, header = TRUE)
flu_combined_data <- do.call(rbind, flu_data)
# 读取nor组数据
nor_files <- list.files(path = "/data/wuwq/noise/Monod_FLU/phys", pattern = "^B_nor\\d+_Extrinsic\\.csv$", full.names = TRUE)
nor_data <- lapply(nor_files, read.csv, header = TRUE)
nor_combined_data <- do.call(rbind, nor_data)
# 找到covid和nor组的交集基因
common_genes <- intersect(unique(flu_combined_data$X), unique(nor_combined_data$X))
# 准备数据框架存储结果
results <- data.frame(
  Gene = character(),
  Mean_Flu_log10.gamma = numeric(),
  Mean_Nor_log10.gamma = numeric(),
  Diff_log10.gamma = numeric(),
  P_Value = numeric(),
  FDR = numeric(),
  stringsAsFactors = FALSE
)
# 进行秩和检验（Mann-Whitney U检验）
for (gene in common_genes) {
  # 提取covid和nor组对应基因的数据
  flu_gene_data <- flu_combined_data[flu_combined_data$X == gene, "log10.gamma"]
  nor_gene_data <- nor_combined_data[nor_combined_data$X == gene, "log10.gamma"]
  # 确保每组数据不为空
  if (length(flu_gene_data) > 0 && length(nor_gene_data) > 0) {
    # 计算均值
    mean_flu <- mean(flu_gene_data, na.rm = TRUE)
    mean_nor <- mean(nor_gene_data, na.rm = TRUE)
    diff <- mean_flu - mean_nor
    # 进行秩和检验（Mann-Whitney U检验）
    test_result <- wilcox.test(flu_gene_data, nor_gene_data)
    # 提取p值并进行FDR校正
    p_value <- test_result$p.value
    results <- rbind(results, data.frame(
      Gene = gene,
      Mean_Flu_log10.gamma = mean_flu,
      Mean_Nor_log10.gamma = mean_nor,
      Diff_log10.gamma = diff,   # 此处对于组间计算时用减法，因为Diff_log10.b实际为log10(Covid/Nor)
      P_Value = p_value,
      FDR = p.adjust(p_value, method = "fdr")
    ))
  }
}
results <- results[order(results$FDR), ]
write.csv(results, "/data/wuwq/noise/Monod_FLU/phys/phys_Extrinsic_sum/flu_nor_wilcox/B_log10gamma.csv", row.names = FALSE)


### 2.2. 汇总 △log10 b, log10 beta/k , △log10 gamma/k, △log10 k   →  即 Diff_log10.b  Diff_log10.beta  Diff_log10.gamma  Diff_log10.k
# △log10b=(log10 b1)-(log10 b2)  → Diff_log10.b
# △log10 beta/k = (log10 beta1)-(log10 beta2)   Diff_log10.beta
# △log10 gamma/k = (log10 gamma1)-(log10 gamma2)   Diff_log10.gamma
# △log10 k = -1/2 * (△log10 beta/k + △log10 gamma/k)      Diff_log10.k   
# Diff_log10.k = -1/2 * (Diff_log10.beta+ Diff_log10.gamma)

rm(list=ls())
# 读取上一步计算出的三个参数组间变化值
log10b <- read.csv("/data/wuwq/noise/Monod_FLU/phys/phys_Extrinsic_sum/flu_nor_wilcox/B_log10b.csv", header = TRUE)
log10beta <- read.csv("/data/wuwq/noise/Monod_FLU/phys/phys_Extrinsic_sum/flu_nor_wilcox/B_log10beta.csv", header = TRUE)
log10gamma <- read.csv("/data/wuwq/noise/Monod_FLU/phys/phys_Extrinsic_sum/flu_nor_wilcox/B_log10gamma.csv", header = TRUE)
# 合并log10b和log10beta表，基于Gene列
merged_data <- merge(
  log10b[, c("Gene", "Diff_log10.b")], 
  log10beta[, c("Gene", "Diff_log10.beta")], 
  by = "Gene"
)
# 再将log10gamma与前面合并结果结合，基于Gene列
flu_nor_phys <- merge(
  merged_data,
  log10gamma[, c("Gene", "Diff_log10.gamma")],
  by = "Gene"
)
# 新增一列Diff_log10.k
flu_nor_phys$Diff_log10.k <- -1/2 * (flu_nor_phys$Diff_log10.beta + flu_nor_phys$Diff_log10.gamma)
write.csv(flu_nor_phys, "/data/wuwq/noise/Monod_FLU/phys/phys_Extrinsic_sum/flu_nor_wilcox/allphys_B.csv", row.names = FALSE)



##### 3.✅️✅️✅️  covid组和 flu组
### 3.1. log10.b     log10.beta     log10.gamma      先计算检验
rm(list=ls())
# 读取covid组数据
covid_files <- list.files(path = "/data/wuwq/noise/Monod_FLU/phys", pattern = "^B_covid\\d+_Extrinsic\\.csv$", full.names = TRUE)
covid_data <- lapply(covid_files, read.csv, header = TRUE)
covid_combined_data <- do.call(rbind, covid_data)
# 读取flu组数据
flu_files <- list.files(path = "/data/wuwq/noise/Monod_FLU/phys", pattern = "^B_flu\\d+_Extrinsic\\.csv$", full.names = TRUE)
flu_data <- lapply(flu_files, read.csv, header = TRUE)
flu_combined_data <- do.call(rbind, flu_data)
# 找到covid和nor组的交集基因
common_genes <- intersect(unique(covid_combined_data$X), unique(flu_combined_data$X))
# 准备数据框架存储结果
results <- data.frame(
  Gene = character(),
  Mean_Covid_log10.gamma = numeric(),
  Mean_Flu_log10.gamma = numeric(),
  Diff_log10.gamma = numeric(),
  P_Value = numeric(),
  FDR = numeric(),
  stringsAsFactors = FALSE
)
# 进行秩和检验（Mann-Whitney U检验）
for (gene in common_genes) {
  # 提取covid和nor组对应基因的数据
  covid_gene_data <- covid_combined_data[covid_combined_data$X == gene, "log10.gamma"]
  flu_gene_data <- flu_combined_data[flu_combined_data$X == gene, "log10.gamma"]
  # 确保每组数据不为空
  if (length(covid_gene_data) > 0 && length(flu_gene_data) > 0) {
    # 计算均值
    mean_covid <- mean(covid_gene_data, na.rm = TRUE)
    mean_flu <- mean(flu_gene_data, na.rm = TRUE)
    diff <- mean_covid - mean_flu
    # 进行秩和检验（Mann-Whitney U检验）
    test_result <- wilcox.test(covid_gene_data, flu_gene_data)
    # 提取p值并进行FDR校正
    p_value <- test_result$p.value
    results <- rbind(results, data.frame(
      Gene = gene,
      Mean_Covid_log10.gamma = mean_covid,
      Mean_Flu_log10.gamma = mean_flu,
      Diff_log10.gamma = diff,   # 此处对于组间计算时用减法，因为Diff_log10.b实际为log10(Covid/Nor)
      P_Value = p_value,
      FDR = p.adjust(p_value, method = "fdr")
    ))
  }
}
results <- results[order(results$FDR), ]
write.csv(results, "/data/wuwq/noise/Monod_FLU/phys/phys_Extrinsic_sum/covid_flu_wilcox/B_log10gamma.csv", row.names = FALSE)


### 3.2. 汇总 △log10 b, log10 beta/k , △log10 gamma/k, △log10 k   →  即 Diff_log10.b  Diff_log10.beta  Diff_log10.gamma  Diff_log10.k
# △log10b=(log10 b1)-(log10 b2)  → Diff_log10.b
# △log10 beta/k = (log10 beta1)-(log10 beta2)   Diff_log10.beta
# △log10 gamma/k = (log10 gamma1)-(log10 gamma2)   Diff_log10.gamma
# △log10 k = -1/2 * (△log10 beta/k + △log10 gamma/k)      Diff_log10.k   
# Diff_log10.k = -1/2 * (Diff_log10.beta+ Diff_log10.gamma)

rm(list=ls())
# 读取上一步计算出的三个参数组间变化值
log10b <- read.csv("/data/wuwq/noise/Monod_FLU/phys/phys_Extrinsic_sum/covid_flu_wilcox/B_log10b.csv", header = TRUE)
log10beta <- read.csv("/data/wuwq/noise/Monod_FLU/phys/phys_Extrinsic_sum/covid_flu_wilcox/B_log10beta.csv", header = TRUE)
log10gamma <- read.csv("/data/wuwq/noise/Monod_FLU/phys/phys_Extrinsic_sum/covid_flu_wilcox/B_log10gamma.csv", header = TRUE)
# 合并log10b和log10beta表，基于Gene列
merged_data <- merge(
  log10b[, c("Gene", "Diff_log10.b")], 
  log10beta[, c("Gene", "Diff_log10.beta")], 
  by = "Gene"
)
# 再将log10gamma与前面合并结果结合，基于Gene列
covid_flu_phys <- merge(
  merged_data,
  log10gamma[, c("Gene", "Diff_log10.gamma")],
  by = "Gene"
)
# 新增一列Diff_log10.k
covid_flu_phys$Diff_log10.k <- -1/2 * (covid_flu_phys$Diff_log10.beta + covid_flu_phys$Diff_log10.gamma)
write.csv(covid_flu_phys, "/data/wuwq/noise/Monod_FLU/phys/phys_Extrinsic_sum/covid_flu_wilcox/allphys_B.csv", row.names = FALSE)











################ ✅️✅️✅️✅️✅️✅️✅️✅️✅️✅️✅️
################ ✅️✅️✅️✅️✅️✅️✅️✅️✅️✅️✅️     Bursty
################ ✅️✅️✅️✅️✅️✅️✅️✅️✅️✅️✅️
##### （一）读取Bursty的结果，计算每组参数的均值
############ 这个数据集只有5个细胞类型--没有otherT !!!
# ['covid1', 'covid2', 'covid3', 'covid4', 'covid5', 'covid6', 'covid7', 'covid8', 'covid9', 'covid10', 'covid11'] 
rm(list=ls())
covid_files <- list.files(path = "/data/wuwq/noise/Monod_FLU/phys", pattern = "^B_covid\\d+_Bursty\\.csv$", full.names = TRUE) # ^代表以CD4开头
covid_data <- lapply(covid_files, read.csv, header = TRUE)
# 合并所有数据框为一个长格式数据框
combined_data <- do.call(rbind, covid_data)
# 计算每个基因的三个参数的均值
covid_mean <- combined_data %>%
  group_by(X) %>%
  summarise(
    mean_log10_b = mean(log10.b, na.rm = TRUE),
    mean_log10_beta = mean(log10.beta, na.rm = TRUE),
    mean_log10_gamma = mean(log10.gamma, na.rm = TRUE)
  )
write.csv(covid_mean, "/data/wuwq/noise/Monod_FLU/phys/phys_Bursty_sum/B_covid.csv", row.names = FALSE)

# ['flu1', 'flu2', 'flu3', 'flu4', 'flu5']
rm(list=ls())
covid_files <- list.files(path = "/data/wuwq/noise/Monod_FLU/phys", pattern = "^B_flu\\d+_Bursty\\.csv$", full.names = TRUE) # ^代表以CD4开头
covid_data <- lapply(covid_files, read.csv, header = TRUE)
# 合并所有数据框为一个长格式数据框
combined_data <- do.call(rbind, covid_data)
# 计算每个基因的三个参数的均值
covid_mean <- combined_data %>%
  group_by(X) %>%
  summarise(
    mean_log10_b = mean(log10.b, na.rm = TRUE),
    mean_log10_beta = mean(log10.beta, na.rm = TRUE),
    mean_log10_gamma = mean(log10.gamma, na.rm = TRUE)
  )
write.csv(covid_mean, "/data/wuwq/noise/Monod_FLU/phys/phys_Bursty_sum/B_flu.csv", row.names = FALSE)

# ['nor1', 'nor2', 'nor3', 'nor4']
rm(list=ls())
covid_files <- list.files(path = "/data/wuwq/noise/Monod_FLU/phys", pattern = "^B_nor\\d+_Bursty\\.csv$", full.names = TRUE) # ^代表以CD4开头
covid_data <- lapply(covid_files, read.csv, header = TRUE)
# 合并所有数据框为一个长格式数据框
combined_data <- do.call(rbind, covid_data)
# 计算每个基因的三个参数的均值
covid_mean <- combined_data %>%
  group_by(X) %>%
  summarise(
    mean_log10_b = mean(log10.b, na.rm = TRUE),
    mean_log10_beta = mean(log10.beta, na.rm = TRUE),
    mean_log10_gamma = mean(log10.gamma, na.rm = TRUE)
  )
write.csv(covid_mean, "/data/wuwq/noise/Monod_FLU/phys/phys_Bursty_sum/B_nor.csv", row.names = FALSE)




##### （二）读取Bursty的结果，计算每组参数的均值  5种细胞类型--没有otherT
#########################################################################
# △log10b=(log10 b1)-(log10 b2)  → Diff_log10.b
# △log10 beta/k = (log10 beta1)-(log10 beta2)   Diff_log10.beta
# △log10 gamma/k = (log10 gamma1)-(log10 gamma2)   Diff_log10.gamma
# △log10 k = -1/2 * (△log10 beta/k + △log10 gamma/k)      Diff_log10.k   
# Diff_log10.k = -1/2 * (Diff_log10.beta+ Diff_log10.gamma)
##########################################################################

##### 1.✅️✅️✅️  covid组和 nor组    
### 1.1. log10.b     log10.beta     log10.gamma      先计算检验
rm(list=ls())
# 读取covid组数据
covid_files <- list.files(path = "/data/wuwq/noise/Monod_FLU/phys", pattern = "^B_covid\\d+_Bursty\\.csv$", full.names = TRUE)
covid_data <- lapply(covid_files, read.csv, header = TRUE)
covid_combined_data <- do.call(rbind, covid_data)
# 读取nor组数据
nor_files <- list.files(path = "/data/wuwq/noise/Monod_FLU/phys", pattern = "^B_nor\\d+_Bursty\\.csv$", full.names = TRUE)
nor_data <- lapply(nor_files, read.csv, header = TRUE)
nor_combined_data <- do.call(rbind, nor_data)
# 找到covid和nor组的交集基因
common_genes <- intersect(unique(covid_combined_data$X), unique(nor_combined_data$X))
# 准备数据框架存储结果
results <- data.frame(
  Gene = character(),
  Mean_Covid_log10.gamma = numeric(),
  Mean_Nor_log10.gamma = numeric(),
  Diff_log10.gamma = numeric(),
  P_Value = numeric(),
  FDR = numeric(),
  stringsAsFactors = FALSE
)
# 进行秩和检验（Mann-Whitney U检验）
for (gene in common_genes) {
  # 提取covid和nor组对应基因的数据
  covid_gene_data <- covid_combined_data[covid_combined_data$X == gene, "log10.gamma"]
  nor_gene_data <- nor_combined_data[nor_combined_data$X == gene, "log10.gamma"]
  # 确保每组数据不为空
  if (length(covid_gene_data) > 0 && length(nor_gene_data) > 0) {
    # 计算均值
    mean_covid <- mean(covid_gene_data, na.rm = TRUE)
    mean_nor <- mean(nor_gene_data, na.rm = TRUE)
    diff <- mean_covid - mean_nor
    # 进行秩和检验（Mann-Whitney U检验）
    test_result <- wilcox.test(covid_gene_data, nor_gene_data)
    # 提取p值并进行FDR校正
    p_value <- test_result$p.value
    results <- rbind(results, data.frame(
      Gene = gene,
      Mean_Covid_log10.gamma = mean_covid,
      Mean_Nor_log10.gamma = mean_nor,
      Diff_log10.gamma = diff,   # 此处对于组间计算时用减法，因为Diff_log10.b实际为log10(Covid/Nor)
      P_Value = p_value,
      FDR = p.adjust(p_value, method = "fdr")
    ))
  }
}
results <- results[order(results$FDR), ]
write.csv(results, "/data/wuwq/noise/Monod_FLU/phys/phys_Bursty_sum/covid_nor_wilcox/B_log10gamma.csv", row.names = FALSE)


### 1.2. 汇总 △log10 b, log10 beta/k , △log10 gamma/k, △log10 k   →  即 Diff_log10.b  Diff_log10.beta  Diff_log10.gamma  Diff_log10.k
# △log10b=(log10 b1)-(log10 b2)  → Diff_log10.b
# △log10 beta/k = (log10 beta1)-(log10 beta2)   Diff_log10.beta
# △log10 gamma/k = (log10 gamma1)-(log10 gamma2)   Diff_log10.gamma
# △log10 k = -1/2 * (△log10 beta/k + △log10 gamma/k)      Diff_log10.k   
# Diff_log10.k = -1/2 * (Diff_log10.beta+ Diff_log10.gamma)

rm(list=ls())
# 读取上一步计算出的三个参数组间变化值
log10b <- read.csv("/data/wuwq/noise/Monod_FLU/phys/phys_Bursty_sum/covid_nor_wilcox/B_log10b.csv", header = TRUE)
log10beta <- read.csv("/data/wuwq/noise/Monod_FLU/phys/phys_Bursty_sum/covid_nor_wilcox/B_log10beta.csv", header = TRUE)
log10gamma <- read.csv("/data/wuwq/noise/Monod_FLU/phys/phys_Bursty_sum/covid_nor_wilcox/B_log10gamma.csv", header = TRUE)
# 合并log10b和log10beta表，基于Gene列
merged_data <- merge(
  log10b[, c("Gene", "Diff_log10.b")], 
  log10beta[, c("Gene", "Diff_log10.beta")], 
  by = "Gene"
)
# 再将log10gamma与前面合并结果结合，基于Gene列
covid_nor_phys <- merge(
  merged_data,
  log10gamma[, c("Gene", "Diff_log10.gamma")],
  by = "Gene"
)
# 新增一列Diff_log10.k
covid_nor_phys$Diff_log10.k <- -1/2 * (covid_nor_phys$Diff_log10.beta + covid_nor_phys$Diff_log10.gamma)
write.csv(covid_nor_phys, "/data/wuwq/noise/Monod_FLU/phys/phys_Bursty_sum/covid_nor_wilcox/allphys_B.csv", row.names = FALSE)



##### 2.✅️✅️✅         flu组 和 nor组
### 2.1. log10.b     log10.beta     log10.gamma      先计算检验
rm(list=ls())
# 读取flu组数据
flu_files <- list.files(path = "/data/wuwq/noise/Monod_FLU/phys", pattern = "^B_flu\\d+_Bursty\\.csv$", full.names = TRUE)
flu_data <- lapply(flu_files, read.csv, header = TRUE)
flu_combined_data <- do.call(rbind, flu_data)
# 读取nor组数据
nor_files <- list.files(path = "/data/wuwq/noise/Monod_FLU/phys", pattern = "^B_nor\\d+_Bursty\\.csv$", full.names = TRUE)
nor_data <- lapply(nor_files, read.csv, header = TRUE)
nor_combined_data <- do.call(rbind, nor_data)
# 找到covid和nor组的交集基因
common_genes <- intersect(unique(flu_combined_data$X), unique(nor_combined_data$X))
# 准备数据框架存储结果
results <- data.frame(
  Gene = character(),
  Mean_Flu_log10.gamma = numeric(),
  Mean_Nor_log10.gamma = numeric(),
  Diff_log10.gamma = numeric(),
  P_Value = numeric(),
  FDR = numeric(),
  stringsAsFactors = FALSE
)
# 进行秩和检验（Mann-Whitney U检验）
for (gene in common_genes) {
  # 提取covid和nor组对应基因的数据
  flu_gene_data <- flu_combined_data[flu_combined_data$X == gene, "log10.gamma"]
  nor_gene_data <- nor_combined_data[nor_combined_data$X == gene, "log10.gamma"]
  # 确保每组数据不为空
  if (length(flu_gene_data) > 0 && length(nor_gene_data) > 0) {
    # 计算均值
    mean_flu <- mean(flu_gene_data, na.rm = TRUE)
    mean_nor <- mean(nor_gene_data, na.rm = TRUE)
    diff <- mean_flu - mean_nor
    # 进行秩和检验（Mann-Whitney U检验）
    test_result <- wilcox.test(flu_gene_data, nor_gene_data)
    # 提取p值并进行FDR校正
    p_value <- test_result$p.value
    results <- rbind(results, data.frame(
      Gene = gene,
      Mean_Flu_log10.gamma = mean_flu,
      Mean_Nor_log10.gamma = mean_nor,
      Diff_log10.gamma = diff,   # 此处对于组间计算时用减法，因为Diff_log10.b实际为log10(Covid/Nor)
      P_Value = p_value,
      FDR = p.adjust(p_value, method = "fdr")
    ))
  }
}
results <- results[order(results$FDR), ]
write.csv(results, "/data/wuwq/noise/Monod_FLU/phys/phys_Bursty_sum/flu_nor_wilcox/B_log10gamma.csv", row.names = FALSE)


### 2.2. 汇总 △log10 b, log10 beta/k , △log10 gamma/k, △log10 k   →  即 Diff_log10.b  Diff_log10.beta  Diff_log10.gamma  Diff_log10.k
# △log10b=(log10 b1)-(log10 b2)  → Diff_log10.b
# △log10 beta/k = (log10 beta1)-(log10 beta2)   Diff_log10.beta
# △log10 gamma/k = (log10 gamma1)-(log10 gamma2)   Diff_log10.gamma
# △log10 k = -1/2 * (△log10 beta/k + △log10 gamma/k)      Diff_log10.k   
# Diff_log10.k = -1/2 * (Diff_log10.beta+ Diff_log10.gamma)

rm(list=ls())
# 读取上一步计算出的三个参数组间变化值
log10b <- read.csv("/data/wuwq/noise/Monod_FLU/phys/phys_Bursty_sum/flu_nor_wilcox/B_log10b.csv", header = TRUE)
log10beta <- read.csv("/data/wuwq/noise/Monod_FLU/phys/phys_Bursty_sum/flu_nor_wilcox/B_log10beta.csv", header = TRUE)
log10gamma <- read.csv("/data/wuwq/noise/Monod_FLU/phys/phys_Bursty_sum/flu_nor_wilcox/B_log10gamma.csv", header = TRUE)
# 合并log10b和log10beta表，基于Gene列
merged_data <- merge(
  log10b[, c("Gene", "Diff_log10.b")], 
  log10beta[, c("Gene", "Diff_log10.beta")], 
  by = "Gene"
)
# 再将log10gamma与前面合并结果结合，基于Gene列
flu_nor_phys <- merge(
  merged_data,
  log10gamma[, c("Gene", "Diff_log10.gamma")],
  by = "Gene"
)
# 新增一列Diff_log10.k
flu_nor_phys$Diff_log10.k <- -1/2 * (flu_nor_phys$Diff_log10.beta + flu_nor_phys$Diff_log10.gamma)
write.csv(flu_nor_phys, "/data/wuwq/noise/Monod_FLU/phys/phys_Bursty_sum/flu_nor_wilcox/allphys_B.csv", row.names = FALSE)



##### 3.✅️✅️✅️  covid组和 flu组
### 3.1. log10.b     log10.beta     log10.gamma      先计算检验
rm(list=ls())
# 读取covid组数据
covid_files <- list.files(path = "/data/wuwq/noise/Monod_FLU/phys", pattern = "^B_covid\\d+_Bursty\\.csv$", full.names = TRUE)
covid_data <- lapply(covid_files, read.csv, header = TRUE)
covid_combined_data <- do.call(rbind, covid_data)
# 读取flu组数据
flu_files <- list.files(path = "/data/wuwq/noise/Monod_FLU/phys", pattern = "^B_flu\\d+_Bursty\\.csv$", full.names = TRUE)
flu_data <- lapply(flu_files, read.csv, header = TRUE)
flu_combined_data <- do.call(rbind, flu_data)
# 找到covid和nor组的交集基因
common_genes <- intersect(unique(covid_combined_data$X), unique(flu_combined_data$X))
# 准备数据框架存储结果
results <- data.frame(
  Gene = character(),
  Mean_Covid_log10.gamma = numeric(),
  Mean_Flu_log10.gamma = numeric(),
  Diff_log10.gamma = numeric(),
  P_Value = numeric(),
  FDR = numeric(),
  stringsAsFactors = FALSE
)
# 进行秩和检验（Mann-Whitney U检验）
for (gene in common_genes) {
  # 提取covid和nor组对应基因的数据
  covid_gene_data <- covid_combined_data[covid_combined_data$X == gene, "log10.gamma"]
  flu_gene_data <- flu_combined_data[flu_combined_data$X == gene, "log10.gamma"]
  # 确保每组数据不为空
  if (length(covid_gene_data) > 0 && length(flu_gene_data) > 0) {
    # 计算均值
    mean_covid <- mean(covid_gene_data, na.rm = TRUE)
    mean_flu <- mean(flu_gene_data, na.rm = TRUE)
    diff <- mean_covid - mean_flu
    # 进行秩和检验（Mann-Whitney U检验）
    test_result <- wilcox.test(covid_gene_data, flu_gene_data)
    # 提取p值并进行FDR校正
    p_value <- test_result$p.value
    results <- rbind(results, data.frame(
      Gene = gene,
      Mean_Covid_log10.gamma = mean_covid,
      Mean_Flu_log10.gamma = mean_flu,
      Diff_log10.gamma = diff,   # 此处对于组间计算时用减法，因为Diff_log10.b实际为log10(Covid/Nor)
      P_Value = p_value,
      FDR = p.adjust(p_value, method = "fdr")
    ))
  }
}
results <- results[order(results$FDR), ]
write.csv(results, "/data/wuwq/noise/Monod_FLU/phys/phys_Bursty_sum/covid_flu_wilcox/B_log10gamma.csv", row.names = FALSE)


### 3.2. 汇总 △log10 b, log10 beta/k , △log10 gamma/k, △log10 k   →  即 Diff_log10.b  Diff_log10.beta  Diff_log10.gamma  Diff_log10.k
# △log10b=(log10 b1)-(log10 b2)  → Diff_log10.b
# △log10 beta/k = (log10 beta1)-(log10 beta2)   Diff_log10.beta
# △log10 gamma/k = (log10 gamma1)-(log10 gamma2)   Diff_log10.gamma
# △log10 k = -1/2 * (△log10 beta/k + △log10 gamma/k)      Diff_log10.k   
# Diff_log10.k = -1/2 * (Diff_log10.beta+ Diff_log10.gamma)

rm(list=ls())
# 读取上一步计算出的三个参数组间变化值
log10b <- read.csv("/data/wuwq/noise/Monod_FLU/phys/phys_Bursty_sum/covid_flu_wilcox/B_log10b.csv", header = TRUE)
log10beta <- read.csv("/data/wuwq/noise/Monod_FLU/phys/phys_Bursty_sum/covid_flu_wilcox/B_log10beta.csv", header = TRUE)
log10gamma <- read.csv("/data/wuwq/noise/Monod_FLU/phys/phys_Bursty_sum/covid_flu_wilcox/B_log10gamma.csv", header = TRUE)
# 合并log10b和log10beta表，基于Gene列
merged_data <- merge(
  log10b[, c("Gene", "Diff_log10.b")], 
  log10beta[, c("Gene", "Diff_log10.beta")], 
  by = "Gene"
)
# 再将log10gamma与前面合并结果结合，基于Gene列
covid_flu_phys <- merge(
  merged_data,
  log10gamma[, c("Gene", "Diff_log10.gamma")],
  by = "Gene"
)
# 新增一列Diff_log10.k
covid_flu_phys$Diff_log10.k <- -1/2 * (covid_flu_phys$Diff_log10.beta + covid_flu_phys$Diff_log10.gamma)
write.csv(covid_flu_phys, "/data/wuwq/noise/Monod_FLU/phys/phys_Bursty_sum/covid_flu_wilcox/allphys_B.csv", row.names = FALSE)




rm(list=ls())
###########################################################
################################################### 坐标系转换的概率密度函数
###########################################################
### covid_nor
raw_axis_covid <- read.delim("/data/wuwq/noise/Monod_FLU/phys/phys_Extrinsic_sum/covid_nor_sig_rawaxis.txt")
A_covid <- raw_axis_covid$Diff_log10.b
B_covid <- raw_axis_covid$Diff_log10.k
# 旋转矩阵
rotation_matrix_covid <- matrix(c(1/sqrt(2), -1/sqrt(2), 1/sqrt(2), 1/sqrt(2)), nrow=2)
# 原坐标点转换为矩阵形式
original_coordinates_covid <- cbind(A_covid, B_covid)
# 使用旋转矩阵转换坐标
new_coordinates_covid <- original_coordinates_covid %*% rotation_matrix_covid
# 提取新坐标系中的 x 坐标 (对应 y = -x 的新轴)
new_x_covid <- new_coordinates_covid[, 1]
# 使用 ggplot2 来绘制密度分布函数
library(ggplot2)
# 绘制密度分布
density_plot_covid <- ggplot(data.frame(new_x_covid), aes(x = new_x_covid)) +
  geom_density(fill = "blue", alpha = 0.5) +
  labs(title = "Density Distribution along New x-axis",
       x = "Transformed x-coordinate", y = "Density") +
  theme_minimal()
# 绘图
print(density_plot_covid)
# 使用 ggsave 来保存为 SVG 格式
ggsave("/data/wuwq/noise/Monod_FLU/phys/phys_Extrinsic_sum/density_plot_covid.svg", plot = density_plot_covid, width = 8, height = 6, units = "in")


### flu_nor
raw_axis_flu <- read.delim("/data/wuwq/noise/Monod_FLU/phys/phys_Extrinsic_sum/flu_nor_sig_rawaxis.txt")
A_flu <- raw_axis_flu$Diff_log10.b
B_flu <- raw_axis_flu$Diff_log10.k
# 旋转矩阵
rotation_matrix_flu <- matrix(c(1/sqrt(2), -1/sqrt(2), 1/sqrt(2), 1/sqrt(2)), nrow=2)
# 原坐标点转换为矩阵形式
original_coordinates_flu <- cbind(A_flu, B_flu)
# 使用旋转矩阵转换坐标
new_coordinates_flu <- original_coordinates_flu %*% rotation_matrix_flu
# 提取新坐标系中的 x 坐标 (对应 y = -x 的新轴)
new_x_flu <- new_coordinates_flu[, 1]
# 使用 ggplot2 来绘制密度分布函数
library(ggplot2)
# 绘制密度分布
density_plot_flu <- ggplot(data.frame(new_x_flu), aes(x = new_x_flu)) +
  geom_density(fill = "pink", alpha = 0.5) +
  labs(title = "Density Distribution along New x-axis",
       x = "Transformed x-coordinate", y = "Density") +
  theme_minimal()
# 绘图
print(density_plot_flu)
# 使用 ggsave 来保存为 SVG 格式
ggsave("/data/wuwq/noise/Monod_FLU/phys/phys_Extrinsic_sum/density_plot_flu.svg", plot = density_plot_flu, width = 8, height = 6, units = "in")









##################################################################
############ ------------ result_3_COMBAT -----------#############
############ ------      感染性疾病预后应用       -------#########
############ ------     COMBAT数据集挖掘       -------############ 
##################################################################
library(dplyr)
library(knitr)
library(tidyr)
library(stats)
library(tibble)
library(ggplot2)
library(cowplot)
library(caret)
library(pROC)
install.packages("randomForest")
library(randomForest)


rm(list=ls())
############################################## （一）整合疾病二级基线信息--COMBAT---❗️因基因名为列名❗️️，故 ✅️特殊读取和 ✅️特殊保存
noise_disease <- read.table("/data/wuwq/noise/DISEASE_ALL/merge_disease_health_otherT.csv", sep = ",", row.names=1, header = TRUE, check.names = FALSE, quote = "\"")
noise_disease <- subset(noise_disease, dataset == "COMBAT")
#### 更新二级结局--COMBAT 无误
#-- donor_id已经无重复，因此用sample_id没有问题，COMBAT_baseline中的10个数据含健康个体，虽然冗余但是根据id过滤疾病个体后并不影响结果
baseline <- read.delim("/data/wuwq/noise/DISEASE_ALL/baseline_COMBAT.txt", header = TRUE, row.names = 1)
rownames(baseline) <- baseline$scRNASeq_sample_ID
common_rows <- intersect(rownames(baseline), rownames(noise_disease))
noise_d <- data.frame(row.names = common_rows)
noise_d <- cbind(noise_d, baseline[common_rows, , drop = FALSE],  noise_disease[common_rows, , drop = FALSE] )
### 去除全部为NA的基因列
noise_d <- Filter(function(x) !all(is.na(x)), noise_d)
unique(noise_d$Disease_group) # "COVID"  "SEPSIS" "FLU"
unique(noise_d$Outcome) # 0  1  NA
write.table(noise_d,file = "/data/wuwq/noise/DISEASE_ALL/app_disease/otherT_COMBAT.csv", sep = ",", row.names = TRUE, col.names = NA, quote = TRUE)


rm(list=ls())
### merge细胞类型   特殊读取 ✅️
noise_disease <- read.table("/data/wuwq/noise/DISEASE_ALL/merge_disease_health_CD4.csv", sep = ",", row.names=1, header = TRUE, check.names = FALSE, quote = "\"")
noise_CD4 <- subset(noise_disease, dataset == "COMBAT")
noise_disease <- read.table("/data/wuwq/noise/DISEASE_ALL/merge_disease_health_CD8.csv", sep = ",", row.names=1, header = TRUE, check.names = FALSE, quote = "\"")
noise_CD8 <- subset(noise_disease, dataset == "COMBAT")
noise_disease <- read.table("/data/wuwq/noise/DISEASE_ALL/merge_disease_health_MNP.csv", sep = ",", row.names=1, header = TRUE, check.names = FALSE, quote = "\"")
noise_MNP <- subset(noise_disease, dataset == "COMBAT")
noise_disease <- read.table("/data/wuwq/noise/DISEASE_ALL/merge_disease_health_NK.csv", sep = ",", row.names=1, header = TRUE, check.names = FALSE, quote = "\"")
noise_NK <- subset(noise_disease, dataset == "COMBAT")
noise_disease <- read.table("/data/wuwq/noise/DISEASE_ALL/merge_disease_health_B.csv", sep = ",", row.names=1, header = TRUE, check.names = FALSE, quote = "\"")
noise_B <- subset(noise_disease, dataset == "COMBAT")
noise_disease <- read.table("/data/wuwq/noise/DISEASE_ALL/merge_disease_health_otherT.csv", sep = ",", row.names=1, header = TRUE, check.names = FALSE, quote = "\"")
noise_otherT <- subset(noise_disease, dataset == "COMBAT")
rm(noise_disease)
# 去掉前5列自带基线
noise_CD4 <- noise_CD4[,-c(1:5)]
noise_CD8 <- noise_CD8[,-c(1:5)]
noise_MNP <- noise_MNP[,-c(1:5)]
noise_NK <- noise_NK[,-c(1:5)]
noise_B <- noise_B[,-c(1:5)]
noise_otherT <- noise_otherT[,-c(1:5)]
# 读取并处理每个数据表
list_df <- list(
  CD4 = noise_CD4,
  CD8 = noise_CD8,
  MNP = noise_MNP,
  NK = noise_NK,
  B = noise_B,
  otherT = noise_otherT
)
# 去除列名中的前缀，并将每个数据框的列名后加上后缀
list_df <- setNames(lapply(names(list_df), function(name) {
  df <- list_df[[name]]
  # 获取原始列名
  original_colnames <- colnames(df)
  # 去掉原本的前缀部分
  colnames(df) <- sub(paste0("^", name, "\\."), "", original_colnames)
  # 给列加后缀（避免重复）
  colnames(df) <- paste0(colnames(df), "_", name)
  return(df)
}), names(list_df))
# 获取所有行名并集
all_rows <- Reduce(union, lapply(list_df, rownames))
# 获取所有列名并集
all_cols <- Reduce(union, lapply(list_df, colnames))
# 补齐缺失的行和列
list_df_filled <- lapply(list_df, function(df) {
  # 补齐缺失行
  missing_rows <- setdiff(all_rows, rownames(df))
  if (length(missing_rows) > 0) {
    add_rows <- as.data.frame(matrix(NA, nrow = length(missing_rows), ncol = ncol(df)))
    colnames(add_rows) <- colnames(df)
    rownames(add_rows) <- missing_rows
    df <- rbind(df, add_rows)
  }
  # 补齐缺失列
  missing_cols <- setdiff(all_cols, colnames(df))
  if (length(missing_cols) > 0) {
    add_cols <- as.data.frame(matrix(NA, nrow = nrow(df), ncol = length(missing_cols)))
    colnames(add_cols) <- missing_cols
    rownames(add_cols) <- rownames(df)
    df <- cbind(df, add_cols)
  }
  # 按照并集的列和行名排序
  df <- df[all_rows, all_cols, drop = FALSE]
  # 将第1列及以后的列转换为数值格式
  df[, 1:ncol(df)] <- lapply(df[, 1:ncol(df)], as.numeric)
  return(df)
})
# 合并所有表
noise_all <- do.call(cbind, list_df_filled)
### 去除全部为NA的基因列
noise_all <- Filter(function(x) !all(is.na(x)), noise_all)
# 去掉列名中'细胞类型.' 的前缀
# 定义要去除的前缀列表
prefixes <- c("CD4.", "CD8.", "MNP.", "NK.", "B.", "otherT.")
# 使用正则表达式构造可以匹配这些前缀的模式，采用管道符 `|` 进行 "或" 匹配
pattern <- paste0("^(", paste(prefixes, collapse = "|"), ")")
# 使用 gsub 去除这些前缀
colnames(noi_all) <- gsub(pattern, "", colnames(noi_all))
write.table(noise_all,file = "/data/wuwq/noise/DISEASE_ALL/app_disease/merge_COMBAT.csv", sep = ",", row.names = TRUE, col.names = NA, quote = TRUE)

# 更新二级结局--COMBAT 无误
#-- donor_id已经无重复，因此用sample_id没有问题，COMBAT_baseline中的10个数据含健康个体，虽然冗余但是根据id过滤疾病个体后并不影响结果
baseline <- read.delim("/data/wuwq/noise/DISEASE_ALL/baseline_COMBAT.txt", header = TRUE, row.names = 1)
rownames(baseline) <- baseline$scRNASeq_sample_ID
common_rows <- intersect(rownames(baseline), rownames(noise_all))
noise_d <- data.frame(row.names = common_rows)
noise_d <- cbind(noise_d, baseline[common_rows, , drop = FALSE],  noise_all[common_rows, , drop = FALSE] )
# 去除全部为NA的基因列
noise_d <- Filter(function(x) !all(is.na(x)), noise_d)
unique(noise_d$Disease_group) # "COVID"  "SEPSIS" "FLU"
unique(noise_d$Outcome) # 0  1  NA
# 把corrected_noise列前移
noise_d <- read.table("/data/wuwq/noise/DISEASE_ALL/app_disease/merge_COMBAT.csv", sep = ",", row.names=1, header = TRUE, check.names = FALSE, quote = "\"")
# 找出以 "corrrected_noise_" 开头的列
corrected_noise_cols <- grep("^corrected_noise_", names(noise_d), value = TRUE)
# 找出非 "corrrected_noise_" 的列
other_cols <- setdiff(names(noise_d), corrected_noise_cols)
# 保留前9列的基线列
baseline_cols <- other_cols[1:9]
# 剩余的非 "corrrected_noise_" 列
remaining_cols <- other_cols[10:length(other_cols)]
# 重新排列列的顺序
new_order <- c(baseline_cols, corrected_noise_cols, remaining_cols)
# 重新排列数据框列
noise_all_reordered <- noise_d %>%
  select(all_of(new_order))

write.table(noise_all_reordered,file = "/data/wuwq/noise/DISEASE_ALL/app_disease/merge_COMBAT.csv", sep = ",", row.names = TRUE, col.names = NA, quote = TRUE)


rm(list=ls())
#### merge 基因表达量数据   ✅️特殊读取
exp_disease <- read.table("/data/wuwq/noise/DISEASE_ALL/exp_u/union_onlydisease_CD4.csv", sep = ",", row.names=1, header = TRUE, check.names = FALSE, quote = "\"")
exp_CD4 <- subset(exp_disease, dataset == "COMBAT")
exp_disease <- read.table("/data/wuwq/noise/DISEASE_ALL/exp_u/union_onlydisease_CD8.csv", sep = ",", row.names=1, header = TRUE, check.names = FALSE, quote = "\"")
exp_CD8 <- subset(exp_disease, dataset == "COMBAT")
exp_disease <- read.table("/data/wuwq/noise/DISEASE_ALL/exp_u/union_onlydisease_MNP.csv", sep = ",", row.names=1, header = TRUE, check.names = FALSE, quote = "\"")
exp_MNP <- subset(exp_disease, dataset == "COMBAT")
exp_disease <- read.table("/data/wuwq/noise/DISEASE_ALL/exp_u/union_onlydisease_NK.csv", sep = ",", row.names=1, header = TRUE, check.names = FALSE, quote = "\"")
exp_NK <- subset(exp_disease, dataset == "COMBAT")
exp_disease <- read.table("/data/wuwq/noise/DISEASE_ALL/exp_u/union_onlydisease_B.csv", sep = ",", row.names=1, header = TRUE, check.names = FALSE, quote = "\"")
exp_B <- subset(exp_disease, dataset == "COMBAT")
exp_disease <- read.table("/data/wuwq/noise/DISEASE_ALL/exp_u/union_onlydisease_otherT.csv", sep = ",", row.names=1, header = TRUE, check.names = FALSE, quote = "\"")
exp_otherT <- subset(exp_disease, dataset == "COMBAT")
rm(exp_disease)

# 去掉前4列自带基线
exp_CD4 <- exp_CD4[,-c(1:4)]
exp_CD8 <- exp_CD8[,-c(1:4)]
exp_MNP <- exp_MNP[,-c(1:4)]
exp_NK <- exp_NK[,-c(1:4)]
exp_B <- exp_B[,-c(1:4)]
exp_otherT <- exp_otherT[,-c(1:4)]
# 读取并处理每个数据表
list_df <- list(
  CD4 = exp_CD4,
  CD8 = exp_CD8,
  MNP = exp_MNP,
  NK = exp_NK,
  B = exp_B,
  otherT = exp_otherT
)
# 去除列名中的前缀，并将每个数据框的列名后加上后缀
list_df <- setNames(lapply(names(list_df), function(name) {
  df <- list_df[[name]]
  # 获取原始列名
  original_colnames <- colnames(df)
  # 去掉原本的前缀部分
  colnames(df) <- sub(paste0("^", name, "\\."), "", original_colnames)
  # 给列加后缀（避免重复）
  colnames(df) <- paste0(colnames(df), "_", name)
  return(df)
}), names(list_df))
# 获取所有行名并集
all_rows <- Reduce(union, lapply(list_df, rownames))
# 获取所有列名并集
all_cols <- Reduce(union, lapply(list_df, colnames))
# 补齐缺失的行和列
list_df_filled <- lapply(list_df, function(df) {
  # 补齐缺失行
  missing_rows <- setdiff(all_rows, rownames(df))
  if (length(missing_rows) > 0) {
    add_rows <- as.data.frame(matrix(NA, nrow = length(missing_rows), ncol = ncol(df)))
    colnames(add_rows) <- colnames(df)
    rownames(add_rows) <- missing_rows
    df <- rbind(df, add_rows)
  }
  # 补齐缺失列
  missing_cols <- setdiff(all_cols, colnames(df))
  if (length(missing_cols) > 0) {
    add_cols <- as.data.frame(matrix(NA, nrow = nrow(df), ncol = length(missing_cols)))
    colnames(add_cols) <- missing_cols
    rownames(add_cols) <- rownames(df)
    df <- cbind(df, add_cols)
  }
  # 按照并集的列和行名排序
  df <- df[all_rows, all_cols, drop = FALSE]
  # 将第1列及以后的列转换为数值格式
  df[, 1:ncol(df)] <- lapply(df[, 1:ncol(df)], as.numeric)
  return(df)
})
# 合并所有表
exp_all <- do.call(cbind, list_df_filled)
### 去除全部为NA的基因列
exp_all <- Filter(function(x) !all(is.na(x)), exp_all)
# 去掉列名中'细胞类型.' 的前缀
# 定义要去除的前缀列表
prefixes <- c("CD4.", "CD8.", "MNP.", "NK.", "B.", "otherT.")
# 使用正则表达式构造可以匹配这些前缀的模式，采用管道符 `|` 进行 "或" 匹配
pattern <- paste0("^(", paste(prefixes, collapse = "|"), ")")
# 使用 gsub 去除这些前缀
colnames(exp_all) <- gsub(pattern, "", colnames(exp_all))



# 更新二级结局--COMBAT 无误
#-- donor_id已经无重复，因此用sample_id没有问题，COMBAT_baseline中的10个数据含健康个体，虽然冗余但是根据id过滤疾病个体后并不影响结果
baseline <- read.delim("/data/wuwq/noise/DISEASE_ALL/baseline_COMBAT.txt", header = TRUE, row.names = 1)
rownames(baseline) <- baseline$scRNASeq_sample_ID
common_rows <- intersect(rownames(baseline), rownames(exp_all))
exp_d <- data.frame(row.names = common_rows)
exp_d <- cbind(exp_d, baseline[common_rows, , drop = FALSE],  exp_all[common_rows, , drop = FALSE] )
# 去除全部为NA的基因列
exp_d <- Filter(function(x) !all(is.na(x)), exp_d)
# 去除全部为0的基因列
exp_d <- Filter(function(x) !all(x == 0, na.rm = TRUE), exp_d)
unique(exp_d$Disease_group) # "COVID"  "SEPSIS" "FLU"
unique(exp_d$Outcome) # 0  1  NA
write.table(exp_d,file = "/data/wuwq/noise/DISEASE_ALL/app_disease/merge_exp_COMBAT.csv", sep = ",", row.names = TRUE, col.names = NA, quote = TRUE)



rm(list=ls())
############################################## （二）和出院Outcome_score显著相关的99个独立于exp的noise基因
##############################################❗️因基因名为列名❗️️，故 ✅️特殊读取和 ✅️特殊保存
##### 1.noise  ✅️spearman 为主   &   pearson也做了   ### Spearman_rho 这个列名为了方便未改
noise_all <- read.table("/data/wuwq/noise/DISEASE_ALL/app_disease/merge_COMBAT.csv", sep = ",", row.names=1, header = TRUE, check.names = FALSE, quote = "\"")
# 基线变量
baseline_cols <- c('Outcome_score', 'Age_mid', 'sex')
# 提取基因-细胞类型对的列名
gene_cell_cols <- names(noise_all)[10:ncol(noise_all)]
# 初始化结果存储（用 data.frame 更安全）
results <- data.frame(
  Gene_CellType = gene_cell_cols,
  Spearman_rho = NA,
  p_value = NA,
  FDR = NA,
  stringsAsFactors = FALSE
)
# 先拟合 Outcome_score ~ Age_mid + sex 的残差
outcome_model <- lm(Outcome_score ~ Age_mid + sex, data = noise_all, na.action = na.exclude)
outcome_residuals <- residuals(outcome_model)
# 循环基因列
for (i in seq_along(gene_cell_cols)) {
  gene_col <- gene_cell_cols[i]
  gene_model <- lm(noise_all[[gene_col]] ~ noise_all$Age_mid + noise_all$sex, 
                   na.action = na.exclude)
  gene_residuals <- residuals(gene_model)
  # 去掉 NA
  valid_idx <- which(!is.na(outcome_residuals) & !is.na(gene_residuals))
  if (length(valid_idx) >= 3) {
    rho <- cor(outcome_residuals[valid_idx],
               gene_residuals[valid_idx],
               method = "pearson")
    p_value <- cor.test(outcome_residuals[valid_idx],
                        gene_residuals[valid_idx],
                        method = "pearson")$p.value
  } else {
    rho <- NA
    p_value <- NA
  }
  results$Spearman_rho[i] <- rho
  results$p_value[i] <- p_value
}
results$FDR <- p.adjust(results$p_value, method = "fdr")
results_sorted <- results[order(results$p_value), ]
#write.csv(results_sorted, file="/data/wuwq/noise/DISEASE_ALL/app_disease/COMBAT/Outcome_score_spearman_noi.csv")
write.csv(results_sorted, file="/data/wuwq/noise/DISEASE_ALL/app_disease/COMBAT/Outcome_score_pearson_noi.csv")


rm(list=ls())
##### 2.exp
exp_all <- read.table("/data/wuwq/noise/DISEASE_ALL/app_disease/merge_exp_COMBAT.csv", sep = ",", row.names=1, header = TRUE, check.names = FALSE, quote = "\"")
# 基线变量
baseline_cols <- c('Outcome_score', 'Age_mid', 'sex')
# 提取基因-细胞类型对的列名
gene_cell_cols <- names(exp_all)[10:ncol(exp_all)]
# 初始化结果存储（用 data.frame 更安全）
results <- data.frame(
  Gene_CellType = gene_cell_cols,
  Spearman_rho = NA,
  p_value = NA,
  FDR = NA,
  stringsAsFactors = FALSE
)
# 先拟合 Outcome_score ~ Age_mid + sex 的残差
outcome_model <- lm(Outcome_score ~ Age_mid + sex, data = exp_all, na.action = na.exclude)
outcome_residuals <- residuals(outcome_model)
# 循环基因列
for (i in seq_along(gene_cell_cols)) {
  gene_col <- gene_cell_cols[i]
  gene_model <- lm(exp_all[[gene_col]] ~ exp_all$Age_mid + exp_all$sex, 
                   na.action = na.exclude)
  gene_residuals <- residuals(gene_model)
  # 去掉 NA和0
  valid_idx <- which(!is.na(outcome_residuals) & 
                       !is.na(gene_residuals) &
                       outcome_residuals != 0 & 
                       gene_residuals != 0)
  if (length(valid_idx) >= 3) {
    rho <- cor(outcome_residuals[valid_idx],
               gene_residuals[valid_idx],
               method = "pearson")
    p_value <- cor.test(outcome_residuals[valid_idx],
                        gene_residuals[valid_idx],
                        method = "pearson")$p.value
  } else {
    rho <- NA
    p_value <- NA
  }
  results$Spearman_rho[i] <- rho
  results$p_value[i] <- p_value
}
results$FDR <- p.adjust(results$p_value, method = "fdr")
results_sorted <- results[order(results$p_value), ]
#write.csv(results_sorted, file="/data/wuwq/noise/DISEASE_ALL/app_disease/COMBAT/Outcome_score_spearman_exp.csv")
write.csv(results_sorted, file="/data/wuwq/noise/DISEASE_ALL/app_disease/COMBAT/Outcome_score_pearson_exp.csv")


rm(list=ls())
##### 3.读取noise和exp结果  得到 Outcome_score相关的 noi_change_exp_nochange 的gene_celltype及相关性情况
### spearman & pearson
noi_Outcome_score <- read.table("/data/wuwq/noise/DISEASE_ALL/app_disease/COMBAT/Outcome_score_spearman_noi.csv", sep = ",", row.names=1, header = TRUE, check.names = FALSE, quote = "\"")
noi_sig <- subset(noi_Outcome_score, (FDR < 0.05) & (Spearman_rho > 0.3 | Spearman_rho < -0.3))
noi_sig_gene <- noi_sig$Gene_CellType

exp_Outcome_score <- read.table("/data/wuwq/noise/DISEASE_ALL/app_disease/COMBAT/Outcome_score_spearman_exp.csv", sep = ",", row.names=1, header = TRUE, check.names = FALSE, quote = "\"")
exp_sig <- subset(exp_Outcome_score, (FDR < 0.05) & (Spearman_rho > 0.3 | Spearman_rho < -0.3))
exp_sig_gene <- exp_sig$Gene_CellType

# noi和exp共变
intersect_gene <- intersect(noi_sig_gene,exp_sig_gene)
# exp变而noi不变的
exp_specific <- setdiff(exp_sig_gene, intersect_gene)
# noi变而exp不变的 ✅️
noi_specific <- setdiff(noi_sig_gene, intersect_gene)
noi_specific <- setdiff(noi_sig_gene, intersect_gene)

# 关注noi变而exp不变的 相关性结果表，将2方面的相关性检验结果合并至一个表中 (卡较严格的FDR <0.05 & |rho|>0.3)
noi_specific_corr <- subset(noi_Outcome_score, Gene_CellType %in% noi_specific)
noi_specific_corr_in_exp <- subset(exp_Outcome_score, Gene_CellType %in% noi_specific)
colnames(noi_specific_corr_in_exp)[-1] <- paste0("exp_", colnames(noi_specific_corr_in_exp)[-1])
combined_data <- merge(noi_specific_corr, noi_specific_corr_in_exp, by = "Gene_CellType", all = TRUE)
combined_data <- subset(combined_data, !grepl("^NA_", Gene_CellType)) #删除以 "NA_" 开头的行
combined_data <- subset(combined_data, !grepl("^corrected_", Gene_CellType)) #删除以 "corrected_" 开头的行
write.csv(combined_data, file="/data/wuwq/noise/DISEASE_ALL/app_disease/COMBAT/Outcome_score_noi_spearman_change_exp_nochange.csv")
#write.csv(combined_data, file="/data/wuwq/noise/DISEASE_ALL/app_disease/COMBAT/Outcome_score_noi_pearson_change_exp_nochange.csv")



##### 4.做exp和noi-specific的相关性系数的比较图的数据准备
# 去掉以 "NA_" 开头的和以 "corrected_" 开头的数据
noi_specific <- noi_specific[!grepl("^(NA_|corrected_)", noi_specific)]
# 读取
noi_data <- read.table("/data/wuwq/noise/DISEASE_ALL/app_disease/merge_COMBAT.csv", sep = ",", row.names=1, header = TRUE, check.names = FALSE, quote = "\"")
noi_data <- noi_data[, c("scRNASeq_sample_ID","Outcome_score", noi_specific)]
exp_data <- read.table("/data/wuwq/noise/DISEASE_ALL/app_disease/merge_exp_COMBAT.csv", sep = ",", row.names=1, header = TRUE, check.names = FALSE, quote = "\"")
exp_data <- exp_data[, c("scRNASeq_sample_ID", noi_specific)]
# 为 noi_data 从第3列开始的列名添加前缀 "noi_"
colnames(noi_data)[-c(1,2)] <- paste0("noi_", colnames(noi_data)[-c(1,2)])
# 为 exp_data 从第2列开始的列名添加前缀 "exp_"
colnames(exp_data)[-1] <- paste0("exp_", colnames(exp_data)[-1])
# 合并两个数据框
noi_exp_data <- merge(noi_data, exp_data, by = "scRNASeq_sample_ID", all = TRUE)
col_noi_specific<-paste0("noi_", noi_specific)
col_exp_specific<-paste0("exp_", noi_specific)
noi_exp_data_plot <- noi_exp_data[, c("scRNASeq_sample_ID","Outcome_score",col_noi_specific)]




rm(list=ls())
############################################## （三）和出院Outcome 二分类变量做 逻辑回归 预测AUC
##############################################❗️因基因名为列名❗️️，故 ✅️特殊读取和 ✅️特殊保存
##### 合并成大的matrix
noi_data <- read.table("/data/wuwq/noise/DISEASE_ALL/app_disease/merge_COMBAT.csv", sep = ",", row.names=1, header = TRUE, check.names = FALSE, quote = "\"")
noi_matrix <- noi_data[,-c(2:6,8:15)]
exp_data <- read.table("/data/wuwq/noise/DISEASE_ALL/app_disease/merge_exp_COMBAT.csv", sep = ",", row.names=1, header = TRUE, check.names = FALSE, quote = "\"")
exp_matrix <- exp_data[, -c(2:6,8:9)]
# 为 noi_data 从第3列开始的列名添加前缀 "noi_"
colnames(noi_matrix)[-c(1,2)] <- paste0("noi_", colnames(noi_matrix)[-c(1,2)])
# 为 exp_data 从第3列开始的列名添加前缀 "exp_"
colnames(exp_matrix)[-c(1,2)] <- paste0("exp_", colnames(exp_matrix)[-c(1,2)])
matrix <- merge(noi_matrix, exp_matrix[,-2], by = "scRNASeq_sample_ID", all = TRUE)
# 保留有Outcome值的数据
noi_matrix <- subset(noi_matrix, Outcome != "NA")
exp_matrix <- subset(exp_matrix, Outcome != "NA")
matrix <- subset(matrix, Outcome != "NA")
# 将NA用0填充
for (i in 3:ncol(noi_matrix)) {
  noi_matrix[[i]][is.na(noi_matrix[[i]])] <- 0
}
for (i in 3:ncol(exp_matrix)) {
  exp_matrix[[i]][is.na(exp_matrix[[i]])] <- 0
}
for (i in 3:ncol(matrix)) {
  matrix[[i]][is.na(matrix[[i]])] <- 0
}



rm(list=ls())
# 直接加载准备好的矩阵
loaded_objs <- readRDS("/data/wuwq/noise/DISEASE_ALL/app_disease/COMBAT/noi_exp_matrix.rds")
exp_matrix <- loaded_objs[[1]]
noi_matrix <- loaded_objs[[2]]
matrix <- loaded_objs[[3]]
rm(loaded_objs)


##### Outcome=0 (outcome_score=1~4), Outcome=1 (outcome_score=5~6)
##### 1. 单变量初筛
### noise
# 创建一个空的数据框，用于存储每个变量的p值
p_noi <- data.frame(variable = character(), p_value = numeric(), coeff = numeric(), stringsAsFactors = FALSE)
# 遍历 noi_matrix 中的每一个变量列，从第三列开始
for (var in colnames(noi_matrix)[c(3:ncol(noi_matrix))]) {
  # 选择当前变量的非缺失值数据
  temp_data <- data.frame(Outcome = noi_matrix$Outcome, variable = noi_matrix[, var])
  temp_data <- na.omit(temp_data)  # 去掉含有 NA 的行
  # 在过滤后的数据上构建逻辑回归模型
  if (nrow(temp_data) > 0) {  # 确保模型有数据可计算
    noi_bi_model <- glm(Outcome ~ variable, data = temp_data, family = binomial)
    # 检查模型的系数矩阵是否有足够的维度
    model_summary <- summary(noi_bi_model)$coefficients
    if (nrow(model_summary) > 1) {
      # 提取 p 值和系数
      p_val <- model_summary[2, 4]
      co <- model_summary[2, 1]
      # 将结果添加到 p_noi 数据框中
      p_noi <- rbind(p_noi, data.frame(variable = var, p_value = p_val, coeff = co))
    } else {
      # 如果没有足够的系数，则记录信息以调试
      cat("Variable ", var, " does not have valid coefficients for Outcome.\n")
    }
  }
}
#计算coeff_abs及fdr
p_noi$coeff_abs <- abs(p_noi$coeff)
p_noi$fdr <- p.adjust(p_noi$p_value,method = "fdr")
p_noi <- p_noi %>% filter(p_value < 0.05)

##### noise_单变量回归中p_value<0.05的基因 的 AUC计算  ❗️此处有修改，再用时核对代码 ❗️
name <- as.character(p_noi$variable)
#name <- c("noi_SGK1_MNP","noi_ANXA1_MNP","noi_VIM_MNP","noi_SEC61B_B","noi_S100A9_MNP",
#          "noi_RETN_MNP","noi_ZFP36_otherT","noi_NFKBIA_NK","noi_FTH1_MNP","noi_TUBB_MNP",
#          "noi_HLA-DQA1_MNP","noi_JUN_NK","noi_IER2_NK","noi_S100A8_MNP","noi_PPP1R15A_NK",
#          "noi_TRBC1_CD4","noi_HMGN2_MNP","noi_SLC2A3_MNP","noi_NFKBIA_otherT","noi_CD52_CD4",
#          "noi_IL32_CD4","noi_CXCR4_NK","noi_S100A4_CD4","noi_HLA-DQB1_MNP","noi_CXCR4_CD4",
#          "noi_S100A11_CD4","noi_MIR23AHG_MNP","noi_RHOB_MNP")

set.seed(1)
# 初始化AUC值的列表（包括每个基因的结果）  5折交叉验证，400次重复抽样
auc_results <- data.frame(variable = character(), Mean_AUC = numeric(), stringsAsFactors = FALSE)
auc_values_list <- list()
folds <- createMultiFolds(y = noi_matrix$Outcome, k = 5, times = 400)
# 循环遍历每个基因名进行单因素逻辑回归
for (gene in name) {
  noi_auc_values <- numeric()
  # 进行交叉验证
  for (i in 1:length(folds)) {
    train <- noi_matrix[folds[[i]], c("Outcome", gene)]
    test <- noi_matrix[-folds[[i]], c("Outcome", gene)]
    noi_model <- glm(Outcome ~ ., family = binomial(link = logit), data = train)
    noi_model_pre <- predict(noi_model, type = 'response', newdata = test)
    auc_value <- as.numeric(auc(as.numeric(test[, "Outcome"]), noi_model_pre))
    noi_auc_values <- c(noi_auc_values, auc_value)
  }
  # 计算统计信息
  mean_auc <- mean(noi_auc_values)
  std_auc <- sd(noi_auc_values)
  se_auc <- std_auc / sqrt(length(noi_auc_values))
  # 保存每个基因的AUC结果
  auc_results <- rbind(auc_results, data.frame(variable = gene, Mean_AUC = mean_auc, Std_AUC = std_auc, SE_AUC = se_auc))
  auc_values_list[[gene]] <- noi_auc_values
}
# auc_results按AUC降序
AUC_noi_sorted <- merge(auc_results, p_noi, by = "variable")
AUC_noi_sorted <- AUC_noi_sorted[order(AUC_noi_sorted$Mean_AUC, decreasing = TRUE), ]

write.csv(AUC_noi_sorted, file="/data/wuwq/noise/DISEASE_ALL/app_disease/COMBAT/AUC_of_noi_sig.csv")
write.csv(auc_values_list, file="/data/wuwq/noise/DISEASE_ALL/app_disease/COMBAT/AUC_of_noi_sig_2000.csv")
#write.csv(auc_results, file="/data/wuwq/noise/DISEASE_ALL/app_disease/COMBAT/AUC_of_noi_sig_top28.csv")
#write.csv(auc_values_list, file="/data/wuwq/noise/DISEASE_ALL/app_disease/COMBAT/AUC_of_noi_matchnoi_top28.csv")



### 2.2.expression
# 创建一个空的数据框，用于存储每个变量的p值
p_exp <- data.frame(variable = character(), p_value = numeric(), coeff = numeric(), stringsAsFactors = FALSE)
# 遍历 exp_matrix 中的每一个变量列，从第三列开始
for (var in colnames(exp_matrix)[c(3:ncol(exp_matrix))]) {
  # 选择当前变量的非缺失值数据
  temp_data <- data.frame(Outcome = exp_matrix$Outcome, variable = exp_matrix[, var])
  temp_data <- na.omit(temp_data)  # 去掉含有 NA 的行
  # 在过滤后的数据上构建逻辑回归模型
  if (nrow(temp_data) > 0) {  # 确保模型有数据可计算
    exp_bi_model <- glm(Outcome ~ variable, data = temp_data, family = binomial)
    # 检查模型的系数矩阵是否有足够的维度
    model_summary <- summary(exp_bi_model)$coefficients
    if (nrow(model_summary) > 1) {
      # 提取 p 值和系数
      p_val <- model_summary[2, 4]
      co <- model_summary[2, 1]
      # 将结果添加到 p_exp 数据框中
      p_exp <- rbind(p_exp, data.frame(variable = var, p_value = p_val, coeff = co))
    } else {
      # 如果没有足够的系数，则记录信息以调试
      cat("Variable ", var, " does not have valid coefficients for Outcome.\n")
    }
  }
}
#计算coeff_abs及fdr
p_exp$coeff_abs <- abs(p_exp$coeff)
p_exp$fdr <- p.adjust(p_exp$p_value,method = "fdr")
p_exp <- p_exp %>% filter(p_value < 0.05)


### expression_单变量回归中p_value<0.05的基因 计算AUC  ❗️此处有修改，再用时核对代码 ❗️
#name <- as.character(p_exp$variable)
name <- c("exp_SGK1_MNP","exp_ANXA1_MNP","exp_VIM_MNP","exp_SEC61B_B","exp_S100A9_MNP",
          "exp_RETN_MNP","exp_ZFP36_otherT","exp_NFKBIA_NK","exp_FTH1_MNP","exp_TUBB_MNP",
          "exp_HLA-DQA1_MNP","exp_JUN_NK","exp_IER2_NK","exp_S100A8_MNP","exp_PPP1R15A_NK",
          "exp_TRBC1_CD4","exp_HMGN2_MNP","exp_SLC2A3_MNP","exp_NFKBIA_otherT","exp_CD52_CD4",
          "exp_IL32_CD4","exp_CXCR4_NK","exp_S100A4_CD4","exp_HLA-DQB1_MNP","exp_CXCR4_CD4",
          "exp_S100A11_CD4","exp_MIR23AHG_MNP","exp_RHOB_MNP"
)
name <- "exp_ID2_MNP"
set.seed(1)
# 初始化AUC值的列表（包括每个基因的结果）  5折交叉验证，400次重复抽样
auc_results <- data.frame(variable = character(), Mean_AUC = numeric(), stringsAsFactors = FALSE)
auc_values_list <- list()
folds <- createMultiFolds(y = exp_matrix$Outcome, k = 5, times = 400)
# 循环遍历每个基因名进行单因素逻辑回归
for (gene in name) {
  exp_auc_values <- numeric()
  # 进行交叉验证
  for (i in 1:length(folds)) {
    train <- exp_matrix[folds[[i]], c("Outcome", gene)]
    test <- exp_matrix[-folds[[i]], c("Outcome", gene)]
    exp_model <- glm(Outcome ~ ., family = binomial(link = logit), data = train)
    exp_model_pre <- predict(exp_model, type = 'response', newdata = test)
    auc_value <- as.numeric(auc(as.numeric(test[, "Outcome"]), exp_model_pre))
    exp_auc_values <- c(exp_auc_values, auc_value)
  }
  # 计算统计信息
  mean_auc <- mean(exp_auc_values)
  std_auc <- sd(exp_auc_values)
  se_auc <- std_auc / sqrt(length(exp_auc_values))
  # 保存每个基因的AUC结果
  auc_results <- rbind(auc_results, data.frame(variable = gene, Mean_AUC = mean_auc, Std_AUC = std_auc, SE_AUC = se_auc))
  auc_values_list[[gene]] <- exp_auc_values
}

#auc_results按AUC降序
#AUC_exp_sorted <- merge(auc_results, p_exp, by = "variable")
#AUC_exp_sorted <- AUC_exp_sorted[order(AUC_exp_sorted$Mean_AUC, decreasing = TRUE), ]
write.csv(auc_results, file="/data/wuwq/noise/DISEASE_ALL/app_disease/COMBAT/AUC_of_exp_matchnoi.csv")
write.csv(auc_values_list, file="/data/wuwq/noise/DISEASE_ALL/app_disease/COMBAT/AUC_of_exp_matchnoi_2000.csv")



rm(list=ls())
##### 3.读取和根据AUC降序选top基因进入构建多变量回归模型 5折交叉验证，100次重复   ✅️✅️✅️  逻辑回归模型
# 读取单变量回归得到的逻辑回归p_value<0.05的基因及AUC
AUC_noi_sorted <- read.table("/data/wuwq/noise/DISEASE_ALL/app_disease/COMBAT/AUC_of_noi_sig.csv", sep = ",", row.names=1, header = TRUE, check.names = FALSE, quote = "\"")
AUC_exp_sorted <- read.table("/data/wuwq/noise/DISEASE_ALL/app_disease/COMBAT/AUC_of_exp_sig.csv", sep = ",", row.names=1, header = TRUE, check.names = FALSE, quote = "\"")


### 3.1.noise_选top**个用于多变量回归
noi_significant_vars <- AUC_noi_sorted  %>% head(20)
significant_vars_name <- as.character(noi_significant_vars$variable)

# noi_matrix 
# set.seed 确保重复结果一致
n_repeats <- 100
k_folds   <- 5
fpr_grid  <- seq(0, 1, length.out = 201)

# 存储每次重复的 AUC
auc_values <- numeric(n_repeats)
# 存储每次重复的 ROC 曲线（在统一的 FPR 网格上插值）
fpr_grid   <- seq(0, 1, length.out = 201)
tpr_matrix <- matrix(NA, nrow = n_repeats, ncol = length(fpr_grid))
for (r in 1:n_repeats) {
  fold_ids <- createFolds(noi_matrix$Outcome, k = k_folds, list = TRUE)
  # 每个样本在这一轮恰好被预测一次
  pred_vec  <- numeric(nrow(noi_matrix))
  label_vec <- numeric(nrow(noi_matrix))
  for (f in seq_along(fold_ids)) {
    test_idx  <- fold_ids[[f]]
    train_idx <- setdiff(seq_len(nrow(noi_matrix)), test_idx)
    train <- noi_matrix[train_idx, c("Outcome", significant_vars_name)]
    test  <- noi_matrix[test_idx,  c("Outcome", significant_vars_name)]
    model <- glm(Outcome ~ ., family = binomial, data = train)
    pred_vec[test_idx]  <- predict(model, newdata = test, type = "response")
    label_vec[test_idx] <- as.numeric(test$Outcome)
  }
  # 这一轮的 ROC（每个样本只出现一次！）
  roc_r <- roc(label_vec, pred_vec, quiet = TRUE)
  auc_values[r] <- as.numeric(auc(roc_r))
  # 在统一 FPR 网格上插值 TPR
  tpr_matrix[r, ] <- approx(
    x = 1 - roc_r$specificities,   # FPR
    y = roc_r$sensitivities,       # TPR
    xout = fpr_grid, rule = 2
  )$y
}
# ====== 汇总结果 ======
mean_auc <- mean(auc_values)
ci_auc   <- quantile(auc_values, c(0.025, 0.975))
mean_tpr  <- colMeans(tpr_matrix)
lower_tpr <- apply(tpr_matrix, 2, quantile, 0.025)
upper_tpr <- apply(tpr_matrix, 2, quantile, 0.975)
# ====== 画图 ======
roc_df <- data.frame(
  FPR = fpr_grid,
  TPR = mean_tpr,
  TPR_lower = lower_tpr,
  TPR_upper = upper_tpr
)
library(ggplot2)
ggplot(roc_df, aes(x = FPR, y = TPR)) +
  geom_ribbon(aes(ymin = TPR_lower, ymax = TPR_upper), 
              fill = "grey70", alpha = 0.4) +
  geom_line(color = "blue", size = 1) +
  geom_abline(intercept = 0, slope = 1, 
              linetype = "dashed", color = "red") +
  labs(
    title = sprintf("ROC (Mean AUC = %.4f, 95%% CI: [%.4f, %.4f])",
                    mean_auc, ci_auc[1], ci_auc[2]),
    x = "False Positive Rate",
    y = "True Positive Rate"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

write.csv(roc_df, file="/data/wuwq/noise/DISEASE_ALL/app_disease/COMBAT/logistic_auc.csv")





### 3.2.expression_选top**个用于多变量回归
exp_significant_vars <- AUC_exp_sorted  %>% head(40)
significant_vars_name <- as.character(exp_significant_vars$variable)
# exp_matrix 
# set.seed 确保重复结果一致
set.seed(1)
# ==============================
# 3. 创建 5 折交叉验证，重复 10 次
# ==============================
folds <- createMultiFolds(y = exp_matrix$Outcome, k = 5, times = 100)
# 使用 list 存储预测和标签，最后再 unlist
all_predictions <- vector("list", length(folds))
all_labels <- vector("list", length(folds))
# ==============================
# 4. 循环每个 fold 拟合逻辑回归
# ==============================
for(i in seq_along(folds)) {
  train_idx <- folds[[i]]
  train <- exp_matrix[train_idx, c("Outcome", significant_vars_name)]
  test  <- exp_matrix[-train_idx, c("Outcome", significant_vars_name)]
  # 检查训练集中 Outcome 是否存在单一类
  if(length(unique(train$Outcome)) < 2) {
    cat("警告：第", i, "折训练集只有单一类别，跳过该fold\n")
    next
  }
  # 捕获 glm 错误
  tryCatch({
    exp_model <- glm(Outcome ~ ., family = binomial(link = "logit"), data = train)
    exp_model_pre <- predict(exp_model, type = "response", newdata = test)
    all_predictions[[i]] <- exp_model_pre
    all_labels[[i]] <- as.numeric(test$Outcome)
  }, error = function(e) {
    cat("警告：模型在第", i, "折发生完全分离，跳过该fold\n")
  })
}
# 合并所有 fold 的预测和标签
all_predictions <- unlist(all_predictions)
all_labels <- unlist(all_labels)
# ==============================
# 5. 计算整体 ROC 和 AUC
# ==============================
roc_obj <- roc(all_labels, all_predictions)
auc_val <- auc(roc_obj)
auc_ci <- ci.auc(roc_obj, conf.level = 0.95)
cat("整体AUC:", round(auc_val,4), "\n")
cat("95% CI: [", round(auc_ci[1],4), ",", round(auc_ci[3],4), "]\n")
# ==============================
# 6. 计算 TPR 的 95% CI
# ==============================
roc_ci <- ci.se(roc_obj, specificities = seq(0, 1, length.out = 100), conf.level = 0.95)
roc_df <- data.frame(
  FPR = 1 - seq(0, 1, length.out = 100), 
  TPR = roc_ci[,2],
  TPR_lower = pmax(0, roc_ci[,1]),
  TPR_upper = pmin(1, roc_ci[,3])
)
# ==============================
# 7. 绘制 ROC 曲线 + 95% CI
# ==============================
ggplot(roc_df, aes(x = FPR, y = TPR)) +
  geom_ribbon(aes(ymin = TPR_lower, ymax = TPR_upper), fill = "grey70", alpha = 0.4) +
  geom_line(color = "blue", size = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  labs(
    title = paste("ROC Curve (AUC =", round(auc_val,4),
                  ", 95% CI: [", round(auc_ci[1],4), ", ", round(auc_ci[3],4), "])"),
    x = "False Positive Rate (1 - Specificity)",
    y = "True Positive Rate (Sensitivity)"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))






### 3.3.合并用于模型构建的noise基因和expression基因
noi_significant_vars <- AUC_noi_sorted  %>% head(20)
exp_significant_vars <- AUC_exp_sorted  %>% head(20)
significant_vars_name <- unique(c(as.character(noi_significant_vars$variable), as.character(exp_significant_vars$variable)))
# matrix 
# set.seed 确保重复结果一致
set.seed(1)
# ==============================
# 3. 创建 5 折交叉验证，重复 10 次
# ==============================
folds <- createMultiFolds(y = matrix$Outcome, k = 5, times = 100)
# 使用 list 存储预测和标签，最后再 unlist
all_predictions <- vector("list", length(folds))
all_labels <- vector("list", length(folds))
# ==============================
# 4. 循环每个 fold 拟合逻辑回归
# ==============================
for(i in seq_along(folds)) {
  train_idx <- folds[[i]]
  train <- matrix[train_idx, c("Outcome", significant_vars_name)]
  test  <- matrix[-train_idx, c("Outcome", significant_vars_name)]
  # 检查训练集中 Outcome 是否存在单一类
  if(length(unique(train$Outcome)) < 2) {
    cat("警告：第", i, "折训练集只有单一类别，跳过该fold\n")
    next
  }
  # 捕获 glm 错误
  tryCatch({
    model <- glm(Outcome ~ ., family = binomial(link = "logit"), data = train)
    model_pre <- predict(model, type = "response", newdata = test)
    all_predictions[[i]] <- model_pre
    all_labels[[i]] <- as.numeric(test$Outcome)
  }, error = function(e) {
    cat("警告：模型在第", i, "折发生完全分离，跳过该fold\n")
  })
}
# 合并所有 fold 的预测和标签
all_predictions <- unlist(all_predictions)
all_labels <- unlist(all_labels)
# ==============================
# 5. 计算整体 ROC 和 AUC
# ==============================
roc_obj <- roc(all_labels, all_predictions)
auc_val <- auc(roc_obj)
auc_ci <- ci.auc(roc_obj, conf.level = 0.95)
cat("整体AUC:", round(auc_val,4), "\n")
cat("95% CI: [", round(auc_ci[1],4), ",", round(auc_ci[3],4), "]\n")
# ==============================
# 6. 计算 TPR 的 95% CI
# ==============================
roc_ci <- ci.se(roc_obj, specificities = seq(0, 1, length.out = 100), conf.level = 0.95)
roc_df <- data.frame(
  FPR = 1 - seq(0, 1, length.out = 100), 
  TPR = roc_ci[,2],
  TPR_lower = pmax(0, roc_ci[,1]),
  TPR_upper = pmin(1, roc_ci[,3])
)
# ==============================
# 7. 绘制 ROC 曲线 + 95% CI
# ==============================
ggplot(roc_df, aes(x = FPR, y = TPR)) +
  geom_ribbon(aes(ymin = TPR_lower, ymax = TPR_upper), fill = "grey70", alpha = 0.4) +
  geom_line(color = "blue", size = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  labs(
    title = paste("ROC Curve (AUC =", round(auc_val,4),
                  ", 95% CI: [", round(auc_ci[1],4), ", ", round(auc_ci[3],4), "])"),
    x = "False Positive Rate (1 - Specificity)",
    y = "True Positive Rate (Sensitivity)"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))





rm(list=ls())
##### 4.读取和根据单变量逻辑回归得到的AUC降序top基因进入构建随机森林模型 5折交叉验证，10次重复      ✅️✅️✅️️   随机森林模型

# 直接加载准备好的矩阵 由单变量逻辑回归得到的AUC_top gene
loaded_objs <- readRDS("/data/wuwq/noise/DISEASE_ALL/app_disease/COMBAT/noi_exp_matrix.rds")
exp_matrix <- loaded_objs[[1]]
noi_matrix <- loaded_objs[[2]]
matrix <- loaded_objs[[3]]
rm(loaded_objs)
# 把 Outcome 转换为因子变量，并确保因子水平合法
noi_matrix$Outcome <- as.factor(noi_matrix$Outcome)
levels(noi_matrix$Outcome) <- make.names(levels(noi_matrix$Outcome))
exp_matrix$Outcome <- as.factor(exp_matrix$Outcome)
levels(exp_matrix$Outcome) <- make.names(levels(exp_matrix$Outcome))
matrix$Outcome <- as.factor(matrix$Outcome)
levels(matrix$Outcome) <- make.names(levels(matrix$Outcome))


# 读取单变量逻辑回归得到的p_value<0.05的基因及AUC
AUC_noi_sorted <- read.table("/data/wuwq/noise/DISEASE_ALL/app_disease/COMBAT/AUC_of_noi_sig.csv", sep = ",", row.names=1, header = TRUE, check.names = FALSE, quote = "\"")
AUC_exp_sorted <- read.table("/data/wuwq/noise/DISEASE_ALL/app_disease/COMBAT/AUC_of_exp_sig.csv", sep = ",", row.names=1, header = TRUE, check.names = FALSE, quote = "\"")


###################################### a. noise
### noise_选top**个用于随机森林
noi_significant_vars <- AUC_noi_sorted  %>% head(6)
significant_vars_name <- as.character(noi_significant_vars$variable)
# Set seed for reproducibility
set.seed(1)
# Create folds for 5-fold cross-validation, repeated 400 times, with stratified sampling
cv_control <- trainControl(
  method = "repeatedcv",
  number = 5,
  repeats = 10,
  sampling = "up",
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = TRUE
)
# Fit the Random Forest model using caret
rf_model <- train(
  Outcome ~ ., 
  data = noi_matrix[, c("Outcome", significant_vars_name)],
  method = "rf",
  metric = "ROC",
  trControl = cv_control
)
# Collect AUC values from cross-validation
all_predictions <- as.numeric(rf_model$pred$pred)
all_labels <- rf_model$pred$obs
# Calculate the Mean AUC across all folds and repetitions
levels(all_labels) #需要保证正类在前，否则要翻转
roc_rf <- roc(all_labels, all_predictions, levels = c("X0", "X1"))
mean_auc_rf <- auc(roc_rf)
# Calculate confidence intervals for the ROC curve
# Use bootstrapping method for robust estimation
ci_roc_rf <- ci.se(roc_rf, specificities=seq(0, 1, length.out=100))
# 看一下结构
print(class(ci_roc_rf))
print(dim(ci_roc_rf))
head(ci_roc_rf)
# 组织作图数据
roc_df_rf <- data.frame(
  Specificity = as.numeric(rownames(ci_roc_rf)),
  FPR = 1 - as.numeric(rownames(ci_roc_rf)),
  TPR = ci_roc_rf[, "50%"],
  Lower_CI = ci_roc_rf[, "2.5%"],
  Upper_CI = ci_roc_rf[, "97.5%"]
)
# Plot the ROC curve for Random Forest with 95% CI
ggplot(roc_df_rf, aes(x = FPR, y = TPR)) +
  geom_line(color = "green", size = 1) +
  geom_ribbon(aes(ymin = Lower_CI, ymax = Upper_CI), fill = "gray80", alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  labs(
    title = paste("ROC Curve (Random Forest, AUC =", round(mean_auc_rf, 4), ")"),
    x = "False Positive Rate (1 - Specificity)",
    y = "True Positive Rate (Sensitivity)"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
cat("平均AUC (Random Forest):", round(mean_auc_rf, 4), "\n")



###################################### b. expression
### expression_选top**个用于随机森林
exp_significant_vars <- AUC_exp_sorted  %>% head(6)
significant_vars_name <- as.character(exp_significant_vars$variable)
# Set seed for reproducibility
set.seed(1)
# Create folds for 5-fold cross-validation, repeated 400 times, with stratified sampling
cv_control <- trainControl(
  method = "repeatedcv",
  number = 5,
  repeats = 10,
  sampling = "up",
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = TRUE
)
# Fit the Random Forest model using caret
rf_model <- train(
  Outcome ~ ., 
  data = exp_matrix[, c("Outcome", significant_vars_name)],
  method = "rf",
  metric = "ROC",
  trControl = cv_control
)
# Collect AUC values from cross-validation
all_predictions <- as.numeric(rf_model$pred$pred)
all_labels <- rf_model$pred$obs
# Calculate the Mean AUC across all folds and repetitions
levels(all_labels) #需要保证正类在前，否则要翻转
roc_rf <- roc(all_labels, all_predictions, levels = c("X0", "X1"))
mean_auc_rf <- auc(roc_rf)
# Calculate confidence intervals for the ROC curve
# Use bootstrapping method for robust estimation
ci_roc_rf <- ci.se(roc_rf, specificities=seq(0, 1, length.out=100))
# 看一下结构
print(class(ci_roc_rf))
print(dim(ci_roc_rf))
head(ci_roc_rf)
# 组织作图数据
roc_df_rf <- data.frame(
  Specificity = as.numeric(rownames(ci_roc_rf)),
  FPR = 1 - as.numeric(rownames(ci_roc_rf)),
  TPR = ci_roc_rf[, "50%"],
  Lower_CI = ci_roc_rf[, "2.5%"],
  Upper_CI = ci_roc_rf[, "97.5%"]
)
# Plot the ROC curve for Random Forest with 95% CI
ggplot(roc_df_rf, aes(x = FPR, y = TPR)) +
  geom_line(color = "green", size = 1) +
  geom_ribbon(aes(ymin = Lower_CI, ymax = Upper_CI), fill = "gray80", alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  labs(
    title = paste("ROC Curve (Random Forest, AUC =", round(mean_auc_rf, 4), ")"),
    x = "False Positive Rate (1 - Specificity)",
    y = "True Positive Rate (Sensitivity)"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
cat("平均AUC (Random Forest):", round(mean_auc_rf, 4), "\n")



###################################### c. noise + expression
noi_significant_vars <- AUC_noi_sorted  %>% head(3)
exp_significant_vars <- AUC_exp_sorted  %>% head(3)
significant_vars_name <- unique(c(as.character(noi_significant_vars$variable), as.character(exp_significant_vars$variable)))
# Set seed for reproducibility
set.seed(1)
# Create folds for 5-fold cross-validation, repeated 400 times, with stratified sampling
cv_control <- trainControl(
  method = "repeatedcv",
  number = 5,
  repeats = 10,
  sampling = "up",
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = TRUE
)
# Fit the Random Forest model using caret
rf_model <- train(
  Outcome ~ ., 
  data = matrix[, c("Outcome", significant_vars_name)],
  method = "rf",
  metric = "ROC",
  trControl = cv_control
)
# Collect AUC values from cross-validation
all_predictions <- as.numeric(rf_model$pred$pred)
all_labels <- rf_model$pred$obs
# Calculate the Mean AUC across all folds and repetitions
levels(all_labels) #需要保证正类在前，否则要翻转
roc_rf <- roc(all_labels, all_predictions, levels = c("X0", "X1"))
mean_auc_rf <- auc(roc_rf)
# Calculate confidence intervals for the ROC curve
# Use bootstrapping method for robust estimation
ci_roc_rf <- ci.se(roc_rf, specificities=seq(0, 1, length.out=100))
# 看一下结构
print(class(ci_roc_rf))
print(dim(ci_roc_rf))
head(ci_roc_rf)
# 组织作图数据
roc_df_rf <- data.frame(
  Specificity = as.numeric(rownames(ci_roc_rf)),
  FPR = 1 - as.numeric(rownames(ci_roc_rf)),
  TPR = ci_roc_rf[, "50%"],
  Lower_CI = ci_roc_rf[, "2.5%"],
  Upper_CI = ci_roc_rf[, "97.5%"]
)
# Plot the ROC curve for Random Forest with 95% CI
ggplot(roc_df_rf, aes(x = FPR, y = TPR)) +
  geom_line(color = "green", size = 1) +
  geom_ribbon(aes(ymin = Lower_CI, ymax = Upper_CI), fill = "gray80", alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  labs(
    title = paste("ROC Curve (Random Forest, AUC =", round(mean_auc_rf, 4), ")"),
    x = "False Positive Rate (1 - Specificity)",
    y = "True Positive Rate (Sensitivity)"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
cat("平均AUC (Random Forest):", round(mean_auc_rf, 4), "\n")


write.csv(roc_df_rf, file="/data/wuwq/noise/DISEASE_ALL/app_disease/COMBAT/RF_auc.csv")










rm(list=ls())
############################################## （四）出院Outcome 二分类变量做 随机森林模型
##############################################❗️因基因名为列名❗️️，故 ✅️特殊读取和 ✅️特殊保存
loaded_objs <- readRDS("/data/wuwq/noise/DISEASE_ALL/app_disease/COMBAT/noi_exp_matrix.rds")
exp_matrix <- loaded_objs[[1]]
noi_matrix <- loaded_objs[[2]]
matrix <- loaded_objs[[3]]
rm(loaded_objs)
# 把 Outcome 转换为因子变量，并确保因子水平合法
noi_matrix$Outcome <- as.factor(noi_matrix$Outcome)
levels(noi_matrix$Outcome) <- make.names(levels(noi_matrix$Outcome))
exp_matrix$Outcome <- as.factor(exp_matrix$Outcome)
levels(exp_matrix$Outcome) <- make.names(levels(exp_matrix$Outcome))
matrix$Outcome <- as.factor(matrix$Outcome)
levels(matrix$Outcome) <- make.names(levels(matrix$Outcome))


##### 1. 单变量初筛 -- 随机森林模型
### noise
# 初始化数据框，用于存储每个基因的AUC结果
auc_results <- data.frame(variable = character(), Mean_AUC = numeric(), Std_AUC = numeric(), SE_AUC = numeric(), stringsAsFactors = FALSE)
auc_values_list <- list()
# 分层 5 折交叉验证和 1 次重复
set.seed(1)
folds <- createMultiFolds(y = noi_matrix$Outcome, k = 5, times = 1)
control <- trainControl(method = "repeatedcv", number = 5, repeats = 1, classProbs = TRUE, summaryFunction = twoClassSummary, savePredictions = TRUE)
# 循环遍历每个变量进行单因素随机森林模型
for (var in colnames(noi_matrix)[3:ncol(noi_matrix)]) {
  # 过滤非缺失值数据
  temp_data <- noi_matrix[, c("Outcome", var)]
  temp_data <- na.omit(temp_data)
  # 确认数据不为空
  if (nrow(temp_data) > 0) {
    # 将Outcome转换为因子类型，确保随机森林可以处理
    temp_data$Outcome <- as.factor(temp_data$Outcome)
    # 使用caret的train函数进行交叉验证和模型训练
    rf_model <- train(Outcome ~ ., data = temp_data, method = "rf", trControl = control, metric = "ROC")
    # 获取AUC值
    auc_values <- rf_model$results$ROC
    # 计算统计信息
    mean_auc <- mean(auc_values)
    std_auc <- sd(auc_values)
    se_auc <- std_auc / sqrt(length(auc_values))
    # 保存每个基因的AUC结果
    auc_results <- rbind(auc_results, data.frame(variable = var, Mean_AUC = mean_auc, Std_AUC = std_auc, SE_AUC = se_auc))
    auc_values_list[[var]] <- auc_values
  }
}
# AUC结果按AUC降序排列
AUC_noi_sorted <- auc_results[order(auc_results$Mean_AUC, decreasing = TRUE), ]
write.csv(AUC_noi_sorted, file="/data/wuwq/noise/DISEASE_ALL/app_disease/COMBAT/RF_AUC_noi.csv")

### expression
# 初始化数据框，用于存储每个基因的AUC结果
# 初始化数据框，用于存储每个基因的AUC结果
auc_results <- data.frame(variable = character(), Mean_AUC = numeric(), Std_AUC = numeric(), SE_AUC = numeric(), stringsAsFactors = FALSE)
auc_values_list <- list()
# 分层 5 折交叉验证和 1 次重复
set.seed(1)
folds <- createMultiFolds(y = exp_matrix$Outcome, k = 5, times = 1)
control <- trainControl(method = "repeatedcv", number = 5, repeats = 1, classProbs = TRUE, summaryFunction = twoClassSummary, savePredictions = TRUE)
# 循环遍历每个变量进行单因素随机森林模型
for (var in colnames(exp_matrix)[3:ncol(exp_matrix)]) {
  # 过滤非缺失值数据
  temp_data <- exp_matrix[, c("Outcome", var)]
  temp_data <- na.omit(temp_data)
  # 确认数据不为空
  if (nrow(temp_data) > 0) {
    # 将Outcome转换为因子类型，确保随机森林可以处理
    temp_data$Outcome <- as.factor(temp_data$Outcome)
    # 使用caret的train函数进行交叉验证和模型训练
    rf_model <- train(Outcome ~ ., data = temp_data, method = "rf", trControl = control, metric = "ROC")
    # 获取AUC值
    auc_values <- rf_model$results$ROC
    # 计算统计信息
    mean_auc <- mean(auc_values)
    std_auc <- sd(auc_values)
    se_auc <- std_auc / sqrt(length(auc_values))
    # 保存每个基因的AUC结果
    auc_results <- rbind(auc_results, data.frame(variable = var, Mean_AUC = mean_auc, Std_AUC = std_auc, SE_AUC = se_auc))
    auc_values_list[[var]] <- auc_values
  }
}
# AUC结果按AUC降序排列
AUC_exp_sorted <- auc_results[order(auc_results$Mean_AUC, decreasing = TRUE), ]
write.csv(AUC_exp_sorted, file="/data/wuwq/expse/DISEASE_ALL/app_disease/COMBAT/RF_AUC_exp.csv")



























rm(list=ls())
############################################## （五）出院Outcome 二分类变量做 秩和检验，锁定noise变而exp不变的基因
##############################################❗️因基因名为列名❗️️，故 ✅️特殊读取和 ✅️特殊保存
# 构建noi_matrix, exp_matrix, matrix
noi_data <- read.table("/data/wuwq/noise/DISEASE_ALL/app_disease/merge_COMBAT.csv", sep = ",", row.names=1, header = TRUE, check.names = FALSE, quote = "\"")
noi_matrix <- noi_data[,-c(2:6,8:15)]
exp_data <- read.table("/data/wuwq/noise/DISEASE_ALL/app_disease/merge_exp_COMBAT.csv", sep = ",", row.names=1, header = TRUE, check.names = FALSE, quote = "\"")
exp_matrix <- exp_data[, -c(2:6,8:9)]
colnames(noi_matrix)[-c(1,2)] <- paste0("noi_", colnames(noi_matrix)[-c(1,2)])
colnames(exp_matrix)[-c(1,2)] <- paste0("exp_", colnames(exp_matrix)[-c(1,2)])
matrix <- merge(noi_matrix, exp_matrix[,-2], by = "scRNASeq_sample_ID", all = TRUE)
# 保留有Outcome值的数据
noi_matrix <- subset(noi_matrix, Outcome != "NA")
exp_matrix <- subset(exp_matrix, Outcome != "NA")
matrix <- subset(matrix, Outcome != "NA")
# 将NA用R填充
for (i in 3:ncol(noi_matrix)) {
  noi_matrix[[i]][is.na(noi_matrix[[i]])] <- 0
}
for (i in 3:ncol(exp_matrix)) {
  exp_matrix[[i]][is.na(exp_matrix[[i]])] <- 0
}
for (i in 3:ncol(matrix)) {
  matrix[[i]][is.na(matrix[[i]])] <- 0
}
# 封装保存以上矩阵
objs_to_save <- list(exp_matrix, noi_matrix, matrix)
saveRDS(objs_to_save, file = "/data/wuwq/noise/DISEASE_ALL/app_disease/COMBAT/noi_exp_matrix.rds")


rm(list=ls())
# 直接加载准备好的矩阵
loaded_objs <- readRDS("/data/wuwq/noise/DISEASE_ALL/app_disease/COMBAT/noi_exp_matrix.rds")
exp_matrix <- loaded_objs[[1]]
noi_matrix <- loaded_objs[[2]]
matrix <- loaded_objs[[3]]
rm(loaded_objs)

# 确定从第3列到最后一列为数值型 & 比较组标签是字符串格式
merge_matrix <- matrix
merge_matrix$Outcome <- as.character(merge_matrix$Outcome)

# 自定义两两比较的函数
pairwise_wilcox_gene <- function(data, group_col, gene_start_col = 7) {
  # 基本检查
  if (!group_col %in% colnames(data)) {
    stop("group_col not found in data")
  }
  groups <- unique(data[[group_col]])
  groups <- groups[!is.na(groups)]
  if (length(groups) < 2) {
    stop("Need at least 2 groups to compare")
  }
  gene_cols <- gene_start_col:ncol(data)
  gene_names <- colnames(data)[gene_cols]
  # 预分配 list
  result_list <- vector("list", length(gene_names) * choose(length(groups), 2))
  idx <- 1
  # 为了效率，先把分组索引提前算好（可选，但大数据集有帮助）
  group_indices <- lapply(groups, function(g) which(data[[group_col]] == g))
  names(group_indices) <- groups
  for (gene in gene_names) {
    expr <- data[[gene]]
    for (i in 1:(length(groups)-1)) {
      for (j in (i+1):length(groups)) {
        g1 <- groups[i]
        g2 <- groups[j]
        idx1 <- group_indices[[g1]]
        idx2 <- group_indices[[g2]]
        x_all <- expr[idx1]
        y_all <- expr[idx2]
        # 总样本数（包含 NA）
        total_n1 <- length(x_all)
        total_n2 <- length(y_all)
        # 有效值
        x <- x_all[!is.na(x_all)]
        y <- y_all[!is.na(y_all)]
        n1 <- length(x)
        n2 <- length(y)
        na_count1 <- total_n1 - n1
        na_count2 <- total_n2 - n2
        # 统一结构（新增 4 列）
        res_row <- data.frame(
          gene       = gene,
          group1     = g1,
          group2     = g2,
          mean1      = NA_real_,
          mean2      = NA_real_,
          log2FC     = NA_real_,
          p_value    = NA_real_,
          fdr        = NA_real_,
          n1         = n1,         # 有效样本数
          n2         = n2,   # 有效样本数
          na_count1  = na_count1,
          na_count2  = na_count2,
          total_n1   = total_n1,
          total_n2   = total_n2,
          stringsAsFactors = FALSE
        )
        # 样本量不足 → 打印提示，保留 NA 行
        if (n1 < 2 || n2 < 2) {
          message(sprintf("Skipped %s: %s vs %s  (effective n1=%d, n2=%d | total_n1=%d, total_n2=%d)", 
                          gene, g1, g2, n1, n2, total_n1, total_n2))
          result_list[[idx]] <- res_row
          idx <- idx + 1
          next
        }
        # 正常计算
        mean1 <- mean(x, na.rm = TRUE)
        mean2 <- mean(y, na.rm = TRUE)
        pseudocount <- 0.1
        log2fc <- log2((mean1 + pseudocount) / (mean2 + pseudocount))
        test_result <- tryCatch(
          wilcox.test(x, y, exact = FALSE),
          error = function(e) {
            message("wilcox.test failed for ", gene, " ", g1, " vs ", g2, " - ", e$message)
            list(p.value = NA_real_)
          },
          warning = function(w) list(p.value = NA_real_)
        )
        res_row$mean1   <- mean1
        res_row$mean2   <- mean2
        res_row$log2FC  <- log2fc
        res_row$p_value <- test_result$p.value
        result_list[[idx]] <- res_row
        idx <- idx + 1
      }
    }
  }
  # 合并所有结果
  results <- do.call(rbind, result_list[!sapply(result_list, is.null)])
  # FDR 校正（只对有效 p 值）
  valid <- !is.na(results$p_value)
  if (any(valid)) {
    results$fdr[valid] <- p.adjust(results$p_value[valid], method = "fdr")
  }
  # 排序：fdr 优先，NA 放最后
  results <- results[order(results$fdr, results$p_value, na.last = TRUE), ]
  # 清理 log2FC
  results$log2FC[is.na(results$mean1) | is.na(results$mean2)] <- NA_real_
  rownames(results) <- NULL
  return(results)
}

result <- pairwise_wilcox_gene(merge_matrix, "Outcome", gene_start_col = 3)  
write.csv(result, file="/data/wuwq/noise/DISEASE_ALL/app_disease/COMBAT/Outcome_wilcox.csv", row.names = TRUE)







rm(list=ls())
############################################## （六）其他 秩和检验
##############################################❗️因基因名为列名❗️️，故 ✅️特殊读取和 ✅️特殊保存
# 构建基线完整的矩阵
noi_data <- read.table("/data/wuwq/noise/DISEASE_ALL/app_disease/merge_COMBAT.csv", sep = ",", row.names=1, header = TRUE, check.names = FALSE, quote = "\"")
noi_matrix <- noi_data
exp_data <- read.table("/data/wuwq/noise/DISEASE_ALL/app_disease/merge_exp_COMBAT.csv", sep = ",", row.names=1, header = TRUE, check.names = FALSE, quote = "\"")
exp_matrix <- exp_data
colnames(noi_matrix)[-c(1:15)] <- paste0("noi_", colnames(noi_matrix)[-c(1:15)])
colnames(exp_matrix)[-c(1:9)] <- paste0("exp_", colnames(exp_matrix)[-c(1:9)])
matrix <- merge(noi_matrix, exp_matrix[,-c(2:9)], by = "scRNASeq_sample_ID", all = TRUE)

unique(matrix$Disease_group) # "COVID"  "SEPSIS" "FLU"
unique(matrix$Disease_stage) # "MILD" NA     "CRIT" "SEV"
unique(matrix$Death28) # 0 1
unique(matrix$Outcome) # 0  1 NA

########  1.Disease_group 为秩和检验分组
merge_matrix <- matrix
merge_matrix$Disease_group <- as.character(merge_matrix$Disease_group)

# 自定义两两比较的函数
pairwise_wilcox_gene <- function(data, group_col, gene_start_col = 10) {
  # 基本检查
  if (!group_col %in% colnames(data)) {
    stop("group_col not found in data")
  }
  groups <- unique(data[[group_col]])
  groups <- groups[!is.na(groups)]
  if (length(groups) < 2) {
    stop("Need at least 2 groups to compare")
  }
  gene_cols <- gene_start_col:ncol(data)
  gene_names <- colnames(data)[gene_cols]
  # 预分配 list
  result_list <- vector("list", length(gene_names) * choose(length(groups), 2))
  idx <- 1
  # 为了效率，先把分组索引提前算好（可选，但大数据集有帮助）
  group_indices <- lapply(groups, function(g) which(data[[group_col]] == g))
  names(group_indices) <- groups
  for (gene in gene_names) {
    expr <- data[[gene]]
    for (i in 1:(length(groups)-1)) {
      for (j in (i+1):length(groups)) {
        g1 <- groups[i]
        g2 <- groups[j]
        idx1 <- group_indices[[g1]]
        idx2 <- group_indices[[g2]]
        x_all <- expr[idx1]
        y_all <- expr[idx2]
        # 总样本数（包含 NA）
        total_n1 <- length(x_all)
        total_n2 <- length(y_all)
        # 有效值
        x <- x_all[!is.na(x_all)]
        y <- y_all[!is.na(y_all)]
        n1 <- length(x)
        n2 <- length(y)
        na_count1 <- total_n1 - n1
        na_count2 <- total_n2 - n2
        # 统一结构（新增 4 列）
        res_row <- data.frame(
          gene       = gene,
          group1     = g1,
          group2     = g2,
          mean1      = NA_real_,
          mean2      = NA_real_,
          log2FC     = NA_real_,
          p_value    = NA_real_,
          fdr        = NA_real_,
          n1         = n1,         # 有效样本数
          n2         = n2,   # 有效样本数
          na_count1  = na_count1,
          na_count2  = na_count2,
          total_n1   = total_n1,
          total_n2   = total_n2,
          stringsAsFactors = FALSE
        )
        # 样本量不足 → 打印提示，保留 NA 行
        if (n1 < 2 || n2 < 2) {
          message(sprintf("Skipped %s: %s vs %s  (effective n1=%d, n2=%d | total_n1=%d, total_n2=%d)", 
                          gene, g1, g2, n1, n2, total_n1, total_n2))
          result_list[[idx]] <- res_row
          idx <- idx + 1
          next
        }
        # 正常计算
        mean1 <- mean(x, na.rm = TRUE)
        mean2 <- mean(y, na.rm = TRUE)
        pseudocount <- 0.1
        log2fc <- log2((mean1 + pseudocount) / (mean2 + pseudocount))
        test_result <- tryCatch(
          wilcox.test(x, y, exact = FALSE),
          error = function(e) {
            message("wilcox.test failed for ", gene, " ", g1, " vs ", g2, " - ", e$message)
            list(p.value = NA_real_)
          },
          warning = function(w) list(p.value = NA_real_)
        )
        res_row$mean1   <- mean1
        res_row$mean2   <- mean2
        res_row$log2FC  <- log2fc
        res_row$p_value <- test_result$p.value
        result_list[[idx]] <- res_row
        idx <- idx + 1
      }
    }
  }
  # 合并所有结果
  results <- do.call(rbind, result_list[!sapply(result_list, is.null)])
  # FDR 校正（只对有效 p 值）
  valid <- !is.na(results$p_value)
  if (any(valid)) {
    results$fdr[valid] <- p.adjust(results$p_value[valid], method = "fdr")
  }
  # 排序：fdr 优先，NA 放最后
  results <- results[order(results$fdr, results$p_value, na.last = TRUE), ]
  # 清理 log2FC
  results$log2FC[is.na(results$mean1) | is.na(results$mean2)] <- NA_real_
  rownames(results) <- NULL
  return(results)
}

result <- pairwise_wilcox_gene(merge_matrix, "Disease_group", gene_start_col = 10)  
write.csv(result, file="/data/wuwq/noise/DISEASE_ALL/app_disease/COMBAT/Disease_group_wilcox.csv", row.names = TRUE)




########  2.Disease_stage 为秩和检验分组 (仅COVID有stage,n=77)
merge_matrix <- subset(matrix,Disease_stage != "NA" )
merge_matrix$Disease_stage <- as.character(merge_matrix$Disease_stage)

# 自定义两两比较的函数
pairwise_wilcox_gene <- function(data, group_col, gene_start_col = 6) {
  # 基本检查
  if (!group_col %in% colnames(data)) {
    stop("group_col not found in data")
  }
  groups <- unique(data[[group_col]])
  groups <- groups[!is.na(groups)]
  if (length(groups) < 2) {
    stop("Need at least 2 groups to compare")
  }
  gene_cols <- gene_start_col:ncol(data)
  gene_names <- colnames(data)[gene_cols]
  # 预分配 list
  result_list <- vector("list", length(gene_names) * choose(length(groups), 2))
  idx <- 1
  # 为了效率，先把分组索引提前算好（可选，但大数据集有帮助）
  group_indices <- lapply(groups, function(g) which(data[[group_col]] == g))
  names(group_indices) <- groups
  for (gene in gene_names) {
    expr <- data[[gene]]
    for (i in 1:(length(groups)-1)) {
      for (j in (i+1):length(groups)) {
        g1 <- groups[i]
        g2 <- groups[j]
        idx1 <- group_indices[[g1]]
        idx2 <- group_indices[[g2]]
        x_all <- expr[idx1]
        y_all <- expr[idx2]
        # 总样本数（包含 NA）
        total_n1 <- length(x_all)
        total_n2 <- length(y_all)
        # 有效值
        x <- x_all[!is.na(x_all)]
        y <- y_all[!is.na(y_all)]
        n1 <- length(x)
        n2 <- length(y)
        na_count1 <- total_n1 - n1
        na_count2 <- total_n2 - n2
        # 统一结构（新增 4 列）
        res_row <- data.frame(
          gene       = gene,
          group1     = g1,
          group2     = g2,
          mean1      = NA_real_,
          mean2      = NA_real_,
          log2FC     = NA_real_,
          p_value    = NA_real_,
          fdr        = NA_real_,
          n1         = n1,         # 有效样本数
          n2         = n2,   # 有效样本数
          na_count1  = na_count1,
          na_count2  = na_count2,
          total_n1   = total_n1,
          total_n2   = total_n2,
          stringsAsFactors = FALSE
        )
        # 样本量不足 → 打印提示，保留 NA 行
        if (n1 < 2 || n2 < 2) {
          message(sprintf("Skipped %s: %s vs %s  (effective n1=%d, n2=%d | total_n1=%d, total_n2=%d)", 
                          gene, g1, g2, n1, n2, total_n1, total_n2))
          result_list[[idx]] <- res_row
          idx <- idx + 1
          next
        }
        # 正常计算
        mean1 <- mean(x, na.rm = TRUE)
        mean2 <- mean(y, na.rm = TRUE)
        pseudocount <- 0.1
        log2fc <- log2((mean1 + pseudocount) / (mean2 + pseudocount))
        test_result <- tryCatch(
          wilcox.test(x, y, exact = FALSE),
          error = function(e) {
            message("wilcox.test failed for ", gene, " ", g1, " vs ", g2, " - ", e$message)
            list(p.value = NA_real_)
          },
          warning = function(w) list(p.value = NA_real_)
        )
        res_row$mean1   <- mean1
        res_row$mean2   <- mean2
        res_row$log2FC  <- log2fc
        res_row$p_value <- test_result$p.value
        result_list[[idx]] <- res_row
        idx <- idx + 1
      }
    }
  }
  # 合并所有结果
  results <- do.call(rbind, result_list[!sapply(result_list, is.null)])
  # FDR 校正（只对有效 p 值）
  valid <- !is.na(results$p_value)
  if (any(valid)) {
    results$fdr[valid] <- p.adjust(results$p_value[valid], method = "fdr")
  }
  # 排序：fdr 优先，NA 放最后
  results <- results[order(results$fdr, results$p_value, na.last = TRUE), ]
  # 清理 log2FC
  results$log2FC[is.na(results$mean1) | is.na(results$mean2)] <- NA_real_
  rownames(results) <- NULL
  return(results)
}

result <- pairwise_wilcox_gene(merge_matrix, "Disease_stage", gene_start_col = 6)  
write.csv(result, file="/data/wuwq/noise/DISEASE_ALL/app_disease/COMBAT/Disease_stage_wilcox.csv", row.names = TRUE)




########  3.Disease_stage 为秩和检验分组 (仅COVID有stage,n=77)
merge_matrix <- matrix
merge_matrix$Death28 <- as.character(merge_matrix$Death28)

# 自定义两两比较的函数
pairwise_wilcox_gene <- function(data, group_col, gene_start_col = 10) {
  # 基本检查
  if (!group_col %in% colnames(data)) {
    stop("group_col not found in data")
  }
  groups <- unique(data[[group_col]])
  groups <- groups[!is.na(groups)]
  if (length(groups) < 2) {
    stop("Need at least 2 groups to compare")
  }
  gene_cols <- gene_start_col:ncol(data)
  gene_names <- colnames(data)[gene_cols]
  # 预分配 list
  result_list <- vector("list", length(gene_names) * choose(length(groups), 2))
  idx <- 1
  # 为了效率，先把分组索引提前算好（可选，但大数据集有帮助）
  group_indices <- lapply(groups, function(g) which(data[[group_col]] == g))
  names(group_indices) <- groups
  for (gene in gene_names) {
    expr <- data[[gene]]
    for (i in 1:(length(groups)-1)) {
      for (j in (i+1):length(groups)) {
        g1 <- groups[i]
        g2 <- groups[j]
        idx1 <- group_indices[[g1]]
        idx2 <- group_indices[[g2]]
        x_all <- expr[idx1]
        y_all <- expr[idx2]
        # 总样本数（包含 NA）
        total_n1 <- length(x_all)
        total_n2 <- length(y_all)
        # 有效值
        x <- x_all[!is.na(x_all)]
        y <- y_all[!is.na(y_all)]
        n1 <- length(x)
        n2 <- length(y)
        na_count1 <- total_n1 - n1
        na_count2 <- total_n2 - n2
        # 统一结构（新增 4 列）
        res_row <- data.frame(
          gene       = gene,
          group1     = g1,
          group2     = g2,
          mean1      = NA_real_,
          mean2      = NA_real_,
          log2FC     = NA_real_,
          p_value    = NA_real_,
          fdr        = NA_real_,
          n1         = n1,         # 有效样本数
          n2         = n2,   # 有效样本数
          na_count1  = na_count1,
          na_count2  = na_count2,
          total_n1   = total_n1,
          total_n2   = total_n2,
          stringsAsFactors = FALSE
        )
        # 样本量不足 → 打印提示，保留 NA 行
        if (n1 < 2 || n2 < 2) {
          message(sprintf("Skipped %s: %s vs %s  (effective n1=%d, n2=%d | total_n1=%d, total_n2=%d)", 
                          gene, g1, g2, n1, n2, total_n1, total_n2))
          result_list[[idx]] <- res_row
          idx <- idx + 1
          next
        }
        # 正常计算
        mean1 <- mean(x, na.rm = TRUE)
        mean2 <- mean(y, na.rm = TRUE)
        pseudocount <- 0.1
        log2fc <- log2((mean1 + pseudocount) / (mean2 + pseudocount))
        test_result <- tryCatch(
          wilcox.test(x, y, exact = FALSE),
          error = function(e) {
            message("wilcox.test failed for ", gene, " ", g1, " vs ", g2, " - ", e$message)
            list(p.value = NA_real_)
          },
          warning = function(w) list(p.value = NA_real_)
        )
        res_row$mean1   <- mean1
        res_row$mean2   <- mean2
        res_row$log2FC  <- log2fc
        res_row$p_value <- test_result$p.value
        result_list[[idx]] <- res_row
        idx <- idx + 1
      }
    }
  }
  # 合并所有结果
  results <- do.call(rbind, result_list[!sapply(result_list, is.null)])
  # FDR 校正（只对有效 p 值）
  valid <- !is.na(results$p_value)
  if (any(valid)) {
    results$fdr[valid] <- p.adjust(results$p_value[valid], method = "fdr")
  }
  # 排序：fdr 优先，NA 放最后
  results <- results[order(results$fdr, results$p_value, na.last = TRUE), ]
  # 清理 log2FC
  results$log2FC[is.na(results$mean1) | is.na(results$mean2)] <- NA_real_
  rownames(results) <- NULL
  return(results)
}

result <- pairwise_wilcox_gene(merge_matrix, "Death28", gene_start_col = 10)  
write.csv(result, file="/data/wuwq/noise/DISEASE_ALL/app_disease/COMBAT/Death28_wilcox.csv", row.names = TRUE)




########  4.Outcome_score 为秩和检验分组
merge_matrix <- subset(matrix,Outcome != "NA")
merge_matrix$Outcome <- as.character(merge_matrix$Outcome)

# 自定义两两比较的函数
pairwise_wilcox_gene <- function(data, group_col, gene_start_col = 10) {
  # 基本检查
  if (!group_col %in% colnames(data)) {
    stop("group_col not found in data")
  }
  groups <- unique(data[[group_col]])
  groups <- groups[!is.na(groups)]
  if (length(groups) < 2) {
    stop("Need at least 2 groups to compare")
  }
  gene_cols <- gene_start_col:ncol(data)
  gene_names <- colnames(data)[gene_cols]
  # 预分配 list
  result_list <- vector("list", length(gene_names) * choose(length(groups), 2))
  idx <- 1
  # 为了效率，先把分组索引提前算好（可选，但大数据集有帮助）
  group_indices <- lapply(groups, function(g) which(data[[group_col]] == g))
  names(group_indices) <- groups
  for (gene in gene_names) {
    expr <- data[[gene]]
    for (i in 1:(length(groups)-1)) {
      for (j in (i+1):length(groups)) {
        g1 <- groups[i]
        g2 <- groups[j]
        idx1 <- group_indices[[g1]]
        idx2 <- group_indices[[g2]]
        x_all <- expr[idx1]
        y_all <- expr[idx2]
        # 总样本数（包含 NA）
        total_n1 <- length(x_all)
        total_n2 <- length(y_all)
        # 有效值
        x <- x_all[!is.na(x_all)]
        y <- y_all[!is.na(y_all)]
        n1 <- length(x)
        n2 <- length(y)
        na_count1 <- total_n1 - n1
        na_count2 <- total_n2 - n2
        # 统一结构（新增 4 列）
        res_row <- data.frame(
          gene       = gene,
          group1     = g1,
          group2     = g2,
          mean1      = NA_real_,
          mean2      = NA_real_,
          log2FC     = NA_real_,
          p_value    = NA_real_,
          fdr        = NA_real_,
          n1         = n1,         # 有效样本数
          n2         = n2,   # 有效样本数
          na_count1  = na_count1,
          na_count2  = na_count2,
          total_n1   = total_n1,
          total_n2   = total_n2,
          stringsAsFactors = FALSE
        )
        # 样本量不足 → 打印提示，保留 NA 行
        if (n1 < 2 || n2 < 2) {
          message(sprintf("Skipped %s: %s vs %s  (effective n1=%d, n2=%d | total_n1=%d, total_n2=%d)", 
                          gene, g1, g2, n1, n2, total_n1, total_n2))
          result_list[[idx]] <- res_row
          idx <- idx + 1
          next
        }
        # 正常计算
        mean1 <- mean(x, na.rm = TRUE)
        mean2 <- mean(y, na.rm = TRUE)
        pseudocount <- 0.1
        log2fc <- log2((mean1 + pseudocount) / (mean2 + pseudocount))
        test_result <- tryCatch(
          wilcox.test(x, y, exact = FALSE),
          error = function(e) {
            message("wilcox.test failed for ", gene, " ", g1, " vs ", g2, " - ", e$message)
            list(p.value = NA_real_)
          },
          warning = function(w) list(p.value = NA_real_)
        )
        res_row$mean1   <- mean1
        res_row$mean2   <- mean2
        res_row$log2FC  <- log2fc
        res_row$p_value <- test_result$p.value
        result_list[[idx]] <- res_row
        idx <- idx + 1
      }
    }
  }
  # 合并所有结果
  results <- do.call(rbind, result_list[!sapply(result_list, is.null)])
  # FDR 校正（只对有效 p 值）
  valid <- !is.na(results$p_value)
  if (any(valid)) {
    results$fdr[valid] <- p.adjust(results$p_value[valid], method = "fdr")
  }
  # 排序：fdr 优先，NA 放最后
  results <- results[order(results$fdr, results$p_value, na.last = TRUE), ]
  # 清理 log2FC
  results$log2FC[is.na(results$mean1) | is.na(results$mean2)] <- NA_real_
  rownames(results) <- NULL
  return(results)
}

result <- pairwise_wilcox_gene(merge_matrix, "Outcome", gene_start_col = 10)  
write.csv(result, file="/data/wuwq/noise/DISEASE_ALL/app_disease/COMBAT/Outcome_wilcox.csv", row.names = TRUE)










rm(list=ls())
############################################## （七） 尝试-随机森林在同一个图里
##### 读取和根据单变量逻辑回归得到的AUC降序top基因进入构建随机森林模型 5折交叉验证，100次重复   随机森林模型

# 直接加载准备好的矩阵 由单变量逻辑回归得到的AUC_top gene
loaded_objs <- readRDS("/data/wuwq/noise/DISEASE_ALL/app_disease/COMBAT/noi_exp_matrix.rds")
exp_matrix <- loaded_objs[[1]]
noi_matrix <- loaded_objs[[2]]
matrix <- loaded_objs[[3]]
rm(loaded_objs)
# 把 Outcome 转换为因子变量，并确保因子水平合法
noi_matrix$Outcome <- as.factor(noi_matrix$Outcome)
levels(noi_matrix$Outcome) <- make.names(levels(noi_matrix$Outcome))
exp_matrix$Outcome <- as.factor(exp_matrix$Outcome)
levels(exp_matrix$Outcome) <- make.names(levels(exp_matrix$Outcome))
matrix$Outcome <- as.factor(matrix$Outcome)
levels(matrix$Outcome) <- make.names(levels(matrix$Outcome))

# 读取单变量逻辑回归得到的p_value<0.05的基因及AUC
AUC_noi_sorted <- read.table("/data/wuwq/noise/DISEASE_ALL/app_disease/COMBAT/AUC_of_noi_sig.csv", sep = ",", row.names=1, header = TRUE, check.names = FALSE, quote = "\"")
AUC_exp_sorted <- read.table("/data/wuwq/noise/DISEASE_ALL/app_disease/COMBAT/AUC_of_exp_sig.csv", sep = ",", row.names=1, header = TRUE, check.names = FALSE, quote = "\"")

# ============================================================
# 3. 参数
# ============================================================
n_repeats <- 100
k_folds   <- 5
set.seed(1)
# ============================================================
# 4. 单轮 5-fold CV：返回完整预测向量（每个样本恰好被预测 1 次）
# ============================================================
run_one_cv <- function(data_matrix, var_names, model_type, seed = 1) {
  X <- as.data.frame(data_matrix[, var_names, drop = FALSE])
  y <- data_matrix$Outcome
  if (is.factor(y)) y_numeric <- as.numeric(y) - 1 else y_numeric <- as.numeric(y)
  y_factor <- factor(y_numeric, levels = c(0, 1))
  set.seed(seed)
  fold_ids <- createFolds(y_numeric, k = k_folds, list = TRUE, returnTrain = FALSE)
  pred_vec <- rep(NA_real_, nrow(X))
  for (f in seq_along(fold_ids)) {
    test_idx  <- fold_ids[[f]]
    train_idx <- setdiff(seq_len(nrow(X)), test_idx)
    if (length(unique(y_numeric[train_idx])) < 2) next
    tryCatch({
      if (model_type == "logistic") {
        train_df <- cbind(Outcome = y_numeric[train_idx], X[train_idx, , drop = FALSE])
        model <- glm(Outcome ~ ., family = binomial(link = "logit"), data = train_df)
        pred_vec[test_idx] <- predict(model, newdata = X[test_idx, , drop = FALSE], type = "response")
      } else if (model_type == "rf") {
        model <- randomForest(x = X[train_idx, , drop = FALSE],
                              y = y_factor[train_idx], ntree = 500)
        pred_vec[test_idx] <- predict(model, newdata = X[test_idx, , drop = FALSE], type = "prob")[, "1"]
      }
    }, error = function(e) {
      cat("  [错误] fold", f, ":", conditionMessage(e), "\n")
    })
  }
  list(pred = pred_vec, label = y_numeric)
}
# ============================================================
# 5. 100 次重复：只收集 AUC 值
# ============================================================
get_auc_repeats <- function(data_matrix, var_names, model_type, n_repeats = 100) {
  auc_values <- numeric(n_repeats)
  for (r in 1:n_repeats) {
    res <- run_one_cv(data_matrix, var_names, model_type, seed = r)
    valid <- !is.na(res$pred)
    if (sum(valid) < 10) { auc_values[r] <- NA; next }
    roc_r <- roc(res$label[valid], res$pred[valid], quiet = TRUE)
    auc_values[r] <- as.numeric(auc(roc_r))
    if (r %% 20 == 0) cat("  完成", r, "/", n_repeats, "\n")
  }
  valid_n <- sum(!is.na(auc_values))
  cat("  有效重复:", valid_n, "/", n_repeats, "\n")
  list(mean = mean(auc_values, na.rm = TRUE),
       sd   = sd(auc_values, na.rm = TRUE))
}



# 6.1. 修正后的方法：直接得到平均ROC值
# ============================================================
get_mean_roc_fixed <- function(data_matrix, var_names, model_type, n_repeats = 100) {
  tpr_mat <- NULL
  fpr_grid <- seq(0, 1, length.out = 15) # 统一 FPR 网格
  for (repeat_idx in 1:n_repeats) {
    res <- run_one_cv(data_matrix, var_names, model_type, seed = repeat_idx)
    valid <- !is.na(res$pred)
    if (sum(valid) < 10) next # 如果预测值不足，则跳过
    roc_obj <- roc(res$label[valid], res$pred[valid], quiet = TRUE)
    fpr_raw <- 1 - roc_obj$specificities
    tpr_raw <- roc_obj$sensitivities
    # 插值计算TPR，确保每次迭代的数据映射到统一网格
    tpr_interp <- approx(fpr_raw, tpr_raw, xout = fpr_grid, method = "linear")$y
    # 把插值的 TPR 存储到矩阵
    tpr_mat <- rbind(tpr_mat, tpr_interp)
    # 状态报告
    if (repeat_idx %% 20 == 0) {
      cat("完成", sprintf("%d/%d 次重复\n", repeat_idx, n_repeats))
    }
  }
  # 计算平均和标准差以获得稳健的 ROC 曲线
  mean_tpr <- colMeans(tpr_mat, na.rm = TRUE)
  sd_tpr <- apply(tpr_mat, 2, sd, na.rm = TRUE)
  cat(sprintf("有效重复: %d/%d\n", nrow(tpr_mat), n_repeats))
  return(data.frame(
    FPR = fpr_grid,
    TPR_mean = mean_tpr,
    TPR_sd = sd_tpr
  ))
}

# ============================================================
# 7.1. 运行：获取平均ROC 坐标 + 100 次 AUC 统计
# ============================================================
### noise_选top**个用于多变量回归
noi_significant_vars <- AUC_noi_sorted  %>% head(20)
significant_vars_name <- as.character(noi_significant_vars$variable)
# noi_matrix 
roc_noi_rf <- get_mean_roc_fixed(noi_matrix, significant_vars_name, "rf", 100)
auc_noi_rf <- get_auc_repeats(noi_matrix, significant_vars_name, "rf", n_repeats) #n_repeats

### expression_选top**个用于多变量回归
exp_significant_vars <- AUC_exp_sorted  %>% head(20)
exp_significant_vars_name <- as.character(exp_significant_vars$variable)
# exp_matrix 
roc_exp_rf <- get_mean_roc_fixed(exp_matrix, exp_significant_vars_name, "rf", 100)
auc_exp_rf <- get_auc_repeats(exp_matrix, exp_significant_vars_name, "rf", n_repeats)

### noise + expression_选top**个用于多变量回归
noi_significant_vars <- AUC_noi_sorted  %>% head(10)
exp_significant_vars <- AUC_exp_sorted  %>% head(10)
mix_significant_vars_name <- unique(c(as.character(noi_significant_vars$variable), as.character(exp_significant_vars$variable)))
# exp_matrix 
roc_mix_rf <- get_mean_roc_fixed(matrix, mix_significant_vars_name, "rf", 100)
auc_mix_rf <- get_auc_repeats(matrix, mix_significant_vars_name, "rf", n_repeats)

cat("Random Forest AUC:", 
    sprintf("%.4f", auc_noi_rf$mean), "±", sprintf("%.4f", auc_noi_rf$sd), "\n")
cat("Random Forest AUC:", 
    sprintf("%.4f", auc_exp_rf$mean), "±", sprintf("%.4f", auc_exp_rf$sd), "\n")
cat("Random Forest AUC:", 
    sprintf("%.4f", auc_mix_rf$mean), "±", sprintf("%.4f", auc_mix_rf$sd), "\n")

# ============================================================
# 8. 画图
# ============================================================
plot_roc_with_ci <- function(roc_df1, auc_stat1, label1, color1,
                             roc_df2, auc_stat2, label2, color2,
                             roc_df3, auc_stat3, label3, color3,
                             title_text, save_path = NULL) {
  # 添加模型列，用于区分不同的绘图数据
  roc_df1$Model <- sprintf("%s (AUC = %.4f ± %.4f)", label1, auc_stat1$mean, auc_stat1$sd)
  roc_df1$Lower <- roc_df1$TPR_mean - 1.96 * roc_df1$TPR_sd / sqrt(nrow(roc_df1)) # 下界 (Lower)
  roc_df1$Upper <- roc_df1$TPR_mean + 1.96 * roc_df1$TPR_sd / sqrt(nrow(roc_df1)) # 上界 (Upper)
  
  roc_df2$Model <- sprintf("%s (AUC = %.4f ± %.4f)", label2, auc_stat2$mean, auc_stat2$sd)
  roc_df2$Lower <- roc_df2$TPR_mean - 1.96 * roc_df2$TPR_sd / sqrt(nrow(roc_df2)) # 下界
  roc_df2$Upper <- roc_df2$TPR_mean + 1.96 * roc_df2$TPR_sd / sqrt(nrow(roc_df2)) # 上界

  roc_df3$Model <- sprintf("%s (AUC = %.4f ± %.4f)", label3, auc_stat3$mean, auc_stat3$sd)
  roc_df3$Lower <- roc_df3$TPR_mean - 1.96 * roc_df3$TPR_sd / sqrt(nrow(roc_df3)) # 下界
  roc_df3$Upper <- roc_df3$TPR_mean + 1.96 * roc_df3$TPR_sd / sqrt(nrow(roc_df3)) # 上界
  # 合并三个 ROC 数据框
  plot_df <- rbind(roc_df1, roc_df2, roc_df3)
  # 使用 ggplot 绘制
  p <- ggplot(plot_df, aes(x = FPR, y = TPR_mean, color = Model)) +
    # 绘制平均 ROC 曲线
    geom_line(size = 1.2) +
    # 绘制置信区间的阴影
    geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = Model), alpha = 0.2, color = NA) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
    scale_color_manual(values = setNames(c(color1, color2, color3),
                                         c(unique(roc_df1$Model), unique(roc_df2$Model), unique(roc_df3$Model)))) +
    scale_fill_manual(values = setNames(c(color1, color2, color3),
                                        c(unique(roc_df1$Model), unique(roc_df2$Model), unique(roc_df3$Model)))) +
    labs(
      title = title_text,
      x = "False Positive Rate (1 - Specificity)",
      y = "True Positive Rate (Sensitivity)",
      color = NULL,
      fill = NULL
    ) +
    coord_equal() +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      legend.position = c(0.62, 0.2),
      legend.background = element_rect(fill = "white", color = "grey80"),
      legend.text = element_text(size = 10)
    )
}
# ============================================================
# 调用修改后的函数，将置信区间绘制到图上
# ============================================================
p2 <- plot_roc_with_ci(
  roc_noi_rf, auc_noi_rf, "Noise", "green",
  roc_exp_rf, auc_exp_rf, "Expression", "#2563EB",
  roc_mix_rf, auc_mix_rf, "Mix", "#DC2626",
  "Outcome_1_0: ROC (5-Fold CV x 100 Repeats with 95% CI)"
)
# 显示图形
print(p2)

### 下载表格
write.csv(roc_noi_rf, file = "/data/wuwq/noise/DISEASE_ALL/app_disease/COMBAT/roc_noi_rf.csv")
write.csv(roc_exp_rf, file = "/data/wuwq/noise/DISEASE_ALL/app_disease/COMBAT/roc_exp_rf.csv")
write.csv(roc_mix_rf, file = "/data/wuwq/noise/DISEASE_ALL/app_disease/COMBAT/roc_mix_rf.csv")










rm(list=ls())
############################################## （八）出院Outcome COMBAT correted noise数据核对
##############################################❗️因基因名为列名❗️️，故 ✅️特殊读取和 ✅️特殊保存
# 构建noi_matrix (95个有Outcome信息的个体) 1=47, 0=48
noi_data <- read.table("/data/wuwq/noise/DISEASE_ALL/app_disease/merge_COMBAT.csv", sep = ",", row.names=1, header = TRUE, check.names = FALSE, quote = "\"")
unique(noi_data$Outcome)
# 保留有Outcome值的数据
noi_matrix <- subset(noi_data, Outcome != "NA")

unique(noi_matrix$Disease_group) #sepsis的都为1
unique(noi_matrix$Outcome)
column_index <- which(colnames(noi_matrix) == "ANXA1_MNP")
noi_matrix_f <- noi_matrix[, c(column_index, setdiff(1:ncol(noi_matrix), column_index))]


noi_matrix_f <- noi_matrix_f[,-c(3:ncol(noi_matrix_f))]

### 保存数据
write.csv(noi_matrix_f, file="/data/wuwq/noise/DISEASE_ALL/app_disease/COMBAT/ANXA1_MNP_individual.csv")

# 将NA用0填充
for (i in 3:ncol(noi_matrix)) {
  noi_matrix[[i]][is.na(noi_matrix[[i]])] <- 0
}






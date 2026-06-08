##################################################################
############ ------------ result_4_NSCLC -----------#############
############ ------      疾病免疫事件预后应用       -------#########
############ ------     NSCLC数据集挖掘       -------############ 
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

rm(list=ls())
############################################## （一）整合疾病二级基线信息--NSCLC---❗️因基因名为列名❗️️，故 ✅️特殊读取和 ✅️特殊保存
noise_disease <- read.table("/data/wuwq/noise/DISEASE_ALL/merge_disease_health_otherT.csv", sep = ",", row.names=1, header = TRUE, check.names = FALSE, quote = "\"")
noise_disease <- subset(noise_disease, dataset == "NSCLC")
#### 更新二级结局--NSCLC
baseline <- read.delim("/data/wuwq/noise/DISEASE_ALL/baseline_NSCLC.txt", header = TRUE, row.names = 1)
rownames(baseline) <- baseline$orig.ident.1
common_rows <- intersect(rownames(baseline), rownames(noise_disease))
noise_d <- data.frame(row.names = common_rows)
noise_d <- cbind(noise_d, baseline[common_rows, , drop = FALSE],  noise_disease[common_rows, , drop = FALSE] )
### 去除全部为NA的基因列
noise_d <- Filter(function(x) !all(is.na(x)), noise_d)
noise_d <- noise_d[,-c(5,7)]
write.table(noise_d,file = "/data/wuwq/noise/DISEASE_ALL/app_disease/otherT_NSCLC.csv", sep = ",", row.names = TRUE, col.names = NA, quote = TRUE)


rm(list=ls())
### merge细胞类型   特殊读取 ✅️
noise_disease <- read.table("/data/wuwq/noise/DISEASE_ALL/merge_disease_health_CD4.csv", sep = ",", row.names=1, header = TRUE, check.names = FALSE, quote = "\"")
noise_CD4 <- subset(noise_disease, dataset == "NSCLC")
noise_disease <- read.table("/data/wuwq/noise/DISEASE_ALL/merge_disease_health_CD8.csv", sep = ",", row.names=1, header = TRUE, check.names = FALSE, quote = "\"")
noise_CD8 <- subset(noise_disease, dataset == "NSCLC")
noise_disease <- read.table("/data/wuwq/noise/DISEASE_ALL/merge_disease_health_MNP.csv", sep = ",", row.names=1, header = TRUE, check.names = FALSE, quote = "\"")
noise_MNP <- subset(noise_disease, dataset == "NSCLC")
noise_disease <- read.table("/data/wuwq/noise/DISEASE_ALL/merge_disease_health_NK.csv", sep = ",", row.names=1, header = TRUE, check.names = FALSE, quote = "\"")
noise_NK <- subset(noise_disease, dataset == "NSCLC")
noise_disease <- read.table("/data/wuwq/noise/DISEASE_ALL/merge_disease_health_B.csv", sep = ",", row.names=1, header = TRUE, check.names = FALSE, quote = "\"")
noise_B <- subset(noise_disease, dataset == "NSCLC")
noise_disease <- read.table("/data/wuwq/noise/DISEASE_ALL/merge_disease_health_otherT.csv", sep = ",", row.names=1, header = TRUE, check.names = FALSE, quote = "\"")
noise_otherT <- subset(noise_disease, dataset == "NSCLC")
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
colnames(noise_all) <- gsub(pattern, "", colnames(noise_all))

# 更新二级结局--NSCLC
baseline <- read.delim("/data/wuwq/noise/DISEASE_ALL/baseline_NSCLC.txt", header = TRUE, row.names = 1)
baseline <- baseline[,-c(5,7)]
rownames(baseline) <- baseline$orig.ident.1
noise_disease <- noise_all

common_rows <- intersect(rownames(baseline), rownames(noise_disease))
noise_d <- data.frame(row.names = common_rows)
noise_d <- cbind(noise_d, baseline[common_rows, , drop = FALSE],  noise_disease[common_rows, , drop = FALSE] )
# 去除全部为NA的基因列
noise_d <- Filter(function(x) !all(is.na(x)), noise_d)
unique(noise_d$SubType_new) # "CR"    "AR"    "other" "PriR"
unique(noise_d$irAE_or_not) # "irAE" "not"
# 把corrected_noise列前移
# 找出以 "corrrected_noise_" 开头的列
corrected_noise_cols <- grep("^corrected_noise_", names(noise_d), value = TRUE)
# 找出非 "corrrected_noise_" 的列
other_cols <- setdiff(names(noise_d), corrected_noise_cols)
# 保留前6列的基线列
baseline_cols <- other_cols[1:6]
# 剩余的非 "corrrected_noise_" 列
remaining_cols <- other_cols[7:length(other_cols)]
# 重新排列列的顺序
new_order <- c(baseline_cols, corrected_noise_cols, remaining_cols)
# 重新排列数据框列
noise_all_reordered <- noise_d %>%
  select(all_of(new_order))
write.table(noise_all_reordered,file = "/data/wuwq/noise/DISEASE_ALL/app_disease/merge_NSCLC.csv", sep = ",", row.names = TRUE, col.names = NA, quote = TRUE)


rm(list=ls())
#### merge 基因表达量数据   ✅️特殊读取
exp_disease <- read.table("/data/wuwq/noise/DISEASE_ALL/exp_u/union_onlydisease_CD4.csv", sep = ",", row.names=1, header = TRUE, check.names = FALSE, quote = "\"")
exp_CD4 <- subset(exp_disease, dataset == "NSCLC")
exp_disease <- read.table("/data/wuwq/noise/DISEASE_ALL/exp_u/union_onlydisease_CD8.csv", sep = ",", row.names=1, header = TRUE, check.names = FALSE, quote = "\"")
exp_CD8 <- subset(exp_disease, dataset == "NSCLC")
exp_disease <- read.table("/data/wuwq/noise/DISEASE_ALL/exp_u/union_onlydisease_MNP.csv", sep = ",", row.names=1, header = TRUE, check.names = FALSE, quote = "\"")
exp_MNP <- subset(exp_disease, dataset == "NSCLC")
exp_disease <- read.table("/data/wuwq/noise/DISEASE_ALL/exp_u/union_onlydisease_NK.csv", sep = ",", row.names=1, header = TRUE, check.names = FALSE, quote = "\"")
exp_NK <- subset(exp_disease, dataset == "NSCLC")
exp_disease <- read.table("/data/wuwq/noise/DISEASE_ALL/exp_u/union_onlydisease_B.csv", sep = ",", row.names=1, header = TRUE, check.names = FALSE, quote = "\"")
exp_B <- subset(exp_disease, dataset == "NSCLC")
exp_disease <- read.table("/data/wuwq/noise/DISEASE_ALL/exp_u/union_onlydisease_otherT.csv", sep = ",", row.names=1, header = TRUE, check.names = FALSE, quote = "\"")
exp_otherT <- subset(exp_disease, dataset == "NSCLC")
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
# 去除全部为0的基因列
exp_all <- Filter(function(col) !all(col == 0), exp_all)
# 去掉列名中'细胞类型.' 的前缀
# 定义要去除的前缀列表
prefixes <- c("CD4.", "CD8.", "MNP.", "NK.", "B.", "otherT.")
# 使用正则表达式构造可以匹配这些前缀的模式，采用管道符 `|` 进行 "或" 匹配
pattern <- paste0("^(", paste(prefixes, collapse = "|"), ")")
# 使用 gsub 去除这些前缀
colnames(exp_all) <- gsub(pattern, "", colnames(exp_all))


# 更新二级结局--NSCLC
baseline <- read.delim("/data/wuwq/noise/DISEASE_ALL/baseline_NSCLC.txt", header = TRUE, row.names = 1)
baseline <- baseline[,-c(5,7)]
rownames(baseline) <- baseline$orig.ident.1
common_rows <- intersect(rownames(baseline), rownames(exp_all))
exp_d <- data.frame(row.names = common_rows)
exp_d <- cbind(exp_d, baseline[common_rows, , drop = FALSE],  exp_all[common_rows, , drop = FALSE] )
# 去除全部为NA的基因列
exp_d <- Filter(function(x) !all(is.na(x)), exp_d)
# 去除全部为0的基因列
exp_d <- Filter(function(x) !all(x == 0, na.rm = TRUE), exp_d)
unique(exp_d$SubType_new) # "CR"    "AR"    "other" "PriR" 
unique(exp_d$irAE_or_not) # "irAE" "not"
write.table(exp_d,file = "/data/wuwq/noise/DISEASE_ALL/app_disease/merge_exp_NSCLC.csv", sep = ",", row.names = TRUE, col.names = NA, quote = TRUE)



rm(list=ls())
############################################## （二）irAE_or_not  分组 2折1000次   预后分析
##############################################❗️因基因名为列名❗️️，故 ✅️特殊读取和 ✅️特殊保存
##### 合并成大的matrix
noi_data <- read.table("/data/wuwq/noise/DISEASE_ALL/app_disease/merge_NSCLC.csv", sep = ",", row.names=1, header = TRUE, check.names = FALSE, quote = "\"")
noi_matrix <- noi_data[,-c(1:2)]
exp_data <- read.table("/data/wuwq/noise/DISEASE_ALL/app_disease/merge_exp_NSCLC.csv", sep = ",", row.names=1, header = TRUE, check.names = FALSE, quote = "\"")
exp_matrix <- exp_data[, -c(1:2)]
# 为 noi_data 从第11列开始的列名添加前缀 "noi_"
colnames(noi_matrix)[-c(1:10)] <- paste0("noi_", colnames(noi_matrix)[-c(1:10)])
# 为 exp_data 从第5列开始的列名添加前缀 "exp_"
colnames(exp_matrix)[-c(1:4)] <- paste0("exp_", colnames(exp_matrix)[-c(1:4)])
matrix <- merge(noi_matrix, exp_matrix[,-c(2:4)], by = "orig.ident.1", all = TRUE)
# 将NA用0填充
for (i in 13:ncol(noi_matrix)) {
  noi_matrix[[i]][is.na(noi_matrix[[i]])] <- 0
}
for (i in 5:ncol(exp_matrix)) {
  exp_matrix[[i]][is.na(exp_matrix[[i]])] <- 0
}
for (i in 11:ncol(matrix)) {
  matrix[[i]][is.na(matrix[[i]])] <- 0
}
# 封装保存以上矩阵
objs_to_save <- list(exp_matrix, noi_matrix, matrix)
saveRDS(objs_to_save, file = "/data/wuwq/noise/DISEASE_ALL/app_disease/NSCLC/noi_exp_matrix.rds")





rm(list=ls())
# 直接加载准备好的矩阵
loaded_objs <- readRDS("/data/wuwq/noise/DISEASE_ALL/app_disease/NSCLC/noi_exp_matrix.rds")
exp_matrix <- loaded_objs[[1]]
noi_matrix <- loaded_objs[[2]]
matrix <- loaded_objs[[3]]
rm(loaded_objs)




##### 用irAE做结局 二分类结局  需要确保其是1 and 0
noi_matrix$irAE_or_not <- ifelse(noi_matrix$irAE_or_not == "irAE", 1, 0)
exp_matrix$irAE_or_not <- ifelse(exp_matrix$irAE_or_not == "irAE", 1, 0)
matrix$irAE_or_not <- ifelse(matrix$irAE_or_not == "irAE", 1, 0)
noi_matrix$irAE_or_not <- as.integer(noi_matrix$irAE_or_not)
exp_matrix$irAE_or_not <- as.integer(exp_matrix$irAE_or_not)
matrix$irAE_or_not <- as.integer(matrix$irAE_or_not)



#####  1.单变量初筛
### 1.1.noise
# 创建一个空的数据框，用于存储每个变量的p值
p_noi <- data.frame(variable = character(), p_value = numeric(), coeff = numeric(), stringsAsFactors = FALSE)
# 遍历 noi_matrix 中的每一个变量列，从第三列开始
for (var in colnames(noi_matrix)[c(11:ncol(noi_matrix))]) {
  # 选择当前变量的非缺失值数据
  temp_data <- data.frame(irAE_or_not = noi_matrix$irAE_or_not, variable = noi_matrix[, var])
  temp_data <- na.omit(temp_data)  # 去掉含有 NA 的行
  # 在过滤后的数据上构建逻辑回归模型
  if (nrow(temp_data) > 0) {  # 确保模型有数据可计算
    noi_bi_model <- glm(irAE_or_not ~ variable, data = temp_data, family = binomial)
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
      cat("Variable ", var, " does not have valid coefficients for irAE_or_not.\n")
    }
  }
}
#计算coeff_abs及fdr
p_noi$coeff_abs <- abs(p_noi$coeff)
p_noi$fdr <- p.adjust(p_noi$p_value,method = "fdr")
p_noi <- p_noi %>% filter(p_value < 0.05)


##### noise_单变量回归中p_value<0.05的基因 的 AUC计算
name <- as.character(p_noi$variable)
set.seed(1)
# 初始化AUC值的列表（包括每个基因的结果）  2折交叉验证，1000次重复抽样
auc_results <- data.frame(variable = character(), Mean_AUC = numeric(), stringsAsFactors = FALSE)
auc_values_list <- list()
folds <- createMultiFolds(y = noi_matrix$irAE_or_not, k = 2, times = 1000)
# 循环遍历每个基因名进行单因素逻辑回归
for (gene in name) {
  noi_auc_values <- numeric()
  # 进行交叉验证
  for (i in 1:length(folds)) {
    train <- noi_matrix[folds[[i]], c("irAE_or_not", gene)]
    test <- noi_matrix[-folds[[i]], c("irAE_or_not", gene)]
    noi_model <- glm(irAE_or_not ~ ., family = binomial(link = logit), data = train)
    noi_model_pre <- predict(noi_model, type = 'response', newdata = test)
    auc_value <- as.numeric(auc(as.numeric(test[, "irAE_or_not"]), noi_model_pre))
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
write.csv(AUC_noi_sorted, file="/data/wuwq/noise/DISEASE_ALL/app_disease/NSCLC/AUC_of_noi_sig.csv")
write.csv(auc_values_list, file="/data/wuwq/noise/DISEASE_ALL/app_disease/NSCLC/AUC_of_noi_sig_2000.csv")




### 1.2.expression
# 创建一个空的数据框，用于存储每个变量的p值
p_exp <- data.frame(variable = character(), p_value = numeric(), coeff = numeric(), stringsAsFactors = FALSE)
# 遍历 exp_matrix 中的每一个变量列，从第5列开始
for (var in colnames(exp_matrix)[c(5:ncol(exp_matrix))]) {
  # 选择当前变量的非缺失值数据
  temp_data <- data.frame(irAE_or_not = exp_matrix$irAE_or_not, variable = exp_matrix[, var])
  temp_data <- na.omit(temp_data)  # 去掉含有 NA 的行
  # 在过滤后的数据上构建逻辑回归模型
  if (nrow(temp_data) > 0) {  # 确保模型有数据可计算
    exp_bi_model <- glm(irAE_or_not ~ variable, data = temp_data, family = binomial)
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
      cat("Variable ", var, " does not have valid coefficients for irAE_or_not.\n")
    }
  }
}
#计算coeff_abs及fdr
p_exp$coeff_abs <- abs(p_exp$coeff)
p_exp$fdr <- p.adjust(p_exp$p_value,method = "fdr")
p_exp <- p_exp %>% filter(p_value < 0.05)

### expression_单变量回归中p_value<0.05的基因 计算AUC
name <- as.character(p_exp$variable)

name <- c("exp_C22ORF-34_CD4"
)

set.seed(1)
# 初始化AUC值的列表（包括每个基因的结果）  2折交叉验证，1000次重复抽样
auc_results <- data.frame(variable = character(), Mean_AUC = numeric(), stringsAsFactors = FALSE)
auc_values_list <- list()
folds <- createMultiFolds(y = exp_matrix$irAE_or_not, k = 2, times = 1000)
# 循环遍历每个基因名进行单因素逻辑回归
for (gene in name) {
  exp_auc_values <- numeric()
  # 进行交叉验证
  for (i in 1:length(folds)) {
    train <- exp_matrix[folds[[i]], c("irAE_or_not", gene)]
    test <- exp_matrix[-folds[[i]], c("irAE_or_not", gene)]
    exp_model <- glm(irAE_or_not ~ ., family = binomial(link = logit), data = train)
    exp_model_pre <- predict(exp_model, type = 'response', newdata = test)
    auc_value <- as.numeric(auc(as.numeric(test[, "irAE_or_not"]), exp_model_pre))
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
AUC_exp_sorted <- merge(auc_results, p_exp, by = "variable")
AUC_exp_sorted <- AUC_exp_sorted[order(AUC_exp_sorted$Mean_AUC, decreasing = TRUE), ]
#write.csv(AUC_exp_sorted, file="/data/wuwq/noise/DISEASE_ALL/app_disease/NSCLC/AUC_of_exp_sig.csv")
#write.csv(auc_values_list, file="/data/wuwq/noise/DISEASE_ALL/app_disease/NSCLC/AUC_of_exp_sig_2000.csv")
write.csv(auc_results, file="/data/wuwq/noise/DISEASE_ALL/app_disease/NSCLC/AUC_of_exp_matchnoise_v2.csv")



rm(list=ls())
##### 3.读取和根据AUC降序选top基因进入构建多变量回归模型 5折交叉验证，400次重复 
# 读取单变量回归得到的逻辑回归p_value<0.05的基因及AUC
AUC_noi_sorted <- read.table("/data/wuwq/noise/DISEASE_ALL/app_disease/NSCLC/AUC_of_noi_sig.csv", sep = ",", row.names=1, header = TRUE, check.names = FALSE, quote = "\"")
AUC_exp_sorted <- read.table("/data/wuwq/noise/DISEASE_ALL/app_disease/NSCLC/AUC_of_exp_sig.csv", sep = ",", row.names=1, header = TRUE, check.names = FALSE, quote = "\"")

# 直接加载准备好的矩阵
loaded_objs <- readRDS("/data/wuwq/noise/DISEASE_ALL/app_disease/NSCLC/noi_exp_matrix.rds")
exp_matrix <- loaded_objs[[1]]
noi_matrix <- loaded_objs[[2]]
matrix <- loaded_objs[[3]]
rm(loaded_objs)
##### 用irAE做结局 二分类结局  需要确保其是1 and 0
noi_matrix$irAE_or_not <- ifelse(noi_matrix$irAE_or_not == "irAE", 1, 0)
exp_matrix$irAE_or_not <- ifelse(exp_matrix$irAE_or_not == "irAE", 1, 0)
matrix$irAE_or_not <- ifelse(matrix$irAE_or_not == "irAE", 1, 0)
noi_matrix$irAE_or_not <- as.integer(noi_matrix$irAE_or_not)
exp_matrix$irAE_or_not <- as.integer(exp_matrix$irAE_or_not)
matrix$irAE_or_not <- as.integer(matrix$irAE_or_not)



### 3.1.noise_选top**个用于多变量回归
noi_significant_vars <- AUC_noi_sorted  %>% head(4)
significant_vars_name <- as.character(noi_significant_vars$variable)
# noi_matrix 
set.seed(1)
# Create folds for 2-fold cross-validation, repeated 1000 times
folds <- createMultiFolds(y = noi_matrix$irAE_or_not, k = 2, times = 1000)
# Initialize vectors for storing predictions and actual labels
all_predictions <- numeric()
all_labels <- numeric()
# Initialize list for storing AUC values for each repetition
auc_values_list <- numeric()
# Perform cross-validation 2000 times
for(i in 1:2000) {
  # Define train and test sets based on folds
  train <- noi_matrix[folds[[i]], c("irAE_or_not", significant_vars_name)]
  test <- noi_matrix[-folds[[i]], c("irAE_or_not", significant_vars_name)]
  # Fit the logistic regression model
  noi_model <- glm(irAE_or_not ~ ., family = binomial(link = logit), data = train)
  # Make predictions on the test set
  noi_model_pre <- predict(noi_model, type = 'response', newdata = test)
  # Store the predictions and actual labels
  all_predictions <- c(all_predictions, noi_model_pre)
  all_labels <- c(all_labels, as.numeric(test$irAE_or_not))
  # Calculate AUC for this fold
  roc_obj <- roc(as.numeric(test$irAE_or_not), noi_model_pre)
  auc_value <- auc(roc_obj)
  # Store AUC values
  auc_values_list <- c(auc_values_list, auc_value)
}
# Compute the mean AUC over all repetitions
mean_auc <- mean(auc_values_list)
# Output the average AUC value
cat("平均AUC:", round(mean_auc, 4), "\n")
# Calculate the ROC curve
roc_obj <- roc(all_labels, all_predictions)
# Prepare the ROC data frame for visualization
roc_df <- data.frame(
  FPR = 1 - roc_obj$specificities,
  TPR = roc_obj$sensitivities
)
# Plot the ROC curve
ggplot(roc_df, aes(x = FPR, y = TPR)) +
  geom_line(color = "blue", size = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  labs(
    title = paste("ROC Curve (AUC =", round(mean_auc, 4), ")"),
    x = "False Positive Rate (1 - Specificity)",
    y = "True Positive Rate (Sensitivity)"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))



### 3.2.expression_选top**个用于多变量回归
exp_significant_vars <- AUC_exp_sorted  %>% head(4)
significant_vars_name <- as.character(exp_significant_vars$variable)
# exp_matrix 
set.seed(1)
# Create folds for 2-fold cross-validation, repeated 1000 times
folds <- createMultiFolds(y = exp_matrix$irAE_or_not, k = 2, times = 1000)
# Initialize vectors for storing predictions and actual labels
all_predictions <- numeric()
all_labels <- numeric()
# Initialize list for storing AUC values for each repetition
auc_values_list <- numeric()
# Perform cross-validation 2000 times
for(i in 1:2000) {
  # Define train and test sets based on folds
  train <- exp_matrix[folds[[i]], c("irAE_or_not", significant_vars_name)]
  test <- exp_matrix[-folds[[i]], c("irAE_or_not", significant_vars_name)]
  # Fit the logistic regression model
  exp_model <- glm(irAE_or_not ~ ., family = binomial(link = logit), data = train)
  # Make predictions on the test set
  exp_model_pre <- predict(exp_model, type = 'response', newdata = test)
  # Store the predictions and actual labels
  all_predictions <- c(all_predictions, exp_model_pre)
  all_labels <- c(all_labels, as.numeric(test$irAE_or_not))
  # Calculate AUC for this fold
  roc_obj <- roc(as.numeric(test$irAE_or_not), exp_model_pre)
  auc_value <- auc(roc_obj)
  # Store AUC values
  auc_values_list <- c(auc_values_list, auc_value)
}
# Compute the mean AUC over all repetitions
mean_auc <- mean(auc_values_list)
# Output the average AUC value
cat("平均AUC:", round(mean_auc, 4), "\n")
# Calculate the ROC curve
roc_obj <- roc(all_labels, all_predictions)
# Prepare the ROC data frame for visualization
roc_df <- data.frame(
  FPR = 1 - roc_obj$specificities,
  TPR = roc_obj$sensitivities
)
# Plot the ROC curve
ggplot(roc_df, aes(x = FPR, y = TPR)) +
  geom_line(color = "blue", size = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  labs(
    title = paste("ROC Curve (AUC =", round(mean_auc, 4), ")"),
    x = "False Positive Rate (1 - Specificity)",
    y = "True Positive Rate (Sensitivity)"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
#dev.off()
#dev.new()


### 3.3.合并用于模型构建的noise基因和expression基因
noi_significant_vars <- AUC_noi_sorted  %>% head(2)
exp_significant_vars <- AUC_exp_sorted  %>% head(2)
significant_vars_name <- unique(c(as.character(noi_significant_vars$variable), as.character(exp_significant_vars$variable)))
# matrix 
set.seed(1)
# Create folds for 2-fold cross-validation, repeated 1000 times
folds <- createMultiFolds(y = matrix$irAE_or_not, k = 2, times = 1000)
# Initialize vectors for storing predictions and actual labels
all_predictions <- numeric()
all_labels <- numeric()
# Initialize list for storing AUC values for each repetition
auc_values_list <- numeric()
# Perform cross-validation 2000 times
for(i in 1:2000) {
  # Define train and test sets based on folds
  train <- matrix[folds[[i]], c("irAE_or_not", significant_vars_name)]
  test <- matrix[-folds[[i]], c("irAE_or_not", significant_vars_name)]
  # Fit the logistic regression model
  exp_model <- glm(irAE_or_not ~ ., family = binomial(link = logit), data = train)
  # Make predictions on the test set
  exp_model_pre <- predict(exp_model, type = 'response', newdata = test)
  # Store the predictions and actual labels
  all_predictions <- c(all_predictions, exp_model_pre)
  all_labels <- c(all_labels, as.numeric(test$irAE_or_not))
  # Calculate AUC for this fold
  roc_obj <- roc(as.numeric(test$irAE_or_not), exp_model_pre)
  auc_value <- auc(roc_obj)
  # Store AUC values
  auc_values_list <- c(auc_values_list, auc_value)
}
# Compute the mean AUC over all repetitions
mean_auc <- mean(auc_values_list)
# Output the average AUC value
cat("平均AUC:", round(mean_auc, 4), "\n")
# Calculate the ROC curve
roc_obj <- roc(all_labels, all_predictions)
# Prepare the ROC data frame for visualization
roc_df <- data.frame(
  FPR = 1 - roc_obj$specificities,
  TPR = roc_obj$sensitivities
)
# Plot the ROC curve
ggplot(roc_df, aes(x = FPR, y = TPR)) +
  geom_line(color = "blue", size = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  labs(
    title = paste("ROC Curve (AUC =", round(mean_auc, 4), ")"),
    x = "False Positive Rate (1 - Specificity)",
    y = "True Positive Rate (Sensitivity)"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))






rm(list=ls())
##### 4.读取和根据单变量逻辑回归得到的AUC降序top基因进入构建随机森林模型 5折交叉验证，10次重复      ✅️✅️✅️️   随机森林模型
# 读取单变量回归得到的逻辑回归p_value<0.05的基因及AUC
AUC_noi_sorted <- read.table("/data/wuwq/noise/DISEASE_ALL/app_disease/NSCLC/AUC_of_noi_sig.csv", sep = ",", row.names=1, header = TRUE, check.names = FALSE, quote = "\"")
AUC_exp_sorted <- read.table("/data/wuwq/noise/DISEASE_ALL/app_disease/NSCLC/AUC_of_exp_sig.csv", sep = ",", row.names=1, header = TRUE, check.names = FALSE, quote = "\"")

# 直接加载准备好的矩阵
loaded_objs <- readRDS("/data/wuwq/noise/DISEASE_ALL/app_disease/NSCLC/noi_exp_matrix.rds")
exp_matrix <- loaded_objs[[1]]
noi_matrix <- loaded_objs[[2]]
matrix <- loaded_objs[[3]]
rm(loaded_objs)
##### 用irAE做结局 二分类结局  需要确保其是1 and 0
noi_matrix$irAE_or_not <- ifelse(noi_matrix$irAE_or_not == "irAE", "yes", "no")
exp_matrix$irAE_or_not <- ifelse(exp_matrix$irAE_or_not == "irAE", "yes", "no")
matrix$irAE_or_not <- ifelse(matrix$irAE_or_not == "irAE", "yes", "no")
# 把 irAE_or_not 转换为因子变量，并确保因子水平合法
noi_matrix$irAE_or_not <- as.factor(noi_matrix$irAE_or_not)
levels(noi_matrix$irAE_or_not) <- make.names(levels(noi_matrix$irAE_or_not))
exp_matrix$irAE_or_not <- as.factor(exp_matrix$irAE_or_not)
levels(exp_matrix$irAE_or_not) <- make.names(levels(exp_matrix$irAE_or_not))
matrix$irAE_or_not <- as.factor(matrix$irAE_or_not)
levels(matrix$irAE_or_not) <- make.names(levels(matrix$irAE_or_not))




###################################### a. noise
### noise_选top**个用于随机森林
noi_significant_vars <- AUC_noi_sorted  %>% head(1)
significant_vars_name <- as.character(noi_significant_vars$variable)
# Set seed for reproducibility
set.seed(1)
# Create folds for 5-fold cross-validation, repeated 10 times, with stratified sampling
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
  irAE_or_not ~ ., 
  data = noi_matrix[, c("irAE_or_not", significant_vars_name)],
  method = "rf",
  metric = "ROC",
  trControl = cv_control
)
# Collect AUC values from cross-validation
all_predictions <- as.numeric(rf_model$pred$pred)
all_labels <- rf_model$pred$obs
# Calculate the Mean AUC across all folds and repetitions
levels(all_labels) #需要保证负类在前，否则要翻转
roc_rf <- roc(all_labels, all_predictions, levels = c("no", "yes"))
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
exp_significant_vars <- AUC_exp_sorted  %>% head(40)
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
  irAE_or_not ~ ., 
  data = exp_matrix[, c("irAE_or_not", significant_vars_name)],
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
noi_significant_vars <- AUC_noi_sorted  %>% head(20)
exp_significant_vars <- AUC_exp_sorted  %>% head(20)
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
  irAE_or_not ~ ., 
  data = matrix[, c("irAE_or_not", significant_vars_name)],
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




rm(list=ls())
############################################## （三）irAE or not 二分类变量做 随机森林模型
##############################################❗️因基因名为列名❗️️，故 ✅️特殊读取和 ✅️特殊保存
# 直接加载准备好的矩阵
loaded_objs <- readRDS("/data/wuwq/noise/DISEASE_ALL/app_disease/NSCLC/noi_exp_matrix.rds")
exp_matrix <- loaded_objs[[1]]
noi_matrix <- loaded_objs[[2]]
matrix <- loaded_objs[[3]]
rm(loaded_objs)
##### 用irAE做结局 二分类结局  需要确保其是1 and 0
noi_matrix$irAE_or_not <- ifelse(noi_matrix$irAE_or_not == "irAE", "yes", "no")
exp_matrix$irAE_or_not <- ifelse(exp_matrix$irAE_or_not == "irAE", "yes", "no")
matrix$irAE_or_not <- ifelse(matrix$irAE_or_not == "irAE", "yes", "no")
# 把 irAE_or_not 转换为因子变量，并确保因子水平合法
noi_matrix$irAE_or_not <- as.factor(noi_matrix$irAE_or_not)
levels(noi_matrix$irAE_or_not) <- make.names(levels(noi_matrix$irAE_or_not))
exp_matrix$irAE_or_not <- as.factor(exp_matrix$irAE_or_not)
levels(exp_matrix$irAE_or_not) <- make.names(levels(exp_matrix$irAE_or_not))
matrix$irAE_or_not <- as.factor(matrix$irAE_or_not)
levels(matrix$irAE_or_not) <- make.names(levels(matrix$irAE_or_not))



##### 1. 单变量初筛 -- 随机森林模型
### noise
# 初始化数据框，用于存储每个基因的AUC结果
auc_results <- data.frame(variable = character(), Mean_AUC = numeric(), Std_AUC = numeric(), SE_AUC = numeric(), stringsAsFactors = FALSE)
auc_values_list <- list()
# 分层 2 折交叉验证和 100 次重复
set.seed(1)
folds <- createFolds(y = noi_matrix$irAE_or_not, k = 2)
control <- trainControl(
  method = "repeatedcv", 
  number = 2, 
  repeats = 100, 
  classProbs = TRUE, 
  summaryFunction = twoClassSummary, 
  savePredictions = TRUE, 
  sampling = "up"  # 或者 "up" 进行采样
)
# 循环遍历每个变量进行单因素随机森林模型
for (var in colnames(noi_matrix)[3:ncol(noi_matrix)]) {
  # 过滤非缺失值数据
  temp_data <- noi_matrix[, c("irAE_or_not", var)]
  temp_data <- na.omit(temp_data)
  # 确认数据不为空
  if (nrow(temp_data) > 0) {
    # 将irAE_or_not转换为因子类型，确保随机森林可以处理
    temp_data$irAE_or_not <- as.factor(temp_data$irAE_or_not)
    # 使用caret的train函数进行交叉验证和模型训练
    rf_model <- train(irAE_or_not ~ ., data = temp_data, method = "rf", trControl = control, metric = "ROC")
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
write.csv(AUC_noi_sorted, file="/data/wuwq/noise/DISEASE_ALL/app_disease/NSCLC/RF_AUC_noi.csv")

### expression
# 初始化数据框，用于存储每个基因的AUC结果
# 初始化数据框，用于存储每个基因的AUC结果
auc_results <- data.frame(variable = character(), Mean_AUC = numeric(), Std_AUC = numeric(), SE_AUC = numeric(), stringsAsFactors = FALSE)
auc_values_list <- list()
# 分层 5 折交叉验证和 1 次重复
set.seed(1)
folds <- createMultiFolds(y = exp_matrix$irAE_or_not, k = 5, times = 1)
control <- trainControl(method = "repeatedcv", number = 5, repeats = 1, classProbs = TRUE, summaryFunction = twoClassSummary, savePredictions = TRUE)
# 循环遍历每个变量进行单因素随机森林模型
for (var in colnames(exp_matrix)[3:ncol(exp_matrix)]) {
  # 过滤非缺失值数据
  temp_data <- exp_matrix[, c("irAE_or_not", var)]
  temp_data <- na.omit(temp_data)
  # 确认数据不为空
  if (nrow(temp_data) > 0) {
    # 将irAE_or_not转换为因子类型，确保随机森林可以处理
    temp_data$irAE_or_not <- as.factor(temp_data$irAE_or_not)
    # 使用caret的train函数进行交叉验证和模型训练
    rf_model <- train(irAE_or_not ~ ., data = temp_data, method = "rf", trControl = control, metric = "ROC")
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
write.csv(AUC_exp_sorted, file="/data/wuwq/expse/DISEASE_ALL/app_disease/NSCLC/RF_AUC_exp.csv")





















rm(list=ls())
############################################## （四）irAE_or_not  分组 秩和检验
##############################################❗️因基因名为列名❗️️，故 ✅️特殊读取和 ✅️特殊保存
# 直接加载准备好的矩阵
loaded_objs <- readRDS("/data/wuwq/noise/DISEASE_ALL/app_disease/NSCLC/noi_exp_matrix.rds")
exp_matrix <- loaded_objs[[1]]
noi_matrix <- loaded_objs[[2]]
matrix <- loaded_objs[[3]]
rm(loaded_objs)

# 确定从基因列为数值型 & 比较组标签是字符串格式
merge_matrix <- matrix
unique(merge_matrix$irAE_or_not)

# 自定义两两比较的函数
pairwise_wilcox_gene <- function(data, group_col, gene_start_col = 5) {
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

# irAE_or_not 两两检验
result <- pairwise_wilcox_gene(merge_matrix, "irAE_or_not", gene_start_col = 5)  
write.csv(result, file="/data/wuwq/noise/DISEASE_ALL/app_disease/NSCLC/irAE_or_not_wilcox.csv", row.names = TRUE)

# irAE程度 两两检验
result1 <- pairwise_wilcox_gene(merge_matrix, "irAE", gene_start_col = 5)  
write.csv(result1, file="/data/wuwq/noise/DISEASE_ALL/app_disease/NSCLC/irAE_degree_wilcox.csv", row.names = TRUE)

# 耐药类别 两两检验
result2 <- pairwise_wilcox_gene(merge_matrix, "SubType_new", gene_start_col = 5)  
write.csv(result2, file="/data/wuwq/noise/DISEASE_ALL/app_disease/NSCLC/SubType_drug_degree_wilcox.csv", row.names = TRUE)




rm(list=ls())
############################################## （五）聚类分析
##############################################❗️因基因名为列名❗️️，故 ✅️特殊读取和 ✅️特殊保存
# 直接加载准备好的矩阵
loaded_objs <- readRDS("/data/wuwq/noise/DISEASE_ALL/app_disease/NSCLC/noi_exp_matrix.rds")
exp_matrix <- loaded_objs[[1]]
noi_matrix <- loaded_objs[[2]]
matrix <- loaded_objs[[3]]
rm(loaded_objs)

unique(noi_matrix$SubType_new)
noise_ind <- subset(noi_matrix, SubType_new != "other")
all_info <- colnames(noise_ind[,c(2,5:10)])
# 准备用于PCA的数据
valid_cols <- intersect(all_info, colnames(noise_ind))
pca_data <- noise_ind[, valid_cols, drop = FALSE]
#pca_data <- pca_data[, c("sample_id", setdiff(names(pca_data), "sample_id"))]
#pca_data <- pca_data[, c("Disease_group", setdiff(names(pca_data), "Disease_group"))]
# 将NA替换为0
data <- pca_data
data[is.na(data)] <- 0
data <- data[, colSums(data != 0) > 0]
# scale. = TRUE表示分析前对数据进行归一化；
com1 <- prcomp(data[ ,2:ncol(data)], center = TRUE,scale. = TRUE)
summ<-summary(com1)
df1<-com1$x
head(df1)
df1<-data.frame(df1,data$SubType_new)
df1$data.SubType_new<-as.factor(df1$data.SubType_new)
df1 <- df1[, c("data.SubType_new", setdiff(names(df1), "data.SubType_new"))]
lab1<-paste0("PC1(",round(summ$importance[2,1]*100,2),"%)")
lab2<-paste0("PC2(",round(summ$importance[2,2]*100,2),"%)")

# 创建PCA图形
Fig1a.taxa.pca <- ggplot(df1, aes(PC1, PC2)) +
  geom_point(size = 2, aes(color = data.SubType_new), show.legend = F) +
  scale_color_manual(values = c("#5686C3", "#75C500", "pink", "orange")) +
  stat_ellipse(aes(color = data.SubType_new), fill = "white", geom = "polygon",
               level = 0.95, alpha = 0.01, show.legend = F) +
  labs(x = lab1, y = lab2) +
  theme_classic() +
  theme(axis.line = element_line(colour = "black"),
        axis.title = element_text(color = "black", face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(color = "black", size = 10, face = "bold"))
Fig1a.taxa.pca
# 创建PC1密度图
Fig1a.taxa.pc1.density <- ggplot(df1) +
  geom_density(aes(x = PC1, group = data.SubType_new, fill = data.SubType_new),
               color = "black", alpha = 0.6, position = 'identity') +
  scale_fill_manual(values = c("#5686C3", "#75C500", "pink", "orange")) +
  theme_bw() +
  labs(fill = "")+
  theme(legend.position = "none")  # 隐藏图例
Fig1a.taxa.pc1.density
# 创建PC2密度图
Fig1a.taxa.pc2.density <- ggplot(df1) +
  geom_density(aes(x = PC2, group = data.SubType_new, fill = data.SubType_new),
               color = "black", alpha = 0.6, position = 'identity') +
  scale_fill_manual(values = c("#5686C3", "#75C500", "pink", "orange")) +
  theme_bw() +
  labs(fill = "") +
  coord_flip()+
  theme(axis.text.y = element_blank())
Fig1a.taxa.pc2.density 
# 使用 patchwork 合并图形
p1 <- (Fig1a.taxa.pca / Fig1a.taxa.pc2.density) +  # 使用'/'来定义上下排列
  Fig1a.taxa.pc1.density +  # 使用'+'来加入右侧的PC1密度图
  plot_layout(ncol = 2, heights = c(6/6, 2/6), widths = c(6/4, 2/6))  # 设置图形的布局
p1

#########################  kmeans (上接PCA)  ###############################
df0<-data[ , -1]
df2<-as.matrix(df0)
df2 <- scale(df2)






# 设置随机种子，让结果可以重现
set.seed(1)
# 计算WSS和WSS下降率
wss <- numeric(10)
wss_rate <- numeric(9)  # 计算 WSS 下降率的向量长度为 9，因为从 k=2 开始
for (k in 1:10) {
  kmeans_result <- kmeans(df2, centers = k, nstart = 25)
  wss[k] <- kmeans_result$tot.withinss
  if (k > 1) {
    wss_rate[k-1] <- (wss[k-1] - wss[k]) / wss[k-1]  # 计算WSS下降率
  }
}
# 可视化WSS
wss_df <- data.frame(k = 1:10, WSS = wss)
ggplot(wss_df, aes(x = k, y = WSS)) +
  geom_line() +
  geom_point() +
  ggtitle("WSS vs. k") +
  xlab("Number of Clusters (k)") +
  ylab("Within-cluster Sum of Squares (WSS)") +
  scale_x_continuous(breaks = seq(1, 10, by = 1), limits = c(1, 10))
#write.csv(wss_df,file="/data/wuwq/noise/DISEASE_ALL/app_disease/all_disease/k3_wss.csv")
# 可视化WSS下降率
wss_rate_df <- data.frame(k = 2:10, WSS_Rate = wss_rate)
ggplot(wss_rate_df, aes(x = k, y = WSS_Rate)) +
  geom_col() +  # 使用 geom_col() 绘制柱形图
  ggtitle("WSS下降率 vs. k") +
  xlab("Number of Clusters (k)") +
  ylab("WSS下降率") +
  scale_x_continuous(breaks = seq(2, 10, by = 1))  # 设置x轴间隔为1
#write.csv(wss_rate_df,file="/data/wuwq/noise/DISEASE_ALL/app_disease/all_disease/k3_wss_rate.csv")
# 计算FPC和FPC变化率
fpc <- numeric(10)
fpc_rate <- numeric(9)
# 在计算 FPC 变化率时，确保不处理 FPC 为零的聚类
for (k in 2:10) {
  kmeans_result <- kmeans(df2, centers = k, nstart = 25)
  # 计算轮廓系数作为 FPC
  silhouette_result <- silhouette(kmeans_result$cluster, dist(df2))
  fpc[k-1] <- mean(silhouette_result[, 3], na.rm = TRUE)
  # 计算 FPC 变化率，跳过 FPC 为 0 的情况
  if (k > 2 && fpc[k-1] > 0) {
    fpc_rate[k-2] <- abs((fpc[k-1] - fpc[k-2])/(fpc[k-2]))
  } else {
    fpc_rate[k-2] <- NA  # 如果 FPC 为 0，则不计算变化率
  }
}
# 可视化 FPC
fpc_df <- data.frame(k = 2:10, FPC = fpc[1:9])
# 绘制 FPC 图（跳过为 0 的 FPC）
ggplot(fpc_df, aes(x = k, y = FPC)) +
  geom_line() +
  geom_point() +
  labs(title = "FPC 随 k 值的变化", x = "k", y = "FPC") +
  scale_x_continuous(breaks = seq(2, 10, by = 1), limits = c(2, 10))
#write.csv(fpc_df,file="/data/wuwq/noise/DISEASE_ALL/app_disease/all_disease/k3_fpc.csv")
# 可视化 FPC 变化率
fpc_rate_df <- data.frame(k = 2:10, FPC_Rate = fpc_rate[1:9])
ggplot(fpc_rate_df, aes(x = k, y = FPC_Rate)) +
  geom_col() +
  geom_point()
#write.csv(fpc_rate_df,file="/data/wuwq/noise/DISEASE_ALL/app_disease/all_disease/k3_fpc_rate.csv")









# 设置随机种子，让结果可以重现
set.seed(1)
# 调用kmeans聚类算法 k = 2
km <- kmeans(df2, centers = 2, nstart = 25)
km$totss
### 把聚类cluster附加到PCA矩阵中：
final_data <- cbind(true_group = data$irAE_or_not, kgroup = km$cluster,PC1 = df1$PC1, PC2 = df1$PC2, PC3 = df1$PC3)
write.csv(final_data,file="/data/wuwq/noise/DISEASE_ALL/app_disease/NSCLC/k2_drug.csv")


#########################  基于kmeans cluster 结果的做 PCA 图 (上接kmeans) ###############################
#将kmeans 得到的cluster数据合并进开始的PCA结果；
df3<-data.frame(km$cluster,df1)
df3$km.cluster<-as.factor(df3$km.cluster)
summ<-summary(com1)
lab1<-paste0("PC1(",round(summ$importance[2,1]*100,2),"%)")
lab2<-paste0("PC2(",round(summ$importance[2,2]*100,2),"%)")
# 创建PCA图形
Fig1a.taxa.pca <- ggplot(df3, aes(PC1, PC2)) +
  geom_point(size = 2, aes(color = km.cluster), show.legend = F) +
  scale_color_manual(values = c("#5686C3", "#75C500", "pink", "orange")) +
  stat_ellipse(aes(color = km.cluster), fill = "white", geom = "polygon",
               level = 0.95, alpha = 0.01, show.legend = F) +
  labs(x = lab1, y = lab2) +
  theme_classic() +
  theme(axis.line = element_line(colour = "black"),
        axis.title = element_text(color = "black", face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(color = "black", size = 10, face = "bold"))
Fig1a.taxa.pca
# 创建PC1密度图
Fig1a.taxa.pc1.density <- ggplot(df3) +
  geom_density(aes(x = PC1, group = km.cluster, fill = km.cluster),
               color = "black", alpha = 0.6, position = 'identity') +
  scale_fill_manual(values = c("#5686C3", "#75C500", "pink", "orange")) +
  theme_bw() +
  labs(fill = "")+
  theme(legend.position = "none")  # 隐藏图例
Fig1a.taxa.pc1.density
# 创建PC2密度图
Fig1a.taxa.pc2.density <- ggplot(df3) +
  geom_density(aes(x = PC2, group = km.cluster, fill = km.cluster),
               color = "black", alpha = 0.6, position = 'identity') +
  scale_fill_manual(values = c("#5686C3", "#75C500", "pink", "orange")) +
  theme_bw() +
  labs(fill = "") +
  coord_flip()+
  theme(axis.text.y = element_blank())
Fig1a.taxa.pc2.density 
# 使用 patchwork 合并图形
p1 <- (Fig1a.taxa.pca / Fig1a.taxa.pc2.density) +  # 使用'/'来定义上下排列
  Fig1a.taxa.pc1.density +  # 使用'+'来加入右侧的PC1密度图
  plot_layout(ncol = 2, heights = c(6/6, 2/6), widths = c(6/4, 2/6))  # 设置图形的布局
p1









rm(list=ls())
############################################## （六）尝试 irAE or not 二分类变量做 随机森林模型
##############################################❗️因基因名为列名❗️️，故 ✅️特殊读取和 ✅️特殊保存
# 直接加载准备好的矩阵
loaded_objs <- readRDS("/data/wuwq/noise/DISEASE_ALL/app_disease/NSCLC/noi_exp_matrix.rds")
exp_matrix <- loaded_objs[[1]]
noi_matrix <- loaded_objs[[2]]
matrix <- loaded_objs[[3]]
rm(loaded_objs)
##### 用irAE做结局 二分类结局  需要确保其是1 and 0
noi_matrix$irAE_or_not <- ifelse(noi_matrix$irAE_or_not == "irAE", "yes", "no")
exp_matrix$irAE_or_not <- ifelse(exp_matrix$irAE_or_not == "irAE", "yes", "no")
matrix$irAE_or_not <- ifelse(matrix$irAE_or_not == "irAE", "yes", "no")
# 把 irAE_or_not 转换为因子变量，并确保因子水平合法
noi_matrix$irAE_or_not <- as.factor(noi_matrix$irAE_or_not)
levels(noi_matrix$irAE_or_not) <- make.names(levels(noi_matrix$irAE_or_not))
exp_matrix$irAE_or_not <- as.factor(exp_matrix$irAE_or_not)
levels(exp_matrix$irAE_or_not) <- make.names(levels(exp_matrix$irAE_or_not))
matrix$irAE_or_not <- as.factor(matrix$irAE_or_not)
levels(matrix$irAE_or_not) <- make.names(levels(matrix$irAE_or_not))
# 读取单变量逻辑回归得到的p_value<0.05的基因及AUC
# 读取单变量回归得到的逻辑回归p_value<0.05的基因及AUC
AUC_noi_sorted <- read.table("/data/wuwq/noise/DISEASE_ALL/app_disease/NSCLC/AUC_of_noi_sig.csv", sep = ",", row.names=1, header = TRUE, check.names = FALSE, quote = "\"")
AUC_exp_sorted <- read.table("/data/wuwq/noise/DISEASE_ALL/app_disease/NSCLC/AUC_of_exp_sig.csv", sep = ",", row.names=1, header = TRUE, check.names = FALSE, quote = "\"")

# ============================================================
# 3. 参数
# ============================================================
n_repeats <- 10
k_folds   <- 5
set.seed(1)
# ============================================================
# 4. 单轮 5-fold CV：返回完整预测向量（每个样本恰好被预测 1 次）
# ============================================================
run_one_cv <- function(data_matrix, var_names, model_type, seed = 1) {
  X <- as.data.frame(data_matrix[, var_names, drop = FALSE])
  y <- data_matrix$irAE_or_not
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
        train_df <- cbind(irAE_or_not = y_numeric[train_idx], X[train_idx, , drop = FALSE])
        model <- glm(irAE_or_not ~ ., family = binomial(link = "logit"), data = train_df)
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
get_auc_repeats <- function(data_matrix, var_names, model_type, n_repeats = 10) {
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
# ============================================================
# 6.1. 从一轮 CV 提取原始 ROC 坐标（保留阶梯形状）
# ============================================================
get_raw_roc <- function(data_matrix, var_names, model_type) {
  res <- run_one_cv(data_matrix, var_names, model_type, seed = 123)
  valid <- !is.na(res$pred)
  cat("  有效预测数:", sum(valid), "/", length(valid), "\n")
  roc_obj <- roc(res$label[valid], res$pred[valid], quiet = TRUE)
  # 提取原始坐标（不做任何插值）
  fpr_raw <- 1 - roc_obj$specificities
  tpr_raw <- roc_obj$sensitivities
  ord <- order(fpr_raw)
  data.frame(FPR = fpr_raw[ord], TPR = tpr_raw[ord])
}
# ============================================================
# 7.1. 运行：获取原始 ROC 坐标 + 100 次 AUC 统计
# ============================================================
### noise_选top**个用于多变量回归
noi_significant_vars <- AUC_noi_sorted  %>% head(44)
significant_vars_name <- as.character(noi_significant_vars$variable)
# noi_matrix 
roc_noi_lr <- get_raw_roc(noi_matrix, significant_vars_name, "logistic")
auc_noi_lr <- get_auc_repeats(noi_matrix, significant_vars_name, "logistic", n_repeats)
roc_noi_rf <- get_raw_roc(noi_matrix, significant_vars_name, "rf")
auc_noi_rf <- get_auc_repeats(noi_matrix, significant_vars_name, "rf", n_repeats)

### expression_选top**个用于多变量回归
exp_significant_vars <- AUC_exp_sorted  %>% head(44)
exp_significant_vars_name <- as.character(exp_significant_vars$variable)
# exp_matrix 
roc_exp_lr <- get_raw_roc(exp_matrix, exp_significant_vars_name, "logistic")
auc_exp_lr <- get_auc_repeats(exp_matrix, exp_significant_vars_name, "logistic", n_repeats)
roc_exp_rf <- get_raw_roc(exp_matrix, exp_significant_vars_name, "rf")
auc_exp_rf <- get_auc_repeats(exp_matrix, exp_significant_vars_name, "rf", n_repeats)

### noise + expression_选top**个用于多变量回归
noi_significant_vars <- AUC_noi_sorted  %>% head(22)
exp_significant_vars <- AUC_exp_sorted  %>% head(22)
mix_significant_vars_name <- unique(c(as.character(noi_significant_vars$variable), as.character(exp_significant_vars$variable)))
# exp_matrix 
roc_mix_lr <- get_raw_roc(matrix, mix_significant_vars_name, "logistic")
auc_mix_lr <- get_auc_repeats(matrix, mix_significant_vars_name, "logistic", n_repeats)
roc_mix_rf <- get_raw_roc(matrix, mix_significant_vars_name, "rf")
auc_mix_rf <- get_auc_repeats(matrix, mix_significant_vars_name, "rf", n_repeats)

cat("Random Forest AUC:", 
    sprintf("%.4f", auc_noi_rf$mean), "±", sprintf("%.4f", auc_noi_rf$sd), "\n")
cat("Random Forest AUC:", 
    sprintf("%.4f", auc_exp_rf$mean), "±", sprintf("%.4f", auc_exp_rf$sd), "\n")
cat("Random Forest AUC:", 
    sprintf("%.4f", auc_mix_rf$mean), "±", sprintf("%.4f", auc_mix_rf$sd), "\n")

cat("Logistic Regression AUC:", 
    sprintf("%.4f", auc_noi_lr$mean), "±", sprintf("%.4f", auc_noi_lr$sd), "\n")
cat("Logistic Regression AUC:", 
    sprintf("%.4f", auc_exp_lr$mean), "±", sprintf("%.4f", auc_exp_lr$sd), "\n")
cat("Logistic Regression AUC:", 
    sprintf("%.4f", auc_mix_lr$mean), "±", sprintf("%.4f", auc_mix_lr$sd), "\n")

# ============================================================
# 8. 画图
# ============================================================
plot_roc <- function(roc_df1, auc_stat1, label1, color1,
                     roc_df2, auc_stat2, label2, color2,
                     title_text, save_path = NULL) {
  roc_df1$Model <- sprintf("%s (AUC = %.4f ± %.4f)", label1, auc_stat1$mean, auc_stat1$sd)
  roc_df2$Model <- sprintf("%s (AUC = %.4f ± %.4f)", label2, auc_stat2$mean, auc_stat2$sd)
  plot_df <- rbind(roc_df1, roc_df2)
  p <- ggplot(plot_df, aes(x = FPR, y = TPR, color = Model)) +
    geom_line(size = 1) +               
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
    scale_color_manual(values = setNames(c(color1, color2),
                                         c(unique(roc_df1$Model), unique(roc_df2$Model)))) +
    labs(
      title = title_text,
      x = "False Positive Rate (1 - Specificity)",
      y = "True Positive Rate (Sensitivity)",
      color = NULL
    ) +
    coord_equal() +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      legend.position = c(0.62, 0.15),
      legend.background = element_rect(fill = "white", color = "grey80"),
      legend.text = element_text(size = 10)
    )
  if (!is.null(save_path)) {
    ggsave(save_path, p, width = 7, height = 7, dpi = 300)
    cat("已保存:", save_path, "\n")
  }
  print(p)
  return(p)
}

p1 <- plot_roc(
  roc_noi_lr, auc_noi_lr, "Logistic Regression", "#2563EB",
  roc_noi_rf, auc_noi_rf, "Random Forest",       "#DC2626",
  "Noise Data: ROC (5-Fold CV x 100 Repeats)"
  #"/data/wuwq/noise/DISEASE_ALL/app_disease/NSCLC/ROC_noise_corrected.pdf"
)
p1
p2 <- plot_roc(
  roc_exp_lr, auc_exp_lr, "Logistic Regression", "#2563EB",
  roc_exp_rf, auc_exp_rf, "Random Forest",       "#DC2626",
  "Expression Data: ROC (5-Fold CV x 100 Repeats)"
  #"/data/wuwq/noise/DISEASE_ALL/app_disease/NSCLC/ROC_expression_corrected.pdf"
)
p2
p3 <- plot_roc(
  roc_mix_lr, auc_mix_lr, "Logistic Regression", "#2563EB",
  roc_mix_rf, auc_mix_rf, "Random Forest",       "#DC2626",
  "Expression Data: ROC (5-Fold CV x 100 Repeats)"
  #"/data/wuwq/noise/DISEASE_ALL/app_disease/NSCLC/ROC_mix_corrected.pdf"
)
p3





rm(list=ls())
############################################## （七）终版 irAE or not 二分类变量做 随机森林模型 ✅️✅️✅️
# 直接加载准备好的矩阵
loaded_objs <- readRDS("/data/wuwq/noise/DISEASE_ALL/app_disease/NSCLC/noi_exp_matrix.rds")
exp_matrix <- loaded_objs[[1]]
noi_matrix <- loaded_objs[[2]]
matrix <- loaded_objs[[3]]
rm(loaded_objs)
##### 用irAE做结局 二分类结局  需要确保其是1 and 0
noi_matrix$irAE_or_not <- ifelse(noi_matrix$irAE_or_not == "irAE", "yes", "no")
exp_matrix$irAE_or_not <- ifelse(exp_matrix$irAE_or_not == "irAE", "yes", "no")
matrix$irAE_or_not <- ifelse(matrix$irAE_or_not == "irAE", "yes", "no")
# 把 irAE_or_not 转换为因子变量，并确保因子水平合法
noi_matrix$irAE_or_not <- as.factor(noi_matrix$irAE_or_not)
levels(noi_matrix$irAE_or_not) <- make.names(levels(noi_matrix$irAE_or_not))
exp_matrix$irAE_or_not <- as.factor(exp_matrix$irAE_or_not)
levels(exp_matrix$irAE_or_not) <- make.names(levels(exp_matrix$irAE_or_not))
matrix$irAE_or_not <- as.factor(matrix$irAE_or_not)
levels(matrix$irAE_or_not) <- make.names(levels(matrix$irAE_or_not))
# 读取单变量逻辑回归得到的p_value<0.05的基因及AUC
# 读取单变量回归得到的逻辑回归p_value<0.05的基因及AUC
AUC_noi_sorted <- read.table("/data/wuwq/noise/DISEASE_ALL/app_disease/NSCLC/AUC_of_noi_sig.csv", sep = ",", row.names=1, header = TRUE, check.names = FALSE, quote = "\"")
AUC_exp_sorted <- read.table("/data/wuwq/noise/DISEASE_ALL/app_disease/NSCLC/AUC_of_exp_sig.csv", sep = ",", row.names=1, header = TRUE, check.names = FALSE, quote = "\"")



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
  y <- data_matrix$irAE_or_not
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
        train_df <- cbind(irAE_or_not = y_numeric[train_idx], X[train_idx, , drop = FALSE])
        model <- glm(irAE_or_not ~ ., family = binomial(link = "logit"), data = train_df)
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
  fpr_grid <- seq(0, 1, length.out = 300) # 统一 FPR 网格
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
noi_significant_vars <- AUC_noi_sorted  %>% head(8)
significant_vars_name <- as.character(noi_significant_vars$variable)
# noi_matrix 
roc_noi_rf <- get_mean_roc_fixed(noi_matrix, significant_vars_name, "rf", n_repeats)
auc_noi_rf <- get_auc_repeats(noi_matrix, significant_vars_name, "rf", n_repeats)

### expression_选top**个用于多变量回归
exp_significant_vars <- AUC_exp_sorted  %>% head(8)
exp_significant_vars_name <- as.character(exp_significant_vars$variable)
# exp_matrix 
roc_exp_rf <- get_mean_roc_fixed(exp_matrix, exp_significant_vars_name, "rf", n_repeats)
auc_exp_rf <- get_auc_repeats(exp_matrix, exp_significant_vars_name, "rf", n_repeats)

### noise + expression_选top**个用于多变量回归
noi_significant_vars <- AUC_noi_sorted  %>% head(4)
exp_significant_vars <- AUC_exp_sorted  %>% head(4)
mix_significant_vars_name <- unique(c(as.character(noi_significant_vars$variable), as.character(exp_significant_vars$variable)))
# exp_matrix 
roc_mix_rf <- get_mean_roc_fixed(matrix, mix_significant_vars_name, "rf", n_repeats)
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
p1 <- plot_roc_with_ci(
  roc_noi_rf, auc_noi_rf, "Noise", "green",
  roc_exp_rf, auc_exp_rf, "Expression", "#2563EB",
  roc_mix_rf, auc_mix_rf, "Mix", "#DC2626",
  "irAE_or_not: ROC (5-Fold CV x 100 Repeats with 95% CI)"
)
# 显示图形
print(p1)

### 下载表格
write.csv(roc_noi_rf, file = "/data/wuwq/noise/DISEASE_ALL/app_disease/NSCLC/roc_noi_rf.csv")
write.csv(roc_exp_rf, file = "/data/wuwq/noise/DISEASE_ALL/app_disease/NSCLC/roc_exp_rf.csv")
write.csv(roc_mix_rf, file = "/data/wuwq/noise/DISEASE_ALL/app_disease/NSCLC/roc_mix_rf.csv")




rm(list=ls())
############################################## （八）irAR_ornot correted noise数据核对
##############################################❗️因基因名为列名❗️️，故 ✅️特殊读取和 ✅️特殊保存
# 构建noi_matrix (95个有Outcome信息的个体) 1=47, 0=48
loaded_objs <- readRDS("/data/wuwq/noise/DISEASE_ALL/app_disease/NSCLC/noi_exp_matrix.rds")
noi_matrix <- loaded_objs[[2]]
rm(loaded_objs)

unique(noi_matrix$irAE_or_not) 
noi_matrix <- noi_matrix[,-c(2,4,11:ncol(noi_matrix))]

### 保存数据
write.csv(noi_matrix, file="/data/wuwq/noise/DISEASE_ALL/app_disease/NSCLC/irAE_individual.csv")

# 将NA用0填充
for (i in 3:ncol(noi_matrix)) {
  noi_matrix[[i]][is.na(noi_matrix[[i]])] <- 0
}




rm(list=ls())
############################################## （八）irAR_ornot correted noise数据核对
##############################################❗️因基因名为列名❗️️，故 ✅️特殊读取和 ✅️特殊保存
# 构建noi_matrix (95个有Outcome信息的个体) 1=47, 0=48
loaded_objs <- readRDS("/data/wuwq/noise/DISEASE_ALL/app_disease/NSCLC/noi_exp_matrix.rds")
noi_matrix <- loaded_objs[[2]]
rm(loaded_objs)

unique(noi_matrix$irAE_or_not) 
column_index <- which(colnames(noi_matrix) == "noi_LINC00513_CD4")
noi_matrix_f <- noi_matrix[, c(column_index, setdiff(1:ncol(noi_matrix), column_index))]
noi_matrix_f <- noi_matrix_f[,-c(3:ncol(noi_matrix_f))]


### 保存数据
write.csv(noi_matrix_f, file="/data/wuwq/noise/DISEASE_ALL/app_disease/NSCLC/irAE_noi_LINC00513_CD4.csv")

# 将NA用0填充
for (i in 3:ncol(noi_matrix)) {
  noi_matrix[[i]][is.na(noi_matrix[[i]])] <- 0
}



##################################################################
############ ------------ result_4-2_NSCLC -----------#############
############ ------      疾病耐药预后应用       -------#########
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
# 直接加载准备好的矩阵
loaded_objs <- readRDS("/data/wuwq/noise/DISEASE_ALL/app_disease/NSCLC/noi_exp_matrix.rds")
exp_matrix <- loaded_objs[[1]]
noi_matrix <- loaded_objs[[2]]
matrix <- loaded_objs[[3]]
rm(loaded_objs)




##### 用耐药做结局 二分类结局  需要确保其  PriR & AR = 1 , CR=0
unique(noi_matrix$SubType_new)
noi_matrix$SubType_new <- ifelse(noi_matrix$irAE_or_not == "irAE", 1, 0)
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

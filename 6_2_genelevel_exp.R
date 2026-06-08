##################################################################
### --------------------- 基因表达量 秩和检验 -----------------###  ✅️ 把NA变成0后做检验
##################################################################
# 将所有NA变成0
# exp_u <- exp_u %>% mutate(across(6:last_col(), ~ replace_na(., 0)))
# sum(is.na(exp_u[6:ncol(exp_u)]))

# 将所有0变成NA
# exp_u <- exp_u %>% mutate(across(6:last_col(), ~ na_if(., 0)))
# sum(is.na(exp_u[6:ncol(exp_u)]))

library(dplyr)
library(knitr)
library(tidyr)
library(stats)
library(tibble)



rm(list=ls())
##### （一）control vs disease ✅️疾病与自带对照
# 1. 读取merge好的全数据集  # 将所有NA变成0
exp_u <- read.table("/data/wuwq/noise/DISEASE_ALL/exp_u/union_controldisease_NK.csv", sep = ",", row.names=1, header = TRUE, check.names = FALSE, quote = "\"")
sum(is.na(exp_u[6:ncol(exp_u)]))
exp_u <- exp_u %>% mutate(across(6:last_col(), ~ replace_na(., 0)))
sum(is.na(exp_u[6:ncol(exp_u)]))

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

# 使用自定义的函数进行两两比较
# 1.秩和检验-Disease_group
result_disease <- pairwise_wilcox_gene(exp_u, "Disease_group", gene_start_col = 6)  
write.csv(result_disease, file="/data/wuwq/noise/DISEASE_ALL/exp_u_wilcox/NK_controldisease_wilcox.csv", row.names = TRUE)


rm(list=ls())
##### （二）health vs disease  ✅️疾病与图谱
# 1. 读取merge好的全数据集  # 将所有NA变成0
exp_u <- read.table("/data/wuwq/noise/DISEASE_ALL/exp_u/union_healthdisease_NK.csv", sep = ",", row.names=1, header = TRUE, check.names = FALSE, quote = "\"")
sum(is.na(exp_u[5:ncol(exp_u)]))
exp_u <- exp_u %>% mutate(across(5:last_col(), ~ replace_na(., 0)))
sum(is.na(exp_u[5:ncol(exp_u)]))

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

# 使用自定义的函数进行两两比较
# 1.秩和检验-Disease_group
result_disease <- pairwise_wilcox_gene(exp_u, "Disease_group", gene_start_col = 5)  
write.csv(result_disease, file="/data/wuwq/noise/DISEASE_ALL/exp_u_wilcox/NK_healthdisease_wilcox.csv", row.names = TRUE)




rm(list=ls())
##### （三）health vs disease ✅️ 读取秩和检验的结果，得到疾病不变的表达基因
# 规则为：表达量不变的定义(数量收紧，尽量少)：-1<log2foldchange<1 或 fdr>0.05
# 处理过程中去掉在其中一组中均值为0的基因，即分析基因需要满足在每组中均有表达
### 读取秩和检验表格
exp_wilcox <- read.csv("/data/wuwq/noise/DISEASE_ALL/exp_u_wilcox/otherT_healthdisease_wilcox.csv", header = TRUE,row.names = 1)
unique(exp_wilcox$group1)
unique(exp_wilcox$group2)
# 筛选满足条件的数据
filtered_data <- exp_wilcox[ 
  # group1 中有表达  mean1 > 0
  (exp_wilcox$mean1 > 0) & 
  # group2 中有表达  mean2 > 0
  (exp_wilcox$mean2 > 0) &
  # group1 或 group2 包含 "health"
  (exp_wilcox$group1 == "health" | exp_wilcox$group2 == "health") & 
    # 满足 log2fc 的范围，或 fdr 的条件
  (exp_wilcox$log2FC > -1 & exp_wilcox$log2FC < 1 | exp_wilcox$fdr > 0.05),
]
write.csv(filtered_data, file="/data/wuwq/noise/DISEASE_ALL/exp_u_wilcox/filt_otherT_healthdisease_wilcox.csv", row.names = TRUE)




rm(list=ls())
##### （四）control vs disease ✅️ 读取秩和检验的结果，得到疾病不变的表达基因   control中为“Health”   而上面图谱中为“health”
# 规则为：表达量不变的定义(数量收紧，尽量少)：-1<log2foldchange<1 或 fdr>0.05
# 处理过程中去掉在其中一组中均值为0的基因，即分析基因需要满足在每组中均有表达
### 读取秩和检验表格
exp_wilcox <- read.csv("/data/wuwq/noise/DISEASE_ALL/exp_u_wilcox/otherT_controldisease_wilcox.csv", header = TRUE,row.names = 1)
unique(exp_wilcox$group1)
unique(exp_wilcox$group2)
# 筛选满足条件的数据
filtered_data <- exp_wilcox[ 
  # group1 中有表达  mean1 > 0
  (exp_wilcox$mean1 > 0) & 
    # group2 中有表达  mean2 > 0
    (exp_wilcox$mean2 > 0) &
    # group1 或 group2 包含 "Health"
    (exp_wilcox$group1 == "Health" | exp_wilcox$group2 == "Health") & 
    # 满足 log2fc 的范围，或 fdr 的条件
    (exp_wilcox$log2FC > -1 & exp_wilcox$log2FC < 1 | exp_wilcox$fdr > 0.05),
]
write.csv(filtered_data, file="/data/wuwq/noise/DISEASE_ALL/exp_u_wilcox/filt_otherT_controldisease_wilcox.csv", row.names = TRUE)



rm(list=ls())
##### （五）取 表达量不变的基因中，control 组 和 图谱组  基因名的交集  ✅️ exp水平不变的基因 
# 换细胞类型
health_disease<-read.csv("/data/wuwq/noise/DISEASE_ALL/exp_u_wilcox/filt_otherT_healthdisease_wilcox.csv", header = TRUE,row.names = 1)
control_disease<-read.csv("/data/wuwq/noise/DISEASE_ALL/exp_u_wilcox/filt_otherT_controldisease_wilcox.csv", header = TRUE,row.names = 1)

# 换疾病 "health" "RA" "SLE" "COVID-19" "post-COVID-19 disorder" "SEPSIS" "FLU" "NSCLC"
unique(health_disease$group1)
unique(health_disease$group2)
unique(control_disease$group1)
unique(control_disease$group2)
filt_health_disease <- health_disease[ (health_disease$group1 == "NSCLC") ,]
filt_control_disease <- control_disease[ (control_disease$group2 == "NSCLC") ,]

# 获取基因名交集
gene_health <- unique(filt_health_disease$gene)
gene_control <- unique(filt_control_disease$gene)
common_genes <- intersect(gene_health, gene_control)
# 保存表达量不变的基因名
write.csv(common_genes, file="/data/wuwq/noise/DISEASE_ALL/exp_u_wilcox/nochange_exp_common_gene/otherT_NSCLC.csv", row.names = TRUE)





rm(list=ls())
##### （六）COVID vs postCOVID ✅️ 读取秩和检验的结果，得到疾病不变的表达基因
# 规则为：表达量不变的定义(数量收紧，尽量少)：-1<log2foldchange<1 或 fdr>0.05
# 处理过程中去掉在其中一组中均值为0的基因，即分析基因需要满足在每组中均有表达
### 读取秩和检验表格
exp_wilcox <- read.csv("/data/wuwq/noise/DISEASE_ALL/exp_u_wilcox/otherT_healthdisease_wilcox.csv", header = TRUE,row.names = 1)
unique(exp_wilcox$group1)
unique(exp_wilcox$group2)
# 筛选满足条件的数据
filtered_data <- exp_wilcox[ 
  # group1 中有表达  mean1 > 0
  (exp_wilcox$mean1 > 0) & 
    # group2 中有表达  mean2 > 0
    (exp_wilcox$mean2 > 0) &
    # group1 或 group2 包含 "health"
    (exp_wilcox$group1 == "COVID-19" | exp_wilcox$group2 == "post-COVID-19 disorder") & 
    # 满足 log2fc 的范围，或 fdr 的条件
    (exp_wilcox$log2FC > -1 & exp_wilcox$log2FC < 1 | exp_wilcox$fdr > 0.05),
]
write.csv(filtered_data, file="/data/wuwq/noise/DISEASE_ALL/exp_u_wilcox/nochange_exp_covid_postcovid/otherT_nochange_exp.csv", row.names = TRUE)







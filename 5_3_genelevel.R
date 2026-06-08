
##################################################################
### ----- (一)基因水平noise 卡方检验 描绘nonNA-noise分布------ ### ✅️
##################################################################
### CD4+ T cell
### CD8+ T cell
### Monocyte
### NK cell
### Dendritic cell
### B cell
### Other T cell
library(dplyr)
library(knitr)
library(tidyr)
library(stats)
library(tibble)


rm(list=ls())
# 读取merge好的全数据集
noise <- read.csv("/data/wuwq/noise/DISEASE_ALL/merge_disease_health_otherT.csv", header = FALSE, row.names = 1)
colnames(noise) <- as.character(noise[1, ])
noise <- noise[-1, , drop = FALSE]
# 从第6列到最后一列变为数值型
noise[, 6:ncol(noise)] <- lapply(noise[, 6:ncol(noise)], function(x) {
  x <- as.numeric(x)
  return(x)
})

# 将所有0变成NA
# sum(is.na(noise[6:ncol(noise)]))
# noise <- noise %>% mutate(across(6:last_col(), ~ na_if(., 0)))
# sum(is.na(noise[6:ncol(noise)]))

# 安全的成对比较函数 —— 以 非NA 比例 作为“表达”定义
pairwise_nonNA_safe <- function(data, group_col, gene_cols) {
  # 获取所有分组名称
  groups <- unique(data[[group_col]])
  result_list <- list()
  # 对每个基因逐一处理
  for (gene in gene_cols) {
    # 初始化该基因的结果数据框
    gene_res <- data.frame(
      gene = character(),
      group1 = character(), group2 = character(),
      OR = numeric(), log2_OR = numeric(),
      p_value = numeric(), fdr = numeric(),
      n_group1 = integer(), n_group2 = integer(),
      pct1_nonNA = numeric(), pct2_nonNA = numeric(),
      n_nonNA_group1 = integer(), n_nonNA_group2 = integer(),
      stringsAsFactors = FALSE
    )
    # 用于收集该基因所有比较的原始 p 值（后面统一 FDR）
    gene_p <- numeric()
    # 预先计算每组的统计量（重点关注 non-NA 信息）
    stats <- lapply(groups, function(g) {
      vals <- data[data[[group_col]] == g, gene]
      n <- length(vals)
      if (n == 0) return(NULL)
      n_nonNA <- sum(!is.na(vals))
      list(
        n          = n,                       # 总细胞/样本数
        n_nonNA    = n_nonNA,                 # 非NA数量
        pct_nonNA  = n_nonNA / n,             # 非NA比例（核心指标）
        n_NA       = n - n_nonNA              # NA数量
      )
    })
    names(stats) <- groups
    # 两两比较（只做 i < j，避免重复）
    for (i in seq_along(groups)) {
      for (j in (i+1):length(groups)) {
        g1 <- groups[i]; g2 <- groups[j]
        s1 <- stats[[g1]]; s2 <- stats[[g2]]
        if (is.null(s1) || is.null(s2)) next
        # 构建 2×2 列联表：行=组，列=是否非NA（有值 vs NA）
        tab <- matrix(c(
          s1$n_nonNA,   s1$n_NA,     # group1: 非NA / NA
          s2$n_nonNA,   s2$n_NA      # group2: 非NA / NA
        ), nrow = 2, byrow = TRUE,
        dimnames = list(c(g1, g2), c("nonNA", "NA"))
        )
        # 计算 OR（以 group1 / group2 为方向）
        # 加 0.5 连续性校正，防止 0 或 Inf
        OR_adj <- (tab[1,1] + 0.5) / (tab[1,2] + 0.5) /
          (tab[2,1] + 0.5) / (tab[2,2] + 0.5)
        # 处理极端 log2OR
        log2OR_adj <- if (OR_adj == 0) -Inf else if (is.infinite(OR_adj)) Inf else log2(OR_adj)
        # 决定用 Fisher 还是 chisq
        use_fisher <- (s1$n < 40 || s2$n < 40) || any(tab < 5)
        if (use_fisher) {
          test_res <- tryCatch(fisher.test(tab), error = function(e) list(p.value = NA))
        } else {
          test_res <- tryCatch(chisq.test(tab, correct = FALSE),
                               error = function(e) list(p.value = NA))
        }
        p <- test_res$p.value
        # 存结果
        gene_res <- rbind(gene_res, data.frame(
          gene       = gene,
          group1     = g1,
          group2     = g2,
          OR         = OR_adj,
          log2_OR    = log2OR_adj,
          p_value    = p,
          fdr        = NA_real_,
          pct1_nonNA = s1$pct_nonNA,
          pct2_nonNA = s2$pct_nonNA,
          n_group1   = s1$n,
          n_group2   = s2$n,
          n_nonNA_group1 = s1$n_nonNA,
          n_nonNA_group2 = s2$n_nonNA,
          n_NA_group1    = s1$n_NA,
          n_NA_group2    = s2$n_NA,
          stringsAsFactors = FALSE
        ))
        gene_p <- c(gene_p, p)
      }
    }
    # 对该基因的所有两两比较统一做 FDR 校正
    if (nrow(gene_res) > 0 && length(gene_p) > 0) {
      gene_res$fdr <- p.adjust(gene_p, method = "fdr")
    }
    result_list[[gene]] <- gene_res
  }
  # 合并所有基因的结果
  do.call(rbind, result_list)
}
# 提取基因列的列名
gene_cols <- setdiff(names(noise), c("dataset","group","category","Disease_group","sex","corrected_noise"))

# 1.卡方检验-group
unique(noise$group)
result_1 <- pairwise_nonNA_safe(noise, "group", gene_cols) 
write.csv(result_1, file="/data/wuwq/noise/DISEASE_ALL/noise_chisq/otherT_group_chisq.csv", row.names = TRUE)

# 2.卡方检验-category
result_2 <- pairwise_nonNA_safe(noise, "category", gene_cols) 
write.csv(result_2, file="/data/wuwq/noise/DISEASE_ALL/noise_chisq/otherT_category_chisq.csv", row.names = TRUE)

# 3.卡方检验-Disease_group
result_3 <- pairwise_nonNA_safe(noise, "Disease_group", gene_cols) 
write.csv(result_3, file="/data/wuwq/noise/DISEASE_ALL/noise_chisq/otherT_Disease_group_chisq.csv", row.names = TRUE)





##################################################################
### ---------- (二)基因水平noise 秩和检验 --------------- ###     ✅️
##################################################################
library(dplyr)
library(knitr)
library(tidyr)
library(stats)
library(tibble)
unique(noise$group)

rm(list=ls())
# 读取merge好的全数据集  health + disease
noise <- read.csv("/data/wuwq/noise/DISEASE_ALL/merge_disease_health_B.csv", header = FALSE, row.names = 1)
colnames(noise) <- as.character(noise[1, ])
noise <- noise[-1, , drop = FALSE]
# 从第6列到最后一列变为数值型
noise[, 6:ncol(noise)] <- lapply(noise[, 6:ncol(noise)], function(x) {
  x <- as.numeric(x)
  return(x)
})
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

# 使用自定义的函数进行两两比较
# 2.秩和检验-group
#result_group <- pairwise_wilcox_gene(noise, "group", gene_start_col = 7)  
#write.csv(result_group, file="/data/wuwq/noise/DISEASE_ALL/noise_wilcox/otherT_group_wilcox.csv", row.names = TRUE)
# 2.秩和检验-category
#result_category <- pairwise_wilcox_gene(noise, "category", gene_start_col = 7)  
#write.csv(result_category, file="/data/wuwq/noise/DISEASE_ALL/noise_wilcox/otherT_category_wilcox.csv", row.names = TRUE)
# 3.秩和检验-Disease_group
result_disease <- pairwise_wilcox_gene(noise, "Disease_group", gene_start_col = 7)  
write.csv(result_disease, file="/data/wuwq/noise/DISEASE_ALL/noise_wilcox/B_Disease_group_wilcox.csv", row.names = TRUE)






##################################################################
### -------(三)读取和处理基因noise结果表-高变举例 ------------ ### ❌️
##################################################################
### CD4+ T cell
### CD8+ T cell
### Monocyte
### NK cell
### Dendritic cell
### B cell
### Other T cell

rm(list=ls())
##### 读取融合过后的symbol数据
noise <- read.csv("/data/wuwq/noise/DISEASE_ALL/merge_disease_health_CD4.csv", header = FALSE, row.names = 1)
colnames(noise) <- as.character(noise[1, ])
noise <- noise[-1, , drop = FALSE]
# 从第6列到最后一列变为数值型
noise[, 6:ncol(noise)] <- lapply(noise[, 6:ncol(noise)], function(x) {
  x <- as.numeric(x)
  return(x)
})

# 读取管家基因数据
hk_genes <- read.csv("/data/wuwq/noise/Housekeeping_GenesHuman.csv", header = TRUE)$Gene.name  # Housekeeping Transcript Atlas 网站
hk_genes <- unique(hk_genes[!is.na(hk_genes) & hk_genes != ""])

### 1.所有数据集aidas, eqtl, health, autoimmu, infect, tumor 的 gene noise均值计算 及 高noise基因鉴定
unique(noise$dataset)
filtered_health <- subset(noise, dataset == "aidas") # aidas, eqtl
#unique(noise$category)
#filtered_health <- subset(noise, category == "tumor") # health, autoimmu, infect, tumor

# 只取基因表达列（从第6列开始到最后） # mean_noise降序
gene_cols <- names(filtered_health)[6:ncol(filtered_health)]
result <- filtered_health %>%
  select(all_of(gene_cols)) %>%
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE))) %>%
  pivot_longer(everything(),
               names_to = "gene",
               values_to = "mean_noise") %>%
  mutate(
    rank_asc  = rank(mean_noise, ties.method = "min"),     # 1 = 最低表达
    rank_desc = rank(-mean_noise, ties.method = "min")     # 1 = 最高表达
  ) %>%
  arrange(desc(mean_noise))   # 主排序保持降序
# mean_noise升序
result %>% arrange(mean_noise)
# mean_noise降序
result %>% arrange(desc(mean_noise))

# result 中增加是否管家基因列（升序看最稳定表达的）
result <- result %>%
  mutate(is_hk = gene %in% hk_genes) %>%          # TRUE = 是管家基因
  arrange(mean_noise)                             # 升序：最稳定（noise 最小）排最前面

write.csv(result, file="/data/wuwq/noise/HEALTH_REAL/out/CD4_tumor_gene_noise.csv")




##################################################################
######## -------------(四)读取卡方检验的结果------------- ######## ✅️ noi占比改变的基因
##################################################################

rm(list=ls())
##### 
chisq_gene <- read.csv("/data/wuwq/noise/DISEASE_ALL/noise_chisq/otherT_Disease_group_chisq.csv", header = TRUE, row.names = 1)
unique(chisq_gene$group1)
unique(chisq_gene$group2)

# 筛选满足条件的数据
# 换疾病 "health" "RA" "SLE" "COVID-19" "post-COVID-19 disorder" "SEPSIS" "FLU" "NSCLC"
# 规则为：noi占比改变的定义：卡方检验fdr<0.05
# 处理过程中至少有一个组中该noi_gene鉴定占比>=5%，即分析基因需要满足在在至少一个组中是noise基因
filt_gene <- chisq_gene[ 
    # group1 和 group2 包含 疾病 & "health"
    (chisq_gene$group1 == "NSCLC" & chisq_gene$group2 == "health") & 
    # 至少有一个组中该noi_gene鉴定占比>=5%
    (chisq_gene$pct1_nonNA >= 0.05 | chisq_gene$pct2_nonNA >= 0.05) & 
    # 满足 log2fc 的范围，或 fdr 的条件
    (chisq_gene$fdr < 0.05) ,]

write.csv(filt_gene, file="/data/wuwq/noise/DISEASE_ALL/noise_chisq/change_percent_noi_gene/otherT_NSCLC.csv", row.names = FALSE)





##################################################################
######## -------------(五)读取noise秩和检验的结果------------- ######## ✅️ noi水平改变的基因 → 放松条件，先不关注基因占比
##################################################################

rm(list=ls())
##### 
wilcox_gene <- read.csv("/data/wuwq/noise/DISEASE_ALL/noise_wilcox/otherT_Disease_group_wilcox.csv", header = TRUE, row.names = 1)
unique(wilcox_gene$group1)
unique(wilcox_gene$group2)
### 增加pct1_nonNA列
#summary(wilcox_gene) #n1为有效样本数
wilcox_gene <- wilcox_gene %>%
  mutate(
    pct1_nonNA = (n1) / (n1+na_count1),
    pct2_nonNA = (n2) / (n2+na_count2)
  )
# 筛选满足条件的数据
# 换疾病 "health" "RA" "SLE" "COVID-19" "post-COVID-19 disorder" "SEPSIS" "FLU" "NSCLC"
# 规则为：noi水平改变的定义：fdr<0.05 且 (log2FC > 1 或 log2FC < -1)
filt_gene <- wilcox_gene[ 
  # group1 和 group2 包含 疾病 & "health"
  (wilcox_gene$group1 == "RA" & wilcox_gene$group2 == "health") & 
  #  fdr 的条件
  (wilcox_gene$fdr < 0.05) & 
  # 满足 log2fc 的范围
  (wilcox_gene$log2FC < -1 | wilcox_gene$log2FC > 1) & !is.na(wilcox_gene$log2FC) &
  # 2组中noi_gene 其一 均满足 占比>=5%
  (wilcox_gene$pct1_nonNA >= 0.05 | wilcox_gene$pct2_nonNA >= 0.05)  ,]
write.csv(filt_gene, file="/data/wuwq/noise/DISEASE_ALL/noise_wilcox/change_level_noi_gene/otherT_NSCLC.csv", row.names = FALSE)







##################################################################
######## -----(六)鉴定表达不变 & noise水平改变的基因---#- ######## ✅️ exp水平不变的基因 + noi水平改变的基因 (+ noi占比改变情况 
##################################################################

rm(list=ls())
# "health" "RA" "SLE" "COVID-19" "post-COVID-19 disorder" "SEPSIS" "FLU" "NSCLC"
exp_nochange_level <- read.csv("/data/wuwq/noise/DISEASE_ALL/exp_u_wilcox/nochange_exp_common_gene/CD4_NSCLC.csv", header = TRUE) 
noi_change_level <- read.csv("/data/wuwq/noise/DISEASE_ALL/noise_wilcox/change_level_noi_gene/CD4_NSCLC.csv", header = TRUE)
noi_change_perc <- read.csv("/data/wuwq/noise/DISEASE_ALL/noise_chisq/change_percent_noi_gene/CD4_NSCLC.csv", header = TRUE)
exp_nochange_level_gene <- exp_nochange_level$x   # exp不变的基因
noi_change_level_gene <- noi_change_level$gene    # noi水平改变的基因
noi_change_perc_gene <- noi_change_perc$gene      # noi占比改变的基因
# exp不变，noi水平改变
exp_nochange_noi_level_change <- intersect(exp_nochange_level_gene, noi_change_level_gene)
# exp不变，noi占比改变
exp_nochange_noi_perc_change <- intersect(exp_nochange_level_gene, noi_change_perc_gene)
# exp不变，noi水平和占比同时改变
exp_nochange_noi_level_perc_change <- Reduce(intersect, list(exp_nochange_level_gene, noi_change_level_gene, noi_change_perc_gene))
# exp不变，noi水平改变，noi占比不变
exp_nochange_noi_level_change_only <- setdiff(exp_nochange_noi_level_change, exp_nochange_noi_level_perc_change)
# exp不变，noi水平不变，noi占比改变
exp_nochange_noi_perc_change_only <- setdiff(exp_nochange_noi_perc_change, exp_nochange_noi_level_perc_change)

### 再对应回原表
# exp 不变，noi 水平改变的基因 #秩和检验
wilcox_noi_level_change <- noi_change_level[noi_change_level$gene %in% exp_nochange_noi_level_change, ]
# exp 不变，noi 水平和占比同时改变的基因 #秩和检验
wilcox_noi_level_perc_change <- noi_change_level[noi_change_level$gene %in% exp_nochange_noi_level_perc_change, ]
# exp 不变，noi 水平改变，noi 占比不变的基因 #秩和检验
wilcox_noi_level_change_only <- noi_change_level[noi_change_level$gene %in% exp_nochange_noi_level_change_only, ]
# 查看 exp 不变，noi 占比改变的基因 #卡方检验
chisq_noi_perc_change <- noi_change_perc[noi_change_perc$gene %in% exp_nochange_noi_perc_change, ]
# 查看 exp 不变，noi 水平和占比同时改变的基因 #卡方检验
chisq_noi_level_perc_change <- noi_change_perc[noi_change_perc$gene %in% exp_nochange_noi_level_perc_change, ]
# 查看 exp 不变，noi 水平不变，noi 占比改变的基因 #卡方检验
chisq_noi_perc_change_only <- noi_change_perc[noi_change_perc$gene %in% exp_nochange_noi_perc_change_only, ]

# 保存
new_folder_path <- "/data/wuwq/noise/DISEASE_ALL/analyse_exp_nochange_noi_change/otherT_NSCLC"  
dir.create(new_folder_path)
setwd(new_folder_path)
write.csv(wilcox_noi_level_change, file="wilcox_noi_level_change.csv", row.names = FALSE)
write.csv(wilcox_noi_level_perc_change, file="wilcox_noi_level_perc_change.csv", row.names = FALSE)
write.csv(wilcox_noi_level_change_only, file="wilcox_noi_level_change_only.csv", row.names = FALSE)
write.csv(chisq_noi_perc_change, file="chisq_noi_perc_change.csv", row.names = FALSE)
write.csv(chisq_noi_level_perc_change, file="chisq_noi_level_perc_change.csv", row.names = FALSE)
write.csv(chisq_noi_perc_change_only, file="chisq_noi_perc_change_only.csv", row.names = FALSE)





































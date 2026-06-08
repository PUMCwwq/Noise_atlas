########################################################### 5_2:  个体和基因层面疾病合并
### 基于5_1的结果

### CD4+ T cell
### CD8+ T cell
### Monocyte
### NK cell
### Dendritic cell
### B cell
### Other T cell
library(dplyr)
library(knitr)
library(gt)

##################################################################
### ------- 个体水平noise波动 & 秩和检验 & PCA-------------- ###
##################################################################
### CD4+ T cell
### CD8+ T cell
### Monocyte
### NK cell
### Dendritic cell
### B cell
### Other T cell

rm(list=ls())
##### (一)、总体健康indi_noise波动    -- 特殊保存
noise_health <- read.csv("/data/wuwq/noise/HEALTH_REAL/sym_B_health_noise.csv", header = TRUE,row.names = 1)
colnames(noise_health) <- as.character(noise_health[1, ])
noise_health <- noise_health[-1, , drop = FALSE]
# 从第3列到最后一列变为数值型
noise_health[, 3:ncol(noise_health)] <- lapply(noise_health[, 3:ncol(noise_health)], function(x) {
  x <- as.numeric(x)
  return(x)
})
# 用规范后的数据覆盖原有数据--特殊保存
write.table(noise_health, file = "/data/wuwq/noise/HEALTH_REAL/sym_B_health_noise.csv", sep = ",", row.names = TRUE, col.names = NA, quote = TRUE)



rm(list=ls())
### 读取规范后的健康人symbol数据
noise_health <- read.csv("/data/wuwq/noise/HEALTH_REAL/sym_CD4_health_noise.csv", header = TRUE,row.names = 1)
# 1.全部健康个体的波动 ✅️
result <- noise_health %>%
  summarise(
    n_samples = n(),
    mean_noise = mean(corrected_noise, na.rm = TRUE),
    sd_noise = sd(corrected_noise, na.rm = TRUE),
    se_noise = sd_noise / sqrt(n_samples),  # 新增的SE计算
    var_noise = var(corrected_noise, na.rm = TRUE),
    # 分位数
    min_noise = min(corrected_noise, na.rm = TRUE),
    p20 = quantile(corrected_noise, 0.20, na.rm = TRUE),
    p40 = quantile(corrected_noise, 0.40, na.rm = TRUE),
    p60 = quantile(corrected_noise, 0.60, na.rm = TRUE),
    p80 = quantile(corrected_noise, 0.80, na.rm = TRUE),
    max_noise = max(corrected_noise, na.rm = TRUE)
    # ±1sd, ±2sd, ±3sd（常用于异常值判断或范围参考）
    #mean_minus_3sd = mean_noise - 3 * sd_noise,
    #mean_minus_2sd = mean_noise - 2 * sd_noise,
    #mean_minus_1sd = mean_noise - 1 * sd_noise,
    #mean_plus_1sd  = mean_noise + 1 * sd_noise,
    #mean_plus_2sd  = mean_noise + 2 * sd_noise,
    #mean_plus_3sd  = mean_noise + 3 * sd_noise,
  )
kable(result, format = "markdown")




rm(list=ls())
##### (二)、总体全疾病indi_noise波动    
noise_disease <- read.csv("/data/wuwq/noise/DISEASE_ALL/union_CD4.csv", header = TRUE,row.names = 1)
# 1.全疾病个体的波动 ✅️
result <- noise_disease %>%
  summarise(
    n_samples = n(),
    mean_noise = mean(corrected_noise, na.rm = TRUE),
    sd_noise = sd(corrected_noise, na.rm = TRUE),
    se_noise = sd_noise / sqrt(n_samples),  # 新增的SE计算
    var_noise = var(corrected_noise, na.rm = TRUE),
    # 分位数
    min_noise = min(corrected_noise, na.rm = TRUE),
    p20 = quantile(corrected_noise, 0.20, na.rm = TRUE),
    p40 = quantile(corrected_noise, 0.40, na.rm = TRUE),
    p60 = quantile(corrected_noise, 0.60, na.rm = TRUE),
    p80 = quantile(corrected_noise, 0.80, na.rm = TRUE),
    max_noise = max(corrected_noise, na.rm = TRUE)
  )
kable(result, format = "markdown")



rm(list=ls())
# 2.疾病数据集 波动 ❌️
noise_disease <- read.csv("/data/wuwq/noise/DISEASE_ALL/union_DC.csv", header = TRUE) %>%
  filter(dataset != "dataset", !is.na(dataset), dataset != "")

result <- noise_disease %>%
  group_by(dataset) %>%
  summarise(
    n_samples = n(),
    mean_noise = mean(corrected_noise, na.rm = TRUE),
    sd_noise = sd(corrected_noise, na.rm = TRUE),
    var_noise = var(corrected_noise, na.rm = TRUE),
    # 分位数
    min_noise = min(corrected_noise, na.rm = TRUE),
    p20 = quantile(corrected_noise, 0.20, na.rm = TRUE),
    p40 = quantile(corrected_noise, 0.40, na.rm = TRUE),
    p60 = quantile(corrected_noise, 0.60, na.rm = TRUE),
    p80 = quantile(corrected_noise, 0.80, na.rm = TRUE),
    max_noise = max(corrected_noise, na.rm = TRUE)
  )
kable(result, format = "markdown")




rm(list=ls())
# 3.疾病大类 波动  ✅️
noise_disease <- read.csv("/data/wuwq/noise/DISEASE_ALL/union_B.csv", header = TRUE,row.names = 1) 
unique(noise_disease$Disease_group)
noise_disease$category <- noise_disease$Disease_group %>%
  recode(
    "COVID-19" = "infect",
    "post-COVID-19 disorder" = "infect",
    "SEPSIS" = "infect",
    "FLU" = "infect",
    "NSCLC" = "tumor",
    "RA" = "autoimmu",
    "SLE" = "autoimmu"
  )
noise_disease <- noise_disease %>% relocate(category, .before = 1)
result <- noise_disease %>%
  group_by(category) %>%
  summarise(
    n_samples = n(),
    mean_noise = mean(corrected_noise, na.rm = TRUE),
    sd_noise = sd(corrected_noise, na.rm = TRUE),
    se_noise = sd_noise / sqrt(n_samples),  # 新增的SE计算
    var_noise = var(corrected_noise, na.rm = TRUE),
    # 分位数
    min_noise = min(corrected_noise, na.rm = TRUE),
    p20 = quantile(corrected_noise, 0.20, na.rm = TRUE),
    p40 = quantile(corrected_noise, 0.40, na.rm = TRUE),
    p60 = quantile(corrected_noise, 0.60, na.rm = TRUE),
    p80 = quantile(corrected_noise, 0.80, na.rm = TRUE),
    max_noise = max(corrected_noise, na.rm = TRUE)
  )
kable(result, format = "markdown")



rm(list=ls())
# 4.疾病小类 波动 ✅️
noise_disease <- read.csv("/data/wuwq/noise/DISEASE_ALL/union_B.csv", header = TRUE,row.names = 1)
unique(noise_disease$Disease_group)
result <- noise_disease %>%
  group_by(Disease_group) %>%
  summarise(
    n_samples = n(),
    mean_noise = mean(corrected_noise, na.rm = TRUE),
    sd_noise = sd(corrected_noise, na.rm = TRUE),
    se_noise = sd_noise / sqrt(n_samples),  # 新增的SE计算
    var_noise = var(corrected_noise, na.rm = TRUE),
    # 分位数
    min_noise = min(corrected_noise, na.rm = TRUE),
    p20 = quantile(corrected_noise, 0.20, na.rm = TRUE),
    p40 = quantile(corrected_noise, 0.40, na.rm = TRUE),
    p60 = quantile(corrected_noise, 0.60, na.rm = TRUE),
    p80 = quantile(corrected_noise, 0.80, na.rm = TRUE),
    max_noise = max(corrected_noise, na.rm = TRUE)
  )
kable(result, format = "markdown")










##################################################################
### ------------------- 个体水平noise 秩和检验 --------------- ###
##################################################################
library(dplyr)
library(tidyr)
library(knitr)
### CD4+ T cell
### CD8+ T cell
### Monocyte
### NK cell
### Dendritic cell
### B cell
### Other T cell

rm(list=ls())
##### 疾病和健康数据合并，秩和检验准备---特殊保存  ✅️
### 疾病规整
noise_disease <- read.csv("/data/wuwq/noise/DISEASE_ALL/union_otherT.csv", header = FALSE, row.names = 1)  # header = FALSE
colnames(noise_disease) <- as.character(noise_disease[1, ])
noise_disease <- noise_disease[-1, , drop = FALSE]
noise_disease <- noise_disease[,-3]
col_order <- c("dataset", "Disease_group", "sex", "corrected_noise", setdiff(names(noise_disease), c("dataset", "Disease_group", "sex", "corrected_noise")))
noise_disease <- noise_disease[ , col_order]
### 健康规整
noise_health <- read.csv("/data/wuwq/noise/HEALTH_REAL/sym_otherT_health_noise.csv", header = FALSE, row.names = 1) # header = FALSE
colnames(noise_health) <- as.character(noise_health[1, ])
noise_health <- noise_health[-1, , drop = FALSE]
noise_health <- noise_health[ , -3, drop = FALSE]
unique(noise_health$group)
colnames(noise_health)[colnames(noise_health) == "group"] <- "dataset"
noise_health <- noise_health %>% mutate(Disease_group = "health")
unique(noise_health$Disease_group)
noise_health <- noise_health %>% select(dataset, Disease_group, sex, corrected_noise,  everything())
# 合并全疾病和健康数据 # 把所有 "NA." 开头的列删除
noise <- bind_rows(
  noise_disease %>% select(-starts_with("NA.")),
  noise_health  %>% select(-starts_with("NA."))
) %>%
  mutate(
    across(
      -c(dataset, Disease_group, sex),
      ~ as.numeric(as.character(.x))
    )
  )
### 大类定义
unique(noise$Disease_group)
noise$category <- noise$Disease_group %>%
  recode(
    "COVID-19" = "infect",
    "post-COVID-19 disorder" = "infect",
    "SEPSIS" = "infect",
    "FLU" = "infect",
    "NSCLC" = "tumor",
    "RA" = "autoimmu",
    "SLE" = "autoimmu",
    "health" = "health"
  )
noise <- noise %>% select(dataset, category, Disease_group, sex, corrected_noise,  everything())
### 总疾病定义
unique(noise$category)
noise$group <- noise$category %>%
  recode(
    "infect" = "disease",
    "autoimmu" = "disease",
    "tumor" = "disease",
    "health" = "health"
  )
noise <- noise %>% select(dataset, group, category, Disease_group, sex, corrected_noise,  everything())
category_counts <- noise %>%
  count(category, name = "n_samples", sort = TRUE)
print(category_counts)
write.table(noise,file = "/data/wuwq/noise/DISEASE_ALL/merge_disease_health_otherT.csv", sep = ",", row.names = TRUE, col.names = NA, quote = TRUE)






rm(list=ls())
# 读取merge好的全数据集 ✅️
noise <- read.csv("/data/wuwq/noise/DISEASE_ALL/merge_disease_health_CD4.csv", header = FALSE, row.names = 1)
colnames(noise) <- as.character(noise[1, ])
noise <- noise[-1, , drop = FALSE]

# 自定义两两比较的函数
pairwise_wilcox <- function(data, group_col, value_col) {
  # ─────────────── 改动 1 ───────────────
  # 不再提前 filter(!is.na(group_col))，保留所有行
  # groups 改用 na.omit 后再取唯一值（避免 NA 作为一个“组”参与比较）
  groups <- data[[group_col]] %>% na.omit() %>% unique() %>% sort()
  comparisons <- combn(groups, 2, simplify = FALSE)
  results <- lapply(comparisons, function(pair) {
    g1 <- pair[1]
    g2 <- pair[2]
    # ─────────────── 改动 2 ───────────────
    # 不再用 na.omit() 提前丢掉 value 的 NA
    # 直接取所有对应行的值（包含 NA）
    v1_all <- data %>% filter(.data[[group_col]] == g1) %>% pull(.data[[value_col]])
    v2_all <- data %>% filter(.data[[group_col]] == g2) %>% pull(.data[[value_col]])
    # 只在计算统计量时自然处理 NA
    v1_clean <- v1_all[!is.na(v1_all)]
    v2_clean <- v2_all[!is.na(v2_all)]
    # 如果有效样本量太少（<2），返回 NA 结果
    if (length(v1_clean) < 2 || length(v2_clean) < 2) {
      return(tibble(
        group1 = g1, group2 = g2,
        mean_group1 = NA_real_, mean_group2 = NA_real_,
        mean_diff = NA_real_, fold_change = NA_real_, log2FC = NA_real_,
        p_value = NA_real_
      ))
    }
    m1 <- mean(v1_clean)
    m2 <- mean(v2_clean)
    tibble(
      group1      = g1,
      group2      = g2,
      mean_group1 = m1,
      mean_group2 = m2,
      mean_diff   = m1 - m2,
      fold_change = if (m2 == 0) NA_real_ else m1 / m2,
      log2FC      = if (m2 == 0 || m1/m2 <= 0) NA_real_ else log2(m1/m2),
      p_value     = wilcox.test(v1_clean, v2_clean, exact = FALSE)$p.value
    )
  }) %>% bind_rows()
  results %>%
    mutate(fdr = if_else(!is.na(p_value),
                         p.adjust(p_value, method = "fdr"),
                         NA_real_)) %>%
    relocate(group1, group2, mean_group1, mean_group2, log2FC, .before = mean_diff)
}
# 使用自定义的函数进行两两比较
result_group <- pairwise_wilcox(noise, "group", "corrected_noise")
#result_category <- pairwise_wilcox(noise, "category", "corrected_noise")
#result_disease <- pairwise_wilcox(noise, "Disease_group", "corrected_noise")




rm(list=ls())
# 读取merge好的全数据集后，计算power_noise，加算power_noise的组间秩和检验 ❗️ 
noise <- read.csv("/data/wuwq/noise/DISEASE_ALL/merge_disease_health_otherT.csv", header = FALSE, row.names = 1)
colnames(noise) <- as.character(noise[1, ])
noise <- noise[-1, , drop = FALSE]
# 确保第6列 corrected_noise 是 numeric
noise$corrected_noise <- as.numeric(noise$corrected_noise)
# 将第7列到最后一列转换为 numeric（防止字符类型）
noise[, 7:ncol(noise)] <- lapply(noise[, 7:ncol(noise)], as.numeric)
# 计算每行非NA值的个数（第7列到最后一列）
non_na_counts <- apply(noise[, 7:ncol(noise)], 1, function(x) sum(!is.na(x)))
# 新增 power_noise 列
noise$power_noise <- noise$corrected_noise / non_na_counts
noise <- noise %>%
  relocate(power_noise, .after = corrected_noise)


# 自定义两两比较的函数
pairwise_wilcox <- function(data, group_col, value_col) {
  # ─────────────── 改动 1 ───────────────
  # 不再提前 filter(!is.na(group_col))，保留所有行
  # groups 改用 na.omit 后再取唯一值（避免 NA 作为一个“组”参与比较）
  groups <- data[[group_col]] %>% na.omit() %>% unique() %>% sort()
  comparisons <- combn(groups, 2, simplify = FALSE)
  results <- lapply(comparisons, function(pair) {
    g1 <- pair[1]
    g2 <- pair[2]
    # ─────────────── 改动 2 ───────────────
    # 不再用 na.omit() 提前丢掉 value 的 NA
    # 直接取所有对应行的值（包含 NA）
    v1_all <- data %>% filter(.data[[group_col]] == g1) %>% pull(.data[[value_col]])
    v2_all <- data %>% filter(.data[[group_col]] == g2) %>% pull(.data[[value_col]])
    # 只在计算统计量时自然处理 NA
    v1_clean <- v1_all[!is.na(v1_all)]
    v2_clean <- v2_all[!is.na(v2_all)]
    # 如果有效样本量太少（<2），返回 NA 结果
    if (length(v1_clean) < 2 || length(v2_clean) < 2) {
      return(tibble(
        group1 = g1, group2 = g2,
        mean_group1 = NA_real_, mean_group2 = NA_real_,
        mean_diff = NA_real_, fold_change = NA_real_, log2FC = NA_real_,
        p_value = NA_real_
      ))
    }
    m1 <- mean(v1_clean)
    m2 <- mean(v2_clean)
    tibble(
      group1      = g1,
      group2      = g2,
      mean_group1 = m1,
      mean_group2 = m2,
      mean_diff   = m1 - m2,
      fold_change = if (m2 == 0) NA_real_ else m1 / m2,
      log2FC      = if (m2 == 0 || m1/m2 <= 0) NA_real_ else log2(m1/m2),
      p_value     = wilcox.test(v1_clean, v2_clean, exact = FALSE)$p.value
    )
  }) %>% bind_rows()
  results %>%
    mutate(fdr = if_else(!is.na(p_value),
                         p.adjust(p_value, method = "fdr"),
                         NA_real_)) %>%
    relocate(group1, group2, mean_group1, mean_group2, log2FC, .before = mean_diff)
}
# 使用自定义的函数进行两两比较
result_disease <- pairwise_wilcox(noise, "Disease_group", "power_noise")



rm(list=ls())
# 加算疾病小类 power_noise波动 ❗️
noise <- read.csv("/data/wuwq/noise/DISEASE_ALL/merge_disease_health_otherT.csv", header = FALSE, row.names = 1)
colnames(noise) <- as.character(noise[1, ])
noise <- noise[-1, , drop = FALSE]
# 确保第6列 corrected_noise 是 numeric
noise$corrected_noise <- as.numeric(noise$corrected_noise)
# 将第7列到最后一列转换为 numeric（防止字符类型）
noise[, 7:ncol(noise)] <- lapply(noise[, 7:ncol(noise)], as.numeric)
# 计算每行非NA值的个数（第7列到最后一列）
non_na_counts <- apply(noise[, 7:ncol(noise)], 1, function(x) sum(!is.na(x)))
# 新增 power_noise 列
noise$power_noise <- noise$corrected_noise / non_na_counts
noise <- noise %>%
  relocate(power_noise, .after = corrected_noise)

unique(noise$Disease_group)
result <- noise %>%
  group_by(Disease_group) %>%
  summarise(
    n_samples = n(),
    mean_noise = mean(power_noise, na.rm = TRUE),
    sd_noise = sd(power_noise, na.rm = TRUE),
    se_noise = sd_noise / sqrt(n_samples),  # 新增的SE计算
    var_noise = var(power_noise, na.rm = TRUE),
    # 分位数
    min_noise = min(power_noise, na.rm = TRUE),
    p20 = quantile(power_noise, 0.20, na.rm = TRUE),
    p40 = quantile(power_noise, 0.40, na.rm = TRUE),
    p60 = quantile(power_noise, 0.60, na.rm = TRUE),
    p80 = quantile(power_noise, 0.80, na.rm = TRUE),
    max_noise = max(power_noise, na.rm = TRUE)
  )
kable(result, format = "markdown")


##################################################################
############ -------- result_1_var_noise_gene -------#############
############ ------      disease vs health       -------##########
############ ------         PubMed辅证           -------########## 
##################################################################
library(dplyr)
library(knitr)
library(tidyr)
library(stats)
library(tibble)
library(data.table)
library(ggplot2)
library(ggdist)     # 画half violin
library(ggbeeswarm) # 抖动散点更好看
install.packages(c("ggdist","ggbeeswarm"))
library(vipor)
install.packages("wordcloud2")
library(wordcloud2)
if (!requireNamespace("svglite", quietly = TRUE)) {
  install.packages("svglite")
}
install.packages("Cairo")
library(Cairo)



rm(list=ls())
##### （一） 整体情况总览：各组 top 高变noise基因 和 百分区间分箱  Disease_group + health  ✅️ 普通保存，否则rank列错误
### 1.读取noise表
noise <- read.csv("/data/wuwq/noise/DISEASE_ALL/merge_disease_health_otherT.csv", header = FALSE, row.names = 1)
colnames(noise) <- as.character(noise[1, ])
noise <- noise[-1, , drop = FALSE]
# 从第7列到最后一列变为数值型
noise[, 7:ncol(noise)] <- lapply(noise[, 7:ncol(noise)], function(x) {
  x <- as.numeric(x)
  return(x)
})
# 读取管家基因数据
hk_genes <- read.csv("/data/wuwq/noise/Housekeeping_GenesHuman.csv", header = TRUE)$Gene.name  # Housekeeping Transcript Atlas 网站
hk_genes <- unique(hk_genes[!is.na(hk_genes) & hk_genes != ""])
### 选不同分组 
unique(noise$Disease_group) #"health" "RA" "SLE" "COVID-19" "post-COVID-19 disorder" "SEPSIS" "FLU" "NSCLC"
filtered_noise <- subset(noise, Disease_group == "NSCLC") 
# 只取基因表达列（从第7列开始到最后） # mean_noise降序
gene_cols <- names(filtered_noise)[7:ncol(filtered_noise)]
# 计算总样本数（行数）
total_samples <- nrow(filtered_noise)
# 第一部分：创建基因级别的详细数据
result_detail <- filtered_noise %>%
  select(all_of(gene_cols)) %>%
  summarise(across(everything(), 
                   list(mean = ~mean(.x, na.rm = TRUE),
                        noise_percentage = ~sum(!is.na(.x))/total_samples),
                   .names = "{.col}@{.fn}")) %>%
  pivot_longer(
    everything(),
    names_to = c("gene", ".value"),
    names_sep = "@"
  ) %>%
  rename(mean_noise = mean) %>%
  mutate(
    noise_rank = rank(-mean_noise, ties.method = "min")
  ) %>%
  mutate(is_hk = gene %in% hk_genes) %>%
  arrange(desc(mean_noise))
# 第二部分：创建分箱统计并合并到详细数据
# 定义自定义分箱区间：0单独一箱，然后是(0,0.1), [0.1,0.2), [0.2,0.4), [0.4,0.6), [0.6,0.8), [0.8,1]
result_combined <- result_detail %>%
  mutate(
    noise_bin = case_when(
      noise_percentage == 0 ~ "0",
      noise_percentage > 0 & noise_percentage < 0.05 ~ "(0,0.05)",
      noise_percentage >= 0.05 & noise_percentage < 0.1 ~ "[0.05,0.1)",
      noise_percentage >= 0.1 & noise_percentage < 0.2 ~ "[0.1,0.2)",
      noise_percentage >= 0.2 & noise_percentage < 0.4 ~ "[0.2,0.4)",
      noise_percentage >= 0.4 & noise_percentage < 0.6 ~ "[0.4,0.6)",
      noise_percentage >= 0.6 & noise_percentage < 0.8 ~ "[0.6,0.8)",
      noise_percentage >= 0.8 & noise_percentage < 1 ~ "[0.8,1)",
      noise_percentage == 1 ~ "1"  # 修正：使用 ==
    ),
    noise_bin = factor(
      noise_bin,
      levels = c("0", "(0,0.05)", "[0.05,0.1)", "[0.1,0.2)", "[0.2,0.4)", 
                 "[0.4,0.6)", "[0.6,0.8)", "[0.8,1)", "1")
    )
  ) %>%
  group_by(noise_bin) %>%
  mutate(
    bin_count = n(),
    bin_rank = rank(-mean_noise, ties.method = "min")
  ) %>%
  ungroup()
#### 读取分箱统计
bin_counts <- result_combined %>%
  group_by(noise_bin) %>%
  summarise(gene_count = n()) %>%
  ungroup()
# 定义所需的顺序
desired_order <- c("0", "(0,0.05)", "[0.05,0.1)", "[0.1,0.2)", 
                   "[0.2,0.4)", "[0.4,0.6)", "[0.6,0.8)", "[0.8,1)", "1")
# 转换并排序
bin_counts_ordered <- bin_counts %>%
  mutate(noise_bin = factor(noise_bin, levels = desired_order)) %>%
  arrange(noise_bin)
# 输出标准表格
knitr::kable(bin_counts_ordered,
             format = "simple",
             col.names = c("Noise_Bin", "Gene_Count"),
             caption = "Gene Count by Noise Percentage Bin (Ordered)",
             align = "c")
write.csv(result_combined, file="/data/wuwq/noise/DISEASE_ALL/noise_var/NSCLC_otherT.csv", row.names = TRUE)

### 查看每个组的样本数
group_counts <- noise %>%
  count(Disease_group, name = "Count")
knitr::kable(group_counts, 
             format = "simple",
             col.names = c("Disease_Group", "Count"),
             caption = "Sample Count by Disease Group")





##### (二) 批量查询 PubMed 的解决方案，使用 rentrez 包进行 API 调用
# 安装必要包（如果尚未安装）
if (!require("rentrez")) install.packages("rentrez")
library(rentrez)


rm(list=ls())
# 获取所有细胞类型的gene名并集
set_entrez_key("0cdfdf17cc0aad9bd72ce3a44bfdd347d60a") # NCBI的API号
data1 <- read.csv("/data/wuwq/noise/DISEASE_ALL/exp_u/union_healthdisease_CD4.csv", header = FALSE, row.names = 1)
data2 <- read.csv("/data/wuwq/noise/DISEASE_ALL/exp_u/union_healthdisease_CD8.csv", header = FALSE, row.names = 1)
data3 <- read.csv("/data/wuwq/noise/DISEASE_ALL/exp_u/union_healthdisease_MNP.csv", header = FALSE, row.names = 1)
data4 <- read.csv("/data/wuwq/noise/DISEASE_ALL/exp_u/union_healthdisease_NK.csv", header = FALSE, row.names = 1)
data5 <- read.csv("/data/wuwq/noise/DISEASE_ALL/exp_u/union_healthdisease_DC.csv", header = FALSE, row.names = 1)
data6 <- read.csv("/data/wuwq/noise/DISEASE_ALL/exp_u/union_healthdisease_B.csv", header = FALSE, row.names = 1)
data7 <- read.csv("/data/wuwq/noise/DISEASE_ALL/exp_u/union_healthdisease_otherT.csv", header = FALSE, row.names = 1)
gene1 <- as.character(data1[1, c(5:ncol(data1))])
gene2 <- as.character(data2[1, c(5:ncol(data2))])
gene3 <- as.character(data3[1, c(5:ncol(data3))])
gene4 <- as.character(data4[1, c(5:ncol(data4))])
gene5 <- as.character(data5[1, c(5:ncol(data5))])
gene6 <- as.character(data6[1, c(5:ncol(data6))])
gene7 <- as.character(data7[1, c(5:ncol(data7))])
gene_list <- unique(c(gene1, gene2, gene3, gene4, gene5, gene6, gene7))

# 创建存储结果的空数据框
results <- data.frame(Gene = character(),
                      PubMed_Count = integer(),
                      stringsAsFactors = FALSE)
# 设置查询参数
query_delay <- 0.5  # 每次查询间隔时间（秒），避免请求过频
# 循环查询每个基因
for (gene in gene_list) {
  # 构建查询字符串
  query <- paste0("(", gene, "[Title/Abstract]) AND (\"cell cycle\"[Title/Abstract])")
  # 尝试查询
  tryCatch({
    # 执行PubMed查询
    search_result <- entrez_search(db = "pubmed", 
                                   term = query, 
                                   retmax = 0)  # retmax=0只返回计数
    # 获取文献数量
    count <- search_result$count
    # 添加到结果数据框
    results <- rbind(results, data.frame(Gene = gene, PubMed_Count = count))
    # 打印进度
    message(sprintf("Gene %s: %d publications found", gene, count))
    # 暂停以避免请求过频
    Sys.sleep(query_delay)
  }, error = function(e) {
    # 错误处理
    message(sprintf("Error querying gene %s: %s", gene, e$message))
    results <<- rbind(results, data.frame(Gene = gene, PubMed_Count = NA))
  })
}
# 保存结果
write.csv(results, file="/data/wuwq/noise/DISEASE_ALL/noise_var/pubmed_cellcycle.csv", row.names = TRUE)




##### (三) noi_gene 辅证---- 健康与对照
rm(list=ls())

### celltype & disease: "health" "RA" "SLE" "COVID-19" "post-COVID-19 disorder" "SEPSIS" "FLU" "NSCLC"
noi_gene <- read.csv("/data/wuwq/noise/DISEASE_ALL/noise_var/health_otherT.csv", header = TRUE, row.names = 1)
pubmed_query <- read.delim("/data/wuwq/noise/DISEASE_ALL/noise_var/pubmed_query.txt", header = TRUE)
cyclebase_ref <- read.delim("/data/wuwq/noise/DISEASE_ALL/noise_var/cyclebase_ref.txt", header = TRUE)
merged_data <- merge(noi_gene, pubmed_query, by = "gene")
# 新增一列 noi_gene

merged_data$noi_gene <- ifelse(
  merged_data$noise_percentage == 0, 
  "zero",  # noise_percentage为0时，设置为zero
  ifelse(
    merged_data$noise_percentage >= 0.05, 
    "yes",  # noise_percentage大于等于0.05
    "less0.05"  # noise_percentage在0到0.05之间
  )
)
unique(merged_data$noi_gene)

col_order <- c("gene","noi_gene","pubmed_apoptosis","pubmed_immune", "pubmed_aging","pubmed_sex","pubmed_cellcycle",  setdiff(names(merged_data), c("gene","noi_gene", "pubmed_apoptosis","pubmed_immune", "pubmed_aging","pubmed_sex","pubmed_cellcycle")))
merged_data <- merged_data[ , col_order]
### 选择检验的列
unique(merged_data$noi_gene)
merged_data <- merged_data[ ,c(1:7)]
### 用0替换掉NA
merged_data[is.na(merged_data)] <- 0

### 计算均值并检验
# 提取不同组合
yes_group <- merged_data %>% filter(noi_gene == "yes")
pubmed_query_group <- pubmed_query
pubmed_cellcycle_group <- cyclebase_ref
# 定义要比较的列
columns_to_compare <- c("pubmed_apoptosis", 
                        "pubmed_immune", 
                        "pubmed_aging", 
                        "pubmed_sex", 
                        "pubmed_cellcycle")
# 初始化结果列表
all_results <- data.frame()
# 设置一个很小的值，防止除0
epsilon <- 1e-6

#dir.create("/data/wuwq/noise/DISEASE_ALL/noise_var/raincloud_plots", showWarnings = FALSE)
set.seed(123) # Set seed for reproducibility in `geom_quasirandom`
for (column_name in columns_to_compare) {
  mean_yes <- mean(yes_group[[column_name]], na.rm = TRUE)
  mean_pubmed_query <- mean(pubmed_query_group[[column_name]], na.rm = TRUE)
  mean_cellcycle <- mean(pubmed_cellcycle_group[[column_name]], na.rm = TRUE)
  log2fc_1_2 <- log2((mean_yes + epsilon) / (mean_pubmed_query + epsilon))
  log2fc_1_3 <- log2((mean_yes + epsilon) / (mean_cellcycle + epsilon))
  log2fc_2_3 <- log2((mean_pubmed_query + epsilon) / (mean_cellcycle + epsilon))
  test_1_2 <- wilcox.test(yes_group[[column_name]], pubmed_query_group[[column_name]], exact = FALSE)
  test_1_3 <- wilcox.test(yes_group[[column_name]], pubmed_cellcycle_group[[column_name]], exact = FALSE)
  test_2_3 <- wilcox.test(pubmed_query_group[[column_name]], pubmed_cellcycle_group[[column_name]], exact = FALSE)
  p_values <- c(test_1_2$p.value, test_1_3$p.value, test_2_3$p.value)
  fdr_values <- p.adjust(p_values, method = "fdr")
  temp <- data.frame(
    Variable = column_name,
    Mean_Yes_Group = mean_yes,
    Mean_Pubmed_Query_Group = mean_pubmed_query,
    Mean_Cellcycle_Group = mean_cellcycle,
    Comparison = c("Yes vs PubmedQuery", "Yes vs Cellcycle", "PubmedQuery vs Cellcycle"),
    log2FC = c(log2fc_1_2, log2fc_1_3, log2fc_2_3),
    p_value = p_values,
    FDR = fdr_values
  )
  all_results <- rbind(all_results, temp)
  # Prepare data for plotting
  plot_data <- data.frame(
    value = c(yes_group[[column_name]],
              pubmed_query_group[[column_name]],
              pubmed_cellcycle_group[[column_name]]),
    group = factor(c(
      rep("Yes", length(yes_group[[column_name]])),
      rep("PubmedQuery", length(pubmed_query_group[[column_name]])),
      rep("Cellcycle", length(pubmed_cellcycle_group[[column_name]]))
    ), levels = c("Cellcycle", "Yes", "PubmedQuery" ))
  )
  plot_data <- plot_data %>% filter(!is.na(value))
  # Correcting the erroneous function name
  p <- ggplot(plot_data, aes(x = group, y = value, fill = group)) +
    geom_violin(trim = FALSE, scale = "width", alpha = 0.5) +
    geom_quasirandom(width = 0.2, alpha = 0.8, size = 1.5) +
    geom_boxplot(width = 0.1, position = position_dodge(0.9), 
                 alpha = 0.5, outlier.shape = NA) +
    scale_y_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    theme_classic() +
    labs(
      title = column_name,
      y = "Value",
      x = NULL
    ) +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5)
    )
  # Save the plot
  #ggsave(
  # filename = paste0("/data/wuwq/noise/DISEASE_ALL/noise_var/raincloud_plots/", column_name, "_otherT.pdf"),
  #  plot = p,
  #  width = 6,
  #  height = 5
  #)
}
### 保存统计数据
write.csv(yes_group, file="/data/wuwq/noise/DISEASE_ALL/noise_var/otherT_health_pubmed_count.csv", row.names = TRUE)
#write.csv(all_results, file="/data/wuwq/noise/DISEASE_ALL/noise_var/otherT_health_pubmed_wilcox.csv", row.names = TRUE)






##### (四-1) noi_gene 辅证---- 疾病与健康 初始，现已不用并更新见下方 ❌️
rm(list=ls())

### celltype & disease: "health" "RA" "SLE" "COVID-19" "post-COVID-19 disorder" "SEPSIS" "FLU" "NSCLC"
noi_gene <- read.csv("/data/wuwq/noise/DISEASE_ALL/noise_var/health_otherT.csv", header = TRUE, row.names = 1)
pubmed_query <- read.delim("/data/wuwq/noise/DISEASE_ALL/noise_var/pubmed_query.txt", header = TRUE)
disease_noi <- read.csv("/data/wuwq/noise/DISEASE_ALL/noise_var/NSCLC_otherT.csv", header = TRUE, row.names = 1)

# 健康新增一列 noi_gene
merged_data <- merge(noi_gene, pubmed_query, by = "gene")
merged_data$noi_gene <- ifelse(
  merged_data$noise_percentage == 0, 
  "zero",  # noise_percentage为0时，设置为zero
  ifelse(
    merged_data$noise_percentage >= 0.05, 
    "yes",  # noise_percentage大于等于0.05
    "less0.05"  # noise_percentage在0到0.05之间
  )
)
unique(merged_data$noi_gene)
col_order <- c("gene","noi_gene","pubmed_cellcycle","pubmed_apoptosis","pubmed_immune", "pubmed_aging","pubmed_sex",  setdiff(names(merged_data), c("gene","noi_gene","pubmed_cellcycle", "pubmed_apoptosis","pubmed_immune", "pubmed_aging","pubmed_sex")))
merged_data <- merged_data[ , col_order]
### 选择待检验的列
merged_data <- merged_data[ ,c(1:7)]
### 用0替换掉NA
merged_data[is.na(merged_data)] <- 0

# 疾病新增一列 noi_gene
merged_data_disease <- merge(disease_noi, pubmed_query, by = "gene")
merged_data_disease$noi_gene <- ifelse(
  merged_data_disease$noise_percentage == 0, 
  "zero",  # noise_percentage为0时，设置为zero
  ifelse(
    merged_data_disease$noise_percentage >= 0.05, 
    "yes",  # noise_percentage大于等于0.05
    "less0.05"  # noise_percentage在0到0.05之间
  )
)
unique(merged_data_disease$noi_gene)
col_order <- c("gene","noi_gene","pubmed_cellcycle","pubmed_apoptosis","pubmed_immune", "pubmed_aging","pubmed_sex",  setdiff(names(merged_data_disease), c("gene","noi_gene","pubmed_cellcycle", "pubmed_apoptosis","pubmed_immune", "pubmed_aging","pubmed_sex")))
merged_data_disease <- merged_data_disease[ , col_order]
### 选择待检验的列
merged_data_disease <- merged_data_disease[ ,c(1:7)]
### 用0替换掉NA
merged_data_disease[is.na(merged_data_disease)] <- 0

### 计算均值并检验
# 提取不同组合
yes_group <- merged_data %>% filter(noi_gene == "yes")
pubmed_query_group <- pubmed_query
yes_group_disease <- merged_data_disease %>% filter(noi_gene == "yes")
# 定义要比较的列
columns_to_compare <- c("pubmed_apoptosis", 
                        "pubmed_immune", 
                        "pubmed_aging", 
                        "pubmed_sex", 
                        "pubmed_cellcycle")
# 初始化结果列表
all_results <- data.frame()
# 设置一个很小的值，防止除0
epsilon <- 1e-6
#dir.create("/data/wuwq/noise/DISEASE_ALL/noise_var/raincloud_plots", showWarnings = FALSE)
set.seed(123) # Set seed for reproducibility in `geom_quasirandom`
for (column_name in columns_to_compare) {
  mean_health <- mean(yes_group[[column_name]], na.rm = TRUE)
  mean_pubmed_query <- mean(pubmed_query_group[[column_name]], na.rm = TRUE)
  mean_disease <- mean(yes_group_disease[[column_name]], na.rm = TRUE)
  log2fc_1_2 <- log2((mean_health + epsilon) / (mean_pubmed_query + epsilon))
  log2fc_1_3 <- log2((mean_health + epsilon) / (mean_disease + epsilon))
  log2fc_2_3 <- log2((mean_pubmed_query + epsilon) / (mean_disease + epsilon))
  test_1_2 <- wilcox.test(yes_group[[column_name]], pubmed_query_group[[column_name]], exact = FALSE)
  test_1_3 <- wilcox.test(yes_group[[column_name]], yes_group_disease[[column_name]], exact = FALSE)
  test_2_3 <- wilcox.test(pubmed_query_group[[column_name]], yes_group_disease[[column_name]], exact = FALSE)
  p_values <- c(test_1_2$p.value, test_1_3$p.value, test_2_3$p.value)
  fdr_values <- p.adjust(p_values, method = "fdr")
  temp <- data.frame(
    Variable = column_name,
    Mean_health_Group = mean_health,
    Mean_Pubmed_Query_Group = mean_pubmed_query,
    Mean_disease_Group = mean_disease,
    Comparison = c("health vs control", "health vs disease", "control vs disease"),
    log2FC = c(log2fc_1_2, log2fc_1_3, log2fc_2_3),
    p_value = p_values,
    FDR = fdr_values
  )
  all_results <- rbind(all_results, temp)
  # Prepare data for plotting
  plot_data <- data.frame(
    value = c(yes_group[[column_name]],
              pubmed_query_group[[column_name]],
              yes_group_disease[[column_name]]),
    group = factor(c(
      rep("health", length(yes_group[[column_name]])),
      rep("control", length(pubmed_query_group[[column_name]])),
      rep("disease", length(yes_group_disease[[column_name]]))
    ), levels = c("disease", "health", "control" ))
  )
  plot_data <- plot_data %>% filter(!is.na(value))
  # Correcting the erroneous function name
  p <- ggplot(plot_data, aes(x = group, y = value, fill = group)) +
    geom_violin(trim = FALSE, scale = "width", alpha = 0.5) +
    geom_quasirandom(width = 0.2, alpha = 0.8, size = 1.5) +
    geom_boxplot(width = 0.1, position = position_dodge(0.9), 
                 alpha = 0.5, outlier.shape = NA) +
    scale_y_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    theme_classic() +
    labs(
      title = column_name,
      y = "Value",
      x = NULL
    ) +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5)
    )
  # Save the plot
  #ggsave(
  # filename = paste0("/data/wuwq/noise/DISEASE_ALL/noise_var/raincloud_plots/", column_name, "_otherT.pdf"),
  #  plot = p,
  #  width = 6,
  #  height = 5
  #)
}
### 保存统计数据
#write.csv(yes_group_disease, file="/data/wuwq/noise/DISEASE_ALL/noise_var/otherT_NSCLC_pubmed_count.csv", row.names = TRUE)
write.csv(all_results, file="/data/wuwq/noise/DISEASE_ALL/noise_var/otherT_NSCLC_pubmed_wilcox.csv", row.names = TRUE)





##### (四-2) noi_gene 辅证---- 疾病与健康  定义生物学过程  ✅️ 这个是初始尝试，优化见下方
rm(list=ls())

### celltype & disease: "health" "RA" "SLE" "COVID-19" "post-COVID-19 disorder" "SEPSIS" "FLU" "NSCLC"
noi_gene <- read.csv("/data/wuwq/noise/DISEASE_ALL/noise_var/health_otherT.csv", header = TRUE, row.names = 1) #读取健康图谱
disease_noi <- read.csv("/data/wuwq/noise/DISEASE_ALL/noise_var/NSCLC_otherT.csv", header = TRUE, row.names = 1) #读取疾病
pubmed_query <- read.delim("/data/wuwq/noise/DISEASE_ALL/noise_var/pubmed_query.txt", header = TRUE)

# 健康新增一列 noi_gene
merged_data <- merge(noi_gene, pubmed_query, by = "gene")
merged_data$noi_gene <- ifelse(
  merged_data$noise_percentage == 0, 
  "zero",  # noise_percentage为0时，设置为zero
  ifelse(
    merged_data$noise_percentage >= 0.05, 
    "yes",  # noise_percentage大于等于0.05
    "less0.05"  # noise_percentage在0到0.05之间
  )
)
unique(merged_data$noi_gene)
col_order <- c("gene","noi_gene","pubmed_cellcycle","pubmed_apoptosis","pubmed_immune", "pubmed_aging","pubmed_sex",  
               "pubmed_cellsignal_cellcommunication","pubmed_metabolism","pubmed_cellmetabolism","pubmed_drugmetabolism",
               "pubmed_HumoralImmunity","pubmed_CellMediatedImmunity","pubmed_InnateImmunity","pubmed_AdaptiveImmunity","pubmed_ImmuneTolerance",
               setdiff(names(merged_data), c("gene","noi_gene","pubmed_cellcycle", "pubmed_apoptosis","pubmed_immune", "pubmed_aging","pubmed_sex", "pubmed_cellsignal_cellcommunication","pubmed_metabolism","pubmed_cellmetabolism","pubmed_drugmetabolism","pubmed_HumoralImmunity","pubmed_CellMediatedImmunity","pubmed_InnateImmunity","pubmed_AdaptiveImmunity","pubmed_ImmuneTolerance")))
merged_data <- merged_data[ , col_order]
### 选择待检验的列
merged_data <- merged_data[ ,c(1:16)]
### 用0替换掉NA
merged_data[is.na(merged_data)] <- 0


# 疾病新增一列 noi_gene
merged_data_disease <- merge(disease_noi, pubmed_query, by = "gene")
merged_data_disease$noi_gene <- ifelse(
  merged_data_disease$noise_percentage == 0, 
  "zero",  # noise_percentage为0时，设置为zero
  ifelse(
    merged_data_disease$noise_percentage >= 0.05, 
    "yes",  # noise_percentage大于等于0.05
    "less0.05"  # noise_percentage在0到0.05之间
  )
)
unique(merged_data_disease$noi_gene)
col_order <- c("gene","noi_gene","pubmed_cellcycle","pubmed_apoptosis","pubmed_immune", "pubmed_aging","pubmed_sex",  
               "pubmed_cellsignal_cellcommunication","pubmed_metabolism","pubmed_cellmetabolism","pubmed_drugmetabolism",
               "pubmed_HumoralImmunity","pubmed_CellMediatedImmunity","pubmed_InnateImmunity","pubmed_AdaptiveImmunity","pubmed_ImmuneTolerance",
               setdiff(names(merged_data), c("gene","noi_gene","pubmed_cellcycle", "pubmed_apoptosis","pubmed_immune", "pubmed_aging","pubmed_sex", "pubmed_cellsignal_cellcommunication","pubmed_metabolism","pubmed_cellmetabolism","pubmed_drugmetabolism","pubmed_HumoralImmunity","pubmed_CellMediatedImmunity","pubmed_InnateImmunity","pubmed_AdaptiveImmunity","pubmed_ImmuneTolerance")))
merged_data_disease <- merged_data_disease[ , col_order]
### 选择待检验的列
merged_data_disease <- merged_data_disease[ ,c(1:16)]
### 用0替换掉NA
merged_data_disease[is.na(merged_data_disease)] <- 0

### 计算均值并检验
# 提取不同组合
yes_group <- merged_data %>% filter(noi_gene == "yes")
pubmed_query_group <- pubmed_query
yes_group_disease <- merged_data_disease %>% filter(noi_gene == "yes")
# 定义要比较的列
columns_to_compare <- c("pubmed_cellcycle", 
                        "pubmed_apoptosis", 
                        "pubmed_immune", 
                        "pubmed_aging", 
                        "pubmed_sex",
                        "pubmed_cellsignal_cellcommunication",
                        "pubmed_metabolism",
                        "pubmed_cellmetabolism",
                        "pubmed_drugmetabolism",
                        "pubmed_HumoralImmunity",
                        "pubmed_CellMediatedImmunity",
                        "pubmed_InnateImmunity",
                        "pubmed_AdaptiveImmunity",
                        "pubmed_ImmuneTolerance"
                        )
# 初始化结果列表
all_results <- data.frame()
# 设置一个很小的值，防止除0
epsilon <- 1e-6
#dir.create("/data/wuwq/noise/DISEASE_ALL/noise_var/raincloud_plots", showWarnings = FALSE)
set.seed(123) # Set seed for reproducibility in `geom_quasirandom`
for (column_name in columns_to_compare) {
  # 计算均值
  mean_health <- mean(yes_group[[column_name]], na.rm = TRUE)
  mean_pubmed_query <- mean(pubmed_query_group[[column_name]], na.rm = TRUE)
  mean_disease <- mean(yes_group_disease[[column_name]], na.rm = TRUE)
  # 计算标准误
  se_health <- sd(yes_group[[column_name]], na.rm = TRUE) / sqrt(sum(!is.na(yes_group[[column_name]])))
  se_pubmed_query <- sd(pubmed_query_group[[column_name]], na.rm = TRUE) / sqrt(sum(!is.na(pubmed_query_group[[column_name]])))
  se_disease <- sd(yes_group_disease[[column_name]], na.rm = TRUE) / sqrt(sum(!is.na(yes_group_disease[[column_name]])))
  # 计算log2 fold change
  log2fc_1_2 <- log2((mean_health + epsilon) / (mean_pubmed_query + epsilon))
  log2fc_3_1 <- log2((mean_disease + epsilon) / (mean_health + epsilon))
  log2fc_3_2 <- log2((mean_disease + epsilon) / (mean_pubmed_query + epsilon))
  # wilcox检验
  test_1_2 <- wilcox.test(yes_group[[column_name]], pubmed_query_group[[column_name]], exact = FALSE)
  test_3_1 <- wilcox.test(yes_group_disease[[column_name]], yes_group[[column_name]], exact = FALSE)
  test_3_2 <- wilcox.test(yes_group_disease[[column_name]], pubmed_query_group[[column_name]], exact = FALSE)
  p_values <- c(test_1_2$p.value, test_3_1$p.value, test_3_2$p.value)
  fdr_values <- p.adjust(p_values, method = "fdr")
  # 将均值和SE都添加到结果中
  temp <- data.frame(
    Variable = column_name,
    Mean_health_Group = mean_health,
    SE_health_Group = se_health,
    Mean_Pubmed_Query_Group = mean_pubmed_query,
    SE_Pubmed_Query_Group = se_pubmed_query,
    Mean_disease_Group = mean_disease,
    SE_disease_Group = se_disease,
    Comparison = c("health vs control", "disease vs health", "disease vs control"),
    log2FC = c(log2fc_1_2, log2fc_3_1, log2fc_3_2),
    p_value = p_values,
    FDR = fdr_values,
    logtrans_FDR = -log10(fdr_values)
  )
  all_results <- rbind(all_results, temp)
  # 准备绘图数据
  plot_data <- data.frame(
    value = c(yes_group[[column_name]],
              pubmed_query_group[[column_name]],
              yes_group_disease[[column_name]]),
    group = factor(c(
      rep("health", length(yes_group[[column_name]])),
      rep("control", length(pubmed_query_group[[column_name]])),
      rep("disease", length(yes_group_disease[[column_name]]))
    ), levels = c("disease", "health", "control"))
  )
  
  plot_data <- plot_data %>% filter(!is.na(value))
  # 绘图
  p <- ggplot(plot_data, aes(x = group, y = value, fill = group)) +
    geom_violin(trim = FALSE, scale = "width", alpha = 0.5) +
    geom_quasirandom(width = 0.2, alpha = 0.8, size = 1.5) +
    geom_boxplot(width = 0.1, position = position_dodge(0.9), alpha = 0.5, outlier.shape = NA) +
    scale_y_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    theme_classic() +
    labs(title = column_name, y = "Value", x = NULL) +
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
  # 保存绘图（可选）
  # ggsave(
  #   filename = paste0("/data/wuwq/noise/DISEASE_ALL/noise_var/raincloud_plots/", column_name, "_otherT.pdf"),
  #   plot = p,
  #   width = 6,
  #   height = 5
  # )
}
### 保存统计数据  yes_group <- merged_data %>% filter(noi_gene == "yes")
#write.csv(yes_group, file="/data/wuwq/noise/DISEASE_ALL/noise_var/otherT_health_pubmed_count.csv", row.names = TRUE)
#write.csv(yes_group_disease, file="/data/wuwq/noise/DISEASE_ALL/noise_var/otherT_NSCLC_pubmed_count.csv", row.names = TRUE)
write.csv(all_results, file="/data/wuwq/noise/DISEASE_ALL/noise_var/otherT_NSCLC_pubmed_wilcox.csv", row.names = TRUE)






rm(list=ls())
##### (五-1) noi_gene 辅证---- 疾病与健康  定义生物学过程  ✅️✅️ 这个是优化尝试v2
### celltype & disease: "health" "RA" "SLE" "COVID-19" "post-COVID-19 disorder" "SEPSIS" "FLU" "NSCLC"

noi_gene <- read.csv("/data/wuwq/noise/DISEASE_ALL/noise_var/health_otherT.csv", header = TRUE, row.names = 1) #读取健康图谱
disease_noi <- read.csv("/data/wuwq/noise/DISEASE_ALL/noise_var/NSCLC_otherT.csv", header = TRUE, row.names = 1) #读取疾病
pubmed_query <- read.delim("/data/wuwq/noise/DISEASE_ALL/noise_var/pubmed_query.txt", header = TRUE)

### 将pubmed基因文献数匹配至健康group
merged_data <- merge(noi_gene, pubmed_query, by = "gene")
col_order <- c("gene","mean_noise","noise_percentage","pubmed_cellcycle","pubmed_apoptosis","pubmed_immune", "pubmed_aging","pubmed_sex",  
               "pubmed_cellsignal_cellcommunication","pubmed_metabolism","pubmed_cellmetabolism","pubmed_drugmetabolism",
               "pubmed_HumoralImmunity","pubmed_CellMediatedImmunity","pubmed_InnateImmunity","pubmed_AdaptiveImmunity","pubmed_ImmuneTolerance",
               setdiff(names(merged_data), c("gene","mean_noise","noise_percentage","pubmed_cellcycle", "pubmed_apoptosis","pubmed_immune", "pubmed_aging","pubmed_sex", "pubmed_cellsignal_cellcommunication","pubmed_metabolism","pubmed_cellmetabolism","pubmed_drugmetabolism","pubmed_HumoralImmunity","pubmed_CellMediatedImmunity","pubmed_InnateImmunity","pubmed_AdaptiveImmunity","pubmed_ImmuneTolerance")))
merged_data <- merged_data[ , col_order]
merged_data <- merged_data[ ,c(1:17)]
# 用0替换掉NA
merged_data[is.na(merged_data)] <- 0

### 将pubmed基因文献数匹配至疾病group
merged_data_disease <- merge(disease_noi, pubmed_query, by = "gene")
col_order <- c("gene","mean_noise","noise_percentage","pubmed_cellcycle","pubmed_apoptosis","pubmed_immune", "pubmed_aging","pubmed_sex",  
               "pubmed_cellsignal_cellcommunication","pubmed_metabolism","pubmed_cellmetabolism","pubmed_drugmetabolism",
               "pubmed_HumoralImmunity","pubmed_CellMediatedImmunity","pubmed_InnateImmunity","pubmed_AdaptiveImmunity","pubmed_ImmuneTolerance",
               setdiff(names(merged_data_disease), c("gene","mean_noise","noise_percentage","pubmed_cellcycle", "pubmed_apoptosis","pubmed_immune", "pubmed_aging","pubmed_sex", "pubmed_cellsignal_cellcommunication","pubmed_metabolism","pubmed_cellmetabolism","pubmed_drugmetabolism","pubmed_HumoralImmunity","pubmed_CellMediatedImmunity","pubmed_InnateImmunity","pubmed_AdaptiveImmunity","pubmed_ImmuneTolerance")))
merged_data_disease <- merged_data_disease[ , col_order]
merged_data_disease <- merged_data_disease[ ,c(1:17)]
# 用0替换掉NA
merged_data_disease[is.na(merged_data_disease)] <- 0

### 计算均值并检验
# 提取不同筛选标准组合 
### v2:>5% & noise水平降序top100
### v3:占比>20% & noise水平降序top500 (限制最大数目500个)
# pubmed基线全部要
pubmed_query_group <- pubmed_query
# 健康组
health_group <- merged_data[merged_data$noise_percentage > 0.2, ]
yes_group <- health_group %>% arrange(desc(mean_noise)) %>% head(500)
# 疾病组
disease_group <- merged_data_disease[merged_data_disease$noise_percentage > 0.2, ]
yes_group_disease <- disease_group %>% arrange(desc(mean_noise)) %>% head(500)
# 定义要比较的列
columns_to_compare <- c("pubmed_cellmetabolism",
                        "pubmed_drugmetabolism",
                        "pubmed_HumoralImmunity",
                        "pubmed_CellMediatedImmunity",
                        "pubmed_InnateImmunity",
                        "pubmed_AdaptiveImmunity",
                        "pubmed_ImmuneTolerance",
                        "pubmed_cellsignal_cellcommunication",
                        "pubmed_immune",
                        "pubmed_metabolism"
)
# 初始化结果列表
all_results <- data.frame()
# 设置一个很小的值，防止除以0
epsilon <- 1e-6
set.seed(123) # Set seed for reproducibility in `geom_quasirandom`
for (column_name in columns_to_compare) {
  # 计算均值
  mean_health <- mean(yes_group[[column_name]], na.rm = TRUE)
  mean_pubmed_query <- mean(pubmed_query_group[[column_name]], na.rm = TRUE)
  mean_disease <- mean(yes_group_disease[[column_name]], na.rm = TRUE)
  # 计算标准误
  se_health <- sd(yes_group[[column_name]], na.rm = TRUE) / sqrt(sum(!is.na(yes_group[[column_name]])))
  se_pubmed_query <- sd(pubmed_query_group[[column_name]], na.rm = TRUE) / sqrt(sum(!is.na(pubmed_query_group[[column_name]])))
  se_disease <- sd(yes_group_disease[[column_name]], na.rm = TRUE) / sqrt(sum(!is.na(yes_group_disease[[column_name]])))
  # 计算log2 fold change
  log2fc_1_2 <- log2((mean_health + epsilon) / (mean_pubmed_query + epsilon))
  log2fc_3_1 <- log2((mean_disease + epsilon) / (mean_health + epsilon))
  log2fc_3_2 <- log2((mean_disease + epsilon) / (mean_pubmed_query + epsilon))
  # wilcox检验
  test_1_2 <- wilcox.test(yes_group[[column_name]], pubmed_query_group[[column_name]], exact = FALSE)
  test_3_1 <- wilcox.test(yes_group_disease[[column_name]], yes_group[[column_name]], exact = FALSE)
  test_3_2 <- wilcox.test(yes_group_disease[[column_name]], pubmed_query_group[[column_name]], exact = FALSE)
  p_values <- c(test_1_2$p.value, test_3_1$p.value, test_3_2$p.value)
  fdr_values <- p.adjust(p_values, method = "fdr")
  # 将均值和SE都添加到结果中
  temp <- data.frame(
    Variable = column_name,
    Mean_health_Group = mean_health,
    SE_health_Group = se_health,
    Mean_Pubmed_Query_Group = mean_pubmed_query,
    SE_Pubmed_Query_Group = se_pubmed_query,
    Mean_disease_Group = mean_disease,
    SE_disease_Group = se_disease,
    Comparison = c("health vs control", "disease vs health", "disease vs control"),
    log2FC = c(log2fc_1_2, log2fc_3_1, log2fc_3_2),
    p_value = p_values,
    FDR = fdr_values,
    logtrans_FDR = -log10(fdr_values)
  )
  all_results <- rbind(all_results, temp)
}
unique(all_results$Comparison)
all_results <- all_results[all_results$Comparison == "disease vs health", ]

### 保存统计数据  
# write.csv(all_results, file="/data/wuwq/noise/DISEASE_ALL/noise_var/pubmed_top100wilcox_v2/NSCLC_otherT.csv", row.names = TRUE)
write.csv(all_results, file="/data/wuwq/noise/DISEASE_ALL/noise_var/pubmed_per20top500wilcox_v3/NSCLC_otherT.csv", row.names = TRUE)






rm(list=ls())
##### (五-2) noi_gene 辅证---- 疾病与健康  定义生物学过程  ✅️   可视化
### 示例
library(tidyverse)
library(viridis)
# 模拟数据集
data <- data.frame(
  individual=paste( "Mister ", seq(1,60), sep=""),
  group=c( rep('A', 10), rep('B', 30), rep('C', 14), rep('D', 6)) ,
  value1=sample( seq(10,100), 60, replace=T),
  value2=sample( seq(10,100), 60, replace=T),
  value3=sample( seq(10,100), 60, replace=T)
)

# 把数据转化成统一的格式
data <- data %>% gather(key = "observation", value="value", -c(1,2)) 

# 再每组末尾添加空白距离
empty_bar <- 2
nObsType <- nlevels(as.factor(data$observation))
to_add <- data.frame( matrix(NA, empty_bar*nlevels(data$group)*nObsType, ncol(data)) )
colnames(to_add) <- colnames(data)
to_add$group <- rep(levels(data$group), each=empty_bar*nObsType )
data <- rbind(data, to_add)
data <- data %>% arrange(group, individual)
data$id <- rep( seq(1, nrow(data)/nObsType) , each=nObsType)

#获取每个标签的名称和y轴的位置
label_data <- data %>% group_by(id, individual) %>% summarize(tot=sum(value))
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)

# 为基线准备一个数据帧
base_data <- data %>% 
  group_by(group) %>% 
  summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))

# 为grid准备一个数据帧
grid_data <- base_data
grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start <- grid_data$start - 1
grid_data <- grid_data[-1,]

# 绘图
p <- ggplot(data) +      
  # 添加堆叠的条形图
  geom_bar(aes(x=as.factor(id), y=value, fill=observation), stat="identity", alpha=0.5) +
  scale_fill_viridis(discrete=TRUE) +
  # Add a val=100/75/50/25 lines.
  geom_segment(data=grid_data, aes(x = end, y = 0, xend = start, yend = 0), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 100, xend = start, yend = 100), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 200, xend = start, yend = 200), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 300, xend = start, yend = 300), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 400, xend = start, yend = 400), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  # 添加显示100/75/50/25的文本
  ggplot2::annotate("text", x = rep(max(data$id),5), y = c(0, 100, 200, 300, 400), label = c("0", "100", "200", "300", "400") , color="grey", size=10 , angle=0, fontface="bold", hjust=1) +
  ylim(-150,max(label_data$tot, na.rm=T)) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm") 
  ) +
  coord_polar() +
  # 在每个条形图顶部添加一个标签
  geom_text(data=label_data, aes(x=id, y=tot+10, label=individual, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=5, angle= label_data$angle, inherit.aes = FALSE ) +
  # 添加基线信息
  geom_segment(data=base_data, aes(x = start, y = -5, xend = end, yend = -5), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE )  +
  geom_text(data=base_data, aes(x = title, y = -18, label=group), hjust=c(1,1,0,0), colour = "black", alpha=0.8, size=4, fontface="bold", inherit.aes = FALSE)
p

#输出图片
#ggsave(p, file="output.png", width=10, height=10)
#write.csv(data, file="/data/wuwq/noise/DISEASE_ALL/noise_var/pubmed_per20top500wilcox_v3/plot/example_data.csv")



rm(list=ls())
######## 真实数据尝试
# 读取数据
library(dplyr)
library(ggplot2)

rm(list=ls())
################################################################################################################# 1.RA
my_data <- read.delim("/data/wuwq/noise/DISEASE_ALL/noise_var/pubmed_per20top500wilcox_v3/plot/data_RA.txt")

# 添加每组末尾空白
my_empty_bar <- 2
my_nObsType <- nlevels(as.factor(my_data$observation))
my_to_add <- data.frame(matrix(NA, my_empty_bar * nlevels(my_data$group) * my_nObsType, ncol = ncol(my_data)))
colnames(my_to_add) <- colnames(my_data)
my_to_add$group <- rep(levels(my_data$group), each = my_empty_bar * my_nObsType)
my_data <- rbind(my_data, my_to_add)
my_data <- my_data %>% arrange(group, Variable)
my_data$id <- rep(seq(1, nrow(my_data) / my_nObsType), each = my_nObsType)

# 计算每个标签的总和
my_label_data <- my_data %>% 
  group_by(id, Variable) %>% 
  summarize(tot = sum(value, na.rm = TRUE), .groups = "drop")

my_number_of_bar <- nrow(my_label_data)
my_angle <- 90 - 360 * (my_label_data$id - 0.5) / my_number_of_bar
my_label_data <- my_label_data %>%
  mutate(hjust = ifelse(my_angle < -90, 1, 0),
         angle = ifelse(my_angle < -90, my_angle + 180, my_angle))

# 基线数据
my_base_data <- my_data %>% 
  group_by(group) %>% 
  summarize(start = min(id), end = max(id) - my_empty_bar, .groups = "drop") %>% 
  rowwise() %>% 
  mutate(title = mean(c(start, end)))

# grid 数据
my_grid_data <- my_base_data
my_grid_data$end <- my_grid_data$end[c(nrow(my_grid_data), 1:(nrow(my_grid_data) - 1))] + 1
my_grid_data$start <- my_grid_data$start - 1
my_grid_data <- my_grid_data[-1,]

# 计算每组基线标签的水平对齐
my_base_data <- my_base_data %>%
  mutate(hjust = (start + end) / 2 / max(my_data$id))

# -----------------------------
# 设置 observation 的堆叠顺序和颜色
# 内圈从 CD4 到外圈 otherT
my_data$observation <- factor(my_data$observation,
                              levels = c("CD4", "CD8", "MNP", "NK", "B", "otherT"))

obs_colors <- c(
  otherT    = "#FDE725FF",
  B    = "#7AD151FF",
  NK    = "#22A884FF",
  MNP     = "#2A788EFF",
  CD8      = "#414487FF",
  CD4 = "#440154FF"
)

# -----------------------------
# 绘图
my_p <- ggplot(my_data) +
  geom_bar(aes(x = as.factor(id), y = value*100, fill = observation), 
           stat = "identity", alpha = 0.5) +
  scale_fill_manual(values = obs_colors) +  # 使用自定义颜色
  # grid lines
  geom_segment(data = my_grid_data, aes(x = end, y = 0, xend = start, yend = 0), colour = "grey", size = 0.3, inherit.aes = FALSE) +
  geom_segment(data = my_grid_data, aes(x = end, y = 200, xend = start, yend = 200), colour = "grey", size = 0.3, inherit.aes = FALSE) +
  geom_segment(data = my_grid_data, aes(x = end, y = 400, xend = start, yend = 400), colour = "grey", size = 0.3, inherit.aes = FALSE) +
  geom_segment(data = my_grid_data, aes(x = end, y = 600, xend = start, yend = 600), colour = "grey", size = 0.3, inherit.aes = FALSE) +
  geom_segment(data = my_grid_data, aes(x = end, y = 800, xend = start, yend = 800), colour = "grey", size = 0.3, inherit.aes = FALSE) +
  # y轴文本
  annotate("text", x = rep(max(my_data$id), 5), y = c(0, 200, 400, 600, 800), 
           label = c("0", "200", "400", "600", "800"), color = "grey", size = 8, fontface = "bold", hjust = 1) +
  # 设置y轴范围
  expand_limits(y = max(my_label_data$tot, na.rm = TRUE) + 50) +
  theme_minimal() +
  theme(
    legend.position = "right",  # 将图例放在右侧
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1, 4), "cm")
  ) +
  coord_polar() +
  # 基线信息
  geom_segment(data = my_base_data, aes(x = start, y = -5, xend = end, yend = -5),
               colour = "black", size = 0.6, inherit.aes = FALSE) +
  geom_text(data = my_base_data, aes(x = title, y = -18, label = group, hjust = hjust),
            colour = "black", size = 4, fontface = "bold", inherit.aes = FALSE)

# 显示图形
my_p
ggsave(my_p, file="/data/wuwq/noise/DISEASE_ALL/noise_var/pubmed_per20top500wilcox_v3/plot/plot_RA_2.svg", width=10, height=10)


rm(list=ls())
################################################################################################################# 2.SLE
my_data <- read.delim("/data/wuwq/noise/DISEASE_ALL/noise_var/pubmed_per20top500wilcox_v3/plot/data_SLE.txt")

# 添加每组末尾空白
my_empty_bar <- 2
my_nObsType <- nlevels(as.factor(my_data$observation))
my_to_add <- data.frame(matrix(NA, my_empty_bar * nlevels(my_data$group) * my_nObsType, ncol = ncol(my_data)))
colnames(my_to_add) <- colnames(my_data)
my_to_add$group <- rep(levels(my_data$group), each = my_empty_bar * my_nObsType)
my_data <- rbind(my_data, my_to_add)
my_data <- my_data %>% arrange(group, Variable)
my_data$id <- rep(seq(1, nrow(my_data) / my_nObsType), each = my_nObsType)

# 计算每个标签的总和
my_label_data <- my_data %>% 
  group_by(id, Variable) %>% 
  summarize(tot = sum(value, na.rm = TRUE), .groups = "drop")

my_number_of_bar <- nrow(my_label_data)
my_angle <- 90 - 360 * (my_label_data$id - 0.5) / my_number_of_bar
my_label_data <- my_label_data %>%
  mutate(hjust = ifelse(my_angle < -90, 1, 0),
         angle = ifelse(my_angle < -90, my_angle + 180, my_angle))

# 基线数据
my_base_data <- my_data %>% 
  group_by(group) %>% 
  summarize(start = min(id), end = max(id) - my_empty_bar, .groups = "drop") %>% 
  rowwise() %>% 
  mutate(title = mean(c(start, end)))

# grid 数据
my_grid_data <- my_base_data
my_grid_data$end <- my_grid_data$end[c(nrow(my_grid_data), 1:(nrow(my_grid_data) - 1))] + 1
my_grid_data$start <- my_grid_data$start - 1
my_grid_data <- my_grid_data[-1,]

# 计算每组基线标签的水平对齐
my_base_data <- my_base_data %>%
  mutate(hjust = (start + end) / 2 / max(my_data$id))

# -----------------------------
# 设置 observation 的堆叠顺序和颜色
# 内圈从 CD4 到外圈 otherT
my_data$observation <- factor(my_data$observation,
                              levels = c("CD4", "CD8", "MNP", "NK", "B", "otherT"))

obs_colors <- c(
  otherT    = "#FDE725FF",
  B    = "#7AD151FF",
  NK    = "#22A884FF",
  MNP     = "#2A788EFF",
  CD8      = "#414487FF",
  CD4 = "#440154FF"
)

# -----------------------------
# 绘图
my_p <- ggplot(my_data) +
  geom_bar(aes(x = as.factor(id), y = value*100, fill = observation), 
           stat = "identity", alpha = 0.5) +
  scale_fill_manual(values = obs_colors) +  # 使用自定义颜色
  # grid lines
  geom_segment(data = my_grid_data, aes(x = end, y = 0, xend = start, yend = 0), colour = "grey", size = 0.3, inherit.aes = FALSE) +
  geom_segment(data = my_grid_data, aes(x = end, y = 55, xend = start, yend = 55), colour = "grey", size = 0.3, inherit.aes = FALSE) +
  geom_segment(data = my_grid_data, aes(x = end, y = 110, xend = start, yend = 110), colour = "grey", size = 0.3, inherit.aes = FALSE) +
  geom_segment(data = my_grid_data, aes(x = end, y = 165, xend = start, yend = 165), colour = "grey", size = 0.3, inherit.aes = FALSE) +
  geom_segment(data = my_grid_data, aes(x = end, y = 220, xend = start, yend = 220), colour = "grey", size = 0.3, inherit.aes = FALSE) +
  # y轴文本
  annotate("text", x = rep(max(my_data$id), 5), y = c(0, 55, 110, 165, 220), 
           label = c("0", "55", "110", "165", "220"), color = "grey", size = 8, fontface = "bold", hjust = 1) +
  # 设置y轴范围
  expand_limits(y = max(my_label_data$tot, na.rm = TRUE) + 50) +
  theme_minimal() +
  theme(
    legend.position = "right",  # 将图例放在右侧
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1, 4), "cm")
  ) +
  coord_polar() +
  # 基线信息
  geom_segment(data = my_base_data, aes(x = start, y = -5, xend = end, yend = -5),
               colour = "black", size = 0.6, inherit.aes = FALSE) +
  geom_text(data = my_base_data, aes(x = title, y = -18, label = group, hjust = hjust),
            colour = "black", size = 4, fontface = "bold", inherit.aes = FALSE)

# 显示图形
my_p
ggsave(my_p, file="/data/wuwq/noise/DISEASE_ALL/noise_var/pubmed_per20top500wilcox_v3/plot/plot_SLE_2.svg", width=10, height=10)


rm(list=ls())
################################################################################################################# 3.COVID
my_data <- read.delim("/data/wuwq/noise/DISEASE_ALL/noise_var/pubmed_per20top500wilcox_v3/plot/data_COVID.txt")

# 添加每组末尾空白
my_empty_bar <- 2
my_nObsType <- nlevels(as.factor(my_data$observation))
my_to_add <- data.frame(matrix(NA, my_empty_bar * nlevels(my_data$group) * my_nObsType, ncol = ncol(my_data)))
colnames(my_to_add) <- colnames(my_data)
my_to_add$group <- rep(levels(my_data$group), each = my_empty_bar * my_nObsType)
my_data <- rbind(my_data, my_to_add)
my_data <- my_data %>% arrange(group, Variable)
my_data$id <- rep(seq(1, nrow(my_data) / my_nObsType), each = my_nObsType)

# 计算每个标签的总和
my_label_data <- my_data %>% 
  group_by(id, Variable) %>% 
  summarize(tot = sum(value, na.rm = TRUE), .groups = "drop")

my_number_of_bar <- nrow(my_label_data)
my_angle <- 90 - 360 * (my_label_data$id - 0.5) / my_number_of_bar
my_label_data <- my_label_data %>%
  mutate(hjust = ifelse(my_angle < -90, 1, 0),
         angle = ifelse(my_angle < -90, my_angle + 180, my_angle))

# 基线数据
my_base_data <- my_data %>% 
  group_by(group) %>% 
  summarize(start = min(id), end = max(id) - my_empty_bar, .groups = "drop") %>% 
  rowwise() %>% 
  mutate(title = mean(c(start, end)))

# grid 数据
my_grid_data <- my_base_data
my_grid_data$end <- my_grid_data$end[c(nrow(my_grid_data), 1:(nrow(my_grid_data) - 1))] + 1
my_grid_data$start <- my_grid_data$start - 1
my_grid_data <- my_grid_data[-1,]

# 计算每组基线标签的水平对齐
my_base_data <- my_base_data %>%
  mutate(hjust = (start + end) / 2 / max(my_data$id))

# -----------------------------
# 设置 observation 的堆叠顺序和颜色
# 内圈从 CD4 到外圈 otherT
my_data$observation <- factor(my_data$observation,
                              levels = c("CD4", "CD8", "MNP", "NK", "B", "otherT"))

obs_colors <- c(
  otherT    = "#FDE725FF",
  B    = "#7AD151FF",
  NK    = "#22A884FF",
  MNP     = "#2A788EFF",
  CD8      = "#414487FF",
  CD4 = "#440154FF"
)


# -----------------------------
# 绘图
my_p <- ggplot(my_data) +
  geom_bar(aes(x = as.factor(id), y = value*100, fill = observation), 
           stat = "identity", alpha = 0.5) +
  scale_fill_manual(values = obs_colors) +  # 使用自定义颜色
  # grid lines
  geom_segment(data = my_grid_data, aes(x = end, y = 0, xend = start, yend = 0), colour = "grey", size = 0.3, inherit.aes = FALSE) +
  geom_segment(data = my_grid_data, aes(x = end, y = 60, xend = start, yend = 60), colour = "grey", size = 0.3, inherit.aes = FALSE) +
  geom_segment(data = my_grid_data, aes(x = end, y = 120, xend = start, yend = 120), colour = "grey", size = 0.3, inherit.aes = FALSE) +
  geom_segment(data = my_grid_data, aes(x = end, y = 180, xend = start, yend = 180), colour = "grey", size = 0.3, inherit.aes = FALSE) +
  geom_segment(data = my_grid_data, aes(x = end, y = 240, xend = start, yend = 240), colour = "grey", size = 0.3, inherit.aes = FALSE) +
  # y轴文本
  annotate("text", x = rep(max(my_data$id), 5), y = c(0, 60, 120, 180, 240), 
           label = c("0", "60", "120", "180", "240"), color = "grey", size = 8, fontface = "bold", hjust = 1) +
  # 设置y轴范围
  expand_limits(y = max(my_label_data$tot, na.rm = TRUE) + 50) +
  theme_minimal() +
  theme(
    legend.position = "right",  # 将图例放在右侧
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1, 4), "cm")
  ) +
  coord_polar() +
  # 基线信息
  geom_segment(data = my_base_data, aes(x = start, y = -5, xend = end, yend = -5),
               colour = "black", size = 0.6, inherit.aes = FALSE) +
  geom_text(data = my_base_data, aes(x = title, y = -18, label = group, hjust = hjust),
            colour = "black", size = 4, fontface = "bold", inherit.aes = FALSE)
# 显示图形
my_p

ggsave(my_p, file="/data/wuwq/noise/DISEASE_ALL/noise_var/pubmed_per20top500wilcox_v3/plot/plot_COVID_2.svg", width=10, height=10)


rm(list=ls())
################################################################################################################# 4.postCOVID
my_data <- read.delim("/data/wuwq/noise/DISEASE_ALL/noise_var/pubmed_per20top500wilcox_v3/plot/data_postCOVID.txt")

# 添加每组末尾空白
my_empty_bar <- 2
my_nObsType <- nlevels(as.factor(my_data$observation))
my_to_add <- data.frame(matrix(NA, my_empty_bar * nlevels(my_data$group) * my_nObsType, ncol = ncol(my_data)))
colnames(my_to_add) <- colnames(my_data)
my_to_add$group <- rep(levels(my_data$group), each = my_empty_bar * my_nObsType)
my_data <- rbind(my_data, my_to_add)
my_data <- my_data %>% arrange(group, Variable)
my_data$id <- rep(seq(1, nrow(my_data) / my_nObsType), each = my_nObsType)

# 计算每个标签的总和
my_label_data <- my_data %>% 
  group_by(id, Variable) %>% 
  summarize(tot = sum(value, na.rm = TRUE), .groups = "drop")

my_number_of_bar <- nrow(my_label_data)
my_angle <- 90 - 360 * (my_label_data$id - 0.5) / my_number_of_bar
my_label_data <- my_label_data %>%
  mutate(hjust = ifelse(my_angle < -90, 1, 0),
         angle = ifelse(my_angle < -90, my_angle + 180, my_angle))

# 基线数据
my_base_data <- my_data %>% 
  group_by(group) %>% 
  summarize(start = min(id), end = max(id) - my_empty_bar, .groups = "drop") %>% 
  rowwise() %>% 
  mutate(title = mean(c(start, end)))

# grid 数据
my_grid_data <- my_base_data
my_grid_data$end <- my_grid_data$end[c(nrow(my_grid_data), 1:(nrow(my_grid_data) - 1))] + 1
my_grid_data$start <- my_grid_data$start - 1
my_grid_data <- my_grid_data[-1,]

# 计算每组基线标签的水平对齐
my_base_data <- my_base_data %>%
  mutate(hjust = (start + end) / 2 / max(my_data$id))

# -----------------------------
# 设置 observation 的堆叠顺序和颜色
# 内圈从 CD4 到外圈 otherT
my_data$observation <- factor(my_data$observation,
                              levels = c("CD4", "CD8", "MNP", "NK", "B", "otherT"))

obs_colors <- c(
  otherT    = "#FDE725FF",
  B    = "#7AD151FF",
  NK    = "#22A884FF",
  MNP     = "#2A788EFF",
  CD8      = "#414487FF",
  CD4 = "#440154FF"
)


# -----------------------------
# 绘图
my_p <- ggplot(my_data) +
  geom_bar(aes(x = as.factor(id), y = value*100, fill = observation), 
           stat = "identity", alpha = 0.5) +
  scale_fill_manual(values = obs_colors) +  # 使用自定义颜色
  # grid lines
  geom_segment(data = my_grid_data, aes(x = end, y = 0, xend = start, yend = 0), colour = "grey", size = 0.3, inherit.aes = FALSE) +
  geom_segment(data = my_grid_data, aes(x = end, y = 90, xend = start, yend = 90), colour = "grey", size = 0.3, inherit.aes = FALSE) +
  geom_segment(data = my_grid_data, aes(x = end, y = 180, xend = start, yend = 180), colour = "grey", size = 0.3, inherit.aes = FALSE) +
  geom_segment(data = my_grid_data, aes(x = end, y = 270, xend = start, yend = 270), colour = "grey", size = 0.3, inherit.aes = FALSE) +
  geom_segment(data = my_grid_data, aes(x = end, y = 360, xend = start, yend = 360), colour = "grey", size = 0.3, inherit.aes = FALSE) +
  # y轴文本
  annotate("text", x = rep(max(my_data$id), 5), y = c(0, 90, 180, 270, 360), 
           label = c("0", "90", "180", "270", "360"), color = "grey", size = 8, fontface = "bold", hjust = 1) +
  # 设置y轴范围
  expand_limits(y = max(my_label_data$tot, na.rm = TRUE) + 50) +
  theme_minimal() +
  theme(
    legend.position = "right",  # 将图例放在右侧
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1, 4), "cm")
  ) +
  coord_polar() +
  # 基线信息
  geom_segment(data = my_base_data, aes(x = start, y = -5, xend = end, yend = -5),
               colour = "black", size = 0.6, inherit.aes = FALSE) +
  geom_text(data = my_base_data, aes(x = title, y = -18, label = group, hjust = hjust),
            colour = "black", size = 4, fontface = "bold", inherit.aes = FALSE)
# 显示图形
my_p

ggsave(my_p, file="/data/wuwq/noise/DISEASE_ALL/noise_var/pubmed_per20top500wilcox_v3/plot/plot_postCOVID_2.svg", width=10, height=10)


rm(list=ls())
################################################################################################################# 5.SEPSIS
my_data <- read.delim("/data/wuwq/noise/DISEASE_ALL/noise_var/pubmed_per20top500wilcox_v3/plot/data_SEPSIS.txt")

# 添加每组末尾空白
my_empty_bar <- 2
my_nObsType <- nlevels(as.factor(my_data$observation))
my_to_add <- data.frame(matrix(NA, my_empty_bar * nlevels(my_data$group) * my_nObsType, ncol = ncol(my_data)))
colnames(my_to_add) <- colnames(my_data)
my_to_add$group <- rep(levels(my_data$group), each = my_empty_bar * my_nObsType)
my_data <- rbind(my_data, my_to_add)
my_data <- my_data %>% arrange(group, Variable)
my_data$id <- rep(seq(1, nrow(my_data) / my_nObsType), each = my_nObsType)

# 计算每个标签的总和
my_label_data <- my_data %>% 
  group_by(id, Variable) %>% 
  summarize(tot = sum(value, na.rm = TRUE), .groups = "drop")

my_number_of_bar <- nrow(my_label_data)
my_angle <- 90 - 360 * (my_label_data$id - 0.5) / my_number_of_bar
my_label_data <- my_label_data %>%
  mutate(hjust = ifelse(my_angle < -90, 1, 0),
         angle = ifelse(my_angle < -90, my_angle + 180, my_angle))

# 基线数据
my_base_data <- my_data %>% 
  group_by(group) %>% 
  summarize(start = min(id), end = max(id) - my_empty_bar, .groups = "drop") %>% 
  rowwise() %>% 
  mutate(title = mean(c(start, end)))

# grid 数据
my_grid_data <- my_base_data
my_grid_data$end <- my_grid_data$end[c(nrow(my_grid_data), 1:(nrow(my_grid_data) - 1))] + 1
my_grid_data$start <- my_grid_data$start - 1
my_grid_data <- my_grid_data[-1,]

# 计算每组基线标签的水平对齐
my_base_data <- my_base_data %>%
  mutate(hjust = (start + end) / 2 / max(my_data$id))

# -----------------------------
# 设置 observation 的堆叠顺序和颜色
# 内圈从 CD4 到外圈 otherT
my_data$observation <- factor(my_data$observation,
                              levels = c("CD4", "CD8", "MNP", "NK", "B", "otherT"))

obs_colors <- c(
  otherT    = "#FDE725FF",
  B    = "#7AD151FF",
  NK    = "#22A884FF",
  MNP     = "#2A788EFF",
  CD8      = "#414487FF",
  CD4 = "#440154FF"
)

# -----------------------------
# 绘图
my_p <- ggplot(my_data) +
  geom_bar(aes(x = as.factor(id), y = value*100, fill = observation), 
           stat = "identity", alpha = 0.5) +
  scale_fill_manual(values = obs_colors) +  # 使用自定义颜色
  # grid lines
  geom_segment(data = my_grid_data, aes(x = end, y = 0, xend = start, yend = 0), colour = "grey", size = 0.3, inherit.aes = FALSE) +
  geom_segment(data = my_grid_data, aes(x = end, y = 85, xend = start, yend = 85), colour = "grey", size = 0.3, inherit.aes = FALSE) +
  geom_segment(data = my_grid_data, aes(x = end, y = 170, xend = start, yend = 170), colour = "grey", size = 0.3, inherit.aes = FALSE) +
  geom_segment(data = my_grid_data, aes(x = end, y = 255, xend = start, yend = 255), colour = "grey", size = 0.3, inherit.aes = FALSE) +
  geom_segment(data = my_grid_data, aes(x = end, y = 340, xend = start, yend = 340), colour = "grey", size = 0.3, inherit.aes = FALSE) +
  # y轴文本
  annotate("text", x = rep(max(my_data$id), 5), y = c(0, 85, 170, 255, 340), 
           label = c("0", "85", "170", "255", "340"), color = "grey", size = 8, fontface = "bold", hjust = 1) +
  # 设置y轴范围
  expand_limits(y = max(my_label_data$tot, na.rm = TRUE) + 50) +
  theme_minimal() +
  theme(
    legend.position = "right",  # 将图例放在右侧
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1, 4), "cm")
  ) +
  coord_polar() +
  # 基线信息
  geom_segment(data = my_base_data, aes(x = start, y = -5, xend = end, yend = -5),
               colour = "black", size = 0.6, inherit.aes = FALSE) +
  geom_text(data = my_base_data, aes(x = title, y = -18, label = group, hjust = hjust),
            colour = "black", size = 4, fontface = "bold", inherit.aes = FALSE)
# 显示图形
my_p

ggsave(my_p, file="/data/wuwq/noise/DISEASE_ALL/noise_var/pubmed_per20top500wilcox_v3/plot/plot_SEPSIS_2.svg", width=10, height=10)


rm(list=ls())
################################################################################################################# 6.FLU
my_data <- read.delim("/data/wuwq/noise/DISEASE_ALL/noise_var/pubmed_per20top500wilcox_v3/plot/data_FLU.txt")

# 添加每组末尾空白
my_empty_bar <- 2
my_nObsType <- nlevels(as.factor(my_data$observation))
my_to_add <- data.frame(matrix(NA, my_empty_bar * nlevels(my_data$group) * my_nObsType, ncol = ncol(my_data)))
colnames(my_to_add) <- colnames(my_data)
my_to_add$group <- rep(levels(my_data$group), each = my_empty_bar * my_nObsType)
my_data <- rbind(my_data, my_to_add)
my_data <- my_data %>% arrange(group, Variable)
my_data$id <- rep(seq(1, nrow(my_data) / my_nObsType), each = my_nObsType)

# 计算每个标签的总和
my_label_data <- my_data %>% 
  group_by(id, Variable) %>% 
  summarize(tot = sum(value, na.rm = TRUE), .groups = "drop")

my_number_of_bar <- nrow(my_label_data)
my_angle <- 90 - 360 * (my_label_data$id - 0.5) / my_number_of_bar
my_label_data <- my_label_data %>%
  mutate(hjust = ifelse(my_angle < -90, 1, 0),
         angle = ifelse(my_angle < -90, my_angle + 180, my_angle))

# 基线数据
my_base_data <- my_data %>% 
  group_by(group) %>% 
  summarize(start = min(id), end = max(id) - my_empty_bar, .groups = "drop") %>% 
  rowwise() %>% 
  mutate(title = mean(c(start, end)))

# grid 数据
my_grid_data <- my_base_data
my_grid_data$end <- my_grid_data$end[c(nrow(my_grid_data), 1:(nrow(my_grid_data) - 1))] + 1
my_grid_data$start <- my_grid_data$start - 1
my_grid_data <- my_grid_data[-1,]

# 计算每组基线标签的水平对齐
my_base_data <- my_base_data %>%
  mutate(hjust = (start + end) / 2 / max(my_data$id))

# -----------------------------
# 设置 observation 的堆叠顺序和颜色
# 内圈从 CD4 到外圈 otherT
my_data$observation <- factor(my_data$observation,
                              levels = c("CD4", "CD8", "MNP", "NK", "B", "otherT"))

obs_colors <- c(
  otherT    = "#FDE725FF",
  B    = "#7AD151FF",
  NK    = "#22A884FF",
  MNP     = "#2A788EFF",
  CD8      = "#414487FF",
  CD4 = "#440154FF"
)

# -----------------------------
# 绘图
my_p <- ggplot(my_data) +
  geom_bar(aes(x = as.factor(id), y = value*100, fill = observation), 
           stat = "identity", alpha = 0.5) +
  scale_fill_manual(values = obs_colors) +  # 使用自定义颜色
  # grid lines
  geom_segment(data = my_grid_data, aes(x = end, y = 0, xend = start, yend = 0), colour = "grey", size = 0.3, inherit.aes = FALSE) +
  geom_segment(data = my_grid_data, aes(x = end, y = 75, xend = start, yend = 75), colour = "grey", size = 0.3, inherit.aes = FALSE) +
  geom_segment(data = my_grid_data, aes(x = end, y = 150, xend = start, yend = 150), colour = "grey", size = 0.3, inherit.aes = FALSE) +
  geom_segment(data = my_grid_data, aes(x = end, y = 225, xend = start, yend = 225), colour = "grey", size = 0.3, inherit.aes = FALSE) +
  geom_segment(data = my_grid_data, aes(x = end, y = 300, xend = start, yend = 300), colour = "grey", size = 0.3, inherit.aes = FALSE) +
  # y轴文本
  annotate("text", x = rep(max(my_data$id), 5), y = c(0, 75, 150, 225, 300), 
           label = c("0", "75", "150", "225", "300"), color = "grey", size = 8, fontface = "bold", hjust = 1) +
  # 设置y轴范围
  expand_limits(y = max(my_label_data$tot, na.rm = TRUE) + 50) +
  theme_minimal() +
  theme(
    legend.position = "right",  # 将图例放在右侧
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1, 4), "cm")
  ) +
  coord_polar() +
  # 基线信息
  geom_segment(data = my_base_data, aes(x = start, y = -5, xend = end, yend = -5),
               colour = "black", size = 0.6, inherit.aes = FALSE) +
  geom_text(data = my_base_data, aes(x = title, y = -18, label = group, hjust = hjust),
            colour = "black", size = 4, fontface = "bold", inherit.aes = FALSE)
# 显示图形
my_p

ggsave(my_p, file="/data/wuwq/noise/DISEASE_ALL/noise_var/pubmed_per20top500wilcox_v3/plot/plot_FLU_2.svg", width=10, height=10)


rm(list=ls())
################################################################################################################# 7.NSCLC
my_data <- read.delim("/data/wuwq/noise/DISEASE_ALL/noise_var/pubmed_per20top500wilcox_v3/plot/data_NSCLC.txt")

# 添加每组末尾空白
my_empty_bar <- 2
my_nObsType <- nlevels(as.factor(my_data$observation))
my_to_add <- data.frame(matrix(NA, my_empty_bar * nlevels(my_data$group) * my_nObsType, ncol = ncol(my_data)))
colnames(my_to_add) <- colnames(my_data)
my_to_add$group <- rep(levels(my_data$group), each = my_empty_bar * my_nObsType)
my_data <- rbind(my_data, my_to_add)
my_data <- my_data %>% arrange(group, Variable)
my_data$id <- rep(seq(1, nrow(my_data) / my_nObsType), each = my_nObsType)

# 计算每个标签的总和
my_label_data <- my_data %>% 
  group_by(id, Variable) %>% 
  summarize(tot = sum(value, na.rm = TRUE), .groups = "drop")

my_number_of_bar <- nrow(my_label_data)
my_angle <- 90 - 360 * (my_label_data$id - 0.5) / my_number_of_bar
my_label_data <- my_label_data %>%
  mutate(hjust = ifelse(my_angle < -90, 1, 0),
         angle = ifelse(my_angle < -90, my_angle + 180, my_angle))

# 基线数据
my_base_data <- my_data %>% 
  group_by(group) %>% 
  summarize(start = min(id), end = max(id) - my_empty_bar, .groups = "drop") %>% 
  rowwise() %>% 
  mutate(title = mean(c(start, end)))

# grid 数据
my_grid_data <- my_base_data
my_grid_data$end <- my_grid_data$end[c(nrow(my_grid_data), 1:(nrow(my_grid_data) - 1))] + 1
my_grid_data$start <- my_grid_data$start - 1
my_grid_data <- my_grid_data[-1,]

# 计算每组基线标签的水平对齐
my_base_data <- my_base_data %>%
  mutate(hjust = (start + end) / 2 / max(my_data$id))

# -----------------------------
# 设置 observation 的堆叠顺序和颜色
# 内圈从 CD4 到外圈 otherT
my_data$observation <- factor(my_data$observation,
                              levels = c("CD4", "CD8", "MNP", "NK", "B", "otherT"))

obs_colors <- c(
  otherT    = "#FDE725FF",
  B    = "#7AD151FF",
  NK    = "#22A884FF",
  MNP     = "#2A788EFF",
  CD8      = "#414487FF",
  CD4 = "#440154FF"
)


# -----------------------------
# 绘图
my_p <- ggplot(my_data) +
  geom_bar(aes(x = as.factor(id), y = value*100, fill = observation), 
           stat = "identity", alpha = 0.5) +
  scale_fill_manual(values = obs_colors) +  # 使用自定义颜色
  # grid lines
  geom_segment(data = my_grid_data, aes(x = end, y = 0, xend = start, yend = 0), colour = "grey", size = 0.3, inherit.aes = FALSE) +
  geom_segment(data = my_grid_data, aes(x = end, y = 190, xend = start, yend = 190), colour = "grey", size = 0.3, inherit.aes = FALSE) +
  geom_segment(data = my_grid_data, aes(x = end, y = 380, xend = start, yend = 380), colour = "grey", size = 0.3, inherit.aes = FALSE) +
  geom_segment(data = my_grid_data, aes(x = end, y = 570, xend = start, yend = 570), colour = "grey", size = 0.3, inherit.aes = FALSE) +
  geom_segment(data = my_grid_data, aes(x = end, y = 760, xend = start, yend = 760), colour = "grey", size = 0.3, inherit.aes = FALSE) +
  # y轴文本
  annotate("text", x = rep(max(my_data$id), 5), y = c(0, 190, 380, 570, 760), 
           label = c("0", "190", "380", "570", "760"), color = "grey", size = 8, fontface = "bold", hjust = 1) +
  # 设置y轴范围
  expand_limits(y = max(my_label_data$tot, na.rm = TRUE) + 50) +
  theme_minimal() +
  theme(
    legend.position = "right",  # 将图例放在右侧
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1, 4), "cm")
  ) +
  coord_polar() +
  # 基线信息
  geom_segment(data = my_base_data, aes(x = start, y = -5, xend = end, yend = -5),
               colour = "black", size = 0.6, inherit.aes = FALSE) +
  geom_text(data = my_base_data, aes(x = title, y = -18, label = group, hjust = hjust),
            colour = "black", size = 4, fontface = "bold", inherit.aes = FALSE)
# 显示图形
my_p

ggsave(my_p, file="/data/wuwq/noise/DISEASE_ALL/noise_var/pubmed_per20top500wilcox_v3/plot/plot_NSCLC_2.svg", width=10, height=10)



##################################################################
######## ----- (六) 疾病特异性noise_gene鉴定及分析 ----- #########  
##################################################################
##### "health" "RA" "SLE" "COVID-19" "post-COVID-19 disorder" "SEPSIS" "FLU" "NSCLC"

rm(list=ls())
### 读取分箱表 每个组均进行noi_percentage>=0.05的筛选；之后鉴定疾病有 & 健康无的基因表
# CD4
CD4_health <- read.csv("/data/wuwq/noise/DISEASE_ALL/noise_var/health_CD4.csv", header = TRUE, row.names = 1)
CD4_health <- subset(CD4_health, noise_percentage >= 0.05)
CD4_RA <- read.csv("/data/wuwq/noise/DISEASE_ALL/noise_var/RA_CD4.csv", header = TRUE, row.names = 1)
CD4_RA <- subset(CD4_RA, noise_percentage >= 0.05)
CD4_SLE <- read.csv("/data/wuwq/noise/DISEASE_ALL/noise_var/SLE_CD4.csv", header = TRUE, row.names = 1)
CD4_SLE <- subset(CD4_SLE, noise_percentage >= 0.05)
CD4_COVID19 <- read.csv("/data/wuwq/noise/DISEASE_ALL/noise_var/COVID-19_CD4.csv", header = TRUE, row.names = 1)
CD4_COVID19 <- subset(CD4_COVID19, noise_percentage >= 0.05)
CD4_postCOVID19 <- read.csv("/data/wuwq/noise/DISEASE_ALL/noise_var/post-COVID-19 disorder_CD4.csv", header = TRUE, row.names = 1)
CD4_postCOVID19 <- subset(CD4_postCOVID19, noise_percentage >= 0.05)
CD4_SEPSIS <- read.csv("/data/wuwq/noise/DISEASE_ALL/noise_var/SEPSIS_CD4.csv", header = TRUE, row.names = 1)
CD4_SEPSIS <- subset(CD4_SEPSIS, noise_percentage >= 0.05)
CD4_FLU <- read.csv("/data/wuwq/noise/DISEASE_ALL/noise_var/FLU_CD4.csv", header = TRUE, row.names = 1)
CD4_FLU <- subset(CD4_FLU, noise_percentage >= 0.05)
CD4_NSCLC <- read.csv("/data/wuwq/noise/DISEASE_ALL/noise_var/NSCLC_CD4.csv", header = TRUE, row.names = 1)
CD4_NSCLC <- subset(CD4_NSCLC, noise_percentage >= 0.05)
gene_CD4_health <- CD4_health$gene
gene_CD4_RA <- CD4_RA$gene
gene_CD4_SLE <- CD4_SLE$gene
gene_CD4_COVID19 <- CD4_COVID19$gene
gene_CD4_postCOVID19 <- CD4_postCOVID19$gene
gene_CD4_SEPSIS <- CD4_SEPSIS$gene
gene_CD4_FLU <- CD4_FLU$gene
gene_CD4_NSCLC <- CD4_NSCLC$gene
specific_CD4_RA <- setdiff(gene_CD4_RA, gene_CD4_health)
specific_CD4_SLE <- setdiff(gene_CD4_SLE, gene_CD4_health)
specific_CD4_COVID19 <- setdiff(gene_CD4_COVID19, gene_CD4_health)
specific_CD4_postCOVID19 <- setdiff(gene_CD4_postCOVID19, gene_CD4_health)
specific_CD4_SEPSIS <- setdiff(gene_CD4_SEPSIS, gene_CD4_health)
specific_CD4_FLU <- setdiff(gene_CD4_FLU, gene_CD4_health)
specific_CD4_NSCLC <- setdiff(gene_CD4_NSCLC, gene_CD4_health)
# CD8
CD8_health <- read.csv("/data/wuwq/noise/DISEASE_ALL/noise_var/health_CD8.csv", header = TRUE, row.names = 1)
CD8_health <- subset(CD8_health, noise_percentage >= 0.05)
CD8_RA <- read.csv("/data/wuwq/noise/DISEASE_ALL/noise_var/RA_CD8.csv", header = TRUE, row.names = 1)
CD8_RA <- subset(CD8_RA, noise_percentage >= 0.05)
CD8_SLE <- read.csv("/data/wuwq/noise/DISEASE_ALL/noise_var/SLE_CD8.csv", header = TRUE, row.names = 1)
CD8_SLE <- subset(CD8_SLE, noise_percentage >= 0.05)
CD8_COVID19 <- read.csv("/data/wuwq/noise/DISEASE_ALL/noise_var/COVID-19_CD8.csv", header = TRUE, row.names = 1)
CD8_COVID19 <- subset(CD8_COVID19, noise_percentage >= 0.05)
CD8_postCOVID19 <- read.csv("/data/wuwq/noise/DISEASE_ALL/noise_var/post-COVID-19 disorder_CD8.csv", header = TRUE, row.names = 1)
CD8_postCOVID19 <- subset(CD8_postCOVID19, noise_percentage >= 0.05)
CD8_SEPSIS <- read.csv("/data/wuwq/noise/DISEASE_ALL/noise_var/SEPSIS_CD8.csv", header = TRUE, row.names = 1)
CD8_SEPSIS <- subset(CD8_SEPSIS, noise_percentage >= 0.05)
CD8_FLU <- read.csv("/data/wuwq/noise/DISEASE_ALL/noise_var/FLU_CD8.csv", header = TRUE, row.names = 1)
CD8_FLU <- subset(CD8_FLU, noise_percentage >= 0.05)
CD8_NSCLC <- read.csv("/data/wuwq/noise/DISEASE_ALL/noise_var/NSCLC_CD8.csv", header = TRUE, row.names = 1)
CD8_NSCLC <- subset(CD8_NSCLC, noise_percentage >= 0.05)
gene_CD8_health <- CD8_health$gene
gene_CD8_RA <- CD8_RA$gene
gene_CD8_SLE <- CD8_SLE$gene
gene_CD8_COVID19 <- CD8_COVID19$gene
gene_CD8_postCOVID19 <- CD8_postCOVID19$gene
gene_CD8_SEPSIS <- CD8_SEPSIS$gene
gene_CD8_FLU <- CD8_FLU$gene
gene_CD8_NSCLC <- CD8_NSCLC$gene
specific_CD8_RA <- setdiff(gene_CD8_RA, gene_CD8_health)
specific_CD8_SLE <- setdiff(gene_CD8_SLE, gene_CD8_health)
specific_CD8_COVID19 <- setdiff(gene_CD8_COVID19, gene_CD8_health)
specific_CD8_postCOVID19 <- setdiff(gene_CD8_postCOVID19, gene_CD8_health)
specific_CD8_SEPSIS <- setdiff(gene_CD8_SEPSIS, gene_CD8_health)
specific_CD8_FLU <- setdiff(gene_CD8_FLU, gene_CD8_health)
specific_CD8_NSCLC <- setdiff(gene_CD8_NSCLC, gene_CD8_health)
# MNP
MNP_health <- read.csv("/data/wuwq/noise/DISEASE_ALL/noise_var/health_MNP.csv", header = TRUE, row.names = 1)
MNP_health <- subset(MNP_health, noise_percentage >= 0.05)
MNP_RA <- read.csv("/data/wuwq/noise/DISEASE_ALL/noise_var/RA_MNP.csv", header = TRUE, row.names = 1)
MNP_RA <- subset(MNP_RA, noise_percentage >= 0.05)
MNP_SLE <- read.csv("/data/wuwq/noise/DISEASE_ALL/noise_var/SLE_MNP.csv", header = TRUE, row.names = 1)
MNP_SLE <- subset(MNP_SLE, noise_percentage >= 0.05)
MNP_COVID19 <- read.csv("/data/wuwq/noise/DISEASE_ALL/noise_var/COVID-19_MNP.csv", header = TRUE, row.names = 1)
MNP_COVID19 <- subset(MNP_COVID19, noise_percentage >= 0.05)
MNP_postCOVID19 <- read.csv("/data/wuwq/noise/DISEASE_ALL/noise_var/post-COVID-19 disorder_MNP.csv", header = TRUE, row.names = 1)
MNP_postCOVID19 <- subset(MNP_postCOVID19, noise_percentage >= 0.05)
MNP_SEPSIS <- read.csv("/data/wuwq/noise/DISEASE_ALL/noise_var/SEPSIS_MNP.csv", header = TRUE, row.names = 1)
MNP_SEPSIS <- subset(MNP_SEPSIS, noise_percentage >= 0.05)
MNP_FLU <- read.csv("/data/wuwq/noise/DISEASE_ALL/noise_var/FLU_MNP.csv", header = TRUE, row.names = 1)
MNP_FLU <- subset(MNP_FLU, noise_percentage >= 0.05)
MNP_NSCLC <- read.csv("/data/wuwq/noise/DISEASE_ALL/noise_var/NSCLC_MNP.csv", header = TRUE, row.names = 1)
MNP_NSCLC <- subset(MNP_NSCLC, noise_percentage >= 0.05)
gene_MNP_health <- MNP_health$gene
gene_MNP_RA <- MNP_RA$gene
gene_MNP_SLE <- MNP_SLE$gene
gene_MNP_COVID19 <- MNP_COVID19$gene
gene_MNP_postCOVID19 <- MNP_postCOVID19$gene
gene_MNP_SEPSIS <- MNP_SEPSIS$gene
gene_MNP_FLU <- MNP_FLU$gene
gene_MNP_NSCLC <- MNP_NSCLC$gene
specific_MNP_RA <- setdiff(gene_MNP_RA, gene_MNP_health)
specific_MNP_SLE <- setdiff(gene_MNP_SLE, gene_MNP_health)
specific_MNP_COVID19 <- setdiff(gene_MNP_COVID19, gene_MNP_health)
specific_MNP_postCOVID19 <- setdiff(gene_MNP_postCOVID19, gene_MNP_health)
specific_MNP_SEPSIS <- setdiff(gene_MNP_SEPSIS, gene_MNP_health)
specific_MNP_FLU <- setdiff(gene_MNP_FLU, gene_MNP_health)
specific_MNP_NSCLC <- setdiff(gene_MNP_NSCLC, gene_MNP_health)
# NK
NK_health <- read.csv("/data/wuwq/noise/DISEASE_ALL/noise_var/health_NK.csv", header = TRUE, row.names = 1)
NK_health <- subset(NK_health, noise_percentage >= 0.05)
NK_RA <- read.csv("/data/wuwq/noise/DISEASE_ALL/noise_var/RA_NK.csv", header = TRUE, row.names = 1)
NK_RA <- subset(NK_RA, noise_percentage >= 0.05)
NK_SLE <- read.csv("/data/wuwq/noise/DISEASE_ALL/noise_var/SLE_NK.csv", header = TRUE, row.names = 1)
NK_SLE <- subset(NK_SLE, noise_percentage >= 0.05)
NK_COVID19 <- read.csv("/data/wuwq/noise/DISEASE_ALL/noise_var/COVID-19_NK.csv", header = TRUE, row.names = 1)
NK_COVID19 <- subset(NK_COVID19, noise_percentage >= 0.05)
NK_postCOVID19 <- read.csv("/data/wuwq/noise/DISEASE_ALL/noise_var/post-COVID-19 disorder_NK.csv", header = TRUE, row.names = 1)
NK_postCOVID19 <- subset(NK_postCOVID19, noise_percentage >= 0.05)
NK_SEPSIS <- read.csv("/data/wuwq/noise/DISEASE_ALL/noise_var/SEPSIS_NK.csv", header = TRUE, row.names = 1)
NK_SEPSIS <- subset(NK_SEPSIS, noise_percentage >= 0.05)
NK_FLU <- read.csv("/data/wuwq/noise/DISEASE_ALL/noise_var/FLU_NK.csv", header = TRUE, row.names = 1)
NK_FLU <- subset(NK_FLU, noise_percentage >= 0.05)
NK_NSCLC <- read.csv("/data/wuwq/noise/DISEASE_ALL/noise_var/NSCLC_NK.csv", header = TRUE, row.names = 1)
NK_NSCLC <- subset(NK_NSCLC, noise_percentage >= 0.05)
gene_NK_health <- NK_health$gene
gene_NK_RA <- NK_RA$gene
gene_NK_SLE <- NK_SLE$gene
gene_NK_COVID19 <- NK_COVID19$gene
gene_NK_postCOVID19 <- NK_postCOVID19$gene
gene_NK_SEPSIS <- NK_SEPSIS$gene
gene_NK_FLU <- NK_FLU$gene
gene_NK_NSCLC <- NK_NSCLC$gene
specific_NK_RA <- setdiff(gene_NK_RA, gene_NK_health)
specific_NK_SLE <- setdiff(gene_NK_SLE, gene_NK_health)
specific_NK_COVID19 <- setdiff(gene_NK_COVID19, gene_NK_health)
specific_NK_postCOVID19 <- setdiff(gene_NK_postCOVID19, gene_NK_health)
specific_NK_SEPSIS <- setdiff(gene_NK_SEPSIS, gene_NK_health)
specific_NK_FLU <- setdiff(gene_NK_FLU, gene_NK_health)
specific_NK_NSCLC <- setdiff(gene_NK_NSCLC, gene_NK_health)
# B
B_health <- read.csv("/data/wuwq/noise/DISEASE_ALL/noise_var/health_B.csv", header = TRUE, row.names = 1)
B_health <- subset(B_health, noise_percentage >= 0.05)
B_RA <- read.csv("/data/wuwq/noise/DISEASE_ALL/noise_var/RA_B.csv", header = TRUE, row.names = 1)
B_RA <- subset(B_RA, noise_percentage >= 0.05)
B_SLE <- read.csv("/data/wuwq/noise/DISEASE_ALL/noise_var/SLE_B.csv", header = TRUE, row.names = 1)
B_SLE <- subset(B_SLE, noise_percentage >= 0.05)
B_COVID19 <- read.csv("/data/wuwq/noise/DISEASE_ALL/noise_var/COVID-19_B.csv", header = TRUE, row.names = 1)
B_COVID19 <- subset(B_COVID19, noise_percentage >= 0.05)
B_postCOVID19 <- read.csv("/data/wuwq/noise/DISEASE_ALL/noise_var/post-COVID-19 disorder_B.csv", header = TRUE, row.names = 1)
B_postCOVID19 <- subset(B_postCOVID19, noise_percentage >= 0.05)
B_SEPSIS <- read.csv("/data/wuwq/noise/DISEASE_ALL/noise_var/SEPSIS_B.csv", header = TRUE, row.names = 1)
B_SEPSIS <- subset(B_SEPSIS, noise_percentage >= 0.05)
B_FLU <- read.csv("/data/wuwq/noise/DISEASE_ALL/noise_var/FLU_B.csv", header = TRUE, row.names = 1)
B_FLU <- subset(B_FLU, noise_percentage >= 0.05)
B_NSCLC <- read.csv("/data/wuwq/noise/DISEASE_ALL/noise_var/NSCLC_B.csv", header = TRUE, row.names = 1)
B_NSCLC <- subset(B_NSCLC, noise_percentage >= 0.05)
gene_B_health <- B_health$gene
gene_B_RA <- B_RA$gene
gene_B_SLE <- B_SLE$gene
gene_B_COVID19 <- B_COVID19$gene
gene_B_postCOVID19 <- B_postCOVID19$gene
gene_B_SEPSIS <- B_SEPSIS$gene
gene_B_FLU <- B_FLU$gene
gene_B_NSCLC <- B_NSCLC$gene
specific_B_RA <- setdiff(gene_B_RA, gene_B_health)
specific_B_SLE <- setdiff(gene_B_SLE, gene_B_health)
specific_B_COVID19 <- setdiff(gene_B_COVID19, gene_B_health)
specific_B_postCOVID19 <- setdiff(gene_B_postCOVID19, gene_B_health)
specific_B_SEPSIS <- setdiff(gene_B_SEPSIS, gene_B_health)
specific_B_FLU <- setdiff(gene_B_FLU, gene_B_health)
specific_B_NSCLC <- setdiff(gene_B_NSCLC, gene_B_health)
# otherT
otherT_health <- read.csv("/data/wuwq/noise/DISEASE_ALL/noise_var/health_otherT.csv", header = TRUE, row.names = 1)
otherT_health <- subset(otherT_health, noise_percentage >= 0.05)
otherT_RA <- read.csv("/data/wuwq/noise/DISEASE_ALL/noise_var/RA_otherT.csv", header = TRUE, row.names = 1)
otherT_RA <- subset(otherT_RA, noise_percentage >= 0.05)
otherT_COVID19 <- read.csv("/data/wuwq/noise/DISEASE_ALL/noise_var/COVID-19_otherT.csv", header = TRUE, row.names = 1)
otherT_COVID19 <- subset(otherT_COVID19, noise_percentage >= 0.05)
otherT_postCOVID19 <- read.csv("/data/wuwq/noise/DISEASE_ALL/noise_var/post-COVID-19 disorder_otherT.csv", header = TRUE, row.names = 1)
otherT_postCOVID19 <- subset(otherT_postCOVID19, noise_percentage >= 0.05)
otherT_SEPSIS <- read.csv("/data/wuwq/noise/DISEASE_ALL/noise_var/SEPSIS_otherT.csv", header = TRUE, row.names = 1)
otherT_SEPSIS <- subset(otherT_SEPSIS, noise_percentage >= 0.05)
otherT_FLU <- read.csv("/data/wuwq/noise/DISEASE_ALL/noise_var/FLU_otherT.csv", header = TRUE, row.names = 1)
otherT_FLU <- subset(otherT_FLU, noise_percentage >= 0.05)
otherT_NSCLC <- read.csv("/data/wuwq/noise/DISEASE_ALL/noise_var/NSCLC_otherT.csv", header = TRUE, row.names = 1)
otherT_NSCLC <- subset(otherT_NSCLC, noise_percentage >= 0.05)
gene_otherT_health <- otherT_health$gene
gene_otherT_RA <- otherT_RA$gene
gene_otherT_COVID19 <- otherT_COVID19$gene
gene_otherT_postCOVID19 <- otherT_postCOVID19$gene
gene_otherT_SEPSIS <- otherT_SEPSIS$gene
gene_otherT_FLU <- otherT_FLU$gene
gene_otherT_NSCLC <- otherT_NSCLC$gene
specific_otherT_RA <- setdiff(gene_otherT_RA, gene_otherT_health)
specific_otherT_COVID19 <- setdiff(gene_otherT_COVID19, gene_otherT_health)
specific_otherT_postCOVID19 <- setdiff(gene_otherT_postCOVID19, gene_otherT_health)
specific_otherT_SEPSIS <- setdiff(gene_otherT_SEPSIS, gene_otherT_health)
specific_otherT_FLU <- setdiff(gene_otherT_FLU, gene_otherT_health)
specific_otherT_NSCLC <- setdiff(gene_otherT_NSCLC, gene_otherT_health)

##############################################################
# 清除环境中除了以"specific_"开头变量外的其他变量
all_vars <- ls()
# 过滤出所有以"specific_"开头的变量名
specific_vars <- grep("^specific_", all_vars, value = TRUE)
# 找到需要清除的变量名，即不包括以"specific_"开头的变量
vars_to_remove <- setdiff(all_vars, specific_vars)
# 删除这些不需要的变量
rm(list = vars_to_remove)
rm(all_vars)
rm(vars_to_remove)
rm(specific_vars)

##############################################################
##### （除SLE是>=3）在大于等于4种细胞类型中均被鉴定到的specific_noise_gene,作为disease-specific-noise_gene
### RA
gene_lists <- list(
  gene_CD4 = specific_CD4_RA,
  gene_CD8 = specific_CD8_RA,
  gene_MNP = specific_MNP_RA,
  gene_NK = specific_NK_RA,
  gene_B = specific_B_RA,
  gene_otherT = specific_otherT_RA
)
gene_counts <- table(unlist(gene_lists))
specific_RA_share_genes <- names(gene_counts[gene_counts >= 4])

### SLE
gene_lists <- list(
  gene_CD4 = specific_CD4_SLE,
  gene_CD8 = specific_CD8_SLE,
  gene_MNP = specific_MNP_SLE,
  gene_NK = specific_NK_SLE,
  gene_B = specific_B_SLE
)
gene_counts <- table(unlist(gene_lists))
specific_SLE_share_genes <- names(gene_counts[gene_counts >= 3])

### COVID19
gene_lists <- list(
  gene_CD4 = specific_CD4_COVID19,
  gene_CD8 = specific_CD8_COVID19,
  gene_MNP = specific_MNP_COVID19,
  gene_NK = specific_NK_COVID19,
  gene_B = specific_B_COVID19,
  gene_otherT = specific_otherT_COVID19
)
gene_counts <- table(unlist(gene_lists))
specific_COVID19_share_genes <- names(gene_counts[gene_counts >= 4])

### postCOVID19
gene_lists <- list(
  gene_CD4 = specific_CD4_postCOVID19,
  gene_CD8 = specific_CD8_postCOVID19,
  gene_MNP = specific_MNP_postCOVID19,
  gene_NK = specific_NK_postCOVID19,
  gene_B = specific_B_postCOVID19,
  gene_otherT = specific_otherT_postCOVID19
)
gene_counts <- table(unlist(gene_lists))
specific_postCOVID19_share_genes <- names(gene_counts[gene_counts >= 4])

### SEPSIS
gene_lists <- list(
  gene_CD4 = specific_CD4_SEPSIS,
  gene_CD8 = specific_CD8_SEPSIS,
  gene_MNP = specific_MNP_SEPSIS,
  gene_NK = specific_NK_SEPSIS,
  gene_B = specific_B_SEPSIS,
  gene_otherT = specific_otherT_SEPSIS
)
gene_counts <- table(unlist(gene_lists))
specific_SEPSIS_share_genes <- names(gene_counts[gene_counts >= 4])

### FLU
gene_lists <- list(
  gene_CD4 = specific_CD4_FLU,
  gene_CD8 = specific_CD8_FLU,
  gene_MNP = specific_MNP_FLU,
  gene_NK = specific_NK_FLU,
  gene_B = specific_B_FLU,
  gene_otherT = specific_otherT_FLU
)
gene_counts <- table(unlist(gene_lists))
specific_FLU_share_genes <- names(gene_counts[gene_counts >= 4])

### NSCLC
gene_lists <- list(
  gene_CD4 = specific_CD4_NSCLC,
  gene_CD8 = specific_CD8_NSCLC,
  gene_MNP = specific_MNP_NSCLC,
  gene_NK = specific_NK_NSCLC,
  gene_B = specific_B_NSCLC,
  gene_otherT = specific_otherT_NSCLC
)
gene_counts <- table(unlist(gene_lists))
specific_NSCLC_share_genes <- names(gene_counts[gene_counts >= 4])

##### 保存疾病特异性noise_gene
write.csv(specific_SLE_share_genes, file="/data/wuwq/noise/DISEASE_ALL/noise_var/SLE_specific_noi_genes.csv", row.names = TRUE)
write.csv(specific_COVID19_share_genes, file="/data/wuwq/noise/DISEASE_ALL/noise_var/COVID19_specific_noi_genes.csv", row.names = TRUE)
write.csv(specific_postCOVID19_share_genes, file="/data/wuwq/noise/DISEASE_ALL/noise_var/postCOVID19_specific_noi_genes.csv", row.names = TRUE)
write.csv(specific_SEPSIS_share_genes, file="/data/wuwq/noise/DISEASE_ALL/noise_var/SEPSIS_specific_noi_genes.csv", row.names = TRUE)
write.csv(specific_FLU_share_genes, file="/data/wuwq/noise/DISEASE_ALL/noise_var/FLU_specific_noi_genes.csv", row.names = TRUE)
write.csv(specific_NSCLC_share_genes, file="/data/wuwq/noise/DISEASE_ALL/noise_var/NSCLC_specific_noi_genes.csv", row.names = TRUE)



##############################################################
# 清除环境中除了以"specific_"开头变量外的其他变量
all_vars <- ls()
# 过滤出所有以"specific_"开头的变量名
specific_vars <- grep("^specific_", all_vars, value = TRUE)
# 找到需要清除的变量名，即不包括以"specific_"开头的变量
vars_to_remove <- setdiff(all_vars, specific_vars)
# 删除这些不需要的变量
rm(list = vars_to_remove)
rm(all_vars)
rm(vars_to_remove)
rm(specific_vars)

##############################################################
##### 可视化 在大于等于2种细胞类型中均被鉴定到的specific_noise_gene频率(disease-specific-noise_gene)
### RA
gene_lists <- list(
  gene_CD4 = specific_CD4_RA,
  gene_CD8 = specific_CD8_RA,
  gene_MNP = specific_MNP_RA,
  gene_NK = specific_NK_RA,
  gene_B = specific_B_RA,
  gene_otherT = specific_otherT_RA
)
gene_counts <- table(unlist(gene_lists))
# 将结果转换为数据框以供绘制使用
gene_counts_df <- as.data.frame(gene_counts)
colnames(gene_counts_df) <- c("word", "freq")

# 筛选频率>=1的基因
gene_counts_df <- subset(gene_counts_df,freq>=1)
# 使用ggplot和ggwordcloud绘制不规则气泡图
p <- ggplot(gene_counts_df, aes(label = word, size = freq, color = factor(freq))) +
  geom_text_wordcloud_area(eccentricity = 1, rm_outside = TRUE) +  
  scale_size_area(max_size = 20) +  
  scale_color_manual(values = c("1" = "grey","2" = "grey","3" = "grey", "4" = "#BD9AAD", "5" = "#9193B4", "6" = "#972D36")) + 
  coord_fixed() +  # 确保坐标比例固定为1:1，从而使画布呈正方形
  theme_void()  # 取消轴和网格线
CairoSVG(file = "/data/wuwq/noise/DISEASE_ALL/noise_var/wordcloud_plot/RA_genecloud.svg", width = 10, height = 10)
# 绘制图形到 SVG 设备
print(p)
# 关闭设备以完成文件保存
dev.off()


### SLE
gene_lists <- list(
  gene_CD4 = specific_CD4_SLE,
  gene_CD8 = specific_CD8_SLE,
  gene_MNP = specific_MNP_SLE,
  gene_NK = specific_NK_SLE,
  gene_B = specific_B_SLE
)
gene_counts <- table(unlist(gene_lists))
# 将结果转换为数据框以供绘制使用
gene_counts_df <- as.data.frame(gene_counts)
colnames(gene_counts_df) <- c("word", "freq")

# 筛选频率>=2的基因
gene_counts_df <- subset(gene_counts_df,freq>=2)
# 使用ggplot和ggwordcloud绘制不规则气泡图
p <- ggplot(gene_counts_df, aes(label = word, size = freq, color = factor(freq))) +
  geom_text_wordcloud_area(eccentricity = 1, rm_outside = TRUE) +  
  scale_size_area(max_size = 20) +  
  scale_color_manual(values = c("1" = "grey","2" = "grey", "3" = "#BD9AAD", "4" = "#9193B4", "5" = "#972D36")) + 
  coord_fixed() +  # 确保坐标比例固定为1:1，从而使画布呈正方形
  theme_void()  # 取消轴和网格线
CairoSVG(file = "/data/wuwq/noise/DISEASE_ALL/noise_var/wordcloud_plot/SLE_genecloud.svg", width = 10, height = 10)
# 绘制图形到 SVG 设备
print(p)
# 关闭设备以完成文件保存
dev.off()


### COVID19
gene_lists <- list(
  gene_CD4 = specific_CD4_COVID19,
  gene_CD8 = specific_CD8_COVID19,
  gene_MNP = specific_MNP_COVID19,
  gene_NK = specific_NK_COVID19,
  gene_B = specific_B_COVID19,
  gene_otherT = specific_otherT_COVID19
)
gene_counts <- table(unlist(gene_lists))
# 将结果转换为数据框以供绘制使用
gene_counts_df <- as.data.frame(gene_counts)
colnames(gene_counts_df) <- c("word", "freq")

# 筛选频率>=4的基因
gene_counts_df <- subset(gene_counts_df,freq>=4)
# 使用ggplot和ggwordcloud绘制不规则气泡图
p <- ggplot(gene_counts_df, aes(label = word, size = freq, color = factor(freq))) +
  geom_text_wordcloud_area(eccentricity = 1, rm_outside = TRUE) +  
  scale_size_area(max_size = 20) +  
  scale_color_manual(values = c("1" = "grey","2" = "grey","3" = "grey", "4" = "#BD9AAD", "5" = "#9193B4", "6" = "#972D36")) + 
  coord_fixed() +  # 确保坐标比例固定为1:1，从而使画布呈正方形
  theme_void()  # 取消轴和网格线
CairoSVG(file = "/data/wuwq/noise/DISEASE_ALL/noise_var/wordcloud_plot/COVID19_genecloud.svg", width = 10, height = 10)
# 绘制图形到 SVG 设备
print(p)
# 关闭设备以完成文件保存
dev.off()


### postCOVID19
gene_lists <- list(
  gene_CD4 = specific_CD4_postCOVID19,
  gene_CD8 = specific_CD8_postCOVID19,
  gene_MNP = specific_MNP_postCOVID19,
  gene_NK = specific_NK_postCOVID19,
  gene_B = specific_B_postCOVID19,
  gene_otherT = specific_otherT_postCOVID19
)
gene_counts <- table(unlist(gene_lists))
# 将结果转换为数据框以供绘制使用
gene_counts_df <- as.data.frame(gene_counts)
colnames(gene_counts_df) <- c("word", "freq")

# 筛选频率>=4的基因
gene_counts_df <- subset(gene_counts_df,freq>=4)
# 使用ggplot和ggwordcloud绘制不规则气泡图
p <- ggplot(gene_counts_df, aes(label = word, size = freq, color = factor(freq))) +
  geom_text_wordcloud_area(eccentricity = 1, rm_outside = TRUE) +  
  scale_size_area(max_size = 20) +  
  scale_color_manual(values = c("1" = "grey","2" = "grey","3" = "grey", "4" = "#BD9AAD", "5" = "#9193B4", "6" = "#972D36")) + 
  coord_fixed() +  # 确保坐标比例固定为1:1，从而使画布呈正方形
  theme_void()  # 取消轴和网格线
CairoSVG(file = "/data/wuwq/noise/DISEASE_ALL/noise_var/wordcloud_plot/postCOVID19_genecloud.svg", width = 10, height = 10)
# 绘制图形到 SVG 设备
print(p)
# 关闭设备以完成文件保存
dev.off()


### SEPSIS
gene_lists <- list(
  gene_CD4 = specific_CD4_SEPSIS,
  gene_CD8 = specific_CD8_SEPSIS,
  gene_MNP = specific_MNP_SEPSIS,
  gene_NK = specific_NK_SEPSIS,
  gene_B = specific_B_SEPSIS,
  gene_otherT = specific_otherT_SEPSIS
)
gene_counts <- table(unlist(gene_lists))
# 将结果转换为数据框以供绘制使用
gene_counts_df <- as.data.frame(gene_counts)
colnames(gene_counts_df) <- c("word", "freq")

# 筛选频率>=4的基因
gene_counts_df <- subset(gene_counts_df,freq>=4)
# 使用ggplot和ggwordcloud绘制不规则气泡图
p <- ggplot(gene_counts_df, aes(label = word, size = freq, color = factor(freq))) +
  geom_text_wordcloud_area(eccentricity = 1, rm_outside = TRUE) +  
  scale_size_area(max_size = 20) +  
  scale_color_manual(values = c("1" = "grey","2" = "grey","3" = "grey", "4" = "#BD9AAD", "5" = "#9193B4", "6" = "#972D36")) + 
  coord_fixed() +  # 确保坐标比例固定为1:1，从而使画布呈正方形
  theme_void()  # 取消轴和网格线
CairoSVG(file = "/data/wuwq/noise/DISEASE_ALL/noise_var/wordcloud_plot/SEPSIS_genecloud.svg", width = 10, height = 10)
# 绘制图形到 SVG 设备
print(p)
# 关闭设备以完成文件保存
dev.off()


### FLU
gene_lists <- list(
  gene_CD4 = specific_CD4_FLU,
  gene_CD8 = specific_CD8_FLU,
  gene_MNP = specific_MNP_FLU,
  gene_NK = specific_NK_FLU,
  gene_B = specific_B_FLU,
  gene_otherT = specific_otherT_FLU
)
gene_counts <- table(unlist(gene_lists))
# 将结果转换为数据框以供绘制使用
gene_counts_df <- as.data.frame(gene_counts)
colnames(gene_counts_df) <- c("word", "freq")

# 筛选频率>=4的基因
gene_counts_df <- subset(gene_counts_df,freq>=4)
# 使用ggplot和ggwordcloud绘制不规则气泡图
p <- ggplot(gene_counts_df, aes(label = word, size = freq, color = factor(freq))) +
  geom_text_wordcloud_area(eccentricity = 1, rm_outside = TRUE) +  
  scale_size_area(max_size = 20) +  
  scale_color_manual(values = c("1" = "grey","2" = "grey","3" = "grey", "4" = "#BD9AAD", "5" = "#9193B4", "6" = "#972D36")) + 
  coord_fixed() +  # 确保坐标比例固定为1:1，从而使画布呈正方形
  theme_void()  # 取消轴和网格线
CairoSVG(file = "/data/wuwq/noise/DISEASE_ALL/noise_var/wordcloud_plot/FLU_genecloud.svg", width = 10, height = 10)
# 绘制图形到 SVG 设备
print(p)
# 关闭设备以完成文件保存
dev.off()


### NSCLC
gene_lists <- list(
  gene_CD4 = specific_CD4_NSCLC,
  gene_CD8 = specific_CD8_NSCLC,
  gene_MNP = specific_MNP_NSCLC,
  gene_NK = specific_NK_NSCLC,
  gene_B = specific_B_NSCLC,
  gene_otherT = specific_otherT_NSCLC
)
gene_counts <- table(unlist(gene_lists))
# 将结果转换为数据框以供绘制使用
gene_counts_df <- as.data.frame(gene_counts)
colnames(gene_counts_df) <- c("word", "freq")

# 筛选频率>=4的基因
gene_counts_df <- subset(gene_counts_df,freq>=4)
# 使用ggplot和ggwordcloud绘制不规则气泡图
p <- ggplot(gene_counts_df, aes(label = word, size = freq, color = factor(freq))) +
  geom_text_wordcloud_area(eccentricity = 1, rm_outside = TRUE) +  
  scale_size_area(max_size = 20) +  
  scale_color_manual(values = c("1" = "grey","2" = "grey","3" = "grey", "4" = "#BD9AAD", "5" = "#9193B4", "6" = "#972D36")) + 
  coord_fixed() +  # 确保坐标比例固定为1:1，从而使画布呈正方形
  theme_void()  # 取消轴和网格线
CairoSVG(file = "/data/wuwq/noise/DISEASE_ALL/noise_var/wordcloud_plot/NSCLC_genecloud.svg", width = 10, height = 10)
# 绘制图形到 SVG 设备
print(p)
# 关闭设备以完成文件保存
dev.off()






########################################################### 整合基因表达量结果
###   CD4_fit, CD8_fit, MNP_fit, NK_fit, DC_fit, B_fit, otherT_fit
library(dplyr)
library(knitr)
library(tidyr)
library(stats)
library(tibble)
##################################################################
################## ----- (一) 表达量提取------ ################### ✅️  ### 保存时有特殊符号的列名要保护
##################################################################


##### (一)、COMBAT
rm(list=ls())
# 设置文件夹路径
folder_path <- "/data/wuwq/noise/COMBAT/output_root_dir/CD4_fit"
# 获取所有子文件夹路径
subfolders <- list.files(folder_path, full.names = TRUE, include.dirs = TRUE)
# 创建一个空的列表来存储所有读取的数据框
data_list <- list()
for (subfolder in subfolders) {
  # 获取该子文件夹下所有 CSV 文件
  csv_files <- list.files(subfolder, pattern = "_gene_analysis\\.csv$", full.names = TRUE)
  # 过滤不包含"origin"的文件
  csv_files <- csv_files[!grepl("origin", basename(csv_files))]
  # 对每个 CSV 文件进行读取
  for (csv_file in csv_files) {
    # 读取 CSV 文件
    data <- read.csv(csv_file)
    # 提取文件名（去掉路径）
    file_base <- basename(csv_file)
    file_name <- gsub("_gene_analysis\\.csv$", "", file_base)
    # 将数据存储到列表中
    data_list[[file_name]] <- data
  }
}
# 合并所有个体的mean结果，得到总的汇总表(无mean赋值为0)
# 获取所有数据框中的 gene 列的基因 #并集
all_genes <- unique(unlist(lapply(data_list, function(df) df$g)))  # 需要用data_list中的基因并集 (因data_list_filter会使未入基因强制为0)
# 创建一个空的数据框，用于存储合并后的结果
mean_data <- data.frame(gene = all_genes)
# 循环遍历每个 data_list 中的元素，并合并 u 列
for (file_name in names(data_list)) {
  # 获取当前数据框
  df <- data_list[[file_name]]
  # 创建一个空的 u 列，默认值为 0
  u_column <- rep(0, length(all_genes))
  # 使用 match 向量化匹配，更高效
  match_idx <- match(all_genes, df$g)
  valid <- !is.na(match_idx)
  u_column[valid] <- df$u[match_idx[valid]]
  mean_data[[file_name]] <- u_column
}
##### 增加匹配的gene_symbol列
annotation_df <-read.delim("/data/wuwq/noise/mapping_geneid2symbol.txt",header = T)
# 使用 merge 左连接，保留 mean_data 中所有 gene
mean_data <- merge(mean_data, 
                   annotation_df, 
                   by.x = "gene",          # mean_data 中的 gene_id 列
                   by.y = "Gene_ID",       # 注释表中的 gene_id 列
                   all.x = TRUE)           # 保留所有 gene，即使没有匹配
# 把Gene_name列前移
col_order <- c("gene", "Gene_name", setdiff(names(mean_data), c("gene", "Gene_name")))
mean_data <- mean_data[, col_order]
write.table(mean_data,file = "/data/wuwq/noise/COMBAT/expression/CD4_u.csv", sep = ",", row.names = TRUE, col.names = NA, quote = TRUE)


##### (二)、FLU
rm(list=ls())
# 设置文件夹路径
folder_path <- "/data/wuwq/noise/FLU/output_root_dir/otherT_fit"
# 获取所有子文件夹路径
subfolders <- list.files(folder_path, full.names = TRUE, include.dirs = TRUE)
# 创建一个空的列表来存储所有读取的数据框
data_list <- list()
for (subfolder in subfolders) {
  # 获取该子文件夹下所有 CSV 文件
  csv_files <- list.files(subfolder, pattern = "_gene_analysis\\.csv$", full.names = TRUE)
  # 过滤不包含"origin"的文件
  csv_files <- csv_files[!grepl("origin", basename(csv_files))]
  # 对每个 CSV 文件进行读取
  for (csv_file in csv_files) {
    # 读取 CSV 文件
    data <- read.csv(csv_file)
    # 提取文件名（去掉路径）
    file_base <- basename(csv_file)
    file_name <- gsub("_gene_analysis\\.csv$", "", file_base)
    # 将数据存储到列表中
    data_list[[file_name]] <- data
  }
}
# 合并所有个体的mean结果，得到总的汇总表(无mean赋值为0)
# 获取所有数据框中的 gene 列的基因 #并集
all_genes <- unique(unlist(lapply(data_list, function(df) df$g)))  # 需要用data_list中的基因并集 (因data_list_filter会使未入基因强制为0)
# 创建一个空的数据框，用于存储合并后的结果
mean_data <- data.frame(gene = all_genes)
# 循环遍历每个 data_list 中的元素，并合并 u 列
for (file_name in names(data_list)) {
  # 获取当前数据框
  df <- data_list[[file_name]]
  # 创建一个空的 u 列，默认值为 0
  u_column <- rep(0, length(all_genes))
  # 使用 match 向量化匹配，更高效
  match_idx <- match(all_genes, df$g)
  valid <- !is.na(match_idx)
  u_column[valid] <- df$u[match_idx[valid]]
  mean_data[[file_name]] <- u_column
}
##### 增加匹配的gene_symbol列
annotation_df <-read.delim("/data/wuwq/noise/mapping_geneid2symbol.txt",header = T)
# 使用 merge 左连接，保留 mean_data 中所有 gene
mean_data <- merge(mean_data, 
                   annotation_df, 
                   by.x = "gene",          # mean_data 中的 gene_id 列
                   by.y = "Gene_ID",       # 注释表中的 gene_id 列
                   all.x = TRUE)           # 保留所有 gene，即使没有匹配
# 把Gene_name列前移
col_order <- c("gene", "Gene_name", setdiff(names(mean_data), c("gene", "Gene_name")))
mean_data <- mean_data[, col_order]
write.table(mean_data,file = "/data/wuwq/noise/FLU/expression/otherT_u.csv", sep = ",", row.names = TRUE, col.names = NA, quote = TRUE)



##### (三)、SLE
rm(list=ls())
# 设置文件夹路径
folder_path <- "/data/wuwq/noise/SLE/output_root_dir/otherT_fit"
# 获取所有子文件夹路径
subfolders <- list.files(folder_path, full.names = TRUE, include.dirs = TRUE)
# 创建一个空的列表来存储所有读取的数据框
data_list <- list()
for (subfolder in subfolders) {
  # 获取该子文件夹下所有 CSV 文件
  csv_files <- list.files(subfolder, pattern = "_gene_analysis\\.csv$", full.names = TRUE)
  # 过滤不包含"origin"的文件
  csv_files <- csv_files[!grepl("origin", basename(csv_files))]
  # 对每个 CSV 文件进行读取
  for (csv_file in csv_files) {
    # 读取 CSV 文件
    data <- read.csv(csv_file)
    # 提取文件名（去掉路径）
    file_base <- basename(csv_file)
    file_name <- gsub("_gene_analysis\\.csv$", "", file_base)
    # 将数据存储到列表中
    data_list[[file_name]] <- data
  }
}
# 合并所有个体的mean结果，得到总的汇总表(无mean赋值为0)
# 获取所有数据框中的 gene 列的基因 #并集
all_genes <- unique(unlist(lapply(data_list, function(df) df$g)))  # 需要用data_list中的基因并集 (因data_list_filter会使未入基因强制为0)
# 创建一个空的数据框，用于存储合并后的结果
mean_data <- data.frame(gene = all_genes)
# 循环遍历每个 data_list 中的元素，并合并 u 列
for (file_name in names(data_list)) {
  # 获取当前数据框
  df <- data_list[[file_name]]
  # 创建一个空的 u 列，默认值为 0
  u_column <- rep(0, length(all_genes))
  # 使用 match 向量化匹配，更高效
  match_idx <- match(all_genes, df$g)
  valid <- !is.na(match_idx)
  u_column[valid] <- df$u[match_idx[valid]]
  mean_data[[file_name]] <- u_column
}
##### 增加匹配的gene_symbol列
annotation_df <-read.delim("/data/wuwq/noise/mapping_geneid2symbol.txt",header = T)
# 使用 merge 左连接，保留 mean_data 中所有 gene
mean_data <- merge(mean_data, 
                   annotation_df, 
                   by.x = "gene",          # mean_data 中的 gene_id 列
                   by.y = "Gene_ID",       # 注释表中的 gene_id 列
                   all.x = TRUE)           # 保留所有 gene，即使没有匹配
# 把Gene_name列前移
col_order <- c("gene", "Gene_name", setdiff(names(mean_data), c("gene", "Gene_name")))
mean_data <- mean_data[, col_order]
write.table(mean_data,file = "/data/wuwq/noise/SLE/expression/otherT_u.csv", sep = ",", row.names = TRUE, col.names = NA, quote = TRUE)




##### (四)、RA
rm(list=ls())
# 设置文件夹路径
folder_path <- "/data/wuwq/noise/RA/output_root_dir/otherT_fit"
# 获取所有子文件夹路径
subfolders <- list.files(folder_path, full.names = TRUE, include.dirs = TRUE)
# 创建一个空的列表来存储所有读取的数据框
data_list <- list()
for (subfolder in subfolders) {
  # 获取该子文件夹下所有 CSV 文件
  csv_files <- list.files(subfolder, pattern = "_gene_analysis\\.csv$", full.names = TRUE)
  # 过滤不包含"origin"的文件
  csv_files <- csv_files[!grepl("origin", basename(csv_files))]
  # 对每个 CSV 文件进行读取
  for (csv_file in csv_files) {
    # 读取 CSV 文件
    data <- read.csv(csv_file)
    # 提取文件名（去掉路径）
    file_base <- basename(csv_file)
    file_name <- gsub("_gene_analysis\\.csv$", "", file_base)
    # 将数据存储到列表中
    data_list[[file_name]] <- data
  }
}
# 合并所有个体的mean结果，得到总的汇总表(无mean赋值为0)
# 获取所有数据框中的 gene 列的基因 #并集
all_genes <- unique(unlist(lapply(data_list, function(df) df$g)))  # 需要用data_list中的基因并集 (因data_list_filter会使未入基因强制为0)
# 创建一个空的数据框，用于存储合并后的结果
mean_data <- data.frame(gene = all_genes)
# 循环遍历每个 data_list 中的元素，并合并 u 列
for (file_name in names(data_list)) {
  # 获取当前数据框
  df <- data_list[[file_name]]
  # 创建一个空的 u 列，默认值为 0
  u_column <- rep(0, length(all_genes))
  # 使用 match 向量化匹配，更高效
  match_idx <- match(all_genes, df$g)
  valid <- !is.na(match_idx)
  u_column[valid] <- df$u[match_idx[valid]]
  mean_data[[file_name]] <- u_column
}
##### 增加匹配的gene_symbol列
annotation_df <-read.delim("/data/wuwq/noise/mapping_geneid2symbol.txt",header = T)
# 使用 merge 左连接，保留 mean_data 中所有 gene
mean_data <- merge(mean_data, 
                   annotation_df, 
                   by.x = "gene",          # mean_data 中的 gene_id 列
                   by.y = "Gene_ID",       # 注释表中的 gene_id 列
                   all.x = TRUE)           # 保留所有 gene，即使没有匹配
# 把Gene_name列前移
col_order <- c("gene", "Gene_name", setdiff(names(mean_data), c("gene", "Gene_name")))
mean_data <- mean_data[, col_order]
write.table(mean_data,file = "/data/wuwq/noise/RA/expression/otherT_u.csv", sep = ",", row.names = TRUE, col.names = NA, quote = TRUE)




##### (五)、longcovid
rm(list=ls())
# 设置文件夹路径
folder_path <- "/data/wuwq/noise/longcovid/output_root_dir/otherT_fit"
# 获取所有子文件夹路径
subfolders <- list.files(folder_path, full.names = TRUE, include.dirs = TRUE)
# 创建一个空的列表来存储所有读取的数据框
data_list <- list()
for (subfolder in subfolders) {
  # 获取该子文件夹下所有 CSV 文件
  csv_files <- list.files(subfolder, pattern = "_gene_analysis\\.csv$", full.names = TRUE)
  # 过滤不包含"origin"的文件
  csv_files <- csv_files[!grepl("origin", basename(csv_files))]
  # 对每个 CSV 文件进行读取
  for (csv_file in csv_files) {
    # 读取 CSV 文件
    data <- read.csv(csv_file)
    # 提取文件名（去掉路径）
    file_base <- basename(csv_file)
    file_name <- gsub("_gene_analysis\\.csv$", "", file_base)
    # 将数据存储到列表中
    data_list[[file_name]] <- data
  }
}
# 合并所有个体的mean结果，得到总的汇总表(无mean赋值为0)
# 获取所有数据框中的 gene 列的基因 #并集
all_genes <- unique(unlist(lapply(data_list, function(df) df$g)))  # 需要用data_list中的基因并集 (因data_list_filter会使未入基因强制为0)
# 创建一个空的数据框，用于存储合并后的结果
mean_data <- data.frame(gene = all_genes)
# 循环遍历每个 data_list 中的元素，并合并 u 列
for (file_name in names(data_list)) {
  # 获取当前数据框
  df <- data_list[[file_name]]
  # 创建一个空的 u 列，默认值为 0
  u_column <- rep(0, length(all_genes))
  # 使用 match 向量化匹配，更高效
  match_idx <- match(all_genes, df$g)
  valid <- !is.na(match_idx)
  u_column[valid] <- df$u[match_idx[valid]]
  mean_data[[file_name]] <- u_column
}
##### 增加匹配的gene_symbol列
annotation_df <-read.delim("/data/wuwq/noise/mapping_geneid2symbol.txt",header = T)
# 使用 merge 左连接，保留 mean_data 中所有 gene
mean_data <- merge(mean_data, 
                   annotation_df, 
                   by.x = "gene",          # mean_data 中的 gene_id 列
                   by.y = "Gene_ID",       # 注释表中的 gene_id 列
                   all.x = TRUE)           # 保留所有 gene，即使没有匹配
# 把Gene_name列前移
col_order <- c("gene", "Gene_name", setdiff(names(mean_data), c("gene", "Gene_name")))
mean_data <- mean_data[, col_order]
write.table(mean_data,file = "/data/wuwq/noise/longcovid/expression/otherT_u.csv", sep = ",", row.names = TRUE, col.names = NA, quote = TRUE)




##### (六)、NSCLC
rm(list=ls())
# 设置文件夹路径
folder_path <- "/data/wuwq/noise/NSCLC/output_root_dir/otherT_fit"
# 获取所有子文件夹路径
subfolders <- list.files(folder_path, full.names = TRUE, include.dirs = TRUE)
# 创建一个空的列表来存储所有读取的数据框
data_list <- list()
for (subfolder in subfolders) {
  # 获取该子文件夹下所有 CSV 文件
  csv_files <- list.files(subfolder, pattern = "_gene_analysis\\.csv$", full.names = TRUE)
  # 过滤不包含"origin"的文件
  csv_files <- csv_files[!grepl("origin", basename(csv_files))]
  # 对每个 CSV 文件进行读取
  for (csv_file in csv_files) {
    # 读取 CSV 文件
    data <- read.csv(csv_file)
    # 提取文件名（去掉路径）
    file_base <- basename(csv_file)
    file_name <- gsub("_gene_analysis\\.csv$", "", file_base)
    # 将数据存储到列表中
    data_list[[file_name]] <- data
  }
}
# 合并所有个体的mean结果，得到总的汇总表(无mean赋值为0)
# 获取所有数据框中的 gene 列的基因 #并集
all_genes <- unique(unlist(lapply(data_list, function(df) df$g)))  # 需要用data_list中的基因并集 (因data_list_filter会使未入基因强制为0)
# 创建一个空的数据框，用于存储合并后的结果
mean_data <- data.frame(gene = all_genes)
# 循环遍历每个 data_list 中的元素，并合并 u 列
for (file_name in names(data_list)) {
  # 获取当前数据框
  df <- data_list[[file_name]]
  # 创建一个空的 u 列，默认值为 0
  u_column <- rep(0, length(all_genes))
  # 使用 match 向量化匹配，更高效
  match_idx <- match(all_genes, df$g)
  valid <- !is.na(match_idx)
  u_column[valid] <- df$u[match_idx[valid]]
  mean_data[[file_name]] <- u_column
}
##### 增加匹配的gene_symbol列
annotation_df <-read.delim("/data/wuwq/noise/mapping_geneid2symbol.txt",header = T)
# 使用 merge 左连接，保留 mean_data 中所有 gene
mean_data <- merge(mean_data, 
                   annotation_df, 
                   by.x = "gene",          # mean_data 中的 gene_id 列
                   by.y = "Gene_ID",       # 注释表中的 gene_id 列
                   all.x = TRUE)           # 保留所有 gene，即使没有匹配
# 把Gene_name列前移
col_order <- c("gene", "Gene_name", setdiff(names(mean_data), c("gene", "Gene_name")))
mean_data <- mean_data[, col_order]
write.table(mean_data,file = "/data/wuwq/noise/NSCLC/expression/otherT_u.csv", sep = ",", row.names = TRUE, col.names = NA, quote = TRUE)






##################################################################
################## ----- (二) 与基线合并------ ################### ✅️ #保留了第一次出现的 ENSG00000168255   POLR2J3; ENSG00000285437   POLR2J3
##################################################################

rm(list=ls())
### 1.读取 对齐后大类基线信息、mean_u表  COMBAT/FLU/SLE/RA/longcovid 
baseline <- read.delim("/data/wuwq/noise/baseline/baseline_NSCLC_2.txt", header = TRUE, row.names = 1,sep = "\t")
row.names(baseline) <- baseline$sample_id
exp_u <- read.table("/data/wuwq/noise/NSCLC/otherT_u.csv", sep = ",", header = TRUE, check.names = FALSE, quote = "\"")
exp_u <- exp_u[ , -1]
exp_u <- exp_u %>%
  mutate(
    Gene_name = if_else(
      is.na(Gene_name) | Gene_name == "",
      gene,
      Gene_name,
      missing = gene   # 可选：明确处理缺失值
    )
  )
# 保留第一次出现的（通常是主要的那个） #ENSG00000168255   POLR2J3; ENSG00000285437   POLR2J3
exp_u <- exp_u[!duplicated(exp_u$Gene_name), ]
rownames(exp_u) <- exp_u$Gene_name
exp_u <- exp_u[ , -c(1:2), drop = FALSE]
exp_u <- as.data.frame(t(exp_u))
# 合并基线和noise表
common_rows <- intersect(rownames(baseline), rownames(exp_u))
exp_base <- data.frame(row.names = common_rows)
exp_base <- cbind(exp_base, baseline[common_rows, , drop = FALSE],  exp_u[common_rows, , drop = FALSE] )
write.table(exp_base,file = "/data/wuwq/noise/NSCLC/expression/otherT_base_u_2.csv", sep = ",", row.names = TRUE, col.names = NA, quote = TRUE)


rm(list=ls())
### 1.2.读取 对齐基线信息、mean_u表     NSCLC
baseline <- read.delim("/data/wuwq/noise/baseline/baseline_NSCLC_2.txt", header = TRUE, row.names = 1,sep = "\t")
row.names(baseline) <- baseline$sample_id
exp_u <- read.table("/data/wuwq/noise/NSCLC/otherT_u.csv", sep = ",", header = TRUE, check.names = FALSE, quote = "\"")
exp_u <- exp_u[ , -c(1,3)]
rownames(exp_u) <- exp_u$gene
exp_u <- exp_u[ , -1, drop = FALSE]
exp_u <- as.data.frame(t(exp_u))
# 合并基线和noise表
common_rows <- intersect(rownames(baseline), rownames(exp_u))
exp_base <- data.frame(row.names = common_rows)
exp_base <- cbind(exp_base, baseline[common_rows, , drop = FALSE],  exp_u[common_rows, , drop = FALSE] )
write.table(exp_base,file = "/data/wuwq/noise/NSCLC/expression/otherT_base_u_2.csv", sep = ",", row.names = TRUE, col.names = NA, quote = TRUE)



rm(list=ls())
### 2.1.读取 小基线信息、mean_u表     COMBAT / FLU / SLE / RA / longcovid
baseline <- read.delim("/data/wuwq/noise/baseline/baseline_COMBAT.txt", header = TRUE, row.names = 1,sep = "\t")
row.names(baseline) <- baseline$orig.ident.1
exp_u <- read.table("/data/wuwq/noise/COMBAT/CD4_u.csv", sep = ",", header = TRUE, check.names = FALSE, quote = "\"")
exp_u <- exp_u[ , -1]
exp_u <- exp_u %>%
  mutate(
    Gene_name = if_else(
      is.na(Gene_name) | Gene_name == "",
      gene,
      Gene_name,
      missing = gene   # 可选：明确处理缺失值
    )
  )
# 保留第一次出现的（通常是主要的那个） #ENSG00000168255   POLR2J3; ENSG00000285437   POLR2J3
exp_u <- exp_u[!duplicated(exp_u$Gene_name), ]
rownames(exp_u) <- exp_u$Gene_name
exp_u <- exp_u[ , -c(1:2), drop = FALSE]
exp_u <- as.data.frame(t(exp_u))
# 合并基线和noise表
common_rows <- intersect(rownames(baseline), rownames(exp_u))
exp_base <- data.frame(row.names = common_rows)
exp_base <- cbind(exp_base, baseline[common_rows, , drop = FALSE],  exp_u[common_rows, , drop = FALSE] )
write.table(exp_base,file = "/data/wuwq/noise/COMBAT/expression/CD4_base_u.csv", sep = ",", row.names = TRUE, col.names = NA, quote = TRUE)

rm(list=ls())
### 2.2.读取 小基线信息、mean_u表     NSCLC
baseline <- read.delim("/data/wuwq/noise/baseline/baseline_NSCLC.txt", header = TRUE, row.names = 1,sep = "\t")
row.names(baseline) <- baseline$orig.ident.1
exp_u <- read.table("/data/wuwq/noise/NSCLC/otherT_u.csv", sep = ",", header = TRUE, check.names = FALSE, quote = "\"")
exp_u <- exp_u[ , -c(1,3)]
rownames(exp_u) <- exp_u$gene
exp_u <- exp_u[ , -1, drop = FALSE]
exp_u <- as.data.frame(t(exp_u))
# 合并基线和noise表
common_rows <- intersect(rownames(baseline), rownames(exp_u))
exp_base <- data.frame(row.names = common_rows)
exp_base <- cbind(exp_base, baseline[common_rows, , drop = FALSE],  exp_u[common_rows, , drop = FALSE] )
write.table(exp_base,file = "/data/wuwq/noise/NSCLC/expression/otherT_base_u.csv", sep = ",", row.names = TRUE, col.names = NA, quote = TRUE)




#exp_base <- read.table("/data/wuwq/noise/COMBAT/CD4_base_u.csv", sep = ",", row.names=1, header = TRUE, check.names = FALSE, quote = "\"")
##################################################################
############## ----- (三) 全部纯疾病expr_u整合------ ############# 
##################################################################

rm(list=ls())
# 读取每个疾病的mean列--- 统一关键列名的_2.csv文件
COMBAT <- read.table("/data/wuwq/noise/COMBAT/expression/otherT_base_u_2.csv", sep = ",", row.names=1, header = TRUE, check.names = FALSE, quote = "\"")
FLU <- read.table("/data/wuwq/noise/FLU/expression/otherT_base_u_2.csv", sep = ",", row.names=1, header = TRUE, check.names = FALSE, quote = "\"")
SLE <- read.table("/data/wuwq/noise/SLE/expression/otherT_base_u_2.csv", sep = ",", row.names=1, header = TRUE, check.names = FALSE, quote = "\"")
RA <- read.table("/data/wuwq/noise/RA/expression/otherT_base_u_2.csv", sep = ",", row.names=1, header = TRUE, check.names = FALSE, quote = "\"")
longcovid <- read.table("/data/wuwq/noise/longcovid/expression/otherT_base_u_2.csv", sep = ",", row.names=1, header = TRUE, check.names = FALSE, quote = "\"")
NSCLC <- read.table("/data/wuwq/noise/NSCLC/expression/otherT_base_u_2.csv", sep = ",", row.names=1, header = TRUE, check.names = FALSE, quote = "\"")
# 去掉健康人的部分，留下纯疾病数据 -- 还是得用subset函数
unique(COMBAT$Disease_group)
COMBAT_p <- subset(COMBAT, Disease_group != "Health")
FLU_p <- subset(FLU, Disease_group != "Health")
SLE_p <- subset(SLE, Disease_group != "Health")
RA_p <- subset(RA, Disease_group != "Health")
longcovid_p <- subset(longcovid, Disease_group != "Health")
NSCLC_p <- subset(NSCLC, Disease_group != "Health")
# 转置
noise_1 <- as.data.frame(t(COMBAT_p))
noise_2 <- as.data.frame(t(FLU_p))
noise_3 <- as.data.frame(t(SLE_p))
noise_4 <- as.data.frame(t(RA_p))
noise_5 <- as.data.frame(t(longcovid_p))
noise_6 <- as.data.frame(t(NSCLC_p))
# 取$并集$
union_rows <- Reduce(union, list(rownames(noise_1), #rownames(noise_2), rownames(noise_3), 
                                 rownames(noise_4), rownames(noise_5), rownames(noise_6)))
# 取$交集$
intersect_rows <- Reduce(intersect, list(rownames(noise_1), #rownames(noise_2), rownames(noise_3), 
                                         rownames(noise_4), rownames(noise_5), rownames(noise_6)))
# 并集合并
merged <- data.frame(gene = union_rows)
for (i in c(1,4:6)) {
  dataset <- get(paste0("noise_", i))
  merged <- merge(merged, dataset, by.x = "gene", by.y = "row.names", all.x = TRUE, all.y = FALSE)
}
merged <- as.data.frame(t(merged))
#设置基因名为列名
colnames(merged) <- as.character(merged[1, ])
merged <- merged[-1, ]
# 把个体信息前移,并保留其余列
merged <- merged %>% select(dataset, sample_id, Disease_group, sex,   everything())
merged[, 5:ncol(merged)] <- lapply(merged[, 5:ncol(merged)], function(x) {
  x <- as.numeric(x)
  return(x)
})
write.table(merged,file = "/data/wuwq/noise/DISEASE_ALL/exp_u/union_onlydisease_otherT.csv", sep = ",", row.names = TRUE, col.names = NA, quote = TRUE)


# 交集合并 → 辅
merged <- data.frame(gene = intersect_rows)
for (i in c(1,4:6)) {
  dataset <- get(paste0("noise_", i))
  merged <- merge(merged, dataset, by.x = "gene", by.y = "row.names", all.x = TRUE, all.y = FALSE)
}
merged <- as.data.frame(t(merged))
#设置基因名为列名
colnames(merged) <- as.character(merged[1, ])
merged <- merged[-1, ]
# 把个体信息前移,并保留其余列
merged <- merged %>% select(dataset, sample_id, Disease_group, sex,   everything())
merged[, 5:ncol(merged)] <- lapply(merged[, 5:ncol(merged)], function(x) {
  x <- as.numeric(x)
  return(x)
})
write.table(merged,file = "/data/wuwq/noise/DISEASE_ALL/exp_u/intersect_onlydisease_otherT.csv", sep = ",", row.names = TRUE, col.names = NA, quote = TRUE)




##################################################################
############ ----- (四) 全部疾病+自带对照expr_u整合------ ########  直接接上面（三）
##################################################################
# 转置
noise_1 <- as.data.frame(t(COMBAT))
noise_2 <- as.data.frame(t(FLU))
noise_3 <- as.data.frame(t(SLE))
noise_4 <- as.data.frame(t(RA))
noise_5 <- as.data.frame(t(longcovid))
noise_6 <- as.data.frame(t(NSCLC))
# 取$并集$
union_rows <- Reduce(union, list(rownames(noise_1), #rownames(noise_2), rownames(noise_3), 
                                 rownames(noise_4), rownames(noise_5), rownames(noise_6)))
# 取$交集$
intersect_rows <- Reduce(intersect, list(rownames(noise_1), #rownames(noise_2), rownames(noise_3), 
                                         rownames(noise_4), rownames(noise_5), rownames(noise_6)))
# 并集合并
merged <- data.frame(gene = union_rows)
for (i in c(1,4:6)) {
  dataset <- get(paste0("noise_", i))
  merged <- merge(merged, dataset, by.x = "gene", by.y = "row.names", all.x = TRUE, all.y = FALSE)
}
merged <- as.data.frame(t(merged))
#设置基因名为列名
colnames(merged) <- as.character(merged[1, ])
merged <- merged[-1, ]
# 把个体信息前移,并保留其余列
merged <- merged %>% select(dataset, sample_id, Disease_group, sex,   everything())
merged[, 5:ncol(merged)] <- lapply(merged[, 5:ncol(merged)], function(x) {
  x <- as.numeric(x)
  return(x)
})
write.table(merged,file = "/data/wuwq/noise/DISEASE_ALL/exp_u/union_controldisease_otherT.csv", sep = ",", row.names = TRUE, col.names = NA, quote = TRUE)

# 交集合并 → 辅
merged <- data.frame(gene = intersect_rows)
for (i in c(1,4:6)) {
  dataset <- get(paste0("noise_", i))
  merged <- merge(merged, dataset, by.x = "gene", by.y = "row.names", all.x = TRUE, all.y = FALSE)
}
merged <- as.data.frame(t(merged))
#设置基因名为列名
colnames(merged) <- as.character(merged[1, ])
merged <- merged[-1, ]
# 把个体信息前移,并保留其余列
merged <- merged %>% select(dataset, sample_id, Disease_group, sex,   everything())
merged[, 5:ncol(merged)] <- lapply(merged[, 5:ncol(merged)], function(x) {
  x <- as.numeric(x)
  return(x)
})
write.table(merged,file = "/data/wuwq/noise/DISEASE_ALL/exp_u/intersect_controldisease_otherT.csv", sep = ",", row.names = TRUE, col.names = NA, quote = TRUE)




##################################################################
####################  ------- (五) 合category-------- ################   
##################################################################

rm(list=ls())
### 批量将 controldisease  合大类 (交集和并集 都)
exp_u <- read.table("/data/wuwq/noise/DISEASE_ALL/exp_u/intersect_controldisease_otherT.csv", sep = ",", row.names=1, header = TRUE, check.names = FALSE, quote = "\"")
unique(exp_u$Disease_group)
exp_u$category <- exp_u$Disease_group %>%
  recode(
    "COVID-19" = "infect",
    "post-COVID-19 disorder" = "infect",
    "SEPSIS" = "infect",
    "FLU" = "infect",
    "NSCLC" = "tumor",
    "RA" = "autoimmu",
    "SLE" = "autoimmu",
    "Health" = "health"
  )
exp_u <- exp_u %>% relocate(category, .before = 1)
write.table(exp_u,file = "/data/wuwq/noise/DISEASE_ALL/exp_u/intersect_controldisease_otherT.csv", sep = ",", row.names = TRUE, col.names = NA, quote = TRUE)




##################################################################
###########  ------- (六) 健康图谱的表达量合并 ------- ###########   
##################################################################
library(dplyr)
library(knitr)
library(tidyr)
library(stats)
library(tibble)
library(stringr)
rm(list=ls())
### 读取健康人的表达量表
aidas_u <- read.table("/data/wuwq/noise/HEALTH_REAL/expression/CD4_u_aidas.csv", sep = ",", row.names=1, header = TRUE, check.names = FALSE, quote = "\"")
eqtl_u <- read.table("/data/wuwq/noise/HEALTH_REAL/expression/CD4_u_eqtl.csv", sep = ",", row.names=1, header = TRUE, check.names = FALSE, quote = "\"")
row.names(aidas_u) <- aidas_u$gene
row.names(eqtl_u) <- eqtl_u$gene
# 提取基因名和基因符号
aidas_genes <- data.frame(
  gene = rownames(aidas_u),
  Gene_name = aidas_u$Gene_name,
  stringsAsFactors = FALSE
)
eqtl_genes <- data.frame(
  gene = rownames(eqtl_u),
  Gene_name = eqtl_u$Gene_name,
  stringsAsFactors = FALSE
)
# 获取所有唯一基因（并集）
all_genes <- unique(c(aidas_genes$gene, eqtl_genes$gene))
# 创建合并后的数据框
merged_data <- data.frame(
  gene = all_genes,
  stringsAsFactors = FALSE
)
# 合并Gene_name（优先使用aidas_u中的符号）
merged_data <- merge(merged_data, aidas_genes, by = "gene", all.x = TRUE)
merged_data <- merge(merged_data, eqtl_genes, by = "gene", all.x = TRUE, suffixes = c("", ".y"))
# 处理Gene_name列：优先使用aidas的符号，缺失时使用eqtl的
merged_data$Gene_name <- ifelse(
  is.na(merged_data$Gene_name),
  merged_data$Gene_name.y,
  merged_data$Gene_name
)
# 移除临时列
merged_data$Gene_name.y <- NULL
# 添加表达量数据
for (gene in all_genes) {
  # 处理aidas_u数据
  if (gene %in% rownames(aidas_u)) {
    merged_data[merged_data$gene == gene, colnames(aidas_u)[-1]] <- aidas_u[gene, -1]
  }
  # 处理eqtl_u数据
  if (gene %in% rownames(eqtl_u)) {
    merged_data[merged_data$gene == gene, colnames(eqtl_u)[-1]] <- eqtl_u[gene, -1]
  }
}
# 将NA值替换为0
merged_data[is.na(merged_data)] <- 0
# 重新排列列：gene, Gene_name, 然后是样本列
sample_cols <- setdiff(colnames(merged_data), c("gene", "Gene_name"))
merged_data <- merged_data[, c("gene", "Gene_name", sample_cols)]
# 去除从第3列开始的列名中的细胞类型后缀
merged_data <- merged_data %>%
  rename_with(~ sub("_CD4T$", "", .x), .cols = 3:ncol(.))

write.table(merged_data,file = "/data/wuwq/noise/HEALTH_REAL/expression/health_u_CD4.csv", sep = ",", row.names = TRUE, col.names = NA, quote = TRUE)




##################################################################
###########  ------- (七) 健康图谱和疾病的表达量合并 ------- ###########   
##################################################################

rm(list=ls())
### 疾病规整
exp_disease <- read.csv("/data/wuwq/noise/DISEASE_ALL/exp_u/union_onlydisease_otherT.csv", header = FALSE, row.names = 1)  # header = FALSE
colnames(exp_disease) <- as.character(exp_disease[1, ])
exp_disease <- exp_disease[-1, , drop = FALSE]
exp_disease <- exp_disease[,-c(2,4)]
col_order <- c("dataset", "Disease_group", setdiff(names(exp_disease), c("dataset", "Disease_group")))
exp_disease <- exp_disease[ , col_order]
unique(exp_disease$Disease_group)
### 健康规整
exp_health <- read.csv("/data/wuwq/noise/HEALTH_REAL/expression/health_u_otherT.csv", header = TRUE, row.names = 1) # header = FALSE
exp_health <- as.data.frame(t(exp_health))
### 将gene_sym定为列名，空值用第一行gene_id填充
exp_health <- exp_health %>%
  # 将列名设置为临时名称
  setNames(paste0("col", 1:ncol(.))) %>%
  # 创建新列名
  mutate(across(everything(), as.character)) %>%
  summarise(across(everything(), ~ ifelse(!is.na(.[2]) & .[2] != "", .[2], .[1]))) %>%
  unlist() %>%
  # 应用新列名并删除前两行
  {setNames(exp_health[-c(1, 2), ], .)} 
exp_health$dataset <- "health"
exp_health$Disease_group <- "health"
col_order <- c("dataset", "Disease_group", setdiff(names(exp_health), c("dataset", "Disease_group")))
exp_health <- exp_health[ , col_order]
# 合并全疾病和健康数据 # 把所有 "NA." 开头的列删除
exp <- bind_rows(
  exp_disease %>% select(-starts_with("NA.")),
  exp_health  %>% select(-starts_with("NA."))
) %>%
  mutate(
    across(
      -c(dataset, Disease_group),
      ~ as.numeric(as.character(.x))
    )
  )
### 大类定义
unique(exp$Disease_group)
exp$category <- exp$Disease_group %>%
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
exp <- exp %>% select(dataset, category, Disease_group,  everything())
### 总疾病定义
unique(exp$category)
exp$group <- exp$category %>%
  recode(
    "infect" = "disease",
    "autoimmu" = "disease",
    "tumor" = "disease",
    "health" = "health"
  )
exp <- exp %>% select(dataset, group, category, Disease_group,  everything())
category_counts <- exp %>%
  count(category, name = "n_samples", sort = TRUE)
print(category_counts)
write.table(exp,file = "/data/wuwq/noise/DISEASE_ALL/exp_u/union_healthdisease_otherT.csv", sep = ",", row.names = TRUE, col.names = NA, quote = TRUE)
















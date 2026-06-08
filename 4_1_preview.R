########################################################### 4_1:批量读取log2 fitratio结果 --- 临时搭框架用
### CD4+ T cell
### CD8+ T cell
### Monocyte
### NK cell
### Dendritic cell
### B cell
### Other T cell

##### COMBAT, FLU, longcovid, NSCLC, SLE, RA
# 基因名统一为gene_symbol

rm(list=ls())
gc()
##### (一)、FLU
###   CD4_fit, CD8_fit, MNP_fit, NK_fit, DC_fit, B_fit, otherT_fit
# 设置文件夹路径
folder_path <- "/data/wuwq/noise/FLU/output_root_dir/B_fit"
# 获取所有子文件夹路径
subfolders <- list.files(folder_path, full.names = TRUE, include.dirs = TRUE)

# 1.创建一个空的列表来存储所有读取的数据框
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
# 2.只保留fitratio>1且fdr<0.05的数据，存储在data_list_filter中
data_list_filter <- list()
data_list_filter <- lapply(data_list, function(df) {
  # 筛选出 fitratio > 1 且 fdr < 0.05 的行
  df_filtered <- df[df$fitratio > 1 & df$fdr < 0.05, ]
  # 返回筛选后的数据框
  return(df_filtered)
})


# 3.合并所有个体的log2fitratio结果，得到总的汇总表(无noise赋值为0)
# 获取所有数据框中的 gene 列的基因 #并集
all_genes <- unique(unlist(lapply(data_list_filter, function(df) df$g)))
# 创建一个空的数据框，用于存储合并后的结果
merged_data <- data.frame(gene = all_genes)
# 循环遍历每个 data_list_filter 中的元素，并合并 log2fitratio 列
for (file_name in names(data_list_filter)) {
  # 获取当前数据框
  df <- data_list_filter[[file_name]]
  # 创建一个空的 fitratio 列，默认值为 0
  log2fitratio_column <- rep(0, length(all_genes))
  # 使用 match 向量化匹配，更高效
  match_idx <- match(all_genes, df$g)
  valid <- !is.na(match_idx)
  log2fitratio_column[valid] <- df$log2fitratio[match_idx[valid]]
  merged_data[[file_name]] <- log2fitratio_column
}
##### 增加匹配的gene_symbol列
annotation_df <-read.delim("/data/wuwq/noise/mapping_geneid2symbol.txt",header = T)
# 使用 merge 左连接，保留 merged_data 中所有 gene
merged_data <- merge(merged_data, 
                     annotation_df, 
                     by.x = "gene",          # merged_data 中的 gene_id 列
                     by.y = "Gene_ID",       # 注释表中的 gene_id 列
                     all.x = TRUE)           # 保留所有 gene，即使没有匹配
# 把Gene_name列前移
col_order <- c("gene", "Gene_name", setdiff(names(merged_data), c("gene", "Gene_name")))
merged_data <- merged_data[, col_order]
write.csv(merged_data, file="/data/wuwq/noise/FLU/preview_result/B_log2fitratio.csv")



# 4.合并所有个体的mean结果，得到总的汇总表(无mean赋值为0)
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
write.csv(mean_data, file="/data/wuwq/noise/FLU/preview_result/B_u.csv")



# 5.获得个体 Adjusted_R_squared、cor_u_cv2.estimate、cor_u_log2fitratio.estimate、fitting_gene_count、细胞数m
target_cols <- c("Adjusted_R_squared", 
                 "cor_u_cv2.estimate", 
                 "cor_u_log2fitratio.estimate", 
                 "fitting_gene_count",
                 "m")
# 创建一个列表，用于存放每个指标的最终结果（个体水平的单值数据框）
result_list <- vector("list", length(target_cols))
names(result_list) <- target_cols
# 逐个指标处理
for (col_name in target_cols) {
  # 创建一个空的数据框，用于存储合并后的结果
  merged_info <- data.frame(gene = all_genes)
  # 循环遍历每个 data_list 中的元素
  for (file_name in names(data_list)) {
    # 获取当前数据框
    df <- data_list[[file_name]]
    # 创建一个空列，默认值为 0
    column <- rep(0, length(all_genes))
    for (i in 1:length(all_genes)) {
      gene_name <- all_genes[i]
      # 如果该基因在当前数据框中存在，则更新对应指标的值
      if (gene_name %in% df$g) {
        column[i] <- df[[col_name]][df$g == gene_name]
      }
    }
    # 将该列加入到 merged_info 数据框中
    merged_info[[file_name]] <- column
  }
  # 设置合并后的数据框的行名为 gene
  rownames(merged_info) <- merged_info$gene
  merged_info <- merged_info[, -1]  # 删除 gene 列
  # ===== 修改的地方：取每一列的第一个非零值 =====
  selected_row <- apply(merged_info, 2, function(col) {
    idx <- which(col != 0)[1]
    if (is.na(idx)) 0 else col[idx]
  })
  # ===== 关键修改：构建数据框时显式设置行名 =====
  merged_info <- data.frame(value = selected_row,               # 先创建一个带名字的向量
                            row.names = names(selected_row),   # 显式指定行名（样本名）
                            stringsAsFactors = FALSE)
  colnames(merged_info) <- col_name                            # 列名设为当前指标名
  # 保存到结果列表（保持为数据框）
  result_list[[col_name]] <- merged_info
}
# 最终合并所有指标
summary_df <- do.call(cbind, result_list)
write.csv(summary_df, "/data/wuwq/noise/FLU/preview_result/B_info.csv", row.names = TRUE, quote = FALSE)





###########################################################################################################################################
rm(list=ls())
##### (二)、COMBAT
###   CD4_fit, CD8_fit, MNP_fit, NK_fit, DC_fit, B_fit, otherT_fit
# 设置文件夹路径
folder_path <- "/data/wuwq/noise/COMBAT/output_root_dir/B_fit"
# 获取所有子文件夹路径
subfolders <- list.files(folder_path, full.names = TRUE, include.dirs = TRUE)

# 1.创建一个空的列表来存储所有读取的数据框
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
# 2.只保留fitratio>1且fdr<0.05的数据，存储在data_list_filter中
data_list_filter <- list()
data_list_filter <- lapply(data_list, function(df) {
  # 筛选出 fitratio > 1 且 fdr < 0.05 的行
  df_filtered <- df[df$fitratio > 1 & df$fdr < 0.05, ]
  # 返回筛选后的数据框
  return(df_filtered)
})

# 3.合并所有个体的log2fitratio结果，得到总的汇总表(无noise赋值为0)
# 获取所有数据框中的 gene 列的基因 #并集
all_genes <- unique(unlist(lapply(data_list_filter, function(df) df$g)))
# 创建一个空的数据框，用于存储合并后的结果
merged_data <- data.frame(gene = all_genes)
# 循环遍历每个 data_list_filter 中的元素，并合并 log2fitratio 列
for (file_name in names(data_list_filter)) {
  # 获取当前数据框
  df <- data_list_filter[[file_name]]
  # 创建一个空的 fitratio 列，默认值为 0
  log2fitratio_column <- rep(0, length(all_genes))
  # 使用 match 向量化匹配，更高效
  match_idx <- match(all_genes, df$g)
  valid <- !is.na(match_idx)
  log2fitratio_column[valid] <- df$log2fitratio[match_idx[valid]]
  merged_data[[file_name]] <- log2fitratio_column
}
##### 增加匹配的gene_symbol列
annotation_df <-read.delim("/data/wuwq/noise/mapping_geneid2symbol.txt",header = T)
# 使用 merge 左连接，保留 merged_data 中所有 gene
merged_data <- merge(merged_data, 
                     annotation_df, 
                     by.x = "gene",          # merged_data 中的 gene_id 列
                     by.y = "Gene_ID",       # 注释表中的 gene_id 列
                     all.x = TRUE)           # 保留所有 gene，即使没有匹配
# 把Gene_name列前移
col_order <- c("gene", "Gene_name", setdiff(names(merged_data), c("gene", "Gene_name")))
merged_data <- merged_data[, col_order]
write.csv(merged_data, file="/data/wuwq/noise/COMBAT/preview_result/B_log2fitratio.csv")

# 4.合并所有个体的mean结果，得到总的汇总表(无mean赋值为0)
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
write.csv(mean_data, file="/data/wuwq/noise/COMBAT/preview_result/B_u.csv")



# 5.获得个体 Adjusted_R_squared、cor_u_cv2.estimate、cor_u_log2fitratio.estimate、fitting_gene_count、细胞数m
target_cols <- c("Adjusted_R_squared", 
                 "cor_u_cv2.estimate", 
                 "cor_u_log2fitratio.estimate", 
                 "fitting_gene_count",
                 "m")
# 创建一个列表，用于存放每个指标的最终结果（个体水平的单值数据框）
result_list <- vector("list", length(target_cols))
names(result_list) <- target_cols
# 逐个指标处理
for (col_name in target_cols) {
  # 创建一个空的数据框，用于存储合并后的结果
  merged_info <- data.frame(gene = all_genes)
  # 循环遍历每个 data_list 中的元素
  for (file_name in names(data_list)) {
    # 获取当前数据框
    df <- data_list[[file_name]]
    # 创建一个空列，默认值为 0
    column <- rep(0, length(all_genes))
    for (i in 1:length(all_genes)) {
      gene_name <- all_genes[i]
      # 如果该基因在当前数据框中存在，则更新对应指标的值
      if (gene_name %in% df$g) {
        column[i] <- df[[col_name]][df$g == gene_name]
      }
    }
    # 将该列加入到 merged_info 数据框中
    merged_info[[file_name]] <- column
  }
  # 设置合并后的数据框的行名为 gene
  rownames(merged_info) <- merged_info$gene
  merged_info <- merged_info[, -1]  # 删除 gene 列
  # ===== 修改的地方：取每一列的第一个非零值 =====
  selected_row <- apply(merged_info, 2, function(col) {
    idx <- which(col != 0)[1]
    if (is.na(idx)) 0 else col[idx]
  })
  # ===== 关键修改：构建数据框时显式设置行名 =====
  merged_info <- data.frame(value = selected_row,               # 先创建一个带名字的向量
                            row.names = names(selected_row),   # 显式指定行名（样本名）
                            stringsAsFactors = FALSE)
  colnames(merged_info) <- col_name                            # 列名设为当前指标名
  # 保存到结果列表（保持为数据框）
  result_list[[col_name]] <- merged_info
}
# 最终合并所有指标
summary_df <- do.call(cbind, result_list)
write.csv(summary_df, "/data/wuwq/noise/COMBAT/preview_result/B_info.csv", row.names = TRUE, quote = FALSE)






###########################################################################################################################################
rm(list=ls())
##### (三)、SLE
###   CD4_fit, CD8_fit, MNP_fit, NK_fit, DC_fit, B_fit, otherT_fit
# 设置文件夹路径
folder_path <- "/data/wuwq/noise/SLE/output_root_dir/CD4_fit"
# 获取所有子文件夹路径
subfolders <- list.files(folder_path, full.names = TRUE, include.dirs = TRUE)

# 1.创建一个空的列表来存储所有读取的数据框
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
# 2.只保留fitratio>1且fdr<0.05的数据，存储在data_list_filter中
data_list_filter <- list()
data_list_filter <- lapply(data_list, function(df) {
  # 筛选出 fitratio > 1 且 fdr < 0.05 的行
  df_filtered <- df[df$fitratio > 1 & df$fdr < 0.05, ]
  # 返回筛选后的数据框
  return(df_filtered)
})

# 3.合并所有个体的log2fitratio结果，得到总的汇总表(无noise赋值为0)
# 获取所有数据框中的 gene 列的基因 #并集
all_genes <- unique(unlist(lapply(data_list_filter, function(df) df$g)))
# 创建一个空的数据框，用于存储合并后的结果
merged_data <- data.frame(gene = all_genes)
# 循环遍历每个 data_list_filter 中的元素，并合并 log2fitratio 列
for (file_name in names(data_list_filter)) {
  # 获取当前数据框
  df <- data_list_filter[[file_name]]
  # 创建一个空的 fitratio 列，默认值为 0
  log2fitratio_column <- rep(0, length(all_genes))
  # 使用 match 向量化匹配，更高效
  match_idx <- match(all_genes, df$g)
  valid <- !is.na(match_idx)
  log2fitratio_column[valid] <- df$log2fitratio[match_idx[valid]]
  merged_data[[file_name]] <- log2fitratio_column
}
##### 增加匹配的gene_symbol列
annotation_df <-read.delim("/data/wuwq/noise/mapping_geneid2symbol.txt",header = T)
# 使用 merge 左连接，保留 merged_data 中所有 gene
merged_data <- merge(merged_data, 
                     annotation_df, 
                     by.x = "gene",          # merged_data 中的 gene_id 列
                     by.y = "Gene_ID",       # 注释表中的 gene_id 列
                     all.x = TRUE)           # 保留所有 gene，即使没有匹配
# 把Gene_name列前移
col_order <- c("gene", "Gene_name", setdiff(names(merged_data), c("gene", "Gene_name")))
merged_data <- merged_data[, col_order]
write.csv(merged_data, file="/data/wuwq/noise/SLE/preview_result/CD4_log2fitratio.csv")

# 4.合并所有个体的mean结果，得到总的汇总表(无mean赋值为0)
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
write.csv(mean_data, file="/data/wuwq/noise/SLE/preview_result/CD4_u.csv")



# 5.获得个体 Adjusted_R_squared、cor_u_cv2.estimate、cor_u_log2fitratio.estimate、fitting_gene_count、细胞数m
target_cols <- c("Adjusted_R_squared", 
                 "cor_u_cv2.estimate", 
                 "cor_u_log2fitratio.estimate", 
                 "fitting_gene_count",
                 "m")
# 创建一个列表，用于存放每个指标的最终结果（个体水平的单值数据框）
result_list <- vector("list", length(target_cols))
names(result_list) <- target_cols
# 逐个指标处理
for (col_name in target_cols) {
  # 创建一个空的数据框，用于存储合并后的结果
  merged_info <- data.frame(gene = all_genes)
  # 循环遍历每个 data_list 中的元素
  for (file_name in names(data_list)) {
    # 获取当前数据框
    df <- data_list[[file_name]]
    # 创建一个空列，默认值为 0
    column <- rep(0, length(all_genes))
    for (i in 1:length(all_genes)) {
      gene_name <- all_genes[i]
      # 如果该基因在当前数据框中存在，则更新对应指标的值
      if (gene_name %in% df$g) {
        column[i] <- df[[col_name]][df$g == gene_name]
      }
    }
    # 将该列加入到 merged_info 数据框中
    merged_info[[file_name]] <- column
  }
  # 设置合并后的数据框的行名为 gene
  rownames(merged_info) <- merged_info$gene
  merged_info <- merged_info[, -1]  # 删除 gene 列
  # ===== 修改的地方：取每一列的第一个非零值 =====
  selected_row <- apply(merged_info, 2, function(col) {
    idx <- which(col != 0)[1]
    if (is.na(idx)) 0 else col[idx]
  })
  # ===== 关键修改：构建数据框时显式设置行名 =====
  merged_info <- data.frame(value = selected_row,               # 先创建一个带名字的向量
                            row.names = names(selected_row),   # 显式指定行名（样本名）
                            stringsAsFactors = FALSE)
  colnames(merged_info) <- col_name                            # 列名设为当前指标名
  # 保存到结果列表（保持为数据框）
  result_list[[col_name]] <- merged_info
}
# 最终合并所有指标
summary_df <- do.call(cbind, result_list)
write.csv(summary_df, "/data/wuwq/noise/SLE/preview_result/CD4_info.csv", row.names = TRUE, quote = FALSE)




###########################################################################################################################################
rm(list=ls())
##### (四)、NSCLC
###   CD4_fit, CD8_fit, MNP_fit, NK_fit, DC_fit, B_fit, otherT_fit
# 设置文件夹路径
folder_path <- "/data/wuwq/noise/NSCLC/output_root_dir/CD4_fit"
# 获取所有子文件夹路径
subfolders <- list.files(folder_path, full.names = TRUE, include.dirs = TRUE)

# 1.创建一个空的列表来存储所有读取的数据框
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
# 2.只保留fitratio>1且fdr<0.05的数据，存储在data_list_filter中
data_list_filter <- list()
data_list_filter <- lapply(data_list, function(df) {
  # 筛选出 fitratio > 1 且 fdr < 0.05 的行
  df_filtered <- df[df$fitratio > 1 & df$fdr < 0.05, ]
  # 返回筛选后的数据框
  return(df_filtered)
})

# 3.合并所有个体的log2fitratio结果，得到总的汇总表(无noise赋值为0)
# 获取所有数据框中的 gene 列的基因 #并集
all_genes <- unique(unlist(lapply(data_list_filter, function(df) df$g)))
# 创建一个空的数据框，用于存储合并后的结果
merged_data <- data.frame(gene = all_genes)
# 循环遍历每个 data_list_filter 中的元素，并合并 log2fitratio 列
for (file_name in names(data_list_filter)) {
  # 获取当前数据框
  df <- data_list_filter[[file_name]]
  # 创建一个空的 fitratio 列，默认值为 0
  log2fitratio_column <- rep(0, length(all_genes))
  # 使用 match 向量化匹配，更高效
  match_idx <- match(all_genes, df$g)
  valid <- !is.na(match_idx)
  log2fitratio_column[valid] <- df$log2fitratio[match_idx[valid]]
  merged_data[[file_name]] <- log2fitratio_column
}
##### 增加匹配的gene_symbol列
annotation_df <-read.delim("/data/wuwq/noise/mapping_geneid2symbol.txt",header = T)
# 使用 merge 左连接，保留 merged_data 中所有 gene
merged_data <- merge(merged_data, 
                     annotation_df, 
                     by.x = "gene",          # merged_data 中的 gene_id 列
                     by.y = "Gene_ID",       # 注释表中的 gene_id 列
                     all.x = TRUE)           # 保留所有 gene，即使没有匹配
# 把Gene_name列前移
col_order <- c("gene", "Gene_name", setdiff(names(merged_data), c("gene", "Gene_name")))
merged_data <- merged_data[, col_order]
write.csv(merged_data, file="/data/wuwq/noise/NSCLC/preview_result/CD4_log2fitratio.csv")

# 4.合并所有个体的mean结果，得到总的汇总表(无mean赋值为0)
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
write.csv(mean_data, file="/data/wuwq/noise/NSCLC/preview_result/CD4_u.csv")



# 5.获得个体 Adjusted_R_squared、cor_u_cv2.estimate、cor_u_log2fitratio.estimate、fitting_gene_count、细胞数m
target_cols <- c("Adjusted_R_squared", 
                 "cor_u_cv2.estimate", 
                 "cor_u_log2fitratio.estimate", 
                 "fitting_gene_count",
                 "m")
# 创建一个列表，用于存放每个指标的最终结果（个体水平的单值数据框）
result_list <- vector("list", length(target_cols))
names(result_list) <- target_cols
# 逐个指标处理
for (col_name in target_cols) {
  # 创建一个空的数据框，用于存储合并后的结果
  merged_info <- data.frame(gene = all_genes)
  # 循环遍历每个 data_list 中的元素
  for (file_name in names(data_list)) {
    # 获取当前数据框
    df <- data_list[[file_name]]
    # 创建一个空列，默认值为 0
    column <- rep(0, length(all_genes))
    for (i in 1:length(all_genes)) {
      gene_name <- all_genes[i]
      # 如果该基因在当前数据框中存在，则更新对应指标的值
      if (gene_name %in% df$g) {
        column[i] <- df[[col_name]][df$g == gene_name]
      }
    }
    # 将该列加入到 merged_info 数据框中
    merged_info[[file_name]] <- column
  }
  # 设置合并后的数据框的行名为 gene
  rownames(merged_info) <- merged_info$gene
  merged_info <- merged_info[, -1]  # 删除 gene 列
  # ===== 修改的地方：取每一列的第一个非零值 =====
  selected_row <- apply(merged_info, 2, function(col) {
    idx <- which(col != 0)[1]
    if (is.na(idx)) 0 else col[idx]
  })
  # ===== 关键修改：构建数据框时显式设置行名 =====
  merged_info <- data.frame(value = selected_row,               # 先创建一个带名字的向量
                            row.names = names(selected_row),   # 显式指定行名（样本名）
                            stringsAsFactors = FALSE)
  colnames(merged_info) <- col_name                            # 列名设为当前指标名
  # 保存到结果列表（保持为数据框）
  result_list[[col_name]] <- merged_info
}
# 最终合并所有指标
summary_df <- do.call(cbind, result_list)
write.csv(summary_df, "/data/wuwq/noise/NSCLC/preview_result/CD4_info.csv", row.names = TRUE, quote = FALSE)



###########################################################################################################################################
rm(list=ls())
##### (五)、RA
###   CD4_fit, CD8_fit, MNP_fit, NK_fit, DC_fit, B_fit, otherT_fit
# 设置文件夹路径
folder_path <- "/data/wuwq/noise/RA/output_root_dir/CD4_fit"
# 获取所有子文件夹路径
subfolders <- list.files(folder_path, full.names = TRUE, include.dirs = TRUE)

# 1.创建一个空的列表来存储所有读取的数据框
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
# 2.只保留fitratio>1且fdr<0.05的数据，存储在data_list_filter中
data_list_filter <- list()
data_list_filter <- lapply(data_list, function(df) {
  # 筛选出 fitratio > 1 且 fdr < 0.05 的行
  df_filtered <- df[df$fitratio > 1 & df$fdr < 0.05, ]
  # 返回筛选后的数据框
  return(df_filtered)
})

# 3.合并所有个体的log2fitratio结果，得到总的汇总表(无noise赋值为0)
# 获取所有数据框中的 gene 列的基因 #并集
all_genes <- unique(unlist(lapply(data_list_filter, function(df) df$g)))
# 创建一个空的数据框，用于存储合并后的结果
merged_data <- data.frame(gene = all_genes)
# 循环遍历每个 data_list_filter 中的元素，并合并 log2fitratio 列
for (file_name in names(data_list_filter)) {
  # 获取当前数据框
  df <- data_list_filter[[file_name]]
  # 创建一个空的 fitratio 列，默认值为 0
  log2fitratio_column <- rep(0, length(all_genes))
  # 使用 match 向量化匹配，更高效
  match_idx <- match(all_genes, df$g)
  valid <- !is.na(match_idx)
  log2fitratio_column[valid] <- df$log2fitratio[match_idx[valid]]
  merged_data[[file_name]] <- log2fitratio_column
}
##### 增加匹配的gene_symbol列
annotation_df <-read.delim("/data/wuwq/noise/mapping_geneid2symbol.txt",header = T)
# 使用 merge 左连接，保留 merged_data 中所有 gene
merged_data <- merge(merged_data, 
                     annotation_df, 
                     by.x = "gene",          # merged_data 中的 gene_id 列
                     by.y = "Gene_ID",       # 注释表中的 gene_id 列
                     all.x = TRUE)           # 保留所有 gene，即使没有匹配
# 把Gene_name列前移
col_order <- c("gene", "Gene_name", setdiff(names(merged_data), c("gene", "Gene_name")))
merged_data <- merged_data[, col_order]
write.csv(merged_data, file="/data/wuwq/noise/RA/preview_result/CD4_log2fitratio.csv")

# 4.合并所有个体的mean结果，得到总的汇总表(无mean赋值为0)
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
write.csv(mean_data, file="/data/wuwq/noise/RA/preview_result/CD4_u.csv")



# 5.获得个体 Adjusted_R_squared、cor_u_cv2.estimate、cor_u_log2fitratio.estimate、fitting_gene_count、细胞数m
target_cols <- c("Adjusted_R_squared", 
                 "cor_u_cv2.estimate", 
                 "cor_u_log2fitratio.estimate", 
                 "fitting_gene_count",
                 "m")
# 创建一个列表，用于存放每个指标的最终结果（个体水平的单值数据框）
result_list <- vector("list", length(target_cols))
names(result_list) <- target_cols
# 逐个指标处理
for (col_name in target_cols) {
  # 创建一个空的数据框，用于存储合并后的结果
  merged_info <- data.frame(gene = all_genes)
  # 循环遍历每个 data_list 中的元素
  for (file_name in names(data_list)) {
    # 获取当前数据框
    df <- data_list[[file_name]]
    # 创建一个空列，默认值为 0
    column <- rep(0, length(all_genes))
    for (i in 1:length(all_genes)) {
      gene_name <- all_genes[i]
      # 如果该基因在当前数据框中存在，则更新对应指标的值
      if (gene_name %in% df$g) {
        column[i] <- df[[col_name]][df$g == gene_name]
      }
    }
    # 将该列加入到 merged_info 数据框中
    merged_info[[file_name]] <- column
  }
  # 设置合并后的数据框的行名为 gene
  rownames(merged_info) <- merged_info$gene
  merged_info <- merged_info[, -1]  # 删除 gene 列
  # ===== 修改的地方：取每一列的第一个非零值 =====
  selected_row <- apply(merged_info, 2, function(col) {
    idx <- which(col != 0)[1]
    if (is.na(idx)) 0 else col[idx]
  })
  # ===== 关键修改：构建数据框时显式设置行名 =====
  merged_info <- data.frame(value = selected_row,               # 先创建一个带名字的向量
                            row.names = names(selected_row),   # 显式指定行名（样本名）
                            stringsAsFactors = FALSE)
  colnames(merged_info) <- col_name                            # 列名设为当前指标名
  # 保存到结果列表（保持为数据框）
  result_list[[col_name]] <- merged_info
}
# 最终合并所有指标
summary_df <- do.call(cbind, result_list)
write.csv(summary_df, "/data/wuwq/noise/RA/preview_result/CD4_info.csv", row.names = TRUE, quote = FALSE)







###########################################################################################################################################
rm(list=ls())
##### (六)、longcovid
###   CD4_fit, CD8_fit, MNP_fit, NK_fit, DC_fit, B_fit, otherT_fit
# 设置文件夹路径
folder_path <- "/data/wuwq/noise/longcovid/output_root_dir/CD4_fit"
# 获取所有子文件夹路径
subfolders <- list.files(folder_path, full.names = TRUE, include.dirs = TRUE)

# 1.创建一个空的列表来存储所有读取的数据框
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
# 2.只保留fitratio>1且fdr<0.05的数据，存储在data_list_filter中
data_list_filter <- list()
data_list_filter <- lapply(data_list, function(df) {
  # 筛选出 fitratio > 1 且 fdr < 0.05 的行
  df_filtered <- df[df$fitratio > 1 & df$fdr < 0.05, ]
  # 返回筛选后的数据框
  return(df_filtered)
})

# 3.合并所有个体的log2fitratio结果，得到总的汇总表(无noise赋值为0)
# 获取所有数据框中的 gene 列的基因 #并集
all_genes <- unique(unlist(lapply(data_list_filter, function(df) df$g)))
# 创建一个空的数据框，用于存储合并后的结果
merged_data <- data.frame(gene = all_genes)
# 循环遍历每个 data_list_filter 中的元素，并合并 log2fitratio 列
for (file_name in names(data_list_filter)) {
  # 获取当前数据框
  df <- data_list_filter[[file_name]]
  # 创建一个空的 fitratio 列，默认值为 0
  log2fitratio_column <- rep(0, length(all_genes))
  # 使用 match 向量化匹配，更高效
  match_idx <- match(all_genes, df$g)
  valid <- !is.na(match_idx)
  log2fitratio_column[valid] <- df$log2fitratio[match_idx[valid]]
  merged_data[[file_name]] <- log2fitratio_column
}
##### 增加匹配的gene_symbol列
annotation_df <-read.delim("/data/wuwq/noise/mapping_geneid2symbol.txt",header = T)
# 使用 merge 左连接，保留 merged_data 中所有 gene
merged_data <- merge(merged_data, 
                     annotation_df, 
                     by.x = "gene",          # merged_data 中的 gene_id 列
                     by.y = "Gene_ID",       # 注释表中的 gene_id 列
                     all.x = TRUE)           # 保留所有 gene，即使没有匹配
# 把Gene_name列前移
col_order <- c("gene", "Gene_name", setdiff(names(merged_data), c("gene", "Gene_name")))
merged_data <- merged_data[, col_order]
write.csv(merged_data, file="/data/wuwq/noise/longcovid/preview_result/CD4_log2fitratio.csv")

# 4.合并所有个体的mean结果，得到总的汇总表(无mean赋值为0)
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
write.csv(mean_data, file="/data/wuwq/noise/longcovid/preview_result/CD4_u.csv")



# 5.获得个体 Adjusted_R_squared、cor_u_cv2.estimate、cor_u_log2fitratio.estimate、fitting_gene_count、细胞数m
target_cols <- c("Adjusted_R_squared", 
                 "cor_u_cv2.estimate", 
                 "cor_u_log2fitratio.estimate", 
                 "fitting_gene_count",
                 "m")
# 创建一个列表，用于存放每个指标的最终结果（个体水平的单值数据框）
result_list <- vector("list", length(target_cols))
names(result_list) <- target_cols
# 逐个指标处理
for (col_name in target_cols) {
  # 创建一个空的数据框，用于存储合并后的结果
  merged_info <- data.frame(gene = all_genes)
  # 循环遍历每个 data_list 中的元素
  for (file_name in names(data_list)) {
    # 获取当前数据框
    df <- data_list[[file_name]]
    # 创建一个空列，默认值为 0
    column <- rep(0, length(all_genes))
    for (i in 1:length(all_genes)) {
      gene_name <- all_genes[i]
      # 如果该基因在当前数据框中存在，则更新对应指标的值
      if (gene_name %in% df$g) {
        column[i] <- df[[col_name]][df$g == gene_name]
      }
    }
    # 将该列加入到 merged_info 数据框中
    merged_info[[file_name]] <- column
  }
  # 设置合并后的数据框的行名为 gene
  rownames(merged_info) <- merged_info$gene
  merged_info <- merged_info[, -1]  # 删除 gene 列
  # ===== 修改的地方：取每一列的第一个非零值 =====
  selected_row <- apply(merged_info, 2, function(col) {
    idx <- which(col != 0)[1]
    if (is.na(idx)) 0 else col[idx]
  })
  # ===== 关键修改：构建数据框时显式设置行名 =====
  merged_info <- data.frame(value = selected_row,               # 先创建一个带名字的向量
                            row.names = names(selected_row),   # 显式指定行名（样本名）
                            stringsAsFactors = FALSE)
  colnames(merged_info) <- col_name                            # 列名设为当前指标名
  # 保存到结果列表（保持为数据框）
  result_list[[col_name]] <- merged_info
}
# 最终合并所有指标
summary_df <- do.call(cbind, result_list)
write.csv(summary_df, "/data/wuwq/noise/longcovid/preview_result/CD4_info.csv", row.names = TRUE, quote = FALSE)





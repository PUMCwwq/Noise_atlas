##################################################
library(Seurat)
library(dplyr)
library(Matrix)
rm(list=ls())
gc()
##################################################
##### COMBAT, FLU, longcovid, NSCLC, SLE, RA

##### 获取细胞类型注释及合并细胞类型
### COMBAT
# 读取SCT后seurat对象
SCT_COMBAT <- readRDS("/data2/wuwq/noise/data_COVID19/round5_COMBAT/SCT_COMBAT.rds")
# 获取 SCT_COMBAT 中的 cell_type 信息
cell_types <- SCT_COMBAT@meta.data$cell_type
unique(cell_types)
merged_cell_types <- recode(cell_types,
                            "classical monocyte" = "Monocyte",
                            "non-classical monocyte" = "Monocyte",
                            "CD4-positive, alpha-beta T cell" = "CD4+ T cell",
                            "CD8-positive, alpha-beta T cell" = "CD8+ T cell",
                            "gamma-delta T cell" = "Other T cell",
                            "mucosal invariant T cell" = "Other T cell",
                            "double negative T regulatory cell" = "Other T cell",
                            "B cell" = "B cell",
                            "plasmablast" = "B cell",
                            "natural killer cell" = "NK cell",
                            "dendritic cell" = "Dendritic cell",
                            "blood cell" = "Other",
                            "double-positive, alpha-beta thymocyte" = "Other",
                            "hematopoietic stem cell" = "Other",
                            "megakaryocyte-erythroid progenitor cell" = "Other",
                            "enucleated reticulocyte" = "Other",
                            "mast cell" = "Other"
)
# 将合并后的细胞类型存储到新的列中
SCT_COMBAT@meta.data$merged_cell_type <- merged_cell_types

# 样本分割准备
merged_cell_types <- SCT_COMBAT@meta.data$merged_cell_type
unique_cell_types <- unique(merged_cell_types)
seurat_list <- list()
for (celltype in unique_cell_types) {
  print(paste("Processing cell type:", celltype))
  # 正确的 subset 条件
  seurat_subset <- subset(SCT_COMBAT, subset = merged_cell_type == celltype)
  # 打印细胞数量和合并后的细胞类型
  print(paste("Number of cells after subsetting for", celltype, ":", length(seurat_subset$merged_cell_type)))
  print(table(seurat_subset$merged_cell_type))
  # 将结果存入列表
  seurat_list[[celltype]] <- seurat_subset
}

# 样本分割
unique(SCT_COMBAT@meta.data$scRNASeq_sample_ID)
# 按照scRNASeq_sample_ID 及  seurat_list中的 merged_cell_type 拆分数据  
output_dir <- "/data/wuwq/noise/COMBAT/rds_files"  # 修改为你想存储结果的文件夹路径
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
# 打印 seurat_list 的长度和基本信息，确保数据被正确加载
print(paste("Length of seurat_list:", length(seurat_list)))  # 打印seurat_list的长度
if (length(seurat_list) == 0) {
  print("Error: seurat_list is empty. Please check the data loading process.")
  stop("seurat_list is empty.")  # 停止执行，如果seurat_list为空
} else {
  print("seurat_list is loaded correctly.")
}
# 外层循环：遍历seurat_list中的每一个细胞类型
for (cell_type in names(seurat_list)) {
  # 获取当前细胞类型的Seurat对象
  seurat_object <- seurat_list[[cell_type]]
  # 检查seurat_object是否有效
  if (is.null(seurat_object)) {
    print(paste("Warning: Seurat object for cell type", cell_type, "is NULL"))
    next
  }
  # 获取counts数据和scRNASeq_sample_ID
  counts_data <- seurat_object@assays$SCT@counts
  sample_ids <- seurat_object@meta.data$scRNASeq_sample_ID
  # 检查counts_data是否为空
  if (is.null(counts_data) || ncol(counts_data) == 0) {
    print(paste("Warning: Counts data is empty for cell type", cell_type))
    next
  }
  # 检查sample_ids是否为空
  if (length(sample_ids) == 0) {
    print(paste("Warning: No sample IDs found for cell type", cell_type))
    next
  }
  # 获取所有唯一的样本ID
  unique_sample_ids <- unique(sample_ids)
  # 创建当前细胞类型的子文件夹来存储分割后的数据
  cell_type_dir <- file.path(output_dir, cell_type)
  if (!dir.exists(cell_type_dir)) {
    dir.create(cell_type_dir)
  }
  # 内层循环：根据样本ID分割数据并保存
  for (sample_id in unique_sample_ids) {
    # 打印调试信息
    print(paste("Processing sample:", sample_id, "for cell type:", cell_type))
    # 获取当前样本ID对应的有效列
    valid_columns <- sample_ids == sample_id
    print(paste("Valid columns:", sum(valid_columns)))  # 打印有效列的数量
    # 打印 valid_columns 的具体情况，确保其正确性
    print(paste("Sample ID:", sample_id))
    print(paste("Number of valid columns:", sum(valid_columns)))
    # 查看 valid_columns 对应的列索引
    valid_column_indices <- which(valid_columns)
    print("Valid column indices:")
    print(valid_column_indices)
    # 如果没有有效列，则跳过当前样本
    if (sum(valid_columns) == 0) {
      print(paste("Warning: No valid data for sample", sample_id, "- Skipping"))
      next
    }
    # 提取数据并保持矩阵结构
    sample_data <- as.matrix(counts_data[, valid_columns, drop = FALSE])
    # 打印 sample_data 的内容，检查数据是否正确
    print("Sample data content:")
    print(sample_data)
    # 检查 sample_data 是否为 NULL 或没有有效数据
    if (is.null(sample_data) || ncol(sample_data) == 0 || all(sample_data == 0)) {
      # 如果样本数据全为零，则跳过该样本的处理
      print(paste("Warning: Sample", sample_id, "has all zero values. Skipping this sample."))
    } else if (any(is.na(sample_data))) {
      # 如果样本数据中包含 NA 值，则跳过该样本的处理
      print(paste("Warning: Sample", sample_id, "contains NA values. Skipping this sample."))
    } else {
      # 如果样本数据有效，处理并保存为 RDS 文件
      print(paste("Valid data for sample", sample_id))
      # 保存 sample_data 为 .rds 文件
      output_file <- file.path(cell_type_dir, paste0(sample_id, "_", cell_type, "_processed.rds"))
      saveRDS(sample_data, file = output_file)
      print(paste("Saved:", output_file))
    }
  }
}


### FLU
rm(list = setdiff(ls(), c("SCT_flu")))
# 读取SCT后seurat对象
SCT_flu <- readRDS("/data2/wuwq/noise/data_COVID19/round7_flu/SCT_flu.rds")
# 获取 SCT_flu 中的 cell_type 信息
cell_types <- SCT_flu@meta.data$cell_type
merged_cell_types <- recode(cell_types,
                            "classical monocyte" = "Monocyte",
                            "intermediate monocyte" = "Monocyte",
                            "non-classical monocyte" = "Monocyte",
                            "CD4-positive helper T cell" = "CD4+ T cell",
                            "effector CD4-positive, alpha-beta T cell" = "CD4+ T cell",
                            "CD8-positive, alpha-beta T cell" = "CD8+ T cell",
                            "effector CD8-positive, alpha-beta T cell" = "CD8+ T cell",
                            "IgG memory B cell" = "B cell",
                            "IgG-negative class switched memory B cell" = "B cell",
                            "natural killer cell" = "NK cell",
                            "dendritic cell" = "Dendritic cell",
                            "platelet" = "Other",
                            "blood cell" = "Other",
                            "erythrocyte" = "Other")
SCT_flu@meta.data$merged_cell_type <- merged_cell_types

# 提取细胞barcode和对应的merged_celltype -- Monod分割用
cell_data <- data.frame(
  barcode = colnames(SCT_flu),  # 获取细胞条形码
  merged_cell_type = SCT_flu@meta.data$merged_cell_type,  # 获取细胞类型注释
  sample_id = SCT_flu@meta.data$Sample.ID
)
unique(colnames(SCT_flu))
write.csv(cell_data, file = "/data2/wuwq/noise/monod_FLU/loom/flu_merged_celltype.csv", row.names = FALSE)


# 分割前准备
merged_cell_types <- SCT_flu@meta.data$merged_cell_type
unique_cell_types <- unique(merged_cell_types)
seurat_list <- list()
for (celltype in unique_cell_types) {
  print(paste("Processing cell type:", celltype))
  seurat_subset <- subset(SCT_flu, subset = merged_cell_type == celltype)
  print(paste("Number of cells after subsetting for", celltype, ":", length(seurat_subset$merged_cell_type)))
  print(table(seurat_subset$merged_cell_type))
  seurat_list[[celltype]] <- seurat_subset
}

# 按照Sample.ID 及  seurat_list中的merged_cell_type 拆分数据   
output_dir <- "/data/wuwq/noise/FLU/rds_files"  # 修改为你想存储结果的文件夹路径
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# 打印 seurat_list 的长度和基本信息，确保数据被正确加载
print(paste("Length of seurat_list:", length(seurat_list)))  # 打印seurat_list的长度
if (length(seurat_list) == 0) {
  print("Error: seurat_list is empty. Please check the data loading process.")
  stop("seurat_list is empty.")  # 停止执行，如果seurat_list为空
} else {
  print("seurat_list is loaded correctly.")
}
# 外层循环：遍历seurat_list中的每一个细胞类型
for (cell_type in names(seurat_list)) {
  # 获取当前细胞类型的Seurat对象
  seurat_object <- seurat_list[[cell_type]]
  
  # 检查seurat_object是否有效
  if (is.null(seurat_object)) {
    print(paste("Warning: Seurat object for cell type", cell_type, "is NULL"))
    next
  }
  # 获取counts数据和scRNASeq_sample_ID
  counts_data <- seurat_object@assays$SCT@counts
  sample_ids <- seurat_object@meta.data$Sample.ID
  # 检查counts_data是否为空
  if (is.null(counts_data) || ncol(counts_data) == 0) {
    print(paste("Warning: Counts data is empty for cell type", cell_type))
    next
  }
  # 检查sample_ids是否为空
  if (length(sample_ids) == 0) {
    print(paste("Warning: No sample IDs found for cell type", cell_type))
    next
  }
  # 获取所有唯一的样本ID
  unique_sample_ids <- unique(sample_ids)
  # 创建当前细胞类型的子文件夹来存储分割后的数据
  cell_type_dir <- file.path(output_dir, cell_type)
  if (!dir.exists(cell_type_dir)) {
    dir.create(cell_type_dir)
  }
  # 内层循环：根据样本ID分割数据并保存
  for (sample_id in unique_sample_ids) {
    # 打印调试信息
    print(paste("Processing sample:", sample_id, "for cell type:", cell_type))
    # 获取当前样本ID对应的有效列
    valid_columns <- sample_ids == sample_id
    print(paste("Valid columns:", sum(valid_columns)))  # 打印有效列的数量
    # 打印 valid_columns 的具体情况，确保其正确性
    print(paste("Sample ID:", sample_id))
    print(paste("Number of valid columns:", sum(valid_columns)))
    # 查看 valid_columns 对应的列索引
    valid_column_indices <- which(valid_columns)
    print("Valid column indices:")
    print(valid_column_indices)
    # 如果没有有效列，则跳过当前样本
    if (sum(valid_columns) == 0) {
      print(paste("Warning: No valid data for sample", sample_id, "- Skipping"))
      next
    }
    # 提取数据并保持矩阵结构
    sample_data <- as.matrix(counts_data[, valid_columns, drop = FALSE])
    # 打印 sample_data 的内容，检查数据是否正确
    print("Sample data content:")
    print(sample_data)
    # 检查 sample_data 是否为 NULL 或没有有效数据
    if (is.null(sample_data) || ncol(sample_data) == 0 || all(sample_data == 0)) {
      # 如果样本数据全为零，则跳过该样本的处理
      print(paste("Warning: Sample", sample_id, "has all zero values. Skipping this sample."))
    } else if (any(is.na(sample_data))) {
      # 如果样本数据中包含 NA 值，则跳过该样本的处理
      print(paste("Warning: Sample", sample_id, "contains NA values. Skipping this sample."))
    } else {
      # 如果样本数据有效，处理并保存为 RDS 文件
      print(paste("Valid data for sample", sample_id))
      # 保存 sample_data 为 .rds 文件
      output_file <- file.path(cell_type_dir, paste0(sample_id, "_", cell_type, "_processed.rds"))
      saveRDS(sample_data, file = output_file)
      print(paste("Saved:", output_file))
    }
  }
}



### longcovid
rm(list = setdiff(ls(), c("SCT_longcovid")))
SCT_longcovid <- readRDS("/data2/wuwq/noise/data_COVID19/round6_long_covid/SCT_longcovid.rds")
cell_types <- SCT_longcovid@meta.data$cell_type
merged_cell_types <- recode(cell_types,
                            "classical monocyte" = "Monocyte",
                            "non-classical monocyte" = "Monocyte",
                            "naive thymus-derived CD4-positive, alpha-beta T cell" = "CD4+ T cell",
                            "CD4-positive helper T cell" = "CD4+ T cell",
                            "CD4-positive, alpha-beta cytotoxic T cell" = "CD4+ T cell",
                            "naive thymus-derived CD8-positive, alpha-beta T cell" = "CD8+ T cell",
                            "CD8-positive, alpha-beta cytotoxic T cell" = "CD8+ T cell",
                            "central memory CD8-positive, alpha-beta T cell" = "CD8+ T cell",
                            "effector memory CD8-positive, alpha-beta T cell, terminally differentiated" = "CD8+ T cell",
                            "effector memory CD8-positive, alpha-beta T cell" = "CD8+ T cell",
                            "gamma-delta T cell" = "Other T cell",
                            "mucosal invariant T cell" = "Other T cell",
                            "regulatory T cell" = "Other T cell",
                            "naive B cell" = "B cell",
                            "class switched memory B cell" = "B cell",
                            "B cell" = "B cell",
                            "plasmablast" = "B cell",
                            "plasma cell" = "B cell",
                            "mature NK T cell" = "NK cell",
                            "natural killer cell" = "NK cell",
                            "CD16-negative, CD56-bright natural killer cell, human" = "NK cell",
                            "plasmacytoid dendritic cell" = "Dendritic cell",
                            "conventional dendritic cell" = "Dendritic cell",
                            "dendritic cell" = "Dendritic cell",
                            "hematopoietic precursor cell" = "Other",
                            "innate lymphoid cell" = "Other",
                            "unknown" = "Other",
                            "platelet" = "Other",
                            "granulocyte" = "Other",
                            "erythrocyte" = "Other")
SCT_longcovid@meta.data$merged_cell_type <- merged_cell_types

# 分割前准备
merged_cell_types <- SCT_longcovid@meta.data$merged_cell_type
unique_cell_types <- unique(merged_cell_types)
seurat_list <- list()
for (celltype in unique_cell_types) {
  print(paste("Processing cell type:", celltype))
  seurat_subset <- subset(SCT_longcovid, subset = merged_cell_type == celltype)
  print(paste("Number of cells after subsetting for", celltype, ":", length(seurat_subset$merged_cell_type)))
  print(table(seurat_subset$merged_cell_type))
  seurat_list[[celltype]] <- seurat_subset
}

# 按照sample_id 及  seurat_list中的merged_cell_type 拆分数据   
output_dir <- "/data/wuwq/noise/longcovid/rds_files"  # 修改为你想存储结果的文件夹路径
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# 打印 seurat_list 的长度和基本信息，确保数据被正确加载
print(paste("Length of seurat_list:", length(seurat_list)))  # 打印seurat_list的长度
if (length(seurat_list) == 0) {
  print("Error: seurat_list is empty. Please check the data loading process.")
  stop("seurat_list is empty.")  # 停止执行，如果seurat_list为空
} else {
  print("seurat_list is loaded correctly.")
}
# 外层循环：遍历seurat_list中的每一个细胞类型
for (cell_type in names(seurat_list)) {
  # 获取当前细胞类型的Seurat对象
  seurat_object <- seurat_list[[cell_type]]
  # 检查seurat_object是否有效
  if (is.null(seurat_object)) {
    print(paste("Warning: Seurat object for cell type", cell_type, "is NULL"))
    next
  }
  # 获取counts数据和scRNASeq_sample_ID
  counts_data <- seurat_object@assays$SCT@counts
  sample_ids <- seurat_object@meta.data$sample_id
  # 检查counts_data是否为空
  if (is.null(counts_data) || ncol(counts_data) == 0) {
    print(paste("Warning: Counts data is empty for cell type", cell_type))
    next
  }
  # 检查sample_ids是否为空
  if (length(sample_ids) == 0) {
    print(paste("Warning: No sample IDs found for cell type", cell_type))
    next
  }
  # 获取所有唯一的样本ID
  unique_sample_ids <- unique(sample_ids)
  # 创建当前细胞类型的子文件夹来存储分割后的数据
  cell_type_dir <- file.path(output_dir, cell_type)
  if (!dir.exists(cell_type_dir)) {
    dir.create(cell_type_dir)
  }
  # 内层循环：根据样本ID分割数据并保存
  for (sample_id in unique_sample_ids) {
    # 打印调试信息
    print(paste("Processing sample:", sample_id, "for cell type:", cell_type))
    # 获取当前样本ID对应的有效列
    valid_columns <- sample_ids == sample_id
    print(paste("Valid columns:", sum(valid_columns)))  # 打印有效列的数量
    # 打印 valid_columns 的具体情况，确保其正确性
    print(paste("Sample ID:", sample_id))
    print(paste("Number of valid columns:", sum(valid_columns)))
    # 查看 valid_columns 对应的列索引
    valid_column_indices <- which(valid_columns)
    print("Valid column indices:")
    print(valid_column_indices)
    # 如果没有有效列，则跳过当前样本
    if (sum(valid_columns) == 0) {
      print(paste("Warning: No valid data for sample", sample_id, "- Skipping"))
      next
    }
    # 提取数据并保持矩阵结构
    sample_data <- as.matrix(counts_data[, valid_columns, drop = FALSE])
    # 打印 sample_data 的内容，检查数据是否正确
    print("Sample data content:")
    print(sample_data)
    # 检查 sample_data 是否为 NULL 或没有有效数据
    if (is.null(sample_data) || ncol(sample_data) == 0 || all(sample_data == 0)) {
      # 如果样本数据全为零，则跳过该样本的处理
      print(paste("Warning: Sample", sample_id, "has all zero values. Skipping this sample."))
    } else if (any(is.na(sample_data))) {
      # 如果样本数据中包含 NA 值，则跳过该样本的处理
      print(paste("Warning: Sample", sample_id, "contains NA values. Skipping this sample."))
    } else {
      # 如果样本数据有效，处理并保存为 RDS 文件
      print(paste("Valid data for sample", sample_id))
      # 保存 sample_data 为 .rds 文件
      output_file <- file.path(cell_type_dir, paste0(sample_id, "_", cell_type, "_processed.rds"))
      saveRDS(sample_data, file = output_file)
      print(paste("Saved:", output_file))
    }
  }
}



### SLE #3107414   nohup Rscript /data2/wuwq/noise/code/nohup_cutrds_SLE.R > cutrds_SLE.log 2>&1 &   
rm(list = setdiff(ls(), c("SCT_SLE")))
SCT_SLE <- readRDS("/data2/wuwq/noise/data_SLE_RA/SCT_SLE.rds")
cell_types <- SCT_SLE@meta.data$cell_type
unique(cell_types)
merged_cell_types <- recode(cell_types,
                            "classical monocyte" = "Monocyte",
                            "non-classical monocyte" = "Monocyte",
                            "CD4-positive, alpha-beta T cell" = "CD4+ T cell",
                            "CD8-positive, alpha-beta T cell" = "CD8+ T cell",
                            "plasmablast" = "B cell",
                            "B cell" = "B cell",
                            "natural killer cell" = "NK cell",
                            "plasmacytoid dendritic cell" = "Dendritic cell",
                            "conventional dendritic cell" = "Dendritic cell",
                            "lymphocyte" = "Other",
                            "progenitor cell" = "Other"
)
# 将合并后的细胞类型存储到新的列中
SCT_SLE@meta.data$merged_cell_type <- merged_cell_types

# 分割准备
merged_cell_types <- SCT_SLE@meta.data$merged_cell_type
unique_cell_types <- unique(merged_cell_types)
seurat_list <- list()
for (celltype in unique_cell_types) {
  print(paste("Processing cell type:", celltype))
  seurat_subset <- subset(SCT_SLE, subset = merged_cell_type == celltype)
  print(paste("Number of cells after subsetting for", celltype, ":", length(seurat_subset$merged_cell_type)))
  print(table(seurat_subset$merged_cell_type))
  seurat_list[[celltype]] <- seurat_subset
}

# 按照sample_uuid 及  seurat_list中的merged_cell_type 拆分数据   
output_dir <- "/data/wuwq/noise/SLE/rds_files"  # 修改为你想存储结果的文件夹路径
# 如果输出文件夹不存在，则创建文件夹
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
# 打印 seurat_list 的长度和基本信息，确保数据被正确加载
print(paste("Length of seurat_list:", length(seurat_list)))  # 打印seurat_list的长度
if (length(seurat_list) == 0) {
  print("Error: seurat_list is empty. Please check the data loading process.")
  stop("seurat_list is empty.")  # 停止执行，如果seurat_list为空
} else {
  print("seurat_list is loaded correctly.")
}
# 外层循环：遍历seurat_list中的每一个细胞类型
for (cell_type in names(seurat_list)) {
  # 获取当前细胞类型的Seurat对象
  seurat_object <- seurat_list[[cell_type]]
  # 检查seurat_object是否有效
  if (is.null(seurat_object)) {
    print(paste("Warning: Seurat object for cell type", cell_type, "is NULL"))
    next
  }
  # 获取counts数据和scRNASeq_sample_ID
  counts_data <- seurat_object@assays$SCT@counts
  sample_ids <- seurat_object@meta.data$sample_uuid
  # 检查counts_data是否为空
  if (is.null(counts_data) || ncol(counts_data) == 0) {
    print(paste("Warning: Counts data is empty for cell type", cell_type))
    next
  }
  # 检查sample_ids是否为空
  if (length(sample_ids) == 0) {
    print(paste("Warning: No sample IDs found for cell type", cell_type))
    next
  }
  # 获取所有唯一的样本ID
  unique_sample_ids <- unique(sample_ids)
  # 创建当前细胞类型的子文件夹来存储分割后的数据
  cell_type_dir <- file.path(output_dir, cell_type)
  if (!dir.exists(cell_type_dir)) {
    dir.create(cell_type_dir)
  }
  # 内层循环：根据样本ID分割数据并保存
  for (sample_id in unique_sample_ids) {
    # 打印调试信息
    print(paste("Processing sample:", sample_id, "for cell type:", cell_type))
    # 获取当前样本ID对应的有效列
    valid_columns <- sample_ids == sample_id
    print(paste("Valid columns:", sum(valid_columns)))  # 打印有效列的数量
    # 打印 valid_columns 的具体情况，确保其正确性
    print(paste("Sample ID:", sample_id))
    print(paste("Number of valid columns:", sum(valid_columns)))
    # 查看 valid_columns 对应的列索引
    valid_column_indices <- which(valid_columns)
    print("Valid column indices:")
    print(valid_column_indices)
    # 如果没有有效列，则跳过当前样本
    if (sum(valid_columns) == 0) {
      print(paste("Warning: No valid data for sample", sample_id, "- Skipping"))
      next
    }
    # 提取数据并保持矩阵结构
    sample_data <- as.matrix(counts_data[, valid_columns, drop = FALSE])
    # 打印 sample_data 的内容，检查数据是否正确
    print("Sample data content:")
    print(sample_data)
    # 检查 sample_data 是否为 NULL 或没有有效数据
    if (is.null(sample_data) || ncol(sample_data) == 0 || all(sample_data == 0)) {
      # 如果样本数据全为零，则跳过该样本的处理
      print(paste("Warning: Sample", sample_id, "has all zero values. Skipping this sample."))
    } else if (any(is.na(sample_data))) {
      # 如果样本数据中包含 NA 值，则跳过该样本的处理
      print(paste("Warning: Sample", sample_id, "contains NA values. Skipping this sample."))
    } else {
      # 如果样本数据有效，处理并保存为 RDS 文件
      print(paste("Valid data for sample", sample_id))
      # 保存 sample_data 为 .rds 文件
      output_file <- file.path(cell_type_dir, paste0(sample_id, "_", cell_type, "_processed.rds"))
      saveRDS(sample_data, file = output_file)
      print(paste("Saved:", output_file))
    }
  }
}












### RA
rm(list = setdiff(ls(), c("SCT_RA")))
SCT_RA <- readRDS("/data2/wuwq/noise/data_SLE_RA/SCT_RA.rds")
cell_types <- SCT_RA@meta.data$cell_type
merged_cell_types <- recode(cell_types,
                            "classical monocyte" = "Monocyte",
                            "non-classical monocyte" = "Monocyte",
                            "naive thymus-derived CD4-positive, alpha-beta T cell" = "CD4+ T cell",
                            "central memory CD4-positive, alpha-beta T cell" = "CD4+ T cell",
                            "effector memory CD4-positive, alpha-beta T cell" = "CD4+ T cell",
                            "CD4-positive, alpha-beta T cell" = "CD4+ T cell",
                            "naive thymus-derived CD8-positive, alpha-beta T cell" = "CD8+ T cell",
                            "CD8-positive, alpha-beta memory T cell" = "CD8+ T cell",
                            "effector memory CD8-positive, alpha-beta T cell, terminally differentiated" = "CD8+ T cell",
                            "gamma-delta T cell" = "Other T cell",
                            "naive B cell" = "B cell",
                            "memory B cell" = "B cell",
                            "plasmablast" = "B cell",
                            "natural killer cell" = "NK cell",
                            "myeloid dendritic cell" = "Dendritic cell"
)
# 将合并后的细胞类型存储到新的列中
SCT_RA@meta.data$merged_cell_type <- merged_cell_types

# 分割前准备
merged_cell_types <- SCT_RA@meta.data$merged_cell_type
unique_cell_types <- unique(merged_cell_types)
seurat_list <- list()
for (celltype in unique_cell_types) {
  print(paste("Processing cell type:", celltype))
  seurat_subset <- subset(SCT_RA, subset = merged_cell_type == celltype)
  print(paste("Number of cells after subsetting for", celltype, ":", length(seurat_subset$merged_cell_type)))
  print(table(seurat_subset$merged_cell_type))
  seurat_list[[celltype]] <- seurat_subset
}

# 按照donor_id 及  seurat_list中的merged_cell_type 拆分数据 
output_dir <- "/data/wuwq/noise/RA/rds_files"  # 修改为你想存储结果的文件夹路径
# 如果输出文件夹不存在，则创建文件夹
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# 打印 seurat_list 的长度和基本信息，确保数据被正确加载
print(paste("Length of seurat_list:", length(seurat_list)))  # 打印seurat_list的长度
if (length(seurat_list) == 0) {
  print("Error: seurat_list is empty. Please check the data loading process.")
  stop("seurat_list is empty.")  # 停止执行，如果seurat_list为空
} else {
  print("seurat_list is loaded correctly.")
}
# 外层循环：遍历seurat_list中的每一个细胞类型
for (cell_type in names(seurat_list)) {
  # 获取当前细胞类型的Seurat对象
  seurat_object <- seurat_list[[cell_type]]
  
  # 检查seurat_object是否有效
  if (is.null(seurat_object)) {
    print(paste("Warning: Seurat object for cell type", cell_type, "is NULL"))
    next
  }
  # 获取counts数据和scRNASeq_sample_ID
  counts_data <- seurat_object@assays$SCT@counts
  sample_ids <- seurat_object@meta.data$donor_id
  # 检查counts_data是否为空
  if (is.null(counts_data) || ncol(counts_data) == 0) {
    print(paste("Warning: Counts data is empty for cell type", cell_type))
    next
  }
  # 检查sample_ids是否为空
  if (length(sample_ids) == 0) {
    print(paste("Warning: No sample IDs found for cell type", cell_type))
    next
  }
  # 获取所有唯一的样本ID
  unique_sample_ids <- unique(sample_ids)
  # 创建当前细胞类型的子文件夹来存储分割后的数据
  cell_type_dir <- file.path(output_dir, cell_type)
  if (!dir.exists(cell_type_dir)) {
    dir.create(cell_type_dir)
  }
  # 内层循环：根据样本ID分割数据并保存
  for (sample_id in unique_sample_ids) {
    # 打印调试信息
    print(paste("Processing sample:", sample_id, "for cell type:", cell_type))
    # 获取当前样本ID对应的有效列
    valid_columns <- sample_ids == sample_id
    print(paste("Valid columns:", sum(valid_columns)))  # 打印有效列的数量
    # 打印 valid_columns 的具体情况，确保其正确性
    print(paste("Sample ID:", sample_id))
    print(paste("Number of valid columns:", sum(valid_columns)))
    # 查看 valid_columns 对应的列索引
    valid_column_indices <- which(valid_columns)
    print("Valid column indices:")
    print(valid_column_indices)
    # 如果没有有效列，则跳过当前样本
    if (sum(valid_columns) == 0) {
      print(paste("Warning: No valid data for sample", sample_id, "- Skipping"))
      next
    }
    # 提取数据并保持矩阵结构
    sample_data <- as.matrix(counts_data[, valid_columns, drop = FALSE])
    # 打印 sample_data 的内容，检查数据是否正确
    print("Sample data content:")
    print(sample_data)
    # 检查 sample_data 是否为 NULL 或没有有效数据
    if (is.null(sample_data) || ncol(sample_data) == 0 || all(sample_data == 0)) {
      # 如果样本数据全为零，则跳过该样本的处理
      print(paste("Warning: Sample", sample_id, "has all zero values. Skipping this sample."))
    } else if (any(is.na(sample_data))) {
      # 如果样本数据中包含 NA 值，则跳过该样本的处理
      print(paste("Warning: Sample", sample_id, "contains NA values. Skipping this sample."))
    } else {
      # 如果样本数据有效，处理并保存为 RDS 文件
      print(paste("Valid data for sample", sample_id))
      # 保存 sample_data 为 .rds 文件
      output_file <- file.path(cell_type_dir, paste0(sample_id, "_", cell_type, "_processed.rds"))
      saveRDS(sample_data, file = output_file)
      print(paste("Saved:", output_file))
    }
  }
}






### NSCLC
# 细胞类型注释复现 → 细胞类型合并
# 1.细胞类型注释复现
SCT_NSCLC <- readRDS("/data2/wuwq/noise/data_NSCLC/SCT_NSCLC.rds")
DefaultAssay(SCT_NSCLC) <- "RNA" 
Art_SCT_NSCLC = NormalizeData(SCT_NSCLC)
Art_SCT_NSCLC <- CellCycleScoring(object = Art_SCT_NSCLC, g2m.features = cc.genes$g2m.genes, 
                                  s.features = cc.genes$s.genes)
Art_SCT_NSCLC <- FindVariableFeatures(Art_SCT_NSCLC, selection.method = "vst", assay = "RNA",
                                      nfeatures = 2000, verbose = FALSE)
Art_SCT_NSCLC <- ScaleData(Art_SCT_NSCLC, vars.to.regress = c("percent.mt", "nCount_RNA","nFeature_RNA"), 
                           assay = "RNA")
DefaultAssay(Art_SCT_NSCLC)
Art_SCT_NSCLC <- RunPCA(Art_SCT_NSCLC, npcs = 50, verbose = F)
Art_SCT_NSCLC <- RunTSNE(Art_SCT_NSCLC, reduction = "pca", dims = 1:30, 
                         perplexity = 30, max_iter = 1000, 
                         theta = 0.5, eta = 200, num_threads = 0)
Art_SCT_NSCLC <- RunUMAP(Art_SCT_NSCLC, reduction = "pca", dims = 1:30, 
                         n.components = 2, n.neighbors = 30, 
                         n.epochs = 200, min.dist = 0.3, 
                         learning.rate = 1, spread = 1)
# 保存降维后的结果
saveRDS(Art_SCT_NSCLC, "/data2/wuwq/noise/data_NSCLC/Art_SCT_NSCLC.rds")
# 使用harmony方法进行批次矫正多样本数据整合
library(harmony)
Art_SCT_NSCLC.harmony <- RunHarmony(Art_SCT_NSCLC, group.by.vars = "orig.ident", reduction.use = "pca", 
                                    dims.use = 1:50, assay = "RNA")
Art_SCT_NSCLC[["harmony"]] <- Art_SCT_NSCLC.harmony[["harmony"]]
# 使用 Harmony 结果中的 PCA 进行聚类
Art_SCT_NSCLC <- RunPCA(Art_SCT_NSCLC, reduction = "harmony", reduction.name = "pca_harmony")
# UMAP
Art_SCT_NSCLC <- RunUMAP(Art_SCT_NSCLC, dims = 1:50, reduction = "harmony", reduction.name = "umap_harmony")
DimPlot(Art_SCT_NSCLC, reduction = "umap_harmony", group.by = "orig.ident") + 
  ggtitle("UMAP Harmony")
# t-SNE
Art_SCT_NSCLC <- RunTSNE(Art_SCT_NSCLC, reduction = "harmony", dims = 1:50, reduction.name = "Tsne_harmony",
                         perplexity = 30, max_iter = 1000, 
                         theta = 0.5, eta = 200, num_threads = 0)
DimPlot(Art_SCT_NSCLC, reduction = "Tsne_harmony", group.by = "orig.ident") + 
  ggtitle("Tsne Harmony")
DefaultAssay(Art_SCT_NSCLC)
# 整合方法结果的可视化
p1 <- DimPlot(Art_SCT_NSCLC, reduction = "umap", group.by = "orig.ident") + ggtitle("UMAP raw_data")
p2 <- DimPlot(Art_SCT_NSCLC, reduction = "umap_harmony", group.by = "orig.ident") + 
  ggtitle("UMAP Harmony")
leg <- get_legend(p1)
gridExtra::grid.arrange(gridExtra::arrangeGrob(p1 + NoLegend() + NoAxes(), p2 + NoLegend() + 
                                                 NoAxes(), nrow = 2), 
                        leg, ncol = 2, widths = c(8, 2))
# 保存数据harmony整合后的结果
saveRDS(Art_SCT_NSCLC, "/data2/wuwq/noise/data_NSCLC/Art_SCT_NSCLC_harmony.rds")
# 读取harmony整合后的结果
Art_SCT_NSCLC_harmony <- readRDS("/data2/wuwq/noise/data_NSCLC/Art_SCT_NSCLC_harmony.rds")
# 计算邻接矩阵（基于 Harmony PCA 结果） # not_21  yes_30个PC #Computing nearest neighbor graph, Computing SNN 
Art_SCT_NSCLC_harmony <- FindNeighbors(Art_SCT_NSCLC_harmony, reduction = "pca_harmony", dims = 1:30, k.param = 60, prune.SNN = 1/15)
names(Art_SCT_NSCLC_harmony@graphs) #"RNA_nn"  "RNA_snn"
pheatmap(Art_SCT_NSCLC_harmony@graphs$RNA_nn[1:200, 1:200], 
         col = c("white", "black"), border_color = "grey90", 
         legend = F, cluster_rows = F, cluster_cols = F, fontsize = 2)
# 进行聚类和可视化（默认使用 Louvain 聚类算法）
for (res in c(0.1, 0.2, 0.4, 1, 1.5, 2)) {
  Art_SCT_NSCLC_harmony <- FindClusters(Art_SCT_NSCLC_harmony, graph.name = "RNA_snn", resolution = res, algorithm = 1)
}
for (res in c(0.15, 0.25, 0.3, 0.35, 0.45, 0.5)) {
  Art_SCT_NSCLC_harmony <- FindClusters(Art_SCT_NSCLC_harmony, graph.name = "RNA_snn", resolution = res, algorithm = 1)
}
Art_SCT_NSCLC_harmony <- RunUMAP(Art_SCT_NSCLC_harmony, dims = 1:50, reduction = "harmony", reduction.name = "umap_harmony")
plot_grid(ncol = 3,
          DimPlot(Art_SCT_NSCLC_harmony, reduction = "umap_harmony", group.by = "RNA_snn_res.0.1", raster = FALSE) + ggtitle("louvain_0.1"),
          DimPlot(Art_SCT_NSCLC_harmony, reduction = "umap_harmony", group.by = "RNA_snn_res.0.15", raster = FALSE) + ggtitle("louvain_0.15"),
          DimPlot(Art_SCT_NSCLC_harmony, reduction = "umap_harmony", group.by = "RNA_snn_res.0.2", raster = FALSE) + ggtitle("louvain_0.2"),
          DimPlot(Art_SCT_NSCLC_harmony, reduction = "umap_harmony", group.by = "RNA_snn_res.0.25", raster = FALSE) + ggtitle("louvain_0.25"),
          DimPlot(Art_SCT_NSCLC_harmony, reduction = "umap_harmony", group.by = "RNA_snn_res.0.3", raster = FALSE) + ggtitle("louvain_0.3"),
          DimPlot(Art_SCT_NSCLC_harmony, reduction = "umap_harmony", group.by = "RNA_snn_res.0.35", raster = FALSE) + ggtitle("louvain_0.35"),
          DimPlot(Art_SCT_NSCLC_harmony, reduction = "umap_harmony", group.by = "RNA_snn_res.0.4", raster = FALSE) + ggtitle("louvain_0.4"),
          DimPlot(Art_SCT_NSCLC_harmony, reduction = "umap_harmony", group.by = "RNA_snn_res.0.45", raster = FALSE) + ggtitle("louvain_0.45"),
          DimPlot(Art_SCT_NSCLC_harmony, reduction = "umap_harmony", group.by = "RNA_snn_res.0.5", raster = FALSE) + ggtitle("louvain_0.5"))
suppressPackageStartupMessages(library(clustree))
clustree(Art_SCT_NSCLC_harmony@meta.data, prefix = "RNA_snn_res.")
# 保存聚类后的结果
saveRDS(Art_SCT_NSCLC_harmony, "/data2/wuwq/noise/data_NSCLC/Art_SCT_NSCLC_har_clu.rds")
suppressPackageStartupMessages({
  library(Seurat)
  library(venn)
  library(dplyr)
  library(cowplot)
  library(ggplot2)
  library(pheatmap)
  library(rafalib)
  library(scPred)
})
unique(Art_SCT_NSCLC_har_clu@meta.data$SubType)
# Set the identity as louvain with resolution 0.4
Art_SCT_NSCLC_har_clu <- SetIdent(Art_SCT_NSCLC_har_clu, value = "RNA_snn_res.0.4")
table(Art_SCT_NSCLC_har_clu@active.ident)
plot_grid(ncol = 3, DimPlot(Art_SCT_NSCLC_har_clu,reduction = "umap_harmony", label = T) + NoAxes(), DimPlot(Art_SCT_NSCLC_har_clu,reduction = "umap_harmony", group.by = "orig.ident") + 
            NoAxes(), DimPlot(Art_SCT_NSCLC_har_clu,reduction = "umap_harmony", group.by = "SubType") + NoAxes())
DefaultAssay(Art_SCT_NSCLC_har_clu)
# 尝试3种方法，最终选了1.2.1.gsea_nes_cell_types的注释结果
# 1.1.scPred包
Art_SCT_NSCLC_har_clu@active.assay = "RNA"
# 提取scPred包中PBMC参考数据集
reference <- scPred::pbmc_1
reference
reference <- reference %>% NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>% 
  RunPCA(verbose = F) %>% RunUMAP(dims = 1:30)
DimPlot(reference, group.by = "cell_type", label = TRUE, repel = TRUE) + NoAxes()
# 对自己数据进行处理
Art_SCT_NSCLC_har_clu <- SetIdent(Art_SCT_NSCLC_har_clu, value = "RNA_snn_res.0.4")
DimPlot(Art_SCT_NSCLC_har_clu, reduction = "umap_harmony", label = TRUE, repel = TRUE) + NoAxes()
# 使用FindTransferAnchors函数寻找query和reference数据集之间的anchors
transfer.anchors <- FindTransferAnchors(reference = reference, query = Art_SCT_NSCLC_har_clu, dims = 1:30)
# 使用TransferData函数进行标签转移的细胞类型预测
predictions <- TransferData(anchorset = transfer.anchors, refdata = reference$cell_type, 
                            dims = 1:30)
# 将细胞类型预测的结果添加到metadata$predicted.id
Art_SCT_NSCLC_har_clu <- AddMetaData(object = Art_SCT_NSCLC_har_clu, metadata = predictions)
DimPlot(Art_SCT_NSCLC_har_clu, group.by = "predicted.id", reduction = "umap_harmony", label = T, repel = T) + NoAxes()
#DimPlot(Art_SCT_NSCLC_har_clu, group.by = "predicted.id", label = T, repel = T) + NoAxes()
# 细胞数可视化
ggplot(Art_SCT_NSCLC_har_clu@meta.data, aes(x = RNA_snn_res.0.4, fill = predicted.id)) + geom_bar() + theme_classic()
# 保存细胞类型注释后的数据集
saveRDS(Art_SCT_NSCLC_har_clu, "/data2/wuwq/noise/data_NSCLC/Final_SCT_NSCLC.rds")

### 1.2.gsea富集
Final_SCT_NSCLC <- readRDS("/data2/wuwq/noise/data_NSCLC/Final_SCT_NSCLC.rds")
Final_SCT_NSCLC <- SetIdent(Final_SCT_NSCLC, value = "RNA_snn_res.0.4")
DGE_table <- FindAllMarkers(Final_SCT_NSCLC, logfc.threshold = 0, test.use = "wilcox", min.pct = 0.1, 
                            min.diff.pct = 0, only.pos = TRUE, max.cells.per.ident = 20, return.thresh = 1, 
                            assay = "RNA")
DGE_list <- split(DGE_table, DGE_table$cluster)
# 读入cellmarker参考基因集列表，并做初步的筛选
markers <- read_xlsx("/data2/wuwq/noise/data_NSCLC/CellMarker_list/Cell_marker_Human.xlsx")
markers <- markers[markers$tissue_class == "Blood", ]
markers <- markers[markers$cancer_type == "Non-small Cell Lung Cancer"| markers$cancer_type == "Normal" , ]
markers <- markers[markers$tissue_type == "Blood" | markers$tissue_type == "Peripheral blood" , ]
celltype_list <- lapply(unique(markers$cell_name), function(x) {
  x <- paste(markers$Symbol[markers$cell_name == x], sep = ",")
  x <- gsub("[[]|[]]| |-", ",", x)
  x <- unlist(strsplit(x, split = ","))
  x <- unique(x[!x %in% c("", "NA", "family")])
  x <- casefold(x, upper = T)
})
names(celltype_list) <- unique(markers$cell_name)
celltype_list <- lapply(celltype_list , function(x) {x[1:min(length(x),50)]} )
celltype_list <- celltype_list[unlist(lapply(celltype_list, length)) < 100]
celltype_list <- celltype_list[unlist(lapply(celltype_list, length)) > 5]
#对已知的基因集进行GSEA富集分析
res <- lapply(DGE_list, function(x) {
  gene_rank <- setNames(x$avg_log2FC, x$gene)
  fgseaRes <- fgsea(pathways = celltype_list, stats = gene_rank, nperm = 10000)
  return(fgseaRes)
})
names(res) <- names(DGE_list)
# You can filter and resort the table based on ES, NES or pvalue
# 筛出p<0.01和size>5的行
res <- lapply(res, function(x) {
  x[x$pval < 0.01, ]
})
res <- lapply(res, function(x) {
  x[x$size > 5, ]
})
# 1.2.1.细胞类型对应规则gsea_nes  
# 按照NES降序排列 √
res <- lapply(res, function(x) {
  x[order(x$NES, decreasing = T), ]
})
# show top 3 for each cluster.
lapply(res, head, 3)
cluster_17 <- Final_SCT_NSCLC@meta.data$RNA_snn_res.0.4
unique(cluster_17)
gsea_nes_cell_types <- recode(cluster_17,
                              "0" = "CD1C+_B dendritic cell",
                              "1" = "Natural killer cell",
                              "2" = "CD8+ T cell",
                              "3" = "Naive CD4+ T cell",
                              "4" = "Naive CD4+ T cell",
                              "5" = "CD4+ T cell",
                              "6" = "B cell",
                              "7" = "CD16+ monocyte",
                              "8" = "B cell",
                              "9" = "Regulatory T(Treg) cell",
                              "10" = "CD1C+_B dendritic cell",
                              "11" = "CD1C+_A dendritic cell",
                              "12" = "Hematopoietic stem cell",
                              "13" = "unknown",
                              "14" = "Myeloid dendritic cell",
                              "15" = "Erythroid cell",
                              "16" = "Megakaryocyte"
)
# 将合并后的细胞类型存储到新的列中
Final_SCT_NSCLC@meta.data$gsea_nes_cell_types <- gsea_nes_cell_types

# 1.2.2.细胞类型对应规则gsea_nes
# 按照p升序排列
res <- lapply(res, function(x) {
  x[order(x$pval), ]
})
# show top 3 for each cluster.
lapply(res, head, 3)
cluster_17 <- Final_SCT_NSCLC@meta.data$RNA_snn_res.0.4
unique(cluster_17)
gsea_p_cell_types <- recode(cluster_17,
                            "0" = "CD1C+_B dendritic cell",
                            "1" = "Effector CD8+ memory T (Tem) cell", #
                            "2" = "Effector CD8+ memory T (Tem) cell", #
                            "3" = "Naive CD8+ T cell", #
                            "4" = "Naive CD4+ T cell",
                            "5" = "CD4+ T cell",
                            "6" = "B cell",
                            "7" = "Myeloid cell", #
                            "8" = "B cell",
                            "9" = "Regulatory T(Treg) cell",
                            "10" = "CD1C+_B dendritic cell",
                            "11" = "CD1C+_A dendritic cell",
                            "12" = "CD1C+_B dendritic cell", #
                            "13" = "unknown",
                            "14" = "Myeloid dendritic cell",
                            "15" = "Erythroid cell",
                            "16" = "Megakaryocyte"
)
# 将合并后的细胞类型存储到新的列中
Final_SCT_NSCLC@meta.data$gsea_p_cell_types <- gsea_p_cell_types

# 2.可视化 #predicted.id  #gsea_nes_cell_types  #gsea_p_cell_types
#基于harmony_umap的比较
cowplot::plot_grid(ncol = 2, DimPlot(Final_SCT_NSCLC, label = T, group.by = "predicted.id",reduction = "umap_harmony") + 
                     NoAxes(), DimPlot(Final_SCT_NSCLC, label = T, group.by = "gsea_nes_cell_types",reduction = "umap_harmony") + NoAxes(),
                   DimPlot(Final_SCT_NSCLC, label = T, group.by = "gsea_p_cell_types",reduction = "umap_harmony"))
#基于原文umap的比较
umap_data <- Final_SCT_NSCLC@meta.data
p1 <- ggplot(umap_data, aes(x = umap1, y = umap2, color = predicted.id)) +
  geom_point(alpha = 0.7) +
  labs(title = "Predicted ID") +
  theme_minimal() +
  theme(legend.position = "right")
print(p1)
p2 <- ggplot(umap_data, aes(x = umap1, y = umap2, color = gsea_nes_cell_types)) +
  geom_point(alpha = 0.7) +
  labs(title = "GSEA NES Cell Types") +
  theme_minimal() +
  theme(legend.position = "right")
print(p2)
p3 <- ggplot(umap_data, aes(x = umap1, y = umap2, color = gsea_p_cell_types)) +
  geom_point(alpha = 0.7) +
  labs(title = "GSEA P Cell Types") +
  theme_minimal() +
  theme(legend.position = "right")
print(p3)
# 使用 cowplot 排列多个图形
plot_grid(p1, p2, p3, ncol = 2)


rm(list=ls())
# 3.细胞类型合并 # 最终选了1.2.1.gsea_nes_cell_types的注释结果
Final_SCT_NSCLC <- readRDS("/data2/wuwq/noise/data_NSCLC/Final_SCT_NSCLC.rds")
cell_types <- Final_SCT_NSCLC@meta.data$gsea_nes_cell_types
unique(cell_types)
# 定义合并规则
merged_cell_types_gsea_nes <- recode(cell_types,
                                     "CD16+ monocyte" = "Monocyte",
                                     "Naive CD4+ T cell" = "CD4+ T cell",
                                     "CD4+ T cell" = "CD4+ T cell",
                                     "CD8+ T cell" = "CD8+ T cell",
                                     "Regulatory T(Treg) cell" = "Other T cell",
                                     "B cell" = "B cell",
                                     "Natural killer cell" = "NK cell",
                                     "CD1C+_A dendritic cell" = "Dendritic cell",
                                     "CD1C+_B dendritic cell" = "Dendritic cell",
                                     "Myeloid dendritic cell" = "Dendritic cell",
                                     "unknown" = "Other",
                                     "Erythroid cell" = "Other",
                                     "Megakaryocyte" = "Other",
                                     "Hematopoietic stem cell" = "Other")
# 将合并后的细胞类型存储到新的列中
Final_SCT_NSCLC@meta.data$merged_cell_types_gsea_nes <- merged_cell_types_gsea_nes
unique_cell_types <- unique(merged_cell_types_gsea_nes)

# 分割准备
seurat_list <- list()
for (celltype in unique_cell_types) {
  print(paste("Processing cell type:", celltype))
  seurat_subset <- subset(Final_SCT_NSCLC, subset = merged_cell_types_gsea_nes == celltype)
  print(paste("Number of cells after subsetting for", celltype, ":", length(seurat_subset$merged_cell_types_gsea_nes)))
  print(table(seurat_subset$merged_cell_types_gsea_nes))
  seurat_list[[celltype]] <- seurat_subset
}

# 按照orig.ident 及  seurat_list中的新celltype 拆分数据
output_dir <- "/data/wuwq/noise/NSCLC/rds_files"  # 修改为你想存储结果的文件夹路径
# 如果输出文件夹不存在，则创建文件夹
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# 打印 seurat_list 的长度和基本信息，确保数据被正确加载
print(paste("Length of seurat_list:", length(seurat_list)))  # 打印seurat_list的长度
if (length(seurat_list) == 0) {
  print("Error: seurat_list is empty. Please check the data loading process.")
  stop("seurat_list is empty.")  # 停止执行，如果seurat_list为空
} else {
  print("seurat_list is loaded correctly.")
}

# 外层循环：遍历seurat_list中的每一个细胞类型
for (cell_type in names(seurat_list)) {
  # 获取当前细胞类型的Seurat对象
  seurat_object <- seurat_list[[cell_type]]
  # 检查seurat_object是否有效
  if (is.null(seurat_object)) {
    print(paste("Warning: Seurat object for cell type", cell_type, "is NULL"))
    next
  }
  # 获取counts数据和scRNASeq_sample_ID
  counts_data <- seurat_object@assays$SCT@counts
  sample_ids <- seurat_object@meta.data$orig.ident
  # 检查counts_data是否为空
  if (is.null(counts_data) || ncol(counts_data) == 0) {
    print(paste("Warning: Counts data is empty for cell type", cell_type))
    next
  }
  # 检查sample_ids是否为空
  if (length(sample_ids) == 0) {
    print(paste("Warning: No sample IDs found for cell type", cell_type))
    next
  }
  # 获取所有唯一的样本ID
  unique_sample_ids <- unique(sample_ids)
  # 创建当前细胞类型的子文件夹来存储分割后的数据
  cell_type_dir <- file.path(output_dir, cell_type)
  if (!dir.exists(cell_type_dir)) {
    dir.create(cell_type_dir)
  }
  # 内层循环：根据样本ID分割数据并保存
  for (sample_id in unique_sample_ids) {
    # 打印调试信息
    print(paste("Processing sample:", sample_id, "for cell type:", cell_type))
    # 获取当前样本ID对应的有效列
    valid_columns <- sample_ids == sample_id
    print(paste("Valid columns:", sum(valid_columns)))  # 打印有效列的数量
    # 打印 valid_columns 的具体情况，确保其正确性
    print(paste("Sample ID:", sample_id))
    print(paste("Number of valid columns:", sum(valid_columns)))
    # 查看 valid_columns 对应的列索引
    valid_column_indices <- which(valid_columns)
    print("Valid column indices:")
    print(valid_column_indices)
    # 如果没有有效列，则跳过当前样本
    if (sum(valid_columns) == 0) {
      print(paste("Warning: No valid data for sample", sample_id, "- Skipping"))
      next
    }
    # 提取数据并保持矩阵结构
    sample_data <- as.matrix(counts_data[, valid_columns, drop = FALSE])
    # 打印 sample_data 的内容，检查数据是否正确
    print("Sample data content:")
    print(sample_data)
    # 检查 sample_data 是否为 NULL 或没有有效数据
    if (is.null(sample_data) || ncol(sample_data) == 0 || all(sample_data == 0)) {
      # 如果样本数据全为零，则跳过该样本的处理
      print(paste("Warning: Sample", sample_id, "has all zero values. Skipping this sample."))
    } else if (any(is.na(sample_data))) {
      # 如果样本数据中包含 NA 值，则跳过该样本的处理
      print(paste("Warning: Sample", sample_id, "contains NA values. Skipping this sample."))
    } else {
      # 如果样本数据有效，处理并保存为 RDS 文件
      print(paste("Valid data for sample", sample_id))
      # 保存 sample_data 为 .rds 文件
      output_file <- file.path(cell_type_dir, paste0(sample_id, "_", cell_type, "_processed.rds"))
      saveRDS(sample_data, file = output_file)
      print(paste("Saved:", output_file))
    }
  }
}


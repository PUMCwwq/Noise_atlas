########################################################### 整合noise矫正的结果
### CD4+ T cell
### CD8+ T cell
### Monocyte
### NK cell
### Dendritic cell
### B cell
### Other T cell


##### (一)、COMBAT
rm(list=ls())
### 1.读取 对齐后大类基线信息、log2fitratio矫正表
baseline <- read.delim("/data/wuwq/noise/baseline/baseline_COMBAT_2.txt", header = TRUE, row.names = 1,sep = "\t")
row.names(baseline) <- baseline$sample_id
log2fitratio <- read.csv("/data/wuwq/noise/COMBAT/CD4/CD4_fixed_baseline_gene_wide_table.csv", header = TRUE, row.names = 1)
#indi_noise <- read.csv("/data/wuwq/noise/COMBAT/CD4/CD4_fixed_baseline_corrected_noise.csv", header = TRUE, row.names = 1)
noise <- as.data.frame(log2fitratio)
noise <- noise[ , -c(1:5,7:10)]
noise[] <- lapply(noise, as.numeric)
# 合并gene_noise和indi_noise
#noise <- merge(noise, 
#               indi_noise["corrected_noise"], 
#               by = "row.names", 
#               all.x = TRUE)
#noise <- noise %>% relocate(corrected_noise, .before = 1)
#row.names(noise) <- noise$Row.names
#noise <- noise[,-2]
#write.csv(noise, "/data/wuwq/noise/COMBAT/cd4_noise.csv", row.names = TRUE)
#write.csv(baseline, "/data/wuwq/noise/COMBAT/cd4_baseline.csv", row.names = TRUE)

# 合并基线和noise表
common_rows <- intersect(rownames(baseline), rownames(noise))
noise_ind <- data.frame(row.names = common_rows)
noise_ind <- cbind(noise_ind, baseline[common_rows, , drop = FALSE],  noise[common_rows, , drop = FALSE] )
write.csv(noise_ind, "/data/wuwq/noise/COMBAT/CD4/COMBAT_2.csv", row.names = TRUE)

rm(list=ls())
### 2.读取 数据集详细基线信息、log2fitratio表
baseline <- read.delim("/data/wuwq/noise/baseline/baseline_COMBAT.txt", header = TRUE, row.names = 1,sep = "\t")
row.names(baseline) <- baseline$scRNASeq_sample_ID
log2fitratio <- read.csv("/data/wuwq/noise/COMBAT/CD4/CD4_fixed_baseline_gene_wide_table.csv", header = TRUE, row.names = 1)
noise <- as.data.frame(log2fitratio)
noise <- noise[ , -c(1:5,7:10)]
noise[] <- lapply(noise, as.numeric)
# 合并基线和noise表
common_rows <- intersect(rownames(baseline), rownames(noise))
noise_ind <- data.frame(row.names = common_rows)
noise_ind <- cbind(noise_ind, baseline[common_rows, , drop = FALSE],  noise[common_rows, , drop = FALSE] )
write.csv(noise_ind, "/data/wuwq/noise/COMBAT/CD4/COMBAT.csv", row.names = TRUE)


##### (二)、FLU
rm(list=ls())
### 1.读取 对齐后大类基线信息、log2fitratio矫正表
baseline <- read.delim("/data/wuwq/noise/baseline/baseline_FLU_2.txt", header = TRUE, row.names = 1,sep = "\t")
row.names(baseline) <- baseline$sample_id
log2fitratio <- read.csv("/data/wuwq/noise/FLU/CD4/CD4_fixed_baseline_gene_wide_table.csv", header = TRUE, row.names = 1)
noise <- as.data.frame(log2fitratio)
noise <- noise[ , -c(1:5,7:10)]
noise[] <- lapply(noise, as.numeric)
# 合并基线和noise表
common_rows <- intersect(rownames(baseline), rownames(noise))
noise_ind <- data.frame(row.names = common_rows)
noise_ind <- cbind(noise_ind, baseline[common_rows, , drop = FALSE],  noise[common_rows, , drop = FALSE] )
write.csv(noise_ind, "/data/wuwq/noise/FLU/CD4/FLU_2.csv", row.names = TRUE)

rm(list=ls())
### 2.读取 数据集详细基线信息、log2fitratio表
baseline <- read.delim("/data/wuwq/noise/baseline/baseline_FLU.txt", header = TRUE, row.names = 1,sep = "\t")
row.names(baseline) <- baseline$Sample.ID
log2fitratio <- read.csv("/data/wuwq/noise/FLU/CD4/CD4_fixed_baseline_gene_wide_table.csv", header = TRUE, row.names = 1)
noise <- as.data.frame(log2fitratio)
noise <- noise[ , -c(1:5,7:10)]
noise[] <- lapply(noise, as.numeric)
# 合并基线和noise表
common_rows <- intersect(rownames(baseline), rownames(noise))
noise_ind <- data.frame(row.names = common_rows)
noise_ind <- cbind(noise_ind, baseline[common_rows, , drop = FALSE],  noise[common_rows, , drop = FALSE] )
write.csv(noise_ind, "/data/wuwq/noise/FLU/CD4/FLU.csv", row.names = TRUE)



##### (三)、SLE
rm(list=ls())
### 1.读取 对齐后大类基线信息、log2fitratio矫正表
baseline <- read.delim("/data/wuwq/noise/baseline/baseline_SLE_2.txt", header = TRUE, row.names = 1,sep = "\t")
row.names(baseline) <- baseline$sample_id
log2fitratio <- read.csv("/data/wuwq/noise/SLE/CD4/CD4_fixed_baseline_gene_wide_table.csv", header = TRUE, row.names = 1)
noise <- as.data.frame(log2fitratio)
noise <- noise[ , -c(1:5,7:10)]
noise[] <- lapply(noise, as.numeric)
# 合并基线和noise表
common_rows <- intersect(rownames(baseline), rownames(noise))
noise_ind <- data.frame(row.names = common_rows)
noise_ind <- cbind(noise_ind, baseline[common_rows, , drop = FALSE],  noise[common_rows, , drop = FALSE] )
write.csv(noise_ind, "/data/wuwq/noise/SLE/CD4/SLE_2.csv", row.names = TRUE)

rm(list=ls())
### 2.读取 数据集详细基线信息、log2fitratio表
baseline <- read.delim("/data/wuwq/noise/baseline/baseline_SLE.txt", header = TRUE, row.names = 1,sep = "\t")
row.names(baseline) <- baseline$sample_id
log2fitratio <- read.csv("/data/wuwq/noise/SLE/CD4/CD4_fixed_baseline_gene_wide_table.csv", header = TRUE, row.names = 1)
noise <- as.data.frame(log2fitratio)
noise <- noise[ , -c(1:5,7:10)]
noise[] <- lapply(noise, as.numeric)
# 合并基线和noise表
common_rows <- intersect(rownames(baseline), rownames(noise))
noise_ind <- data.frame(row.names = common_rows)
noise_ind <- cbind(noise_ind, baseline[common_rows, , drop = FALSE],  noise[common_rows, , drop = FALSE] )
write.csv(noise_ind, "/data/wuwq/noise/SLE/CD4/SLE.csv", row.names = TRUE)



##### (四)、RA
rm(list=ls())
### 1.读取 对齐后大类基线信息、log2fitratio矫正表
baseline <- read.delim("/data/wuwq/noise/baseline/baseline_RA_2.txt", header = TRUE, row.names = 1,sep = "\t")
row.names(baseline) <- baseline$sample_id
log2fitratio <- read.csv("/data/wuwq/noise/RA/CD4/CD4_fixed_baseline_gene_wide_table.csv", header = TRUE, row.names = 1)
noise <- as.data.frame(log2fitratio)
noise <- noise[ , -c(1:5,7:10)]
noise[] <- lapply(noise, as.numeric)
# 合并基线和noise表
common_rows <- intersect(rownames(baseline), rownames(noise))
noise_ind <- data.frame(row.names = common_rows)
noise_ind <- cbind(noise_ind, baseline[common_rows, , drop = FALSE],  noise[common_rows, , drop = FALSE] )
write.csv(noise_ind, "/data/wuwq/noise/RA/CD4/RA_2.csv", row.names = TRUE)

rm(list=ls())
### 2.读取 数据集详细基线信息、log2fitratio表
baseline <- read.delim("/data/wuwq/noise/baseline/baseline_RA.txt", header = TRUE, row.names = 1,sep = "\t")
row.names(baseline) <- baseline$donor_id
log2fitratio <- read.csv("/data/wuwq/noise/RA/CD4/CD4_fixed_baseline_gene_wide_table.csv", header = TRUE, row.names = 1)
noise <- as.data.frame(log2fitratio)
noise <- noise[ , -c(1:5,7:10)]
noise[] <- lapply(noise, as.numeric)
# 合并基线和noise表
common_rows <- intersect(rownames(baseline), rownames(noise))
noise_ind <- data.frame(row.names = common_rows)
noise_ind <- cbind(noise_ind, baseline[common_rows, , drop = FALSE],  noise[common_rows, , drop = FALSE] )
write.csv(noise_ind, "/data/wuwq/noise/RA/CD4/RA.csv", row.names = TRUE)



##### (五)、NSCLC
rm(list=ls())
### 1.读取 对齐后大类基线信息、log2fitratio矫正表
baseline <- read.delim("/data/wuwq/noise/baseline/baseline_NSCLC_2.txt", header = TRUE, row.names = 1,sep = "\t")
row.names(baseline) <- baseline$sample_id
log2fitratio <- read.csv("/data/wuwq/noise/NSCLC/CD4/CD4_fixed_baseline_gene_wide_table.csv", header = TRUE, row.names = 1)
noise <- as.data.frame(log2fitratio)
noise <- noise[ , -c(1:5,7:10)]
noise[] <- lapply(noise, as.numeric)
# 合并基线和noise表
common_rows <- intersect(rownames(baseline), rownames(noise))
noise_ind <- data.frame(row.names = common_rows)
noise_ind <- cbind(noise_ind, baseline[common_rows, , drop = FALSE],  noise[common_rows, , drop = FALSE] )
write.csv(noise_ind, "/data/wuwq/noise/NSCLC/CD4/NSCLC_2.csv", row.names = TRUE)

rm(list=ls())
### 2.读取 数据集详细基线信息、log2fitratio表
baseline <- read.delim("/data/wuwq/noise/baseline/baseline_NSCLC.txt", header = TRUE, row.names = 1,sep = "\t")
row.names(baseline) <- baseline$orig.ident.1
log2fitratio <- read.csv("/data/wuwq/noise/NSCLC/CD4/CD4_fixed_baseline_gene_wide_table.csv", header = TRUE, row.names = 1)
noise <- as.data.frame(log2fitratio)
noise <- noise[ , -c(1:5,7:10)]
noise[] <- lapply(noise, as.numeric)
# 合并基线和noise表
common_rows <- intersect(rownames(baseline), rownames(noise))
noise_ind <- data.frame(row.names = common_rows)
noise_ind <- cbind(noise_ind, baseline[common_rows, , drop = FALSE],  noise[common_rows, , drop = FALSE] )
write.csv(noise_ind, "/data/wuwq/noise/NSCLC/CD4/NSCLC.csv", row.names = TRUE)



##### (六)、longcovid
rm(list=ls())
### 1.读取 对齐后大类基线信息、log2fitratio矫正表
baseline <- read.delim("/data/wuwq/noise/baseline/baseline_longcovid_2.txt", header = TRUE, row.names = 1,sep = "\t")
row.names(baseline) <- baseline$sample_id
log2fitratio <- read.csv("/data/wuwq/noise/longcovid/CD4/CD4_fixed_baseline_gene_wide_table.csv", header = TRUE, row.names = 1)
noise <- as.data.frame(log2fitratio)
noise <- noise[ , -c(1:5,7:10)]
noise[] <- lapply(noise, as.numeric)
# 合并基线和noise表
common_rows <- intersect(rownames(baseline), rownames(noise))
noise_ind <- data.frame(row.names = common_rows)
noise_ind <- cbind(noise_ind, baseline[common_rows, , drop = FALSE],  noise[common_rows, , drop = FALSE] )
write.csv(noise_ind, "/data/wuwq/noise/longcovid/CD4/longcovid_2.csv", row.names = TRUE)

rm(list=ls())
### 2.读取 数据集详细基线信息、log2fitratio表
baseline <- read.delim("/data/wuwq/noise/baseline/baseline_longcovid.txt", header = TRUE, row.names = 1,sep = "\t")
row.names(baseline) <- baseline$sample_id
log2fitratio <- read.csv("/data/wuwq/noise/longcovid/CD4/CD4_fixed_baseline_gene_wide_table.csv", header = TRUE, row.names = 1)
noise <- as.data.frame(log2fitratio)
noise <- noise[ , -c(1:5,7:10)]
noise[] <- lapply(noise, as.numeric)
# 合并基线和noise表
common_rows <- intersect(rownames(baseline), rownames(noise))
noise_ind <- data.frame(row.names = common_rows)
noise_ind <- cbind(noise_ind, baseline[common_rows, , drop = FALSE],  noise[common_rows, , drop = FALSE] )
write.csv(noise_ind, "/data/wuwq/noise/longcovid/CD4/longcovid.csv", row.names = TRUE)






rm(list=ls())
############################################### 健康人增加gene_symbol列并保存，### 也需特殊保存
noise <- read.csv("/data/wuwq/noise/HEALTH_REAL/final/othertmerge.csv",row.names = 1)
noise <- noise %>%
  # 删除正数第4-5列
  select(-c(4:5)) %>%
  # 删除倒数1-4列（最后4列）
  select(-tail(names(.), 4))
noise <- as.data.frame(t(noise))
##### 增加匹配的gene_symbol列
annotation_df <-read.delim("/data/wuwq/noise/mapping_geneid2symbol.txt",header = T)
# sym_noise_1
sym_noise <- merge(noise, 
                     annotation_df, 
                     by.x = 0,               # noise_x 中的 行名
                     by.y = "Gene_ID",       # 注释表中的 gene_id 列
                     all.x = TRUE)           # 保留所有 gene，即使没有匹配
sym_noise <- sym_noise %>% relocate(Gene_name, .before = 1)
row.names(sym_noise)<-sym_noise$Row.names
sym_noise <- sym_noise[,-2]
sym_noise <- as.data.frame(t(sym_noise))
col_order <- c("group", "sex", "age","corrected_noise", setdiff(names(sym_noise), c("group", "sex", "age","corrected_noise")))
sym_noise <- sym_noise[ , col_order]
sym_noise[1, ] <- ifelse(is.na(sym_noise[1, ]), colnames(sym_noise), sym_noise[1, ])
write.table(sym_noise,file = "/data/wuwq/noise/HEALTH_REAL/sym_otherT_health_noise.csv", sep = ",", row.names = TRUE, col.names = NA, quote = TRUE)






rm(list=ls())
############################################### 疾病增加gene_symbol列并保存
###################################### 逐个将noise表的基因id，增加gene_symbol列并保存
####### (1). 7种细胞类型 & 5个数据集(NSCLC处理见下)
noise <- read.csv("/data/wuwq/noise/COMBAT/CD4/COMBAT_2.csv", row.names = 1)
row.names(noise) <- noise$sample_id
noise <- as.data.frame(t(noise))
##### 增加匹配的gene_symbol列
annotation_df <-read.delim("/data/wuwq/noise/mapping_geneid2symbol.txt",header = T)
# sym_noise_1
sym_noise <- merge(noise, 
                   annotation_df, 
                   by.x = 0,               # noise_x 中的 行名
                   by.y = "Gene_ID",       # 注释表中的 gene_id 列
                   all.x = TRUE)           # 保留所有 gene，即使没有匹配
sym_noise <- sym_noise %>% relocate(Gene_name, .before = 1)
row.names(sym_noise)<-sym_noise$Row.names
sym_noise <- sym_noise[,-2]
sym_noise <- as.data.frame(t(sym_noise))
col_order <- c("corrected_noise", "dataset","sample_id","Disease_group", "sex", setdiff(names(sym_noise), c("corrected_noise", "dataset","sample_id","Disease_group", "sex")))
sym_noise <- sym_noise[ , col_order]
sym_noise[1, ] <- ifelse(is.na(sym_noise[1, ]), colnames(sym_noise), sym_noise[1, ])
write.csv(sym_noise, "/data/wuwq/noise/COMBAT/CD4/sym_COMBAT_2.csv", row.names = TRUE)



rm(list=ls())
####### (2).NSCLC是symbol，单独处理
noise <- read.csv("/data/wuwq/noise/NSCLC/CD4/NSCLC_2.csv", row.names = 1) 
row.names(noise) <- noise$sample_id
noise <- as.data.frame(t(noise))
noise <- cbind(Gene_name = rownames(noise), noise)
noise[, 1] <- gsub("\\.", "-", as.character(noise[, 1]))
rownames(noise) <- as.character(noise[, 1])
sym_noise <- as.data.frame(t(noise))
col_order <- c("corrected_noise", "dataset","sample_id","Disease_group", "sex", setdiff(names(sym_noise), c("corrected_noise", "dataset","sample_id","Disease_group", "sex")))
sym_noise <- sym_noise[ , col_order]
write.csv(sym_noise, "/data/wuwq/noise/NSCLC/CD4/sym_NSCLC_2.csv", row.names = TRUE)







####################################### 全部疾病整合表--symbol设为列名后--需特殊保存
rm(list=ls())
##### （一）、individual_noise、gene_noise按照疾病大类合并处理
# 读取每个疾病的sum_noise列--- 统一关键列名的_2.csv文件
COMBAT <- read.csv("/data/wuwq/noise/COMBAT/otherT/sym_COMBAT_2.csv", header = TRUE,row.names = 1)
FLU <- read.csv("/data/wuwq/noise/FLU/otherT/sym_FLU_2.csv", header = TRUE,row.names = 1)
SLE <- read.csv("/data/wuwq/noise/SLE/otherT/sym_SLE_2.csv", header = TRUE,row.names = 1)
RA <- read.csv("/data/wuwq/noise/RA/otherT/sym_RA_2.csv", header = TRUE,row.names = 1)
longcovid <- read.csv("/data/wuwq/noise/longcovid/otherT/sym_longcovid_2.csv", header = TRUE,row.names = 1)
NSCLC <- read.csv("/data/wuwq/noise/NSCLC/otherT/sym_NSCLC_2.csv",header = TRUE,row.names = 1)
colnames(COMBAT) <- as.character(COMBAT[1, ])
colnames(FLU) <- as.character(FLU[1, ])
colnames(SLE) <- as.character(SLE[1, ])
colnames(RA) <- as.character(RA[1, ])
colnames(longcovid) <- as.character(longcovid[1, ])
colnames(NSCLC) <- as.character(NSCLC[1, ])
COMBAT <- COMBAT[-1, , drop = FALSE]
FLU <- FLU[-1, , drop = FALSE]
SLE <- SLE[-1, , drop = FALSE]
RA <- RA[-1, , drop = FALSE]
longcovid <- longcovid[-1, , drop = FALSE]
NSCLC <- NSCLC[-1, , drop = FALSE]
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
# 并集合并 → 主
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
merged <- merged %>% select(corrected_noise, dataset, sample_id, Disease_group, sex,   everything())
merged[, 6:ncol(merged)] <- lapply(merged[, 6:ncol(merged)], function(x) {
  x <- as.numeric(x)
  return(x)
})
write.table(merged,file = "/data/wuwq/noise/DISEASE_ALL/union_otherT.csv", sep = ",", row.names = TRUE, col.names = NA, quote = TRUE)

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
merged <- merged %>% select(corrected_noise, dataset, sample_id, Disease_group, sex,   everything())
merged[, 6:ncol(merged)] <- lapply(merged[, 6:ncol(merged)], function(x) {
  x <- as.numeric(x)
  return(x)
})
write.table(merged,file = "/data/wuwq/noise/DISEASE_ALL/intersect_otherT.csv", sep = ",", row.names = TRUE, col.names = NA, quote = TRUE)







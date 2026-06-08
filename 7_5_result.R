##################################################################
############ ------------ result_5_全疾病 -----------#############
############ ------       基于个体noise的聚类     -------#########
############ --------         无监督聚类       -------############ 
##################################################################
library(dplyr)
library(knitr)
library(tidyr)
library(stats)
library(tibble)
library(ggplot2)
library(cowplot)
library(caret)
library(cluster)



rm(list=ls())
### merge细胞类型   特殊读取 ✅️
noise_CD4 <- read.table("/data/wuwq/noise/DISEASE_ALL/union_CD4.csv", sep = ",", row.names=1, header = TRUE, check.names = FALSE, quote = "\"")
noise_CD8 <- read.table("/data/wuwq/noise/DISEASE_ALL/union_CD8.csv", sep = ",", row.names=1, header = TRUE, check.names = FALSE, quote = "\"")
noise_MNP <- read.table("/data/wuwq/noise/DISEASE_ALL/union_MNP.csv", sep = ",", row.names=1, header = TRUE, check.names = FALSE, quote = "\"")
noise_NK <- read.table("/data/wuwq/noise/DISEASE_ALL/union_NK.csv", sep = ",", row.names=1, header = TRUE, check.names = FALSE, quote = "\"")
noise_B <- read.table("/data/wuwq/noise/DISEASE_ALL/union_B.csv", sep = ",", row.names=1, header = TRUE, check.names = FALSE, quote = "\"")
noise_otherT <- read.table("/data/wuwq/noise/DISEASE_ALL/union_otherT.csv", sep = ",", row.names=1, header = TRUE, check.names = FALSE, quote = "\"")


# 去掉2-5列自带基线
noise_CD4 <- noise_CD4[,-c(2:5)]
noise_CD8 <- noise_CD8[,-c(2:5)]
noise_MNP <- noise_MNP[,-c(2:5)]
noise_NK <- noise_NK[,-c(2:5)]
noise_B <- noise_B[,-c(2:5)]
noise_otherT <- noise_otherT[,-c(2:5)]
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

# 更新总体结局_2
baseline <- read.delim("/data/wuwq/noise/DISEASE_ALL/app_disease/baseline_merge_disease.txt", header = TRUE, row.names = 1)
rownames(baseline) <- baseline$sample_id
unique(baseline$dataset)
baseline <- baseline[,-1]


common_rows <- intersect(rownames(baseline), rownames(noise_all))
noise_d <- data.frame(row.names = common_rows)
noise_d <- cbind(noise_d, baseline[common_rows, , drop = FALSE],  noise_all[common_rows, , drop = FALSE] )
# 去除全部为NA的基因列
noise_d <- Filter(function(x) !all(is.na(x)), noise_d)
unique(noise_d$Disease_group) # "CR"    "AR"    "other" "PriR"
unique(noise_d$sex) # "female" "male"   NA"
# 把corrected_noise列前移
# 找出以 "corrrected_noise_" 开头的列
corrected_noise_cols <- grep("^corrected_noise_", names(noise_d), value = TRUE)
# 找出非 "corrrected_noise_" 的列
other_cols <- setdiff(names(noise_d), corrected_noise_cols)
# 保留前3列的基线列
baseline_cols <- other_cols[1:3]
# 剩余的非 "corrrected_noise_" 列
remaining_cols <- other_cols[4:length(other_cols)]
# 重新排列列的顺序
new_order <- c(baseline_cols, corrected_noise_cols, remaining_cols)
# 重新排列数据框列
noise_all_reordered <- noise_d %>%
  select(all_of(new_order))
write.table(noise_all_reordered,file = "/data/wuwq/noise/DISEASE_ALL/app_disease/merge_all_disease.csv", sep = ",", row.names = TRUE, col.names = NA, quote = TRUE)



rm(list=ls())
noise_all_reordered <- read.table("/data/wuwq/noise/DISEASE_ALL/app_disease/merge_all_disease.csv", sep = ",", row.names=1, header = TRUE, check.names = FALSE, quote = "\"")
unique(noise_all_reordered$Disease_group)
#all_info <- colnames(noise_all_reordered[,c(2,10:ncol(noise_all_reordered))])
all_info <- colnames(noise_all_reordered[,c(2,4:9)])
#all_info <- colnames(noise_all_reordered[,c(2,4:ncol(noise_all_reordered))])
noise_ind <- noise_all_reordered
# 准备用于PCA的数据
valid_cols <- intersect(all_info, colnames(noise_ind))
pca_data <- noise_ind[, valid_cols, drop = FALSE]
#pca_data <- pca_data[, c("sample_id", setdiff(names(pca_data), "sample_id"))]
#pca_data <- pca_data[, c("Disease_group", setdiff(names(pca_data), "Disease_group"))]
# 将NA替换为0
pca_data[is.na(pca_data)] <- 0
pca_data <- pca_data[, colSums(pca_data != 0) > 0]
### 保留占比大于**%的基因列
#non_zero_prop <- colSums(pca_data != 0) / nrow(pca_data)
# 设置非零比例的门槛为**% 并过滤列
#data <- pca_data[, non_zero_prop > 0.8]



# scale. = TRUE表示分析前对数据进行归一化；
com1 <- prcomp(data[ ,2:ncol(data)], center = TRUE,scale. = TRUE)
summ<-summary(com1)
df1<-com1$x
head(df1)
df1<-data.frame(df1,data$Disease_group)
df1$data.Disease_group<-as.factor(df1$data.Disease_group)
df1 <- df1[, c("data.Disease_group", setdiff(names(df1), "data.Disease_group"))]
lab1<-paste0("PC1(",round(summ$importance[2,1]*100,2),"%)")
lab2<-paste0("PC2(",round(summ$importance[2,2]*100,2),"%)")

# 创建PCA图形
Fig1a.taxa.pca <- ggplot(df1, aes(PC1, PC2)) +
  geom_point(size = 2, aes(color = data.Disease_group), show.legend = F) +
  scale_color_manual(values = c("#5686C3", "#75C500", "pink", "orange","red","brown","grey")) +
  stat_ellipse(aes(color = data.Disease_group), fill = "white", geom = "polygon",
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
  geom_density(aes(x = PC1, group = data.Disease_group, fill = data.Disease_group),
               color = "black", alpha = 0.6, position = 'identity') +
  scale_color_manual(values = c("#5686C3", "#75C500", "pink", "orange","red","brown","grey")) +
  theme_bw() +
  labs(fill = "")+
  theme(legend.position = "none")  # 隐藏图例
Fig1a.taxa.pc1.density
# 创建PC2密度图
Fig1a.taxa.pc2.density <- ggplot(df1) +
  geom_density(aes(x = PC2, group = data.Disease_group, fill = data.Disease_group),
               color = "black", alpha = 0.6, position = 'identity') +
  scale_color_manual(values = c("#5686C3", "#75C500", "pink", "orange","red","brown","grey")) +
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
write.csv(wss_df,file="/data/wuwq/noise/DISEASE_ALL/app_disease/all_disease/k3_wss.csv")
# 可视化WSS下降率
wss_rate_df <- data.frame(k = 2:10, WSS_Rate = wss_rate)
ggplot(wss_rate_df, aes(x = k, y = WSS_Rate)) +
  geom_col() +  # 使用 geom_col() 绘制柱形图
  ggtitle("WSS下降率 vs. k") +
  xlab("Number of Clusters (k)") +
  ylab("WSS下降率") +
  scale_x_continuous(breaks = seq(2, 10, by = 1))  # 设置x轴间隔为1
write.csv(wss_rate_df,file="/data/wuwq/noise/DISEASE_ALL/app_disease/all_disease/k3_wss_rate.csv")
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
write.csv(fpc_df,file="/data/wuwq/noise/DISEASE_ALL/app_disease/all_disease/k3_fpc.csv")
# 可视化 FPC 变化率
fpc_rate_df <- data.frame(k = 2:10, FPC_Rate = fpc_rate[1:9])
ggplot(fpc_rate_df, aes(x = k, y = FPC_Rate)) +
  geom_col() +
  geom_point()
write.csv(fpc_rate_df,file="/data/wuwq/noise/DISEASE_ALL/app_disease/all_disease/k3_fpc_rate.csv")





# 设置随机种子，让结果可以重现
set.seed(1)
# 调用kmeans聚类算法 k = 3
km <- kmeans(df2, centers = 4, nstart = 25)
km$totss
### 把聚类cluster附加到PCA矩阵中：
final_data <- cbind(true_group = data$Disease_group, kgroup = km$cluster,PC1 = df1$PC1, PC2 = df1$PC2, PC3 = df1$PC3)
write.csv(final_data,file="/data/wuwq/noise/DISEASE_ALL/app_disease/all_disease/k3.csv")


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




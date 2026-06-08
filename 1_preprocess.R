##################################################
library(Seurat)
library(dplyr)
library(Matrix)
rm(list=ls())
gc()
##################################################
##### COMBAT, FLU, longcovid, NSCLC, SLE, RA


##### (一) 从matrix到seurat对象
### COMBAT
matrix_COMBAT <-readMM("/data2/wuwq/noise/data_COVID19/rds/raw_adata_COMBAT.mtx")
meta<- read.csv("/data2/wuwq/noise/data_COVID19/rds/cellinfo_COMBAT.csv",header = TRUE, row.names = 1)
geneinfo <- read.csv("/data2/wuwq/noise/data_COVID19/rds/geneinfo_COMBAT.csv",header = TRUE, row.names = 1)
mtx <- as.matrix(matrix_COMBAT)
##列是细胞 行为基因
mtx <- t(mtx)
colnames(mtx) <- row.names(meta)
rownames(mtx) <- row.names(geneinfo)
# 核对细胞和基因名信息是否一致
identical(rownames(meta),colnames(mtx)) 
# 拼装seurat对象
sce.all=CreateSeuratObject(counts = mtx ,
                           meta.data = meta)
saveRDS(sce.all, "/data2/wuwq/noise/data_COVID19/rds/S_COMBAT.rds")

### FLU
matrix_flu <-readMM("/data2/wuwq/noise/data_COVID19/round7_flu/raw_adata_flu.mtx")
meta<- read.csv("/data2/wuwq/noise/data_COVID19/round7_flu/cellinfo_flu.csv",header = TRUE, row.names = 1)
geneinfo <- read.csv("/data2/wuwq/noise/data_COVID19/round7_flu/geneinfo_flu.csv",header = TRUE, row.names = 1)
mtx <- as.matrix(matrix_flu)
##列是细胞 行为基因
mtx <- t(mtx)
colnames(mtx) <- row.names(meta)
rownames(mtx) <- row.names(geneinfo)
# 核对细胞和基因名信息是否一致
identical(rownames(meta),colnames(mtx)) 
# 拼装seurat对象
sce.all=CreateSeuratObject(counts = mtx ,
                           meta.data = meta)
saveRDS(sce.all, "/data2/wuwq/noise/data_COVID19/round7_flu/data_flu.rds")

### longcovid
matrix_long_covid <-readMM("/data2/wuwq/noise/data_COVID19/round6_long_covid/raw_adata_longcovid.mtx")
meta<- read.csv("/data2/wuwq/noise/data_COVID19/round6_long_covid/cellinfo_longcovid.csv",header = TRUE, row.names = 1)
geneinfo <- read.csv("/data2/wuwq/noise/data_COVID19/round6_long_covid/geneinfo_longcovid.csv",header = TRUE, row.names = 1)
mtx <- as.matrix(matrix_long_covid)
mtx <- t(mtx)
colnames(mtx) <- row.names(meta)
rownames(mtx) <- row.names(geneinfo)
identical(rownames(meta),colnames(mtx)) 
sce.all=CreateSeuratObject(counts = mtx ,
                           meta.data = meta)
saveRDS(sce.all, "/data2/wuwq/noise/data_COVID19/round6_long_covid/data_longcovid.rds")

### NSCLC
# 指定数据集链接
url <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE285nnn/GSE285888/suppl/GSE285888_RAW.tar"  # 替换为实际链接
# 指定服务器上的保存路径
save_path <- "/data2/wuwq/noise/data_NSCLC/GSE285888_RAW.tar"  # 替换为实际路径
# 下载数据集到指定路径
download.file(url, destfile = save_path)  # 使用curl方法下载
# 解压文件
tar_file <- "/data2/wuwq/noise/data_NSCLC/GSE285888_RAW.tar"
untar(tar_file, exdir = "/data2/wuwq/noise/data_NSCLC")  
# 指定压缩文件的路径
gz_file <- "/data2/wuwq/noise/data_NSCLC/GSM8712033_matrix.txt.gz"
# 解压文件到当前工作目录
gunzip(gz_file, overwrite = TRUE)  # 如果文件已存在，设置 overwrite = TRUE 覆盖
matrix_NSCLC <-read.table("/data2/wuwq/noise/data_NSCLC/GSM8712033_matrix.txt",header = TRUE,row.names = 1,sep = "\t")
meta<- read.csv("/data2/wuwq/noise/data_NSCLC/GSM8712033_metadata.csv",header = TRUE, row.names = 1)
mtx <- as.matrix(matrix_NSCLC)
identical(rownames(meta),colnames(mtx)) 
sce.all=CreateSeuratObject(counts = mtx ,
                           meta.data = cl)
saveRDS(sce.all, "/data2/wuwq/noise/data_NSCLC/data_NSCLC.rds")

### SLE
matrix_SLE <-readMM("/data2/wuwq/noise/data_SLE_RA/raw_adata_SLE.mtx")
meta<- read.csv("/data2/wuwq/noise/data_SLE_RA/SLE/cellinfo_SLE.csv",header = TRUE, row.names = 1)
geneinfo <- read.csv("/data2/wuwq/noise/data_SLE_RA/SLE/geneinfo_SLE.csv",header = TRUE, row.names = 1)
mtx <- as.matrix(matrix_SLE)
mtx <- t(mtx)
colnames(mtx) <- row.names(meta)
rownames(mtx) <- row.names(geneinfo)
identical(rownames(meta),colnames(mtx)) 
sce.all=CreateSeuratObject(counts = mtx ,
                           meta.data = meta)
saveRDS(sce.all, "/data2/wuwq/noise/data_SLE_RA/SLE/data_SLE.rds")

### RA
matrix_RA <-readMM("/data2/wuwq/noise/data_RA/raw_adata_RA.mtx")
meta<- read.csv("/data2/wuwq/noise/data_RA/cellinfo_RA.csv",header = TRUE, row.names = 1)
geneinfo <- read.csv("/data2/wuwq/noise/data_RA/geneinfo_RA.csv",header = TRUE, row.names = 1)
mtx <- as.matrix(mtx)
mtx <- t(mtx)
colnames(mtx) <- row.names(meta)
rownames(mtx) <- row.names(geneinfo)
identical(rownames(meta),colnames(mtx)) 
sce.all=CreateSeuratObject(counts = mtx ,
                           meta.data = meta)
saveRDS(sce.all, "/data2/wuwq/noise/data_RA/data_RA.rds")


##### (二) 质控及SCTransform
### COMBAT
S_COMBAT <- readRDS("/data2/wuwq/noise/data_COVID19/rds/S_COMBAT.rds")
S_COMBAT[["percent.mt"]] <- PercentageFeatureSet(S_COMBAT, pattern = "^MT-")
#根据基因数、线粒体基因比例和 UMI 总数筛选细胞
Q_COMBAT <- subset(S_COMBAT, 
                   subset = nFeature_RNA > 300 &    # 至少检测到300个基因
                     nFeature_RNA < 6000 &           # 最多不超过6000个基因
                     percent.mt < 20 &               # 线粒体基因占比<20%
                     nCount_RNA > 1000               # UMI总数>1000
)
#保留在至少 5 个细胞中表达的基因
Q_COMBAT <- subset(Q_COMBAT, 
                   features = rownames(Q_COMBAT)[
                     Matrix::rowSums(Q_COMBAT@assays$RNA$counts > 0) >= 5]
)
# vars.to.regress = "percent.mt"   
SCT_COMBAT <- SCTransform(Q_COMBAT,vars.to.regress = c("percent.mt"), verbose = TRUE)
saveRDS(SCT_COMBAT, file = "/data2/wuwq/noise/data_COVID19/round5_COMBAT/SCT_COMBAT.rds")

### FLU
data_flu <- readRDS("/data2/wuwq/noise/data_COVID19/round7_flu/data_flu.rds")
data_flu[["percent.mt"]] <- PercentageFeatureSet(data_flu, pattern = "^MT-")
Q_flu <- subset(data_flu, 
                subset = nFeature_RNA > 300 &    # 至少检测到300个基因
                  nFeature_RNA < 6000 &           # 最多不超过6000个基因
                  percent.mt < 20 &               # 线粒体基因占比<20%
                  nCount_RNA > 1000               # UMI总数>1000
)
Q_flu <- subset(Q_flu, 
                features = rownames(Q_flu)[
                  Matrix::rowSums(Q_flu@assays$RNA$counts > 0) >= 5]
)
SCT_flu <- SCTransform(Q_flu,vars.to.regress = c("percent.mt"), verbose = TRUE)
saveRDS(SCT_flu, file = "/data2/wuwq/noise/data_COVID19/round7_flu/SCT_flu.rds")

### longcovid
data_longcovid <- readRDS("/data2/wuwq/noise/data_COVID19/round6_long_covid/data_longcovid.rds")
data_longcovid[["percent.mt"]] <- PercentageFeatureSet(data_longcovid, pattern = "^MT-")
Q_longcovid <- subset(data_longcovid, 
                      subset = nFeature_RNA > 300 &    # 至少检测到300个基因
                        nFeature_RNA < 6000 &           # 最多不超过6000个基因
                        percent.mt < 20 &               # 线粒体基因占比<20%
                        nCount_RNA > 1000               # UMI总数>1000
)
Q_longcovid <- subset(Q_longcovid, 
                      features = rownames(Q_longcovid)[
                        Matrix::rowSums(Q_longcovid@assays$RNA$counts > 0) >= 5]
)
SCT_longcovid <- SCTransform(Q_longcovid,vars.to.regress = c("percent.mt","sequencing_library"), verbose = TRUE)
saveRDS(SCT_longcovid, file = "/data2/wuwq/noise/data_COVID19/round6_long_covid/SCT_longcovid.rds")

### NSCLC
data_NSCLC <- readRDS("/data2/wuwq/noise/data_NSCLC/data_NSCLC.rds")
data_NSCLC[["percent.mt"]] <- PercentageFeatureSet(data_NSCLC, pattern = "^MT-")
Q_NSCLC <- subset(data_NSCLC, 
                  subset = nFeature_RNA > 300 &    # 至少检测到300个基因
                    nFeature_RNA < 6000 &           # 最多不超过6000个基因
                    percent.mt < 20 &               # 线粒体基因占比<20%
                    nCount_RNA > 1000               # UMI总数>1000
)
Q_NSCLC <- subset(Q_NSCLC, 
                  features = rownames(Q_NSCLC)[
                    Matrix::rowSums(Q_NSCLC@assays$RNA$counts > 0) >= 5
                  ]
)
SCT_NSCLC <- SCTransform(Q_NSCLC,vars.to.regress = c("percent.mt"), verbose = TRUE)
saveRDS(SCT_NSCLC, file = "/data2/wuwq/noise/data_NSCLC/SCT_NSCLC.rds")

### SLE
data_SLE <- readRDS("/data2/wuwq/noise/data_SLE_RA/SLE/data_SLE.rds")
data_SLE[["percent.mt"]] <- PercentageFeatureSet(data_SLE, pattern = "^MT-")
Q_SLE <- subset(data_SLE, 
                subset = nFeature_RNA > 300 &    # 至少检测到300个基因
                  nFeature_RNA < 6000 &           # 最多不超过6000个基因
                  percent.mt < 20 &               # 线粒体基因占比<20%
                  nCount_RNA > 1000               # UMI总数>1000
)
Q_SLE <- subset(Q_SLE, 
                features = rownames(Q_SLE)[
                  Matrix::rowSums(Q_SLE@assays$RNA$counts > 0) >= 5
                ]
)
SCT_SLE <- SCTransform(Q_SLE,vars.to.regress = c("percent.mt"), verbose = TRUE)
saveRDS(SCT_SLE, file = "/data2/wuwq/noise/data_SLE_RA/SCT_SLE.rds")

### RA
data_RA <- readRDS("/data2/wuwq/noise/data_SLE_RA/data_RA.rds")
data_RA[["percent.mt"]] <- PercentageFeatureSet(data_RA, pattern = "^MT-")
Q_RA <- subset(data_RA, 
               subset = nFeature_RNA > 300 &    # 至少检测到300个基因
                 nFeature_RNA < 6000 &           # 最多不超过6000个基因
                 percent.mt < 20 &               # 线粒体基因占比<20%
                 nCount_RNA > 1000               # UMI总数>1000
)
Q_RA <- subset(Q_RA, 
               features = rownames(Q_RA)[
                 Matrix::rowSums(Q_RA@assays$RNA$counts > 0) >= 5
               ]
)
SCT_RA <- SCTransform(Q_RA,vars.to.regress = c("percent.mt", "batch"), verbose = TRUE)
saveRDS(SCT_RA, file = "/data2/wuwq/noise/data_SLE_RA/SCT_RA.rds")






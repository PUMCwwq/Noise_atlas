# ------------------------------------------------------------------------------
library(Seurat)
library(matrixStats)   # 用于rowVars
library(stringr)
library(ggplot2)
library(stats)
# ------------------------------------------------------------------------------
# 批量处理函数：输入RDS文件夹，输出每个donor的分析结果

### CD4+ T cell 
### CD8+ T cell
### Other T cell
### Monocyte
### NK cell
### Dendritic cell
### B cell
rm(list=ls())


batch_process_donor_rds <- function(
    rds_dir,          # 存放所有donor的SCT标准化后RDS文件的目录
    output_root_dir,  # 结果保存的根目录（每个donor会生成独立子目录）
    fdr_threshold = 0.05,  # 筛选条件：FDR阈值
    fitratio_threshold = 1,  # 筛选条件：fitratio阈值
    # 新增质控参数（默认值可调整，不影响原始逻辑）
    min_genes_per_cell = 300,  # 细胞质控：最小表达基因数
    max_genes_per_cell = 6000,  # 细胞质控：最大表达基因数
    min_cells_per_gene = 5      # 基因质控：最小表达细胞数
) {
  # 1. 初始化环境与路径检查
  # 创建输出根目录（不存在则新建）
  if (!dir.exists(output_root_dir)) {
    dir.create(output_root_dir, recursive = TRUE)
    message("✅ 创建输出根目录: ", output_root_dir)
  }
  
  # 获取所有RDS文件路径（仅匹配含"SCT"或"processed"的RDS，避免非目标文件）
  rds_files <- list.files(
    path = rds_dir,  
    full.names = TRUE, 
    recursive = FALSE
  )
  
  # 检查是否存在RDS文件
  if (length(rds_files) == 0) {
    stop("❌ 在目录 ", rds_dir, " 中未找到任何RDS文件，请检查路径或文件命名！")
  }
  
  message("\n📊 共发现 ", length(rds_files), " 个donor的RDS文件，开始批量处理...")
  
  # 2. 循环处理每个donor的RDS文件
  for (rds_path in rds_files) {
    # 提取donor名称（从文件名中获取，需根据实际命名调整，示例："donor_001_seurat_processed.rds" → "donor_001"）
    donor_name <- str_remove(basename(rds_path), "_NK cell.*\\.rds$")  # 正则匹配移除后缀
    donor_output_dir <- file.path(output_root_dir, donor_name)  # 每个donor的独立输出目录
    
    # 跳过已处理的donor（避免重复运行）
    if (dir.exists(donor_output_dir)) {
      message("\n⚠️ donor ", donor_name, " 已存在结果目录，跳过处理！")
      next
    }
    
    # 创建当前donor的输出目录
    dir.create(donor_output_dir, recursive = TRUE)
    message("\n=====================================")
    message("🔍 开始处理 donor: ", donor_name)
    message("📁 RDS文件路径: ", rds_path)
    message("💾 结果保存目录: ", donor_output_dir)
    message("=====================================")
    
    tryCatch({
      # ------------------------------------------------------------------------
      # 步骤1：读取Seurat对象并提取SCT矩阵
      # ------------------------------------------------------------------------
      message("1/8: 读取Seurat对象...")
      seurat_obj <- readRDS(rds_path)
      # 新增：检查细胞数，少于10则跳过
      cell_count <- ncol(seurat_obj)  # Seurat对象的列数即细胞数
      if (cell_count < 10) {
        message(paste0("   ⚠️ donor ", donor_name, " 细胞数仅 ", cell_count, " 个（<10），跳过处理！"))
        next  # 直接进入下一个donor的循环
      }
      
      # 提取SCT counts矩阵（转为普通矩阵，避免稀疏矩阵问题）
      X <- as.matrix(seurat_obj)
      message(paste0("   ✅ 原始SCT矩阵: ", nrow(X), " 基因 × ", ncol(X), " 细胞"))
      
      # ------------------------------------------------------------------------
      # 新增步骤：细胞+基因质控（仅筛选，不修改任何后续计算逻辑）
      # ------------------------------------------------------------------------
      message("2/8: 数据质控（细胞+基因）...")
      # 细胞质控：筛选表达基因数在[min_genes_per_cell, max_genes_per_cell]之间的细胞
      genes_per_cell <- colSums(X > 0, na.rm = TRUE)
      quality_cells <- which(genes_per_cell >= min_genes_per_cell & genes_per_cell <= max_genes_per_cell)
      
      # 基因质控：筛选至少在min_cells_per_gene个细胞中表达的基因
      cells_per_gene <- rowSums(X > 0, na.rm = TRUE)
      quality_genes <- which(cells_per_gene >= min_cells_per_gene)
      
      # 提取质控后矩阵（保持原始矩阵结构，仅筛选行和列）
      X <- X[quality_genes, quality_cells, drop = FALSE]
      
      # 输出质控统计信息
      message(paste0("   📊 质控后矩阵: ", nrow(X), " 基因 × ", ncol(X), " 细胞"))
      message(paste0("   📊 保留基因数: ", length(quality_genes), " (≥", min_cells_per_gene, "个细胞表达)"))
      message(paste0("   📊 保留细胞数: ", length(quality_cells), " (表达基因数", min_genes_per_cell, "~", max_genes_per_cell, ")"))
      
      # 质控后校验：避免后续拟合失败
      if (nrow(X) < 100) {
        message(paste0("   ⚠️ donor ", donor_name, " 质控后基因数仅 ", nrow(X), " 个（<100），跳过处理！"))
        next
      }
      if (ncol(X) < 10) {
        message(paste0("   ⚠️ donor ", donor_name, " 质控后细胞数仅 ", ncol(X), " 个（<10），跳过处理！"))
        next
      }
      
      
      # ------------------------------------------------------------------------
      # 步骤3：计算dropout率（基于过滤后的数据）
      # ------------------------------------------------------------------------
      message("3/8: 计算基因dropout率...")
      gene_detection_rate <- apply(X > 0, 1, sum) / ncol(X)
      dropr <- 1 - gene_detection_rate
      message(paste0("   📊 dropout率范围: ", round(min(dropr), 3), " ~ ", round(max(dropr), 3)))
      
      # 步骤3：计算基因表达的均值（u）、方差（vx）、CV²
      # ------------------------------------------------------------------------
      message("3/8: 计算基因均值、方差与CV²...")
      if (any(is.na(X))) {
        message("   ⚠️ 矩阵中存在NA值，计算时将忽略NA...")
        u <- rowMeans(X, na.rm = TRUE)    # 均值（忽略NA）
        vx <- rowVars(X, na.rm = TRUE)    # 方差（忽略NA，需statmod包）
      } else {
        u <- rowMeans(X)                  # 均值
        vx <- rowVars(X)                  # 方差
      }
      # 记录原始基因数量
      original_count <- length(u)
      # 筛选u>0的基因索引
      valid_idx <- which(u > 0)
      # 计算被删除的u=0基因数量
      removed_zero_zero_count <- original_count - length(valid_idx)
      
      # 输出筛选信息
      message(paste0("   📊 原始基因总数: ", original_count))
      message(paste0("   ⚠️ 发现 ", removed_zero_zero_count, " 个u=0的基因，已直接删除"))
      message(paste0("   ✅ 保留有效基因数量: ", length(valid_idx)))
      
      # 3. 仅保留u>0的基因进行后续计算
      gene_names <- rownames(X)
      u <- u[valid_idx]               # 筛选后的均值
      vx <- vx[valid_idx]             # 筛选后的方差
      gene_names <- gene_names[valid_idx]  # 筛选后的基因名称
      X <- X[valid_idx, , drop = FALSE]    # 筛选后的表达矩阵（可选，若后续需用）
      cv2 <- (vx / (u^2))  # CV²（变异系数的平方）
      
      # ------------------------------------------------------------------------
      # 步骤4：准备GLM拟合数据
      # ------------------------------------------------------------------------
      message("4/8: 准备GLM拟合数据...")
      xi <- 1 / u  # 自变量：1/均值
      yi <- cv2    # 因变量：CV²
      m <- ncol(X)  # 细胞数
      df <- m - 1            # 自由度
      
      # The lower limit of average expression was defined
      minMeanForFit <- unname( quantile( u[ which(cv2> .3 ) ], .85 ) )
      useForFit <- xi <= 1/minMeanForFit & cv2<10 & xi != 0 & yi != 0 & xi != Inf & yi != Inf
      # min(100, all_genes) genes were used when more than the lower limit of average expression was less than 100 genes
      if(sum(useForFit) <= 100) {
        minMeanForFit <- unname(sort(u, decreasing = TRUE)[min(100, length(u))])
        useForFit <- xi <= 1/minMeanForFit & xi != 0 & yi != 0 & xi != Inf & yi != Inf
      }
      
      glm_x <- xi[useForFit]
      glm_y <- yi[useForFit]
      
      # Outliers are filtered using 3 x sd
      sd_threshold <- 3
      mean_glm_y <- mean(glm_y)
      sd_glm_y <- sd(glm_y)
      filtered_rows <- abs(glm_y - mean_glm_y) <= sd_threshold * sd_glm_y
      useforglm <- data.frame(glm_x = glm_x, glm_y = glm_y)[filtered_rows, ]
      # ---------------------- 新增：统计符合拟合条件的基因数 ----------------------
      fitting_gene_count <- nrow(useforglm)  # 核心：useforglm的行数=拟合基因数
      message(paste0("   ✅ 符合拟合条件的基因总数: ", fitting_gene_count, " 个"))
      # ---------------------- 新增：统计拟合基因的CV²（glm_y）分布特征 ----------------------
      if (fitting_gene_count > 0) {  # 避免无拟合基因时统计报错
        # 计算CV²的分位数（25%、50%、75%）和最值
        cv2_quantiles <- quantile(useforglm$glm_y, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
        cv2_min <- min(useforglm$glm_y, na.rm = TRUE)
        cv2_max <- max(useforglm$glm_y, na.rm = TRUE)
        
        # 打印统计结果（保留3位小数，格式清晰）
        message(paste0("   📊 拟合基因CV²分布："))
        message(paste0("      - 最小值: ", round(cv2_min, 3)))
        message(paste0("      - 25分位数: ", round(cv2_quantiles[1], 3)))
        message(paste0("      - 50分位数（中位数）: ", round(cv2_quantiles[2], 3)))
        message(paste0("      - 75分位数: ", round(cv2_quantiles[3], 3)))
        message(paste0("      - 最大值: ", round(cv2_max, 3)))
      } else {
        message(paste0("   ⚠️ 无符合条件的拟合基因，无法统计CV²分布"))
      }
      
      # A generalized linear model (GLM) is used to fit the relationship between xi and yi
      message("begin fit")
      start_values <- c(intercept = 0.1, slope = 0.1)
      mdl <- glm(glm_y ~ glm_x, data = useforglm, family = Gamma(link = "identity"), start = start_values)
      cv2fit <- predict(mdl, newdata = data.frame(glm_x = 1/u))
      message(paste0("   📊 原始cv2fit范围: ", round(min(cv2fit), 3), " ~ ", round(max(cv2fit), 3)))
      message("5.5/8: 检查cv2fit有效值范围...")
      
      # 提取非NA值
      valid_cv2fit <- cv2fit[!is.na(cv2fit)]
      
      if (length(valid_cv2fit) == 0) {
        # 无有效值（全为NA）
        message("   ⚠️ cv2fit中无有效值（全为NA）")
      } else {
        # 有有效值，计算真实范围
        valid_min <- round(min(valid_cv2fit), 3)
        valid_max <- round(max(valid_cv2fit), 3)
        valid_count <- length(valid_cv2fit)
        total_count <- length(cv2fit)
        valid_ratio <- round(valid_count / total_count * 100, 2)
        
        message(paste0("   📊 cv2fit有效值数量: ", valid_count, " (占比 ", valid_ratio, "%)"))
        message(paste0("   📊 cv2fit有效值范围: ", valid_min, " ~ ", valid_max))
      }
      # 1. 检测cv2fit中的Inf
      inf_cv2fit_count <- sum(is.infinite(cv2fit))
      if (inf_cv2fit_count > 0) {
        # 2. 计算填充阈值（取cv2fit非Inf值的99.9%分位数，避免主观值）
        valid_cv2fit <- cv2fit[!is.infinite(cv2fit) & !is.na(cv2fit)]
        fill_threshold <- quantile(valid_cv2fit, 0.999, na.rm = TRUE)  # 卡极大阈值
        
        # 3. 用阈值替换Inf（消除无穷大）
        cv2fit[is.infinite(cv2fit)] <- fill_threshold
        
        # 4. 日志记录
        message(paste0("   ⚠️ 发现 ", inf_cv2fit_count, " 个cv2fit=Inf，用99.9%分位数（", round(fill_threshold, 3), "）填充"))
      }
      
      deviance_res <- residuals(mdl, type = "deviance")
      total_deviance <- sum(deviance_res^2)
      null_deviance <- mdl$null.deviance
      cat(sprintf("模型偏差: %.3f (Null偏差: %.3f)\n", 
                  total_deviance, null_deviance))
      # Extract and calculate model coefficients
      zero_number <- sum(u == 0 | cv2fit == 0)
      zero_number
      b <- coef(mdl)
      Intercept <- b[['(Intercept)']]
      slope <- b[['glm_x']]
      glm_x <- useforglm$glm_x
      glm_y <- useforglm$glm_y
      fitting_genes_number <- length(glm_x)
      predicted <- predict(mdl, newdata = data.frame(glm_x = glm_x))
      residuals <- glm_y - predicted
      SST <- sum((glm_y - mean(glm_y))^2)
      SSE <- sum(residuals^2)
      R_squared <- 1 - SSE/SST
      n <- length(glm_y)
      p <- 2
      Adjusted_R_squared <- 1 - (1 - R_squared) * ((n - 1) / (n - p - 1))
      message(paste0("   样本", donor_name, "的adjusted R²: ", round(Adjusted_R_squared, 3)))
      # 5. 过滤adjusted R² < 0.5的样本
      if (Adjusted_R_squared < 0.3) {
        message(paste0("   ⚠️ adjusted R² < 0.3,不保存结果"))
        next
      }
      
      fitratio <- ifelse(u == 0 | cv2fit == 0, NA, cv2 / cv2fit)
      # 提取非NA值
      valid_fitratio <- fitratio[!is.na(fitratio)]
      
      if (length(valid_fitratio) == 0) {
        # 无有效值（全为NA）
        message("   ⚠️ ratio中无有效值（全为NA）")
      } else {
        # 有有效值，计算真实范围
        valid_min <- round(min(valid_fitratio), 3)
        valid_max <- round(max(valid_fitratio), 3)
        valid_count <- length(valid_fitratio)
        total_count <- length(fitratio)
        valid_ratio <- round(valid_count / total_count * 100, 2)
        
        message(paste0("   📊 ratio有效值数量: ", valid_count, " (占比 ", valid_ratio, "%)"))
        message(paste0("   📊 ratio有效值范围: ", valid_min, " ~ ", valid_max))
      }
      log2fitratio <- log2(fitratio)  # Calculate log-transformed: log2(cv2)-log2(cv2fit)
      
      # pval <- pchisq(fitratio * df, df, lower.tail = FALSE)  # Calculating p-value
      pval <- pchisq(fitratio * df, df, lower.tail = FALSE)
      
      # FDR was calculated using the Bonferroni correction method
      fdr <- p.adjust(pval, method = "bonferroni")
      
      T <- data.frame(g = gene_names,u, cv2, cv2fit, fitratio, log2fitratio, pval, fdr,n,SST,SSE,R_squared,Adjusted_R_squared,Intercept,slope )
      T$cv2fit[is.infinite(T$cv2fit)] <- NA # Convert Inf to NA
      T$log2fitratio[is.infinite(T$log2fitratio)] <- NA
      # 保存结果表格（CSV格式，便于后续分析）
      result_csv_path <- file.path(donor_output_dir, paste0(donor_name, "_origin_gene_analysis.csv"))
      write.csv(T, file = result_csv_path, quote = TRUE, row.names = FALSE)
      message("   ✅ 结果表格已保存至: ", result_csv_path)
      
      ### 相关性检验 默认为pearson      
      # 计算u和CV2的相关系数以及p值
      cor_u_cv2 <- cor.test(log10(T$u), log10(T$cv2))
      cor_u_cv2
      
      # 计算u和fitratio的相关系数以及p值
      cor_u_fitratio <- cor.test(log10(T$u), log10(T$fitratio))
      cor_u_fitratio
      
      # 计算u和 log2fitratio 的相关系数以及p值
      cor_u_log2fitratio <- cor.test(log10(T$u), T$log2fitratio)
      cor_u_log2fitratio
      
      # Create result table
      T1 <- data.frame(g = gene_names,u, m,cv2, cv2fit, fitratio, log2fitratio, pval, fdr,
                       n,SST,SSE,R_squared,Adjusted_R_squared,Intercept,slope,fitting_gene_count,
                       cor_u_log2fitratio$estimate, cor_u_cv2$estimate,cor_u_fitratio$estimate)
      # 创建条件逻辑向量：fitratio > 1 且 fdr < 0.05
      condition <- T1$fitratio > fitratio_threshold & T1$fdr < fdr_threshold
      
      # 确定需要处理的数值列（排除基因ID相关列，这里假设前2列是基因ID）
      # 根据实际数据调整需要保留的列索引
      id_cols <- 1:2  # 基因ID相关列（不处理）
      value_cols <- 4:ncol(T1)  # 数值列（需要处理）
      
      # 对不满足条件的行，将数值列赋为0
      T1[!condition, value_cols] <- 0
      
      # 保存结果表格（CSV格式，便于后续分析）
      result_csv_path <- file.path(donor_output_dir, paste0(donor_name, "_gene_analysis.csv"))
      write.csv(T1, file = result_csv_path, quote = TRUE, row.names = FALSE)
      message("   ✅ 结果表格已保存至: ", result_csv_path)
      
      # 绘制CV²拟合曲线
      cv2_plot_file <- file.path(donor_output_dir, paste0(donor_name, "_cv2.pdf"))
      red_points <- T$fdr < 0.05 & T$fitratio > 1
      
      pdf(cv2_plot_file, width = 8, height = 6)
      col <- "#3881B8"
      plot(NULL, xaxt = "n", yaxt = "n",
           log = "xy", xlim = c(1e-1, 100), ylim = c(0.005, 100),
           xlab = "average normalized counts", ylab = "CV²")
      axis(1, 10^(-1:2), c("0.1", "1", "10", "100"))
      axis(2, 10^(-2:2), c("0.01", "0.1", "1", "10", "100"), las = 2)
      abline(h = 10^(-2:2), v = 10^(-1:2), col = "#D0D0D0", lwd = 2)
      
      # 绘制数据点
      points(u, cv2, pch = 20, cex = 1, col = col)
      points(u[red_points], cv2[red_points], pch = 20, cex = 1, col = "pink")
      
      # 绘制拟合曲线
      fitted_curve <- function(x) predict(mdl, newdata = data.frame(glm_x = 1/x))
      xg <- 10^seq(-2, 2, length.out = 1000)
      lines(xg, fitted_curve(xg), col = "#FF000080", lwd = 3)
      
      # 绘制置信区间
      lines(xg, (fitted_curve(xg)) * qchisq(0.975, df) / df,
            col = "#FF000080", lwd = 2, lty = "dashed")
      lines(xg, (fitted_curve(xg)) * qchisq(0.025, df) / df,
            col = "#FF000080", lwd = 2, lty = "dashed")
      
      # 添加图例
      legend("topright", legend = c(donor_name, 
                                    paste0("cor: ", round(cor_u_cv2$estimate, 3))),
             col = col, pch = 20, bty = "n")
      dev.off()
      cat("CV²拟合曲线已保存至:", cv2_plot_file, "\n")
      
      # 绘制fitratio图
      fitratio_plot_file <- file.path(donor_output_dir, paste0(donor_name, "_fitratio.pdf"))
      pdf(fitratio_plot_file, width = 8, height = 6)
      col <- "#3881B8"
      plot(NULL, xaxt = "n", yaxt = "n",
           log = "xy", xlim = c(1e-1, 100), ylim = c(0.005, 100),
           xlab = "average normalized counts", ylab = "fitratio(CV²obs/CV²fit)")
      axis(1, 10^(-1:2), c("0.1", "1", "10", "100"))
      axis(2, 10^(-2:2), c("0.01", "0.1", "1", "10", "100"), las = 2)
      abline(h = 10^(-2:2), v = 10^(-1:2), col = "#D0D0D0", lwd = 2)
      abline(h = 1, col = "black", lty = 2)  # 添加fitratio=1参考线
      
      # 绘制数据点
      points(T$u, T$fitratio, pch = 20, cex = 1.5, col = col)
      points(T$u[red_points], T$fitratio[red_points], 
             pch = 20, cex = 1.5, col = "pink")
      
      # 添加图例
      legend("topright", legend = c(donor_name,
                                    paste0("cor: ", round(cor_u_log2fitratio$estimate, 3))),
             col = col, pch = 20, bty = "n")
      dev.off()
      cat("fitratio图已保存至:", fitratio_plot_file, "\n")
      
    }, error = function(e) {
      cat("处理", donor_name, "时出错:", e$message, "\n")
    })
  }
  
  cat("\n=====================================\n")
  cat("所有样本处理完成!\n")
}

# 运行函数
batch_process_donor_rds(
  rds_dir="/data/wuwq/noise/COMBAT/rds_files/NK cell",          # 存放所有donor的SCT标准化后RDS文件的目录
  output_root_dir="/data/wuwq/noise/COMBAT/output_root_dir/NK_fit",  # 结果保存的根目录（每个donor会生成独立子目录）
  fdr_threshold = 0.05,  # 筛选条件：FDR阈值
  fitratio_threshold = 1  # 筛选条件：fitratio阈值
)

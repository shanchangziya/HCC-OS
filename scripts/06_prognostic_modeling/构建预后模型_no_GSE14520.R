#!/usr/bin/env Rscript
if (dev.cur() > 1) dev.off()
rm(list = ls())
gc()
# ============================================================
# Prognostic Model Building - TCGA-LIHC and ICGC-LIRI Only
# Excluding GSE14520 dataset
# ============================================================

cat("============================================================\n")
cat("  Prognostic Model Building (TCGA + ICGC only)\n")
cat("============================================================\n\n")

# Load required library
library(Mime1)

# Load original data
cat("=== Loading Data ===\n")
load("QWEN0208.Rdata")

# Set dataset names
names(list_train_vali_Data) <- c("TCGA-LIHC", "ICGC-LIRI", "GSE14520")

cat("Original datasets:\n")
for(i in 1:length(list_train_vali_Data)) {
  cat(sprintf("  %s: %d samples, %d features\n",
              names(list_train_vali_Data)[i],
              nrow(list_train_vali_Data[[i]]),
              ncol(list_train_vali_Data[[i]])))
}

# Remove GSE14520
cat("\n=== Removing GSE14520 ===\n")
list_train_vali_Data <- list_train_vali_Data[c("TCGA-LIHC", "ICGC-LIRI")]

cat("Datasets after removal:\n")
for(i in 1:length(list_train_vali_Data)) {
  cat(sprintf("  %s: %d samples, %d features\n",
              names(list_train_vali_Data)[i],
              nrow(list_train_vali_Data[[i]]),
              ncol(list_train_vali_Data[[i]])))
}

# Data validation and cleaning
cat("\n=== Data Validation and Cleaning ===\n")

for(i in 1:length(list_train_vali_Data)) {
  dataset_name <- names(list_train_vali_Data)[i]
  cat(sprintf("\nProcessing %s:\n", dataset_name))

  data <- list_train_vali_Data[[i]]

  # Check for infinite values
  if(any(is.infinite(as.matrix(data[, -c(1:3)])))) {
    cat("  - Replacing infinite values with NA\n")
    data[, -c(1:3)] <- apply(data[, -c(1:3)], 2, function(x) {
      x[is.infinite(x)] <- NA
      return(x)
    })
  }

  # Check for extreme values (cap at 99th and 1st percentile)
  cat("  - Capping extreme values\n")
  gene_cols <- 4:ncol(data)
  for(j in gene_cols) {
    col_data <- as.numeric(data[, j])
    q99 <- quantile(col_data, 0.99, na.rm = TRUE)
    q01 <- quantile(col_data, 0.01, na.rm = TRUE)

    col_data[col_data > q99] <- q99
    col_data[col_data < q01] <- q01

    data[, j] <- col_data
  }

  # Replace NA with column mean
  cat("  - Imputing missing values\n")
  data[, -c(1:3)] <- apply(data[, -c(1:3)], 2, function(x) {
    x[is.na(x)] <- mean(x, na.rm = TRUE)
    return(x)
  })

  list_train_vali_Data[[i]] <- data
  cat(sprintf("  ✓ Cleaned: %d samples, %d features\n", nrow(data), ncol(data)-3))
}

# Validate candidate genes
cat("\n=== Candidate Genes Validation ===\n")
cat(sprintf("Total candidate genes: %d\n", length(g)))

# Remove genes with zero or near-zero variance
train_data <- list_train_vali_Data[[1]]
gene_vars <- apply(train_data[, g], 2, var, na.rm = TRUE)
low_var_genes <- names(gene_vars[gene_vars < 1e-6])

if(length(low_var_genes) > 0) {
  cat(sprintf("Removing %d genes with near-zero variance\n", length(low_var_genes)))
  g <- setdiff(g, low_var_genes)
}

cat(sprintf("Final candidate genes: %d\n", length(g)))

# Save cleaned data
cat("\n=== Saving Cleaned Data ===\n")
save(list_train_vali_Data, g, file = "QWEN0208_no_GSE14520.Rdata")
cat("Cleaned data saved to: QWEN0208_no_GSE14520.Rdata\n")

# Run ML.Dev.Prog.Sig with error handling
cat("\n============================================================\n")
cat("  Building Prognostic Models\n")
cat("============================================================\n\n")

cat("Training dataset: TCGA-LIHC\n")
cat("Validation datasets: ICGC-LIRI\n")
cat("Candidate genes:", length(g), "\n")
cat("Mode: all (comprehensive model testing)\n")
cat("Random seed: 5201314\n\n")

start_time <- Sys.time()
cat("Start time:", format(start_time, "%Y-%m-%d %H:%M:%S"), "\n\n")

tryCatch({
  res <- ML.Dev.Prog.Sig(
    train_data = list_train_vali_Data$`TCGA-LIHC`,
    list_train_vali_Data = list_train_vali_Data,
    unicox.filter.for.candi = TRUE,
    unicox_p_cutoff = 0.01,
    candidate_genes = g,
    mode = 'all',
    nodesize = 10,
    seed = 5201314
  )

  end_time <- Sys.time()
  elapsed_time <- difftime(end_time, start_time, units = "mins")

  cat("\n============================================================\n")
  cat("  Model Building Successful!\n")
  cat("============================================================\n\n")

  cat("End time:", format(end_time, "%Y-%m-%d %H:%M:%S"), "\n")
  cat("Elapsed time:", round(elapsed_time, 2), "minutes\n\n")

  # Save results
  save(res, file = "res_no_GSE14520.Rdata")
  cat("Results saved to: res_no_GSE14520.Rdata\n\n")

  # Print summary
  cat("=== Model Summary ===\n")
  if(!is.null(res$Cindex.res)) {
    cat(sprintf("C-index results: %d entries\n", nrow(res$Cindex.res)))
    cat("\nDatasets in results:\n")
    print(table(res$Cindex.res$ID))
  }

  if(!is.null(res$ml.res)) {
    cat(sprintf("\nNumber of models built: %d\n", length(res$ml.res)))
    cat("\nModel types:\n")
    print(head(names(res$ml.res), 20))
  }

  if(!is.null(res$Sig.genes)) {
    cat(sprintf("\nSignificant genes identified: %d\n", length(res$Sig.genes)))
    cat("Top 10 genes:", paste(head(res$Sig.genes, 10), collapse = ", "), "\n")
  }

}, error = function(e) {
  cat("\n============================================================\n")
  cat("  ERROR OCCURRED\n")
  cat("============================================================\n\n")

  cat("Error message:", e$message, "\n\n")
  cat("Detailed traceback:\n")
  print(traceback())

  # Save error information
  error_info <- list(
    message = e$message,
    call = e$call,
    timestamp = Sys.time(),
    data_info = list(
      n_samples_train = nrow(list_train_vali_Data[[1]]),
      n_samples_valid = nrow(list_train_vali_Data[[2]]),
      n_genes = length(g),
      datasets = names(list_train_vali_Data)
    )
  )
  save(error_info, file = "error_info_no_GSE14520.Rdata")
  cat("\nError information saved to: error_info_no_GSE14520.Rdata\n")
})

cat("\n============================================================\n")
cat("  Script Completed\n")
cat("============================================================\n")

# Clean up
gc()


load("res_no_GSE14520.Rdata")
load("QWEN0208_no_GSE14520.Rdata")
cindex_dis_all(res,
               validate_set = names(list_train_vali_Data)[-1],
               order =names(list_train_vali_Data),width = 0.45)

# os ----------------------------------------------------------------------
all.auc.1y <- cal_AUC_ml_res(res.by.ML.Dev.Prog.Sig = res,train_data = list_train_vali_Data[[1]],
                             inputmatrix.list = list_train_vali_Data,mode = 'all',AUC_time = 1,
                             auc_cal_method="KM")
all.auc.3y <- cal_AUC_ml_res(res.by.ML.Dev.Prog.Sig = res,train_data = list_train_vali_Data[[1]],
                             inputmatrix.list = list_train_vali_Data,mode = 'all',AUC_time = 3,
                             auc_cal_method="KM")
all.auc.5y <- cal_AUC_ml_res(res.by.ML.Dev.Prog.Sig = res,train_data = list_train_vali_Data[[1]],
                             inputmatrix.list = list_train_vali_Data,mode = 'all',AUC_time = 5,
                             auc_cal_method="KM")
auc_dis_select(list(all.auc.1y,all.auc.3y,all.auc.5y),
               model_name="StepCox[forward] + GBM",
               dataset = names(list_train_vali_Data),
               order= names(list_train_vali_Data),
               year=c(1,3,5))
unicox.rs.res <- cal_unicox_ml_res(res.by.ML.Dev.Prog.Sig = res,optimal.model = "StepCox[forward] + GBM",type ='categorical')
metamodel <- cal_unicox_meta_ml_res(input = unicox.rs.res)
meta_unicox_vis(metamodel,
                dataset = names(list_train_vali_Data))


survplot <- vector("list",2) 
for (i in c(1:2)) {
  print(survplot[[i]]<-rs_sur(res, model_name = "StepCox[forward] + GBM",dataset = names(list_train_vali_Data)[i],
                              #color=c("blue","green"),
                              median.line = "hv",
                              cutoff = 0.5,
                              conf.int = T,
                              xlab="Day",pval.coord=c(1000,0.9)))
}
aplot::plot_list(gglist=survplot,ncol=2)



# 雷达图 ---------------------------------------------------------------------
cindex_dis_all_compact <- function(object, # output of ML.Dev.Prog.Sig mode = 'all'
                                   color = NULL, # three color value for cindex and two color value for mean cindex
                                   dataset_col = NULL, # color value for cohort
                                   validate_set, # input validate datasets name (仍保留，用于排序模型)
                                   order = NULL, # cohort order plot
                                   width = NULL, # width of right plot (p3)
                                   height = NULL, # height of top plot (p2)
                                   tile_size = 0.25,      # 热图格子边框线宽(越小越紧凑)
                                   text_size = 1.9,       # 热图数字字号
                                   axis_y_size = 6.5,     # 模型名字号
                                   tile_height = 0.55     # 每行高度压缩(越小越紧凑，0.45~0.7可调)
) {
  suppressPackageStartupMessages({
    library(ggplot2)
    library(aplot)
  })
  
  if (is.null(width))  width  <- 0.28  # 默认比原来更窄
  if (is.null(height)) height <- 0.008 # 默认比原来更薄
  
  if (is.null(color)) {
    color <- c("#0084A7", "#F5FACD", "#E05D00", "#79AF97", "#8491B4")
  }
  if (is.null(dataset_col)) {
    dataset_col <- c(
      "#3182BDFF", "#E6550DFF", "#31A354FF", "#756BB1FF", "#636363FF", "#6BAED6FF", "#FD8D3CFF", "#74C476FF",
      "#9E9AC8FF", "#969696FF", "#9ECAE1FF", "#FDAE6BFF", "#A1D99BFF", "#BCBDDCFF", "#BDBDBDFF", "#C6DBEFFF",
      "#FDD0A2FF", "#C7E9C0FF", "#DADAEBFF", "#D9D9D9FF"
    )
  }
  
  cindex <- object[["Cindex.res"]]
  cindex$Cindex <- as.numeric(sprintf("%.2f", cindex$Cindex))
  
  # --- 计算总体mean（右侧条形图） ---
  mean_cindex <- aggregate(x = cindex$Cindex, by = list(cindex$Model), FUN = mean)
  colnames(mean_cindex) <- c("Model", "mean")
  mean_cindex$mean <- as.numeric(sprintf("%.3f", mean_cindex$mean))
  mean_cindex$Value <- "Mean C-index in all cohorts"
  
  # --- 用验证集均值来排序模型（仍然用validate_set，只是不再画验证集mean板块） ---
  cindex_validate <- cindex[cindex$ID %in% validate_set, ]
  mean_validate <- aggregate(x = cindex_validate$Cindex, by = list(cindex_validate$Model), FUN = mean)
  colnames(mean_validate) <- c("Model", "mean")
  mean_validate$mean <- as.numeric(sprintf("%.3f", mean_validate$mean))
  mean_validate <- mean_validate[order(mean_validate$mean, decreasing = FALSE), ]
  
  labels <- as.data.frame(unique(cindex$ID))
  colnames(labels) <- "Cohort"
  
  # cohort顺序
  if (!is.null(order)) {
    labels$Cohort <- factor(labels$Cohort, levels = order)
    cindex$ID <- factor(cindex$ID, levels = order)
  }
  
  # 模型顺序按验证集均值（从低到高，保证图上更稳定一致）
  cindex$Model <- factor(cindex$Model, levels = mean_validate$Model)
  mean_cindex$Model <- factor(mean_cindex$Model, levels = mean_validate$Model)
  
  # --- p1：C-index 热图（更紧凑） ---
  p1 <- ggplot(cindex, aes(x = ID, y = Model)) +
    geom_tile(aes(fill = Cindex), color = "white", size = tile_size, height = tile_height) +
    geom_text(aes(label = Cindex), vjust = 0.5, color = "black", size = text_size) +
    scale_fill_gradient2(
      low = color[1], mid = color[2], high = color[3],
      midpoint = median(cindex$Cindex, na.rm = TRUE),
      name = "C-index"
    ) +
    theme_minimal(base_size = 9) +
    theme(
      axis.title = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_text(size = axis_y_size),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      plot.margin = margin(2, 2, 2, 2)
    )
  
  # --- p2：顶部 cohort 色条（更薄） ---
  p2 <- ggplot(labels, aes(Cohort, y = 1)) +
    geom_tile(aes(fill = Cohort), color = "white", size = tile_size) +
    scale_fill_manual(values = dataset_col, name = "Cohort") +
    theme_void() +
    theme(plot.margin = margin(0, 2, 0, 2))
  
  # --- p3：右侧总体 mean 条形图（更窄更紧凑） ---
  p3 <- ggplot(mean_cindex, aes(x = Model, y = mean)) +
    geom_col(fill = color[4], width = 0.75) +
    geom_text(aes(label = mean), hjust = 1.05, size = 2) +
    coord_flip() +
    theme_minimal(base_size = 9) +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      plot.margin = margin(2, 2, 2, 2)
    )
  
  # --- 拼图：去掉p4，只拼 p2 + p1 + p3 ---
  print(
    p1 %>%
      insert_top(p2, height = height) %>%
      insert_right(p3, width = width)
  )
}

p <- cindex_dis_all_compact(
  object = res,
  validate_set = c("ICGC-LIRI"),
  order = c("TCGA-LIHC","ICGC-LIRI")
)


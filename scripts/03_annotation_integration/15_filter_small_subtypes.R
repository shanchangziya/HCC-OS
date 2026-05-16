#==== 删除细胞数少于100的亚群 ====
# 日期: 2026-02-01
# 目的: 过滤掉celltype_fine中细胞数<100的亚群

rm(list = ls())
gc()
library(Seurat)
library(dplyr)

#==== 1. 加载数据 ====
cat("=== 加载整合注释对象 ===\n")
load("sce.all.annotated.Rdata")
cat("原始总细胞数:", ncol(sce.all), "\n\n")

#==== 2. 统计各亚群细胞数 ====
cat("=== 各亚群细胞数统计 ===\n")
celltype_counts <- table(sce.all$celltype_fine)
celltype_counts_sorted <- sort(celltype_counts, decreasing = FALSE)
print(celltype_counts_sorted)

#==== 3. 识别细胞数<100的亚群 ====
cat("\n=== 识别需要删除的亚群 ===\n")
low_count_subtypes <- names(celltype_counts[celltype_counts < 100])
cat("细胞数<100的亚群:\n")
for(subtype in low_count_subtypes) {
  cat(sprintf("  - %s: %d个细胞\n", subtype, celltype_counts[subtype]))
}

n_cells_to_remove <- sum(celltype_counts[celltype_counts < 100])
cat(sprintf("\n总计需要删除: %d个细胞 (%.2f%%)\n",
            n_cells_to_remove,
            n_cells_to_remove / ncol(sce.all) * 100))

#==== 4. 删除这些亚群 ====
if(length(low_count_subtypes) > 0) {
  cat("\n=== 执行删除 ===\n")

  # 记录删除前的信息
  before_cells <- ncol(sce.all)
  before_subtypes <- length(unique(sce.all$celltype_fine))

  # 执行删除
  cells_to_keep <- !(sce.all$celltype_fine %in% low_count_subtypes)
  sce.all.filtered <- subset(sce.all, cells = colnames(sce.all)[cells_to_keep])

  # 删除后的统计
  after_cells <- ncol(sce.all.filtered)
  after_subtypes <- length(unique(sce.all.filtered$celltype_fine))

  cat(sprintf("删除前: %d个细胞, %d个亚群\n", before_cells, before_subtypes))
  cat(sprintf("删除后: %d个细胞, %d个亚群\n", after_cells, after_subtypes))
  cat(sprintf("实际删除: %d个细胞, %d个亚群\n",
              before_cells - after_cells,
              before_subtypes - after_subtypes))

  #==== 5. 保存过滤后的对象 ====
  cat("\n=== 保存过滤后的对象 ===\n")
  sce.all <- sce.all.filtered
  save(sce.all, file = "sce.all.annotated.filtered.Rdata")
  cat("已保存: sce.all.annotated.filtered.Rdata\n")

  #==== 6. 更新汇总表 ====
  cat("\n=== 生成新的亚群统计表 ===\n")
  annotation_table <- sce.all@meta.data %>%
    group_by(celltype, celltype_fine) %>%
    summarise(n_cells = n(), .groups = "drop") %>%
    arrange(celltype, desc(n_cells))

  write.csv(annotation_table, "cell_annotation_summary.filtered.csv", row.names = FALSE)
  cat("已保存: cell_annotation_summary.filtered.csv\n")

  #==== 7. 显示过滤后的详细统计 ====
  cat("\n=== 过滤后的各大类细胞亚群统计 ===\n")
  for(ct in unique(sce.all$celltype)) {
    cat(sprintf("\n【%s】 (总计: %d个细胞)\n", ct, sum(sce.all$celltype == ct)))
    subtypes <- annotation_table %>% filter(celltype == ct)
    for(i in 1:nrow(subtypes)) {
      cat(sprintf("  %s: %d\n", subtypes$celltype_fine[i], subtypes$n_cells[i]))
    }
  }

  cat("\n=== 最终汇总 ===\n")
  cat(sprintf("保留的总细胞数: %d\n", ncol(sce.all)))
  cat(sprintf("保留的总亚群数: %d\n", length(unique(sce.all$celltype_fine))))
  cat("\n删除的亚群: ", paste(low_count_subtypes, collapse = ", "), "\n")

} else {
  cat("\n没有需要删除的亚群（所有亚群细胞数均≥100）\n")
}

cat("\n=== 完成！===\n")

#==== Header ====
# 脚本名称: 18_check_all_cell_objects.R
# 目的: 检查所有细胞群对象的情况，为整合做准备
# 日期: 2026-02-02

rm(list = ls())
gc()
library(Seurat)
library(dplyr)

cat("=== 检查所有细胞群对象 ===\n\n")

# 创建汇总表
summary_df <- data.frame(
  File = character(),
  Object_Name = character(),
  N_Cells = numeric(),
  Key_Metadata = character(),
  Has_OS_Score = character(),
  Has_Malignant = character(),
  Notes = character(),
  stringsAsFactors = FALSE
)

#==== 1. 检查肿瘤细胞（氧化应激分组）====
cat("=== 1. 肿瘤细胞（氧化应激分组）===\n")
if(file.exists("05osnmf.Rdata")) {
  load("05osnmf.Rdata")
  cat("05osnmf.Rdata loaded\n")
  cat("Objects in file:", ls(), "\n")

  if(exists("seu")) {
    cat("\n对象名: seu\n")
    cat("细胞数:", ncol(seu), "\n")
    cat("元数据列:\n")
    print(head(colnames(seu@meta.data), 30))

    # 检查氧化应激相关列
    os_cols <- grep("OS|os|OxStress", colnames(seu@meta.data), value = TRUE)
    cat("\n氧化应激相关列:", paste(os_cols, collapse = ", "), "\n")

    # 检查OS_group
    if("OS_group" %in% colnames(seu@meta.data)) {
      cat("\nOS_group分布:\n")
      print(table(seu$OS_group))
    }

    if("OS_group_valley" %in% colnames(seu@meta.data)) {
      cat("\nOS_group_valley分布:\n")
      print(table(seu$OS_group_valley))
    }

    # 检查malignant_status
    malignant_cols <- grep("malignant|cnv|tumor|normal", colnames(seu@meta.data),
                           value = TRUE, ignore.case = TRUE)
    cat("\n肿瘤/正常相关列:", paste(malignant_cols, collapse = ", "), "\n")

    summary_df <- rbind(summary_df, data.frame(
      File = "05osnmf.Rdata",
      Object_Name = "seu",
      N_Cells = ncol(seu),
      Key_Metadata = paste(head(colnames(seu@meta.data), 10), collapse = "; "),
      Has_OS_Score = ifelse(length(os_cols) > 0, "Yes", "No"),
      Has_Malignant = ifelse(length(malignant_cols) > 0, "Yes", "No"),
      Notes = "NMF + OS analysis"
    ))
  }
  rm(list = setdiff(ls(), c("summary_df")))
  gc()
}

cat("\n\n=== 2. 上皮细胞CNV分析 ===\n")
if(file.exists("03epcnv.Rdata")) {
  load("03epcnv.Rdata")
  cat("03epcnv.Rdata loaded\n")
  cat("Objects in file:", ls(), "\n")

  if(exists("sce2_merge")) {
    cat("\n对象名: sce2_merge\n")
    cat("细胞数:", ncol(sce2_merge), "\n")

    # 检查CNV相关
    if("cnv_class" %in% colnames(sce2_merge@meta.data)) {
      cat("\ncnv_class分布:\n")
      print(table(sce2_merge$cnv_class, useNA = "ifany"))
    }

    if("cnv_fraction" %in% colnames(sce2_merge@meta.data)) {
      cat("\ncnv_fraction范围:", range(sce2_merge$cnv_fraction, na.rm = TRUE), "\n")
    }

    summary_df <- rbind(summary_df, data.frame(
      File = "03epcnv.Rdata",
      Object_Name = "sce2_merge",
      N_Cells = ncol(sce2_merge),
      Key_Metadata = paste(head(colnames(sce2_merge@meta.data), 10), collapse = "; "),
      Has_OS_Score = "No",
      Has_Malignant = "Yes (cnv_class)",
      Notes = "CNV analysis, includes NK-T reference"
    ))
  }
  rm(list = setdiff(ls(), c("summary_df")))
  gc()
}

cat("\n\n=== 3. 清洗后的上皮细胞 ===\n")
if(file.exists("sce.epi.clean.Rdata")) {
  load("sce.epi.clean.Rdata")
  cat("sce.epi.clean.Rdata loaded\n")
  cat("Objects in file:", ls(), "\n")

  if(exists("sce.epi.clean")) {
    cat("\n对象名: sce.epi.clean\n")
    cat("细胞数:", ncol(sce.epi.clean), "\n")
    cat("元数据前20列:\n")
    print(head(colnames(sce.epi.clean@meta.data), 20))

    summary_df <- rbind(summary_df, data.frame(
      File = "sce.epi.clean.Rdata",
      Object_Name = "sce.epi.clean",
      N_Cells = ncol(sce.epi.clean),
      Key_Metadata = paste(head(colnames(sce.epi.clean@meta.data), 10), collapse = "; "),
      Has_OS_Score = "Check",
      Has_Malignant = "Check",
      Notes = "Clean epithelial cells"
    ))
  }
  rm(list = setdiff(ls(), c("summary_df")))
  gc()
}

cat("\n\n=== 4. 免疫细胞细分 ===\n")
immune_files <- c(
  "07_NK_annotated.Rdata",
  "07_T_annotated.Rdata",
  "17_TrmCD8_refined_annotated.Rdata",
  "17_CytoCD4_refined_annotated.Rdata",
  "08_Mac_annotated.Rdata",
  "09_B_annotated.Rdata",
  "10_DC_annotated.Rdata",
  "13_Neu_annotated.Rdata"
)

for(file in immune_files) {
  if(file.exists(file)) {
    cat("\n检查:", file, "\n")
    load(file)
    obj_names <- ls()
    obj_names <- obj_names[obj_names != "summary_df"]

    for(obj in obj_names) {
      if(class(get(obj))[1] %in% c("Seurat", "seurat")) {
        obj_data <- get(obj)
        cat("  对象:", obj, "- 细胞数:", ncol(obj_data), "\n")

        # 查找亚型列
        subtype_cols <- grep("subtype|celltype|annotation",
                             colnames(obj_data@meta.data),
                             value = TRUE, ignore.case = TRUE)
        cat("  亚型列:", paste(subtype_cols, collapse = ", "), "\n")

        if(length(subtype_cols) > 0) {
          cat("  亚型分布:\n")
          print(table(obj_data@meta.data[[subtype_cols[length(subtype_cols)]]]))
        }

        summary_df <- rbind(summary_df, data.frame(
          File = file,
          Object_Name = obj,
          N_Cells = ncol(obj_data),
          Key_Metadata = paste(subtype_cols, collapse = "; "),
          Has_OS_Score = "No",
          Has_Malignant = "No",
          Notes = "Immune cell subdivision"
        ))
      }
    }
    rm(list = setdiff(ls(), c("summary_df", "immune_files", "file")))
    gc()
  } else {
    cat(file, "不存在\n")
  }
}

cat("\n\n=== 5. 基质细胞细分 ===\n")
stromal_files <- c(
  "11_CAF_annotated.Rdata",
  "12_Endo_annotated.Rdata"
)

for(file in stromal_files) {
  if(file.exists(file)) {
    cat("\n检查:", file, "\n")
    load(file)
    obj_names <- ls()
    obj_names <- obj_names[obj_names != "summary_df"]

    for(obj in obj_names) {
      if(class(get(obj))[1] %in% c("Seurat", "seurat")) {
        obj_data <- get(obj)
        cat("  对象:", obj, "- 细胞数:", ncol(obj_data), "\n")

        # 查找亚型列
        subtype_cols <- grep("subtype|celltype|annotation",
                             colnames(obj_data@meta.data),
                             value = TRUE, ignore.case = TRUE)
        cat("  亚型列:", paste(subtype_cols, collapse = ", "), "\n")

        if(length(subtype_cols) > 0) {
          cat("  亚型分布:\n")
          print(table(obj_data@meta.data[[subtype_cols[length(subtype_cols)]]]))
        }

        summary_df <- rbind(summary_df, data.frame(
          File = file,
          Object_Name = obj,
          N_Cells = ncol(obj_data),
          Key_Metadata = paste(subtype_cols, collapse = "; "),
          Has_OS_Score = "No",
          Has_Malignant = "No",
          Notes = "Stromal cell subdivision"
        ))
      }
    }
    rm(list = setdiff(ls(), c("summary_df", "stromal_files", "file")))
    gc()
  } else {
    cat(file, "不存在\n")
  }
}

cat("\n\n=== 汇总表 ===\n")
print(summary_df)

write.csv(summary_df, "18_all_objects_summary.csv", row.names = FALSE)
cat("\n汇总已保存到: 18_all_objects_summary.csv\n")

cat("\n=== 检查完成 ===\n")

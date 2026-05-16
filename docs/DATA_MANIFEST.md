# Data Manifest

This repository tracks code only. The following files are expected as external inputs or generated intermediates.

## External Inputs

- `data/03注释.Rdata`: initial annotated scRNA-seq Seurat object.
- `data/02fastcnvST.Rdata`: spatial transcriptomics object used for ST oxidative-stress scoring.
- `visium_ST_processed.RData`: processed Visium object for RCTD.
- `bulk/TCGA-LIHC.Rdata`: TCGA-LIHC bulk expression and survival input.
- `bulk/exp1surv1.rdata`: bulk expression and survival phenotype for Scissor/downstream analysis.
- `bulk/TCGA临床信息.txt`: TCGA clinical metadata.
- `QWEN0208.Rdata`: multi-cohort expression-survival object and candidate gene vector.
- `resnet50_features.csv`: WSI tile-level ResNet50 feature matrix.

## Generated Intermediate Objects

- `sceall.Rdata`
- `01sceall.Rdata`
- `02ep.Rdata`
- `03epcnv.Rdata`
- `03tc.Rdata`
- `tcos.Rdata`
- `04nmfsce.Rdata`
- `05osnmf.Rdata`
- `07_*_annotated.Rdata`
- `08_Mac_annotated.Rdata`
- `09_B_annotated.Rdata`
- `10_DC_annotated.Rdata`
- `11_CAF_annotated.Rdata`
- `12_Endo_annotated.Rdata`
- `13_Neu_annotated.Rdata`
- `17_*_refined_annotated.Rdata`
- `19_sce.all.integrated.Rdata`
- `注释sceall.Rdata`
- `QWEN0208_no_GSE14520.Rdata`
- `res_no_GSE14520.Rdata`
- `Resnet.Rdata`
- `pathdat.Rdata`

Do not commit these files unless the journal explicitly requires small derived examples and data-use permissions allow redistribution.

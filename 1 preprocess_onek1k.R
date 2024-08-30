library(Matrix)
library(readr)
library(Seurat)
library(SeuratObject)
library(data.table)

setwd("")

input_dir <- file.path(getwd(), "input")
output_dir <- file.path(getwd(), "output")

# -------------------------------------------------

file_path_onek1k <- file.path(input_dir, "OneK1K.rds")
onek1k <- readRDS(file_path_onek1k)

#str(onek1k)

onek1k@assays[["RNA"]]@counts <- onek1k@assays[["RNA"]]@data

cat("OneK1K loaded. \n")

# -------------------------------------------------

cat("Filtering out low-quality cells and potential doublets... \n")

cat("Number of cells before filtering:", ncol(onek1k), "\n")

onek1k <- subset(onek1k, subset = predicted.celltype.l2 != "Doublet")

cat("Number of cells after doublet removal:", ncol(onek1k), "\n")

onek1k <- subset(onek1k, subset = (nCount_RNA <= 30000 & nFeature_RNA <= 5000))

cat("Number of cells after quality control:", ncol(onek1k), "\n")

# -------------------------------------------------

cat("Running NormalizeData()...\n")
onek1k <- NormalizeData(onek1k)

cat("Running FindVariableFeatures()...\n")
onek1k <- FindVariableFeatures(onek1k)

cat("Running ScaleData()...\n")
onek1k <- ScaleData(onek1k, vars.to.regress = c("percent.mt", "pool_number"))

cat("Running RunPCA()...\n")
onek1k <- RunPCA(onek1k)

cat("Running FindNeighbors()...\n")
onek1k <- FindNeighbors(onek1k, dims = 1:10)

cat("Running FindClusters()...\n")
onek1k <- FindClusters(onek1k, resolution = 0.5)

cat("Saving OneK1K_preprocessed.rds...\n")
onek1k@assays$RNA@data <- as.matrix(0)
onek1k@assays$RNA@scale.data <- as.matrix(0)

saveRDS(onek1k, file = file.path(input_dir, "OneK1K_preprocessed.rds"))

cat("Saving meta data...\n")
meta_data <- onek1k@meta.data
file_path_meta_data_csv <- file.path(input_dir, "OneK1K_meta_data.csv")
write.csv(meta_data, file_path_meta_data_csv, row.names = TRUE)

cat("Done!")


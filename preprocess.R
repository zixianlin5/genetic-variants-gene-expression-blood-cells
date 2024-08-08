library(Matrix)
library(readr)
library(Seurat)
library(SeuratObject)
library(data.table)

setwd("")

input_dir <- file.path(getwd(), "input")
output_dir <- file.path(getwd(), "output")

cat("Prepared. \n")

# -------------------------------------------------

# read OneK1K
file_path_onek1k <- file.path(input_dir, "OneK1K.rds")
onek1k <- readRDS(file_path_onek1k)

cat("OneK1K1 read.\n")
str(onek1k)

onek1k@assays[["RNA"]]@counts <- onek1k@assays[["RNA"]]@data

# run PCA
cat("Running NormalizeData().\n")
onek1k <- NormalizeData(onek1k)

cat("Running FindVariableFeatures().\n")
onek1k <- FindVariableFeatures(onek1k)

cat("Running ScaleData().\n")
onek1k <- ScaleData(onek1k, vars.to.regress = c("percent.mt", "pool_number"))

cat("Running RunPCA().\n")
onek1k <- RunPCA(onek1k)

cat("Saving OneK1K_modified0_unf.rds.\n")
onek1k1 <- onek1k
onek1k1@assays$RNA@data <- as.matrix(0)
onek1k1@assays$RNA@scale.data <- as.matrix(0)
saveRDS(onek1k1, file = file.path(input_dir, "OneK1K_modified0_unf.rds"))
rm(onek1k1)

cat("Running FindNeighbors().\n")
onek1k <- FindNeighbors(onek1k, dims = 1:10)

cat("Running FindClusters().\n")
onek1k <- FindClusters(onek1k, resolution = 0.5)

cat("Saving OneK1K_modified0.rds.\n")
onek1k@assays$RNA@data <- as.matrix(0)
onek1k@assays$RNA@scale.data <- as.matrix(0)
file.remove(file.path(input_dir, "OneK1K_modified0_unf.rds"))
saveRDS(onek1k, file = file.path(input_dir, "OneK1K_modified0.rds"))

cat("\n========================================================\n")
cat("\nPreprocessing completed.")


 
library(Seurat)
library(SeuratDisk)

setwd("")

cat("Loading data...\n")

input_dir <- file.path(getwd(), "input")
output_dir <- file.path(getwd(), "output")

file_path_onek1k <- file.path(input_dir, "OneK1K_modified0.rds")
onek1k <- readRDS(file_path_onek1k)

# ---------------------------------------------------------

cat("Saving meta data...\n")
meta_data <- onek1k@meta.data
file_path_meta_data_csv <- file.path(input_dir, "OneK1K_meta_data.csv")
write.csv(meta_data, file_path_meta_data_csv, row.names = TRUE)

# ---------------------------------------------------------

cat("All tasks completed successfully.\n")

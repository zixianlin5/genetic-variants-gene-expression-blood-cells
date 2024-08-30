library(dplyr)
library(tibble)
library(Matrix)
library(Seurat)
library(SeuratObject)

setwd("")

input_dir <- file.path(getwd(), "input")
output_dir <- file.path(getwd(), "output")

# ---------------------------------------------------------

file_path_onek1k <- file.path(input_dir, "OneK1K_preprocessed.rds")
onek1k <- readRDS(file_path_onek1k)

variable_genes <- VariableFeatures(onek1k)
individuals <- colnames(onek1k@assays[["RNA"]]@counts)

file_path_esnp <- file.path(input_dir, "esnp_table.tsv.gz")
con <- gzfile(file_path_esnp, "rt", encoding = "UTF-8")
esnp_table <- read.table(con, header = TRUE, sep = "\t")
close(con)

cat("Data loaded. \n")

# ---------------------------------------------------------

filtered_df <- esnp_table %>%
  dplyr::filter(ROUND == 1, GENE_ID %in% variable_genes) %>%
  dplyr::select(SNPID, GENE_ID) %>%
  dplyr::distinct()

# Generate a list of SNPIDs corresponding to each GENE_ID
gene_snp_list <- filtered_df %>%
  dplyr::group_by(GENE_ID) %>%
  dplyr::summarize(SNPIDs = list(SNPID), .groups = 'drop') %>%
  tibble::deframe()

cat("Saving gene_snp_list...\n")

saveRDS(gene_snp_list, file.path(input_dir, "gene_snp_list.rds"))

cat("Done!")


library(monocle3)
library(Matrix)
library(dplyr)
library(tibble)

setwd("")

input_dir <- file.path(getwd(), "input")
output_dir <- file.path(getwd(), "output")

# ---------------------------------------------------------

file_path_sparse <- file.path(input_dir, "sparse_matrix.mtx")
sparse_matrix <- readMM(file_path_sparse)

expres_matrix <- as(sparse_matrix, "CsparseMatrix")

file_path_rows <- file.path(input_dir, "rows.csv")
file_path_cols <- file.path(input_dir, "cols.csv")

rows <- read.csv(file_path_rows, header = FALSE)$V1
cols <- read.csv(file_path_cols, header = FALSE)$V1

rownames(expres_matrix) <- rows
colnames(expres_matrix) <- cols

rm(sparse_matrix)

cat("Expression matrix generated.\n")

# ---------------------------------------------------------

file_path_gene_snp_list <- file.path(input_dir, "gene_snp_list.rds")
gene_snp_list <- readRDS(file_path_gene_snp_list)

file_path_esnp <- file.path(input_dir, "esnp_table.tsv.gz")
con <- gzfile(file_path_esnp, "rt", encoding = "UTF-8")
esnp_table <- read.table(con, header = TRUE, sep = "\t")
close(con)

gene_snp_df <- data.frame(GENE_ID = rep(names(gene_snp_list), lengths(gene_snp_list)),
                          SNPID = unlist(gene_snp_list))

rownames(gene_snp_df) <- NULL

gene_snp_df <- merge(gene_snp_df, 
                     esnp_table[, c("GENE_ID", "GENE", "SNPID", "RSID")], 
                     by = c("GENE_ID", "SNPID"))

gene_snp_df <- gene_snp_df[, c("GENE_ID", "GENE", "SNPID", "RSID")]

file_path_gene_snp_df <- file.path(input_dir, "gene_snp_df.csv")
write.csv(gene_snp_df, file_path_gene_snp_df, row.names = FALSE)

gene_metadata <- gene_snp_df %>%
  select(c("GENE_ID", "GENE")) %>%
  distinct() %>%
  dplyr::rename(gene_short_name = GENE)

gene_metadata1 <- gene_metadata %>%
  dplyr::slice(match(rownames(expres_matrix), GENE_ID))

rownames(gene_metadata) <- gene_metadata$GENE_ID

cat("Gene metadata generated.\n")

# ---------------------------------------------------------

file_path_meta_data_csv <- file.path(input_dir, "OneK1K_meta_data.csv")
cell_metadata <- read.csv(file = file_path_meta_data_csv, row.names = 1)

cell_type_abbreviations <- list(
  "central memory CD4-positive, alpha-beta T cell" = "CD4+ TCM",
  "natural killer cell" = "NK",
  "naive thymus-derived CD4-positive, alpha-beta T cell" = "CD4+ TN",
  "effector memory CD8-positive, alpha-beta T cell" = "CD8+ TEM",
  "naive B cell" = "B naive",
  "central memory CD8-positive, alpha-beta T cell" = "CD8+ TCM",
  "effector memory CD4-positive, alpha-beta T cell" = "CD4+ TEM",
  "memory B cell" = "B mem",
  "erythrocyte" = "Ery",
  "naive thymus-derived CD8-positive, alpha-beta T cell" = "CD8+ TN",
  "CD16-negative, CD56-bright natural killer cell, human" = "CD56bright NK",
  "transitional stage B cell" = "B trans",
  "CD14-low, CD16-positive monocyte" = "CD14low CD16+ Mono",
  "mucosal invariant T cell" = "MAIT",
  "regulatory T cell" = "Treg",
  "CD14-positive monocyte" = "CD14+ Mono",
  "gamma-delta T cell" = "Tγδ",
  "CD4-positive, alpha-beta cytotoxic T cell" = "CD4+ CTL",
  "double negative thymocyte" = "DN Thymo",
  "hematopoietic precursor cell" = "HPC",
  "platelet" = "Plt",
  "plasmablast" = "PB",
  "innate lymphoid cell" = "ILC",
  "conventional dendritic cell" = "cDC",
  "plasmacytoid dendritic cell" = "pDC",
  "CD8-positive, alpha-beta T cell" = "CD8+ T",
  "CD4-positive, alpha-beta T cell" = "CD4+ T",
  "dendritic cell" = "DC"
)

cell_metadata <- cell_metadata %>%
  select(
    nCount_RNA,
    nFeature_RNA,
    percent.mt,
    pool_number,
    donor_id,
    predicted.celltype.l2,
    age,
    cell_type,
    sex
  ) %>%
  rename(celltype_pred = predicted.celltype.l2)

cell_metadata$barcode <- rownames(cell_metadata)

cell_metadata <- cell_metadata %>%
  mutate(cell_type_abbr = recode(cell_type, !!!cell_type_abbreviations))

rownames(cell_metadata) <- cell_metadata$barcode

cell_metadata <- cell_metadata[match(colnames(expres_matrix), 
                                     rownames(cell_metadata)), ]

cat("Cell metadata generated.\n")

file_path_metadata_new <- file.path(input_dir, "cell_metadata_new.rds")
saveRDS(cell_metadata, file_path_metadata_new)

cat("Cell metadata saved.\n")

# ---------------------------------------------------------

cds <- new_cell_data_set(expres_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_metadata)

cat("Cell data set created.\n")

# ---------------------------------------------------------

file_path_cds <- file.path(input_dir, "cds_object.rds")
saveRDS(cds, file = file_path_cds)

cat("CDS object saved to:", file_path_cds, "\n")

cat("Done!")


library(Matrix)
library(Seurat)
library(SeuratObject)
library(ggplot2)
library(patchwork)
library(cowplot)
library(openxlsx)

setwd("")

input_dir <- file.path(getwd(), "input")
output_dir <- file.path(getwd(), "output")

file_path_onek1k <- file.path(input_dir, "OneK1K_preprocessed.rds")
onek1k <- readRDS(file_path_onek1k)

cat("Data loaded.\n")

dim_PCs <- 20

# -------------------------------------------------

# ===================
# Determine the ‘dimensionality’
# ===================

cat("Saving elbow plot... \n")

# visualize the top principal components to understand how much variance
# they capture and their importance in your dataset

elbow_plot <- ElbowPlot(onek1k, ndims = 40) +
  scale_x_continuous(breaks = seq(0, 40, by = 5))

ggsave(
  file.path(output_dir, "elbow_plot.png"),
  plot = elbow_plot,
  width = 8,
  height = 3.8
)

# ===================
# Explore Gene Loadings on PCs
# ===================

cat("Saving top loadings...\n")

loadings <- Loadings(onek1k[["pca"]])
top_loadings <- loadings[, 1:dim_PCs]

# Initialize a list to store data frames for each PC
loadings_list <- vector("list", length = dim_PCs)

# Generate data frames for each PC in the loop
for (i in 1:dim_PCs) {
  top_genes <- head(sort(abs(top_loadings[, i]), decreasing = TRUE), n = 15)
  sorted_genes <- names(top_genes)
  df <- data.frame(Gene = sorted_genes,
                   Loading = loadings[sorted_genes, i],
                   PC = paste("PC", i))
  loadings_list[[i]] <- df
}

# Combine all data frames into one
combined_loadings_df <- do.call(rbind, loadings_list)

# Save the combined data frame to an Excel file
wb <- createWorkbook()
addWorksheet(wb, "Top_Loadings")
writeData(wb, sheet = "Top_Loadings", x = combined_loadings_df)
saveWorkbook(wb, file.path(output_dir, "Top_Loadings.xlsx"), overwrite = TRUE)

# plots will be ploted by Python

# ===================
# Save data for visualization
# ===================

cat("Saving data for visualization... \n")

pca_coordinates <- as.data.frame(Embeddings(onek1k, reduction = "pca")[, 1:25])
umap_coords <- as.data.frame(Embeddings(onek1k, reduction = "umap"))
metadata <- onek1k@meta.data[, c("cell_type", "sex", "age")]

visual_data <- cbind(pca_coordinates, umap_coords, metadata)
visual_data$cell_name <- rownames(visual_data)

file_path_visual_data <- file.path(input_dir, "visual_data.csv")
write.csv(visual_data, file = file_path_visual_data, row.names = FALSE)

cat("Done!")


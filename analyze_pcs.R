library(Matrix, lib.loc = getwd())
library(Seurat)
library(SeuratObject)
library(ggplot2)
library(patchwork)
library(cowplot)
library(openxlsx)

setwd("")

input_dir <- file.path(getwd(), "input")
output_dir <- file.path(getwd(), "output")

dim_PCs <- 20

cat("Preparation complete. \n")

# -------------------------------------------------

# read OneK1K
file_path_onek1k <- file.path(input_dir, "OneK1K_modified.rds")
onek1k <- readRDS(file_path_onek1k)

cat("Data loaded.\n")

# -------------------------------------------------

# ===================
# Determine the ‘dimensionality’
# ===================

cat("Printing DimHeatmap.png.\n")

# visualize the top principal components to understand how much variance 
# they capture and their importance in your dataset

elbow_plot <- ElbowPlot(onek1k, ndims = 25)
ggsave(file.path(output_dir, "elbow_rna.png"), 
       plot = elbow_plot, width = 8, height = 3.8)

#heatmap_plot <- DimHeatmap(onek1k, dims = 1:dim_PCs, cells = 500, balanced = TRUE)
#ggsave(file.path(output_dir, "DimHeatmap.png"), plot = heatmap_plot, width = 20, height = 15, dpi = 300, limitsize = FALSE)

# ===================
# Explore Gene Loadings on PCs
# ===================

cat("Saving plots of top loadings.\n")

loadings <- Loadings(onek1k[["pca"]])
top_loadings <- loadings[, 1:dim_PCs]

# Initialize a list to store data frames for each PC
loadings_list <- vector("list", length = dim_PCs)

# Generate data frames for each PC in the loop
for (i in 1:dim_PCs) {
  top_genes <- head(sort(abs(top_loadings[, i]), decreasing = TRUE), n = 15)
  sorted_genes <- names(top_genes)
  df <- data.frame(
    Gene = sorted_genes,
    Loading = loadings[sorted_genes, i],
    PC = paste("PC", i)
  )
  loadings_list[[i]] <- df
}

# Combine all data frames into one
combined_loadings_df <- do.call(rbind, loadings_list)

# Save the combined data frame to an Excel file
wb <- createWorkbook()
addWorksheet(wb, "Top_Loadings")
writeData(wb, sheet = "Top_Loadings", x = combined_loadings_df)
saveWorkbook(wb, file.path(output_dir, "Top_Loadings.xlsx"), overwrite = TRUE)

# Function to generate plots
generate_top_genes_plot <- function(dims) {
  plot_list <- lapply(dims, function(i) {
    top_genes <- head(sort(abs(top_loadings[, i]), decreasing = TRUE), n = 15)
    sorted_genes <- names(top_genes)
    plot <- ggplot(data.frame(Gene = sorted_genes, 
                              Loading = loadings[sorted_genes, i]), 
                   aes(x = reorder(Gene, abs(Loading)), y = Loading)) +
      geom_col() +
      coord_flip() +
      labs(title = paste("Top genes for PC", i)) +
      theme(plot.title = element_text(hjust = 0.5, size = 10))
    return(plot)
  })
  combined_plot <- wrap_plots(plot_list, ncol = 3)
  return(combined_plot)
}

top_genes_plot_1 <- generate_top_genes_plot(1:12)
top_genes_plot_2 <- generate_top_genes_plot(13:dim_PCs)

ggsave(file.path(output_dir, "Top_genes_PC_1_12.png"), 
       plot = top_genes_plot_1, 
       width = 18, height = 20,
       dpi = 300, limitsize = FALSE)

ggsave(file.path(output_dir, "Top_genes_PC_13_20.png"), 
       plot = top_genes_plot_2, 
       width = 18, height = 15,
       dpi = 300, limitsize = FALSE)

# ===================
# Visualization of PCs
# ===================

cat("Visualizing PCs\n")

plots <- list()
for (i in 2:dim_PCs) {
  
  plot <- DimPlot(onek1k, reduction = "pca", 
                  group.by = "cell_type", dims = c(i-1, i)) +
    labs(title = sprintf("PCA Plot of PC%d and PC%d", i-1, i)) +
    theme(legend.position = "none",
          plot.title = element_text(size = 10),
          axis.title = element_text(size = 8))
  
  plots[[i-1]] <- plot
}


legend_plot <- DimPlot(onek1k, reduction = "pca", 
                       group.by = "cell_type", dims = c(1, 2)) +
  theme(legend.position = "bottom", 
        legend.text = element_text(size = 8),
        legend.title = element_blank(),
        legend.background = element_rect(fill = "white", color = NA),
        legend.spacing.x = unit(1, "cm"))

legend <- get_legend(legend_plot)

combined_plot <- wrap_plots(plots, ncol = 4)

final_plot <- combined_plot + 
  plot_layout(guides = 'collect') & theme(legend.position = 'bottom')

ggsave(file.path(output_dir, "PCA_Plots.png"),
       plot = final_plot, 
       width = 18, height = 20, dpi = 300)

# -------------------------------------------------

cat("\n========================================================\n")

cat("All finished!!!!")


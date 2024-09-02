library(arrow)
library(ggplot2)

setwd("")

input_dir <- file.path(getwd(), "input")
output_dir <- file.path(getwd(), "output")

fitted_df_null <- readRDS(file.path(input_dir, "fitted_value_null.rds"))
fitted_df_full <- readRDS(file.path(input_dir, "fitted_value_full.rds"))
fitted_df_interaction <- readRDS(file.path(input_dir, "fitted_value_interaction.rds"))
model_comparison_df <- readRDS(file.path(input_dir, "model_comparison.rds"))

cat("Data loaded.\n")

# ---------------------------------------------------------

n <- nrow(fitted_df_null)
adjusted_threshold <- 0.05 / n

cat("Number of genes:", n, "\n")
cat("Adjusted threshold:", adjusted_threshold, "\n")

target_df <- matrix(NA, nrow = nrow(fitted_df_null), ncol = ncol(fitted_df_null))
rownames(target_df) <- rownames(fitted_df_null)
colnames(target_df) <- colnames(fitted_df_null)
target_df <- as.data.frame(target_df)

model_log_df <- data.frame(gene = model_comparison_df$gene, used_model = NA)

for (i in 1:n) {
  
  gene_id <- model_comparison_df$gene[i]
  p_val_null_vs_full <- model_comparison_df$p_val_null_vs_full[i]
  p_val_full_vs_interaction <- model_comparison_df$p_val_full_vs_interaction[i]
  
  if (p_val_null_vs_full >= adjusted_threshold &&
      p_val_full_vs_interaction >= adjusted_threshold) {
    
    target_df[gene_id, ] <- fitted_df_null[gene_id, ]
    model_log_df$used_model[i] <- "Model 1: Basic model"
    
  } else if (p_val_null_vs_full < adjusted_threshold &&
             p_val_full_vs_interaction >= adjusted_threshold) {
    
    target_df[gene_id, ] <- fitted_df_full[gene_id, ]
    model_log_df$used_model[i] <- "Model 2: Model with genotypes"
    
  } else if (p_val_full_vs_interaction < adjusted_threshold) {
    
    target_df[gene_id, ] <- fitted_df_interaction[gene_id, ]
    model_log_df$used_model[i] <- "Model 3: Model with interactions"
    
  } else {
    
    # Fallback condition in case none of the above conditions are met
    target_df[gene_id, ] <- fitted_df_full[gene_id, ]  # Default to full model
    model_log_df$used_model[i] <- "Model 2: Model with genotypes (Default)"
    
    cat(
      "Warning: None of the conditions were met for gene ID ",
      gene_id,
      ". Defaulting to the model with genotypes.\n"
    )
  }
  
  if (i %% 10 == 0 || i == n) {
    progress_percentage <- round((i / n) * 100, 2)
    cat("Progress:", progress_percentage, "%\n")
  }
  
}

# ---------------------------------------------------------

cat("Saving data... \n")

saveRDS(target_df, file.path(input_dir, "adjusted_expres_df.rds"))
saveRDS(model_log_df,
        file.path(input_dir, "model_selection_log.rds"))

cat(".rds Data saved.\n")

target_df$rownames <- rownames(target_df)

file_path_feather <- file.path(input_dir, "adjusted_expres_df.feather")
write_feather(target_df, file_path_feather)

cat(".feather Data saved.\n")

# ---------------------------------------------------------

model_count <- as.data.frame(table(model_log_df$used_model))
colnames(model_count) <- c("Model", "Count")

plot <- ggplot(model_count, aes(x = Model, y = Count, fill = Model)) +
  geom_bar(stat = "identity",
           width = 0.4,
           position = position_dodge(width = 0.5)) +
  geom_text(aes(label = Count), vjust = -0.3, size = 4.7) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5,
      size = 11
    ),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    legend.position = "none"
  ) +
  labs(title = "Distribution of Selected Models", 
       x = "Model Type", 
       y = "Number of Genes") +
  scale_fill_brewer(palette = "Set2")

output_path <- file.path(output_dir, "model_selection_distribution.png")
ggsave(
  output_path,
  plot = plot,
  width = 9,
  height = 5,
  dpi = 300
)

cat("Plot saved to:", output_path, "\n")

cat("Done!")


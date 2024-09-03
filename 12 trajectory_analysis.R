library(monocle3)
library(Matrix)
library(dplyr)
library(tibble)
library(ggplot2)

setwd("")

input_dir <- file.path(getwd(), "input")
output_dir <- file.path(getwd(), "output")

file_path_cds <- file.path(input_dir, "cds_object.rds")
cds <- readRDS(file_path_cds)

file_path_metadata <- file.path(input_dir, "cell_metadata_new.rds")
cell_metadata <- readRDS(file_path_metadata)

cat("Data loaded. \n")

# ---------------------------------------------------------

n_repeats <- 100
total_target_cells <- 100000
min_cells_per_type <- 50

control_params <- list(
  nn.k = 20,
  nn.method = "annoy",
  prune_graph = TRUE,
  nn.n_trees = 30
)

# ---------------------------------------------------------

get_earliest_principal_node <- function(cds, cell_type = "hematopoietic precursor cell") {
  cell_ids <- which(colData(cds)[, "cell_type"] == cell_type)
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids, ]))))]
  
  root_pr_nodes
}

all_deg_results <- list()
all_pseudotime_list <- list()
all_cd_downsamples <- list()

for (i in 1:n_repeats) {
  cat("Iteration", i, "of downsampling and analysis... \n")
  
  # Downsampling process
  cell_type_counts <- cell_metadata %>%
    group_by(cell_type) %>%
    summarise(count = n()) %>%
    mutate(
      proportion = count / sum(count),
      target_sample_size = round(proportion * total_target_cells)
    ) %>%
    mutate(final_sample_size = pmin(pmax(
      target_sample_size, min_cells_per_type
    ), count))
  
  final_total_cells <- sum(cell_type_counts$final_sample_size)
  
  if (final_total_cells > total_target_cells) {
    scaling_factor <- total_target_cells / final_total_cells
    cell_type_counts <- cell_type_counts %>%
      mutate(final_sample_size = round(final_sample_size * scaling_factor))
  }
  
  downsampled_cells <- cell_metadata %>%
    inner_join(cell_type_counts, by = "cell_type") %>%
    group_by(cell_type) %>%
    group_modify(~ {
      sample_size <- unique(.x$final_sample_size)
      available_cells <- nrow(.x)
      
      # Ensure sample size does not exceed available cells
      if (sample_size > available_cells) {
        sample_size <- available_cells
      }
      
      .x[sample(seq_len(available_cells),
                size = sample_size,
                replace = FALSE), ]
    }) %>%
    pull(barcode)
  
  cds_downsampled <- cds[, downsampled_cells]
  
  cat("  Data sampled. \n")
  
  # Preprocessing
  cds_downsampled <- preprocess_cds(cds_downsampled, num_dim = 50)
  cds_downsampled <- reduce_dimension(cds_downsampled, reduction_method = 'UMAP')
  cds_downsampled <- cluster_cells(cds_downsampled)
  
  cat("  Data preprocessed. \n")
  
  # Learn graph and order cells
  tryCatch({
    cds_downsampled <- learn_graph(
      cds_downsampled,
      learn_graph_control = control_params,
      use_partition = FALSE,
      close_loop = TRUE
    )
  }, error = function(e) {
    cat("Error in learn_graph: ", e$message, "\n")
    next
  })
  
  cat("  Graph learned. \n")
  
  tryCatch({
    cds_downsampled <- order_cells(
      cds_downsampled, 
      root_pr_nodes = get_earliest_principal_node(cds_downsampled)
      )
  }, error = function(e) {
    cat("Error in order_cells:", e$message, "\n")
    next
  })
  
  cat("  Cells ordered. \n")
  
  pseudotime <- cds_downsampled@principal_graph_aux[['UMAP']]$pseudotime
  all_pseudotime_list[[i]] <- pseudotime
  
  all_cd_downsamples[[i]] <- cds_downsampled
  
  deg_results <- graph_test(cds_downsampled, neighbor_graph = "principal_graph")
  cat("  graph_test completed. \n")
  
  all_deg_results[[i]] <- deg_results
  
  cat("Iteration", i, "completed.\n")
  rm(cds_downsampled, downsampled_cells)
  
}

# ---------------------------------------------------------

pseudotime_matrix <- do.call(cbind, all_pseudotime_list)

# Calculate median pseudotime across iterations for each cell
median_pseudotime <- apply(pseudotime_matrix, 1, median, na.rm = TRUE)

# Calculate deviation from the median for each iteration
deviation_list <- sapply(1:ncol(pseudotime_matrix), function(i) {
  sum(abs(pseudotime_matrix[, i] - median_pseudotime), na.rm = TRUE)
})

# Identify the iteration with the smallest deviation (most representative)
most_representative_iteration <- which.min(deviation_list)

cat("Most representative iteration is:",
    most_representative_iteration,
    "\n")

cds_most_representative <- all_cd_downsamples[[most_representative_iteration]]

file_path_cds_most_representative <- file.path(input_dir, "cds_most_representative.rds")
saveRDS(cds_most_representative, file_path_cds_most_representative)

cat("cds_most_representative saved. \n")

plot_pseudotime <- plot_cells(
  cds_most_representative,
  color_cells_by = "pseudotime",
  label_groups_by_cluster = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE,
  label_cell_groups = FALSE,
  graph_label_size = 1.5
) +
  theme(text = element_text(size = 10))

ggsave(
  filename = file.path(output_dir, "trajectory_pseudotime.png"),
  plot = plot_pseudotime,
  width = 9,
  height = 6,
  dpi = 300
)

plot_celltype <- plot_cells(
  cds_most_representative,
  color_cells_by = "cell_type_abbr",
  label_groups_by_cluster = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE,
  group_label_size = 3,
  label_cell_groups = TRUE,
  graph_label_size = 1.5
) +
  theme(text = element_text(size = 10))

ggsave(
  filename = file.path(output_dir, "trajectory_celltype.png"),
  plot = plot_celltype,
  width = 8.2,
  height = 6,
  dpi = 300
)

# ---------------------------------------------------------

# Calculate the number of cells expressing each gene in the original CDS
expression_matrix <- exprs(cds)
gene_expression_counts <- rowSums(expression_matrix > 0)
gene_expression_df <- data.frame(GENE_ID = rownames(expression_matrix),
                                 num_cells_expressed = gene_expression_counts)

# ---------------------------------------------------------

combined_deg_results <- bind_rows(all_deg_results) %>%
  filter(!is.na(p_value))

# Combine p-values using Fisher's method
combined_deg_summary <- combined_deg_results %>%
  group_by(GENE_ID) %>%
  summarise(combined_p_value = pchisq(-2 * sum(log(p_value)), 
                                      df = 2 * n(), 
                                      lower.tail = FALSE))

# Adjust for multiple comparisons
combined_deg_summary <- combined_deg_summary %>%
  mutate(
    q_value = p.adjust(combined_p_value, method = "BH"),
    significant = ifelse(q_value < 0.05, "Yes", "No")
  )

cat(
  "Number of significantly differentially expressed genes identified:",
  sum(combined_deg_summary$significant == "Yes"),
  "\n"
)

combined_deg_summary <- combined_deg_summary %>%
  left_join(gene_expression_df, by = "GENE_ID")

cat("Saving combined_deg_summary.csv... \n")
write.csv(
  combined_deg_summary,
  file.path(input_dir, "combined_deg_summary.csv"),
  row.names = FALSE
)

cat("Done!\n")


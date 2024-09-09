library(Matrix)
library(Seurat)
library(SeuratObject)
library(lme4)
library(openxlsx)
library(data.table)
library(parallel)

setwd("")

input_dir <- file.path(getwd(), "input")
output_dir <- file.path(getwd(), "output")

# Number of expPCs to use
expPC_num <- 14

# Initialize centralized log file
args <- commandArgs(trailingOnly = TRUE)
log_suffix <- args[1]
log_name <- sprintf("generate_model_coeffs_log_%s.txt", log_suffix)
log_file <- file.path(getwd(), log_name)
write("", file = log_file)

cat("Preparation complete. \n")

# ---------------------------------------------------------

file_path_onek1k <- file.path(input_dir, "OneK1K_preprocessed.rds")
onek1k <- readRDS(file_path_onek1k)

file_path_gene_snp_list <- file.path(input_dir, "gene_snp_list.rds")
gene_snp_list <- readRDS(file_path_gene_snp_list)

file_path_geno_mat <- file.path(input_dir, "genotype_matrix.csv")
genotype_matrix <- fread(file_path_geno_mat)

file_path_model_log <- file.path(input_dir, "model_selection_log.rds")
model_selection_log <- readRDS(file_path_model_log)

cat("Data loaded.\n")

# ---------------------------------------------------------

genes_list <- names(gene_snp_list)
counts_matrix <- onek1k@assays[["RNA"]]@counts

# Calculate update frequency (5% intervals)
n <- length(genes_list)
update_frequency <- ceiling(n * 0.05)

nUMI <- scale(log(colSums(counts_matrix)))
perc_mt <- scale(onek1k@meta.data[["percent.mt"]])
age <- scale(onek1k@meta.data[["age"]])
sex <- as.factor(onek1k@meta.data[["sex"]])
donor_id <- as.factor(onek1k@meta.data[["donor_id"]])
pool_num <- as.factor(onek1k@meta.data[["pool_number"]])
expPCs <- onek1k@reductions[["pca"]]@cell.embeddings

# ---------------------------------------------------------

# Function to convert SNP columns to factors
convert_snp_to_factors <- function(data, snps) {
  for (snp in snps) {
    data[[snp]] <- as.factor(data[[snp]])
  }
  return(data)
}


vif_func <-  function(lm_model) {
  # Extract the design matrix from the model, removing the intercept term
  design_matrix <- model.matrix(lm_model)
  if ("(Intercept)" %in% colnames(design_matrix)) {
    design_matrix <- design_matrix[, -which(colnames(design_matrix) == "(Intercept)")]
  }
  
  # Initialize a vector to store VIF values
  vif_values <- numeric(ncol(design_matrix))
  
  # Calculate VIF for each variable
  for (i in seq_along(vif_values)) {
    # Define the current variable as the dependent variable in the model
    dependent_var <- design_matrix[, i]
    independent_vars <- design_matrix[, -i]
    
    # Fit a linear model for the current variable
    model <- lm(dependent_var ~ independent_vars + 0)
    
    # Calculate the R-squared value
    r_squared <- summary(model)$r.squared
    
    # Calculate VIF
    vif_values[i] <- 1 / (1 - r_squared)
  }
  
  # Associate VIF values with their corresponding variable names
  names(vif_values) <- colnames(design_matrix)
  
  return(vif_values)
}


# Function to check and handle multicollinearity using VIF
check_multicollinearity <- function(data,
                                    fixed_effects,
                                    snp_columns,
                                    gene,
                                    threshold = 10,
                                    log_file) {
  all_effects <- c(fixed_effects, snp_columns)
  
  formula_string <- paste("expres ~", paste(all_effects, collapse = " + "))
  formula_obj <- as.formula(formula_string)
  
  vif_values <- vif_func(lm(formula_obj, data = data))
  
  while (any(vif_values > threshold)) {
    which_max <- which.max(vif_values)
    max_vif <- vif_values[which_max]
    
    if (max_vif <= threshold) {
      break
    }
    
    variable_to_remove <- names(vif_values)[which_max]
    
    # Check if the variable to remove is an SNP element and remove the digit suffix
    if (grepl("^SNP", variable_to_remove)) {
      variable_to_remove <- gsub("\\d+$", "", variable_to_remove)
    }
    
    message <- sprintf(
      "High multicollinearity detected for gene %s. Removing variable: %s (VIF = %.2f)\n",
      gene,
      variable_to_remove,
      max_vif
    )
    write(message, file = log_file, append = TRUE)
    
    # Remove all elements related to the SNP
    all_effects <- all_effects[!grepl(variable_to_remove, all_effects)]
    
    if (length(all_effects) == 0) {
      stop_message <- sprintf(
        "All variables removed due to multicollinearity. Cannot fit model for gene %s\n",
        gene
      )
      write(stop_message, file = log_file, append = TRUE)
      stop(stop_message)
    }
    
    formula_string <- paste("expres ~", paste(all_effects, collapse = " + "))
    formula_obj <- as.formula(formula_string)
    vif_values <- vif_func(lm(formula_obj, data = data))
  }
  
  final_fixed_effects <- setdiff(all_effects, snp_columns)
  final_snp_columns <- intersect(all_effects, snp_columns)
  
  return(
    list(
      final_fixed_effects = final_fixed_effects,
      final_snp_columns = final_snp_columns
    )
  )
}


# Function to generate model formulas (includes multicollinearity check)
generate_formulas <- function(data,
                              fixed_effects_formula,
                              snp_columns,
                              gene,
                              log_file) {
  fixed_effects <- unlist(strsplit(fixed_effects_formula, " \\+ "))
  
  # Check multicollinearity and get the final fixed effects and SNPs
  collinearity_results <- check_multicollinearity(data, fixed_effects, 
                                                  snp_columns, gene, 
                                                  log_file = log_file)
  
  final_fixed_effects <- collinearity_results$final_fixed_effects
  final_snp_columns <- collinearity_results$final_snp_columns
  
  # Construct formulas
  fixed_effects_formula <- paste(final_fixed_effects, collapse = " + ")
  null_formula <- paste("expres ~ (1 | donor_id) + (1 | pool_num) +", fixed_effects_formula)
  null_model_formula <- as.formula(null_formula)
  
  snp_formula <- paste(final_snp_columns, collapse = " + ")
  full_formula <- paste(null_formula, "+", snp_formula)
  full_model_formula <- as.formula(full_formula)
  
  # Interaction model formula
  interaction_terms <- paste(sapply(final_snp_columns, function(snp)
    paste(snp, ":", paste0("expPC", 1:expPC_num), collapse = " + ")), collapse = " + ")
  interaction_formula <- paste(full_formula, "+", interaction_terms)
  interaction_model_formula <- as.formula(interaction_formula)
  
  return(
    list(
      null_model_formula = null_model_formula,
      full_model_formula = full_model_formula,
      interaction_model_formula = interaction_model_formula
    )
  )
}

# ---------------------------------------------------------

cat("Starting gene processing loop... \n")

# Detect the number of cores allocated by SLURM
num_cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))
if (is.na(num_cores))
  num_cores <- detectCores() - 1

# Parallel setup
cl <- makeCluster(num_cores)

clusterEvalQ(cl, {
  library(Matrix)
  library(lme4)
  library(data.table)
})

clusterExport(
  cl,
  c(
    "genes_list",
    "counts_matrix",
    "nUMI",
    "perc_mt",
    "age",
    "sex",
    "expPCs",
    "donor_id",
    "pool_num",
    "genotype_matrix",
    "gene_snp_list",
    "expPC_num",
    "convert_snp_to_factors",
    "generate_formulas",
    "n",
    "output_dir",
    "log_file",
    "check_multicollinearity",
    "vif_func",
    "model_selection_log"
  )
)

tasks_per_core <- ceiling(length(genes_list) / length(cl))
clusterExport(cl, "tasks_per_core")
clusterEvalQ(cl, {
  options(tasks_completed = 0)
})

# Process each gene
process_gene <- function(i, assigned_genes, start_time) {
  this_gene <- genes_list[i]
  expres <- counts_matrix[this_gene, ]
  
  data_list <- list(
    expres = expres,
    nUMI = nUMI,
    perc_mt = perc_mt,
    age = age,
    sex = sex,
    donor_id = donor_id,
    pool_num = pool_num
  )
  
  for (j in 1:expPC_num) {
    data_list[[paste0("expPC", j)]] <- expPCs[, j]
  }
  
  data <- data.frame(data_list)
  
  # Process SNP columns
  snps <- gene_snp_list[[this_gene]]
  
  this_geno_data <- genotype_matrix[, c("IID", snps), with = FALSE]
  
  snp_prefix <- function(name) {
    paste0("SNP", gsub(":", ".", name))
  }
  
  colnames(this_geno_data)[-1] <- sapply(colnames(this_geno_data)[-1], snp_prefix)
  snps <- sapply(snps, snp_prefix)
  colnames(this_geno_data) <- as.character(colnames(this_geno_data))
  snps <- as.character(snps)
  
  data$row_index <- seq_len(nrow(data))
  data <- merge(
    data,
    this_geno_data,
    by.x = "donor_id",
    by.y = "IID",
    all.x = TRUE,
    all.y = FALSE
  )
  data <- data[order(data$row_index), ]
  data$row_index <- NULL
  
  data <- convert_snp_to_factors(data, snps)
  
  missing_rows <- data[!complete.cases(data), ]
  
  if (nrow(missing_rows) > 0) {
    missing_data_log_message <- sprintf("Missing data detected: %d rows affected\n", nrow(missing_rows))
    write(missing_data_log_message,
          file = log_file,
          append = TRUE)
  }
  
  data <- data[complete.cases(data), ]
  
  # ---------------------------------------------------------
  
  # Check multicollinearity and generate model formulas
  fixed_effects <- c("nUMI", "perc_mt", "age", "sex", paste0("expPC", 1:expPC_num))
  fixed_effects_formula <- paste(fixed_effects, collapse = " + ")
  
  snp_columns <- colnames(this_geno_data)[-1]
  
  tryCatch({
    formulas <- generate_formulas(data,
                                  fixed_effects_formula,
                                  snp_columns,
                                  this_gene,
                                  log_file)
  }, error = function(e) {
    message <- sprintf("Error generating formulas for gene %s: %s\n",
                       this_gene,
                       e$message)
    write(message, file = log_file, append = TRUE)
    NULL
  })
  
  if (is.null(formulas)) {
    error_message <- sprintf("Failed to generate formulas for gene %s due to error.\n", this_gene)
    write(error_message, file = log_file, append = TRUE)
  }
  
  # ---------------------------------------------------------
  
  model_type <- as.numeric(gsub(".*([0-9]+):.*", "\\1", model_selection_log[model_selection_log$gene == this_gene, "used_model"]))
  
  # Fit the interaction model and extract results
  result <- tryCatch({
    withCallingHandlers({
      
      if (model_type == 1) {
        model <- lme4::glmer(formula = formulas$null_model_formula,
                             data = data,
                             family = "poisson",
                             nAGQ = 0,
                             control = glmerControl(optimizer = "nloptwrap"))
      } else if (model_type == 2) {
        model <- lme4::glmer(formula = formulas$full_model_formula,
                             data = data,
                             family = "poisson",
                             nAGQ = 0,
                             control = glmerControl(optimizer = "nloptwrap"))
      } else if (model_type == 3) {
        model <- lme4::glmer(formula = formulas$interaction_model_formula,
                             data = data,
                             family = "poisson",
                             nAGQ = 0,
                             control = glmerControl(optimizer = "nloptwrap"))
      } else {
        error_message <- sprintf("Failed to select model for gene %s due to error.\n", this_gene)
        write(error_message, file = log_file, append = TRUE)
      }
      
      coefs <- summary(model)$coefficients
      pvals <- summary(model)$coefficients[, "Pr(>|z|)"]
      
      sig_coefs <- coefs[, "Pr(>|z|)"] < 0.05
      sig_num_sig_coefs <- sum(sig_coefs, na.rm = TRUE)
      
      if (any(sig_coefs)) {
        sig_indices <- which(sig_coefs)
        sig_ests <- coefs[sig_coefs, "Estimate"]
        sig_coefs <- rownames(coefs)[sig_coefs]
        
        sig_expPCs <- sig_coefs[grepl("expPC", sig_coefs) & !grepl(":", sig_coefs)]
        
        sig_SNPs <- sig_coefs[grepl("SNP", sig_coefs) & !grepl(":", sig_coefs)]
        
        if (model_type == 3) {
          sig_interactions <- sig_coefs[grepl(":", sig_coefs)]
          sig_interactions_indices <- sig_indices[grepl(":", sig_coefs)]
          sig_interactions_ests <- sig_ests[grepl(":", sig_coefs)]
          
          if (length(sig_interactions) > 0) {
            sig_interactions_dropdup <- unique(substring(sig_interactions, 1, 
                                                         nchar(sig_interactions) - 1))
          } else {
            sig_interactions_dropdup <- NA
          }
          
        } else {
          sig_interactions <- NA
          sig_interactions_indices <- NA
          sig_interactions_ests <- NA
          sig_interactions_dropdup <- NA
        }
        
      } else {
        sig_indices <- NA
        sig_ests <- NA
        sig_coefs <- NA
        
        sig_expPCs <- NA
        sig_SNPs <- NA
        
        sig_interactions <- NA
        sig_interactions_indices <- NA
        sig_interactions_ests <- NA
        sig_interactions_dropdup <- NA
      }
      
      list(
        gene = this_gene,
        model = model_type,
        
        coefs = coefs,
        pvals = pvals,
        sig_num_sig_coefs = sig_num_sig_coefs,
        
        sig_indices = sig_indices,
        sig_ests = sig_ests,
        sig_coefs = sig_coefs,
        sig_expPCs = sig_expPCs,
        sig_SNPs = sig_SNPs,
        
        sig_interactions = sig_interactions,
        sig_interactions_indices = sig_interactions_indices,
        sig_interactions_ests = sig_interactions_ests,
        sig_interactions_dropdup = sig_interactions_dropdup
      )
      
    }, warning = function(w) {
      message <- sprintf("Warning for gene %s: %s\n", this_gene, w$message)
      write(message, file = log_file, append = TRUE)
      invokeRestart("muffleWarning")
    })
  }, error = function(e) {
    message <- sprintf("Error in gene %s: %s\n", this_gene, e$message)
    write(message, file = log_file, append = TRUE)
    list(gene = this_gene,
         error = e$message,
         missing_rows = missing_rows)
  })
  
  options(tasks_completed = getOption("tasks_completed") + 1)
  tasks_completed <- getOption("tasks_completed")
  core_progress <- (tasks_completed / tasks_per_core) * 100
  progress_message <- sprintf(
    "Core %d: Processed gene %d of %d assigned (%.2f%% completed)\n",
    Sys.getpid(),
    i,
    n,
    core_progress
  )
  write(progress_message, file = log_file, append = TRUE)
  
  return(result)
}

results <- parLapply(cl, seq_along(genes_list), process_gene)

stopCluster(cl)

# ---------------------------------------------------------

cat("Combining results... \n\n")

final_results <- list()

for (i in seq_along(results)) {
  result <- results[[i]]
  this_gene <- result$gene
  
  if (!is.null(result$error)) {
    cat("Error for gene:", this_gene, "\n")
    cat("Error message:", result$error, "\n")
    saveRDS(result$data, file.path(output_dir, paste0("Error_data_", this_gene, ".rds")))
  } else {
    final_results[[this_gene]] <- result
  }
  
  progress_percentage <- round((i / n) * 100, 2)
  cat("Progress:", progress_percentage, "%\n")
}

# ---------------------------------------------------------

cat("\nSaving final results... \n")

save_data_to_rds <- function(data, file_name) {
  rds_file_name <- file.path(input_dir, paste0(file_name, ".rds"))
  saveRDS(data, rds_file_name, compress = "gzip")
  cat(paste("Saved", file_name, "to", rds_file_name, "\n"))
}

save_data_to_rds(final_results, "model_coeffs")

cat("Done!")

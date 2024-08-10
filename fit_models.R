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
log_name <- sprintf("fit_models_log_%s.txt", log_suffix)
log_file <- file.path(getwd(), log_name)
write("", file = log_file)

cat("Preparation complete. \n")

# ---------------------------------------------------------

file_path_onek1k <- file.path(input_dir, "OneK1K_modified0.rds")
onek1k <- readRDS(file_path_onek1k)

file_path_gene_snp_list <- file.path(input_dir, "gene_snp_list.rds")
gene_snp_list <- readRDS(file_path_gene_snp_list)

file_path_geno_mat <- file.path(input_dir, "genotype_matrix.csv")
genotype_matrix <- fread(file_path_geno_mat)

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
check_multicollinearity <- function(data, fixed_effects, snp_columns, 
                                    gene, threshold = 10, log_file) {
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
    
    message <- sprintf("High multicollinearity detected for gene %s. Removing variable: %s (VIF = %.2f)\n", 
                       gene, variable_to_remove, max_vif)
    write(message, file = log_file, append = TRUE)
    
    # Remove all elements related to the SNP
    all_effects <- all_effects[!grepl(variable_to_remove, all_effects)]
    
    if (length(all_effects) == 0) {
      stop_message <- sprintf("All variables removed due to multicollinearity. Cannot fit model for gene %s\n", gene)
      write(stop_message, file = log_file, append = TRUE)
      stop(stop_message)
    }
    
    formula_string <- paste("expres ~", paste(all_effects, collapse = " + "))
    formula_obj <- as.formula(formula_string)
    vif_values <- vif_func(lm(formula_obj, data = data))
  }
  
  final_fixed_effects <- setdiff(all_effects, snp_columns)
  final_snp_columns <- intersect(all_effects, snp_columns)
  
  return(list(final_fixed_effects = final_fixed_effects, final_snp_columns = final_snp_columns))
}


# Function to generate model formulas (includes multicollinearity check)
generate_formulas <- function(data, fixed_effects_formula, snp_columns, gene, log_file) {
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
  interaction_terms <- paste(sapply(final_snp_columns, 
                                    function(snp) paste(snp, ":", paste0("expPC", 1:expPC_num), 
                                                        collapse=" + ")), collapse=" + ")
  interaction_formula <- paste(full_formula, "+", interaction_terms)
  interaction_model_formula <- as.formula(interaction_formula)
  
  return(list(null_model_formula = null_model_formula, 
              full_model_formula = full_model_formula, 
              interaction_model_formula = interaction_model_formula))
}


# ---------------------------------------------------------

cat("Starting gene processing loop... \n")

# Detect the number of cores allocated by SLURM
num_cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))
if (is.na(num_cores)) num_cores <- detectCores() - 1

# Parallel setup
cl <- makeCluster(num_cores)

clusterEvalQ(cl, {
  library(Matrix)
  library(lme4)
  library(data.table)
})

# Export necessary data and functions to the cluster
clusterExport(cl, c("genes_list", "counts_matrix", "nUMI", "perc_mt", "age", 
                    "sex", "expPCs", "donor_id", "pool_num", "genotype_matrix", 
                    "gene_snp_list", "expPC_num", "convert_snp_to_factors", 
                    "generate_formulas", "n", "output_dir", "log_file",
                    "check_multicollinearity", "vif_func", "output_dir"))

tasks_per_core <- ceiling(length(genes_list) / length(cl))
clusterExport(cl, "tasks_per_core")
clusterEvalQ(cl, {
  options(tasks_completed = 0)
})

# Process each gene
process_gene <- function(i, assigned_genes, start_time) {
  this_gene <- genes_list[i]
  expres <- counts_matrix[this_gene,]
  
  data_list <- list(expres = expres, nUMI = nUMI, perc_mt = perc_mt, age = age, 
                    sex = sex, donor_id = donor_id, pool_num = pool_num)
  
  for (j in 1:expPC_num) {
    data_list[[paste0("expPC", j)]] <- expPCs[, j]
  }
  
  data <- data.frame(data_list)
  
  # Process SNP columns
  snps <- gene_snp_list[[this_gene]]
  geno_data <- genotype_matrix[, c("IID", snps), with = FALSE]
  
  snp_prefix <- function(name) {
    paste0("SNP", gsub(":", ".", name))
  }
  
  colnames(geno_data)[-1] <- sapply(colnames(geno_data)[-1], snp_prefix)
  snps <- sapply(snps, snp_prefix)
  colnames(geno_data) <- as.character(colnames(geno_data))
  snps <- as.character(snps)
  
  data <- merge(data, geno_data, 
                by.x = "donor_id", by.y = "IID", 
                all.x = TRUE, all.y = FALSE)
  
  data <- convert_snp_to_factors(data, snps)
  
  missing_rows <- data[!complete.cases(data), ]
  
  if (nrow(missing_rows) > 0) {
    missing_data_log_message <- sprintf("Missing data detected: %d rows affected\n", 
                                        nrow(missing_rows))
    write(missing_data_log_message, file = log_file, append = TRUE)
  }
  
  data <- data[complete.cases(data), ]
  
  # ---------------------------------------------------------
  
  # Generate model formulas (includes collinearity check)
  fixed_effects <- c("nUMI", "perc_mt", "age", "sex", 
                     paste0("expPC", 1:expPC_num))
  fixed_effects_formula <- paste(fixed_effects, collapse = " + ")
  
  snp_columns <- colnames(geno_data)[-1]
  
  tryCatch({
    formulas <- generate_formulas(data, fixed_effects_formula, 
                                  snp_columns, this_gene, log_file)
  }, error = function(e) {
    message <- sprintf("Error generating formulas for gene %s: %s\n", 
                       this_gene, e$message)
    write(message, file = log_file, append = TRUE)
    NULL
  })
  
  if (is.null(formulas)) {
    error_message <- sprintf("Failed to generate formulas for gene %s due to error.\n", 
                             this_gene)
    write(error_message, file = log_file, append = TRUE)
  } 
  
  # ---------------------------------------------------------
  
  # Fit models and extract fitted values
  result <- tryCatch({
    withCallingHandlers({
      null_model <- lme4::glmer(
        formula = formulas$null_model_formula, 
        data = data, family = "poisson", nAGQ = 0, 
        control = glmerControl(optimizer = "nloptwrap")
      )
      
      full_model <- lme4::glmer(
        formula = formulas$full_model_formula, 
        data = data, family = "poisson", nAGQ = 0, 
        control = glmerControl(optimizer = "nloptwrap")
      )
      
      interaction_model <- lme4::glmer(
        formula = formulas$interaction_model_formula, 
        data = data, family = "poisson", nAGQ = 0, 
        control = glmerControl(optimizer = "nloptwrap")
      )
      
      fitted_vals_null <- fitted(null_model)
      fitted_vals_full <- fitted(full_model)
      fitted_vals_interaction <- fitted(interaction_model)
      
      # Model comparison
      null_aic <- AIC(null_model)
      full_aic <- AIC(full_model)
      interaction_aic <- AIC(interaction_model)
      null_bic <- BIC(null_model)
      full_bic <- BIC(full_model)
      interaction_bic <- BIC(interaction_model)
      null_logLik <- logLik(null_model)
      full_logLik <- logLik(full_model)
      interaction_logLik <- logLik(interaction_model)
      
      model_comparison <- data.frame(
        gene = this_gene,
        null_aic = null_aic,
        full_aic = full_aic,
        interaction_aic = interaction_aic,
        null_bic = null_bic,
        full_bic = full_bic,
        interaction_bic = interaction_bic,
        null_logLik = as.numeric(null_logLik),
        full_logLik = as.numeric(full_logLik),
        interaction_logLik = as.numeric(interaction_logLik)
      )
      
      list(
        gene = this_gene,
        fitted_vals_null = fitted_vals_null,
        fitted_vals_full = fitted_vals_full,
        fitted_vals_interaction = fitted_vals_interaction,
        model_comparison = model_comparison,
        missing_rows = missing_rows,
        error = NULL
      )
    }, warning = function(w) {
      message <- sprintf("Warning for gene %s: %s\n", this_gene, w$message)
      write(message, file = log_file, append = TRUE)
      invokeRestart("muffleWarning")
    })
  }, error = function(e) {
    message <- sprintf("Error in gene %s: %s\n", this_gene, e$message)
    write(message, file = log_file, append = TRUE)
    list(
      gene = this_gene,
      error = e$message,
      missing_rows = missing_rows
    )
  })
  
  options(tasks_completed = getOption("tasks_completed") + 1)
  tasks_completed <- getOption("tasks_completed")
  core_progress <- (tasks_completed / tasks_per_core) * 100
  progress_message <- sprintf("Core %d: Processed gene %d of %d assigned (%.2f%% completed)\n", 
                              Sys.getpid(), i, n, core_progress)
  write(progress_message, file = log_file, append = TRUE)
  
  return(result)
}

results <- parLapply(cl, seq_along(genes_list), process_gene)

stopCluster(cl)

# ---------------------------------------------------------

cat("Combining results... \n")

n_genes <- length(genes_list)
n_columns_null <- ncol(results[[1]]$fitted_vals_null)
n_columns_full <- ncol(results[[1]]$fitted_vals_full)
n_columns_interaction <- ncol(results[[1]]$fitted_vals_interaction)

fitted_df_null <- matrix(NA, nrow = n_genes, ncol = n_columns_null)
fitted_df_full <- matrix(NA, nrow = n_genes, ncol = n_columns_full)
fitted_df_interaction <- matrix(NA, nrow = n_genes, ncol = n_columns_interaction)
model_comparison_df <- matrix(NA, nrow = n_genes, 
                              ncol = ncol(results[[1]]$model_comparison))

if (!exists("missing_data_df")) {
  missing_data_df <- data.frame()
}

progress_step <- floor(n_genes * 0.05)

for (i in seq_along(results)) {
  result <- results[[i]]
  this_gene <- result$gene
  
  if (!is.null(result$error)) {
    cat("Error for gene:", this_gene, "\n")
    cat("Error message:", result$error, "\n")
    saveRDS(result$data, 
            file.path(output_dir, paste0("Error_data_", this_gene, ".rds")))
  } else {
    if (nrow(result$missing_rows) > 0) {
      missing_data_df <- rbind(missing_data_df, result$missing_rows)
    }
    
    fitted_df_null[i, ] <- result$fitted_vals_null
    fitted_df_full[i, ] <- result$fitted_vals_full
    fitted_df_interaction[i, ] <- result$fitted_vals_interaction
    model_comparison_df[i, ] <- result$model_comparison
  }
  
  if (i %% progress_step == 0 || i == n_genes) {
    cat("Progress:", round((i / n_genes) * 100), "%\n")
  }
}

fitted_df_null <- as.data.frame(fitted_df_null)
fitted_df_full <- as.data.frame(fitted_df_full)
fitted_df_interaction <- as.data.frame(fitted_df_interaction)
model_comparison_df <- as.data.frame(model_comparison_df)

rownames(fitted_df_null) <- genes_list
rownames(fitted_df_full) <- genes_list
rownames(fitted_df_interaction) <- genes_list
rownames(model_comparison_df) <- genes_list

# ---------------------------------------------------------

cat("Saving data frames... \n")

wb <- createWorkbook()

addWorksheet(wb, "Fitted Values_Null")
writeData(wb, "Fitted Values_Null", fitted_df_null)

addWorksheet(wb, "Fitted Values_Full")
writeData(wb, "Fitted Values_Full", fitted_df_full)

addWorksheet(wb, "Fitted Values_Interaction")
writeData(wb, "Fitted Values_Interaction", fitted_df_interaction)

addWorksheet(wb, "Missing Data")
writeData(wb, "Missing Data", missing_data_df)

addWorksheet(wb, "Model Comparison")
writeData(wb, "Model Comparison", model_comparison_df)

saveWorkbook(wb, file.path(output_dir, "Fitted_Expres_Models.xlsx"), 
             overwrite = TRUE)

cat("\n========================================================\n")

cat("All finished!!!!\n")


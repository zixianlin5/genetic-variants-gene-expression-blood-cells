library(dplyr)
library(tibble)
library(ggplot2)
library(forcats)
library(gtools)
library(ggpattern)
library(data.table)
library(forcats)
library(gtools)
library(monocle3)

setwd("")

input_dir <- file.path(getwd(), "input")
output_dir <- file.path(getwd(), "output")

# ---------------------------------------------------------

file_path_model_selection <- file.path(input_dir, "model_selection_log.rds")
model_selection_log <- readRDS(file_path_model_selection)

model_dict <- split(model_selection_log$gene, model_selection_log$used_model)

file_path_gene_snp_df <- file.path(input_dir, "gene_snp_df.csv")
gene_snp_df <- read.csv(file_path_gene_snp_df) %>% distinct()

file_path_model_coeffs <- file.path(input_dir, "model_coeffs.rds")
model_coeffs <- readRDS(file_path_model_coeffs)

# ---------------------------------------------------------

## State-Dependent Gene Expression

df_sig_epc <- data.frame(
  gene = character(),
  model = integer(),
  num_sig_expPCs = integer(),
  sig_expPCs = character(),
  stringsAsFactors = FALSE
)

expPC_analysis <- list(counts = integer(), genes = list())

for (gene_data in model_coeffs) {
  gene_name <- gene_data$gene
  model_num <- gene_data$model
  
  sig_expPCs <- gene_data$sig_expPCs
  
  if (!is.null(gene_data$sig_interactions_dropdup)) {
    if (length(gene_data$sig_interactions_dropdup) > 0) {
      interaction_expPCs <- unlist(lapply(gene_data$sig_interactions_dropdup, function(interaction) {
        if (!is.na(interaction) &&
            is.character(interaction) && length(interaction) > 0) {
          strsplit(interaction, ":")[[1]][1]
        } else {
          NA  # Return NA if interaction is not valid
        }
      }))
      
      # Remove NAs from the resulting vector
      interaction_expPCs <- interaction_expPCs[!is.na(interaction_expPCs)]
      
      sig_expPCs <- unique(c(sig_expPCs, interaction_expPCs))
    }
  }
  
  sig_expPCs_str <- paste(sig_expPCs, collapse = ",")
  num_sig_expPCs <- length(sig_expPCs)
  
  df_sig_epc <- rbind(
    df_sig_epc,
    data.frame(
      gene = gene_name,
      model = model_num,
      num_sig_expPCs = num_sig_expPCs,
      sig_expPCs = sig_expPCs_str,
      stringsAsFactors = FALSE
    )
  )
  
  for (expPC in sig_expPCs) {
    if (!(expPC %in% names(expPC_analysis$counts))) {
      expPC_analysis$counts[[expPC]] <- 0
      expPC_analysis$genes[[expPC]] <- list()
    }
    expPC_analysis$counts[[expPC]] <- expPC_analysis$counts[[expPC]] + 1
    
    if (!is.null(expPC_analysis$genes[[expPC]])) {
      expPC_analysis$genes[[expPC]] <- c(expPC_analysis$genes[[expPC]], gene_name)
    } else {
      expPC_analysis$genes[[expPC]] <- list(gene_name)
    }
    names(expPC_analysis$genes[[expPC]]) <- expPC_analysis$genes[[expPC]]
  }
}

# ---------------------------------------------------------

# barplot epc_gene_count
expPC_counts_df <- as.data.frame(expPC_analysis$counts)
expPC_counts_df$ePC <- rownames(expPC_counts_df)
colnames(expPC_counts_df) <- c("gene_count", "ePC")

expPC_counts_df$ePC <- gsub("expPC", "ePC", expPC_counts_df$ePC)
expPC_counts_df$ePC <- factor(expPC_counts_df$ePC, levels = paste0("ePC", sort(as.numeric(
  gsub("ePC", "", expPC_counts_df$ePC)
))))


plot_epc_gene_count <- ggplot(expPC_counts_df, aes(x = ePC, y = gene_count)) +
  geom_bar_pattern(
    stat = "identity",
    width = 0.6,
    pattern = "stripe",           
    pattern_angle = 45,           
    pattern_density = 0.05,       
    pattern_spacing = 0.03,       
    fill = "white",               
    color = "black"               
  ) +
  geom_text(aes(label = gene_count), vjust = -0.5, size = 4) +
  xlab("ePC") +
  ylab("Number of Genes") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 12
    ),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    panel.grid.major = element_line(size = 0.8),
    panel.grid.minor = element_blank()
  )

output_path_epc <- file.path(output_dir, "gene_count_dist_across_ePCs.png")
ggsave(
  output_path_epc,
  plot = plot_epc_gene_count,
  width = 9,
  height = 4.8,
  dpi = 300
)


# barplot num_sig_expPCs
plot_num_sig_expPCs <- ggplot(df_sig_epc, aes(x = factor(num_sig_expPCs))) +
  geom_bar_pattern(
    stat = "count",
    pattern = "stripe",          
    pattern_angle = 45,         
    pattern_density = 0.05,      
    pattern_spacing = 0.03,     
    fill = "white",            
    color = "black"              
  ) +
  geom_text(
    stat = 'count',
    aes(label = ..count..),
    vjust = -0.5,
    size = 4
  ) +
  xlab("Number of ePCs Involved") +
  ylab("Count of Genes") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5,
      size = 12
    ),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    panel.grid.major = element_line(size = 0.8),
    panel.grid.minor = element_blank()
  )

output_path_num_sig <- file.path(output_dir, "num_sig_expPCs_distribution.png")
ggsave(
  output_path_num_sig,
  plot = plot_num_sig_expPCs,
  width = 7,
  height = 5,
  dpi = 300
)

# ---------------------------------------------------------

## Cell State-Specific Effects of Variants


file_path_gwas <- file.path(input_dir, "gwas_catalog_associations.tsv")
gwas_data <- fread(file_path_gwas,
                   header = TRUE,
                   sep = "\t",
                   quote = "") %>%
  filter(
    !grepl(
      "Phage display sequencing reveals that genetic, environmental, and intrinsic factors influence variation of human antibody epitope repertoire",
      `STUDY`,
      ignore.case = TRUE
    )
  )

df_sig_interactions <- do.call(rbind, lapply(model_coeffs, function(x) {
  
  if (x$model != 3) {
    return(NULL)
  }
  
  interactions <- x$sig_interactions_dropdup
  
  if (length(interactions) == 0) {
    return(NULL)
  }
  
  snp <- sub(".*:", "", interactions)
  epc <- sub(":.*", "", interactions)
  
  gene_rep <- rep(x$gene, length(interactions))
  
  result <- data.frame(
    gene = gene_rep,
    sig_interactions = interactions,
    SNP = snp,
    ePC = epc,
    stringsAsFactors = FALSE
  )
  
  return(na.omit(result))
}))

length(unique(df_sig_interactions$SNP))

gene_snp_df_filtered <- gene_snp_df %>%
  filter(GENE_ID %in% model_dict$`Model 3: Model with interactions`)
length(unique(gene_snp_df_filtered$SNPID))
rm(gene_snp_df_filtered)

df_sig_interactions$ePC <- gsub("expPC", "ePC", df_sig_interactions$ePC)

snp_epc_count <- df_sig_interactions %>%
  group_by(SNP) %>%
  summarise(ePC_count = n_distinct(ePC)) %>%
  count(ePC_count) %>%
  rename(snp_count = n)

plot_snp_epc_distribution <- ggplot(snp_epc_count, 
                                    aes(x = ePC_count, y = snp_count)) +
  geom_bar_pattern(
    stat = "identity",
    pattern = "circle",   
    pattern_density = 0.08,   
    pattern_spacing = 0.02,   
    fill = "white",          
    color = "black"          
  ) +
  geom_text(aes(label = snp_count), vjust = -0.5, size = 5) +
  xlab("Number of ePCs Each SNP Interacted With") +
  ylab("Count of SNPs") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    panel.grid.major = element_line(size = 0.8),
    panel.grid.minor = element_blank()
  )

output_path_snp <- file.path(output_dir, "distribution_of_SNPs_by_number_of_ePC_interactions.png")
ggsave(
  output_path_snp,
  plot = plot_snp_epc_distribution,
  width = 7,
  height = 6,
  dpi = 300
)

epc_count <- df_sig_interactions %>%
  count(ePC) %>%
  rename(epc_count = n)

epc_count$ePC <- factor(epc_count$ePC, levels = paste0("ePC", 1:max(as.numeric(gsub("ePC", "", epc_count$ePC)))))

plot_epc_count_distribution <- ggplot(epc_count, 
                                      aes(x = ePC, y = epc_count)) +
  geom_bar_pattern(
    stat = "identity",
    pattern = "circle", 
    width = 0.6,
    pattern_density = 0.08,   
    pattern_spacing = 0.02,   
    fill = "white",          
    color = "black"          
  ) +
  geom_text(aes(label = epc_count), vjust = -0.5, size = 5) +
  xlab("ePC") +
  ylab("Number of SNPs") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    panel.grid.major = element_line(size = 0.8),
    panel.grid.minor = element_blank()
  )

output_path_epc <- file.path(output_dir, "distribution_of_ePCs_by_occurrences.png")
ggsave(
  output_path_epc,
  plot = plot_epc_count_distribution,
  width = 9,
  height = 5.5,
  dpi = 300
)

# ---------------------------------------------------------

perform_fisher_test <- function(df_sig_interactions, gene_snp_df, gwas_data, 
                                  blood_traits, specific_ePC = NULL) {
  
  gwas_blood_snps <- gwas_data %>%
    filter(grepl(paste(blood_traits, collapse="|"), `DISEASE/TRAIT`, ignore.case=TRUE) |
             grepl(paste(blood_traits, collapse="|"), `STUDY`, ignore.case=TRUE)) %>%
    dplyr::select(SNPS, STUDY, `DISEASE/TRAIT`) %>%
    distinct()
  
  reference_df <- gwas_blood_snps %>%
    dplyr::select(STUDY, `DISEASE/TRAIT`) %>%
    distinct()
  
  if (!is.null(specific_ePC)) {
    # filter SNPs that interact with any ePC in the list
    interacting_snps <- df_sig_interactions %>%
      filter(ePC %in% specific_ePC) %>%
      dplyr::select(SNP) %>%
      distinct()
  } else {
    interacting_snps <- df_sig_interactions %>%
      dplyr::select(SNP) %>%
      distinct()
  }
  
  interacting_snps <- interacting_snps %>%
    mutate(SNPID = sub("^SNP", "", SNP),
           SNPID = gsub("\\.", ":", SNPID)) %>%
    left_join(gene_snp_df, by = c("SNPID" = "SNPID")) %>%
    dplyr::select(RSID) %>%
    distinct() %>%
    filter(RSID %in% gwas_data$SNPS)
  
  all_snps <- gene_snp_df %>%
    dplyr::select(RSID) %>%
    distinct() %>%
    filter(RSID %in% gwas_data$SNPS)
  
  # Construct the contingency table
  snp_interacts_with_epc_and_is_gwas <- interacting_snps %>%
    filter(RSID %in% gwas_blood_snps$SNPS) %>%
    nrow()
  
  snp_interacts_with_epc_and_is_not_gwas <- interacting_snps %>%
    filter(!RSID %in% gwas_blood_snps$SNPS) %>%
    nrow()
  
  snp_does_not_interact_with_epc_and_is_gwas <- all_snps %>%
    filter(RSID %in% gwas_blood_snps$SNPS & !RSID %in% interacting_snps$RSID) %>%
    nrow()
  
  snp_does_not_interact_with_epc_and_is_not_gwas <- all_snps %>%
    filter(!RSID %in% gwas_blood_snps$SNPS & !RSID %in% interacting_snps$RSID) %>%
    nrow()
  
  # Create contingency table
  contingency_table <- matrix(c(
    snp_interacts_with_epc_and_is_gwas,
    snp_interacts_with_epc_and_is_not_gwas,
    snp_does_not_interact_with_epc_and_is_gwas,
    snp_does_not_interact_with_epc_and_is_not_gwas
  ), nrow = 2, byrow = TRUE)
  
  colnames(contingency_table) <- c("blood_GWAS_SNP", "Not_blood_GWAS_SNP")
  rownames(contingency_table) <- c("Interacts_with_ePC", "Does_not_interact_with_ePC")
  
  # Perform Fisher's Exact Test or Chi-Square Test
  if (any(contingency_table < 5)) {
    test_result <- fisher.test(contingency_table)
  } else {
    test_result <- chisq.test(contingency_table)
  }
  
  list(
    contingency_table = contingency_table,
    test_result = test_result,
    reference_df = reference_df
  )
}

keywords <- c("agranulocytosis", "anemia", "antibodies", "antibody", "antigen", 
              "autoimmune", "b cell", "blood cell", "blood", "bone marrow",  
              "coagulation", "corpuscular", "cytokine", "cytokines", "clotting",
              "dendritic cell", "duffy", "erythrocyte", "eosinophil", 
              "eosinophils", "factor v", "factor viii", "g6pd", "granulocyte", 
              "hematology", "hematopoietic", "hemoglobin", "hemoglobinemic", 
              "hemolysis", "hemolytic", "hemophilia", "immune", 
              "immunity", "immunology", "inflammation", "itp", "leukemia", 
              "leukemias", "leukocyte", "lymphatic", "lymphocyte", 
              "lymphocytes", "lymphoma", "lymphomas", "macrophage", "monocyte", 
              "myeloid", "myeloma", "neutropenia", "neutrophil", "phagocyte",
              "platelet", "platelets", "plasma", "plasmas", "red cell", 
              "red cells", "reticulocyte", "reticulocytes", "serum", "sickle", 
              "sickles", "spherocytosis", "spherocytoses", "spleen", 
              "stem cell", "t cell", "thalassemia", "thalassemias", 
              "thrombocyte", "thrombocythemia", "thrombocythemias", 
              "thrombocytopenia", "thrombocytopenias", "thrombocytopenic", 
              "thrombosis", "thromboses", "von willebrand", "white cell", 
              "white cells")

fisher_result <- perform_fisher_test(df_sig_interactions, gene_snp_df, 
                                       gwas_data, keywords)

# ---------------------------------------------------------

## Inferred Cell Differentiation Pathway

file_path_cds <- file.path(input_dir, "cds_most_representative.rds")
cds <- readRDS(file_path_cds)

file_path_metadata <- file.path(input_dir, "visual_data.csv")
metadata <- read.csv(file_path_metadata)

pseudotime <- cds@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]

pseudotime_df <- data.frame(barcode = names(pseudotime), pseudotime = pseudotime)
metadata_subset <- metadata[, c("cell_name", paste0("PC_", 1:14))]
merged_df <- merge(pseudotime_df, metadata_subset, by.x = "barcode", by.y = "cell_name")
rm(metadata_subset, pseudotime_df)

for (i in 1:14) {
  pc_column <- paste0("PC_", i)
  cds@colData@listData[[pc_column]] <- merged_df[[pc_column]]
}

# calculate correlations between pseudotime and PCs
correlations <- list()

for (i in 1:14) {
  pc_column <- paste0("PC_", i)
  cor_value <- cor(merged_df$pseudotime, merged_df[[pc_column]], method = "pearson")
  correlations[[pc_column]] <- cor_value
}

# linear combine PCs to fit pseudotime
lm_model <- lm(pseudotime ~ PC_1 + PC_2 + PC_3 + PC_4 + PC_5 + PC_6 + PC_7 + 
                 PC_8 + PC_9 + PC_10 + PC_11 + PC_12 + PC_13 + PC_14, data = merged_df)

summary(lm_model)

coefficients <- lm_model$coefficients
print(coefficients)

predicted_pseudotime <- predict(lm_model, newdata = merged_df)

correlation_with_predicted <- cor(predicted_pseudotime, merged_df$pseudotime)
print(paste("Correlation between predicted and actual pseudotime:", correlation_with_predicted))

# ---------------------------------------------------------

library(MASS)

dens <- kde2d(merged_df$pseudotime, predicted_pseudotime, n = 200)

output_file <- file.path(output_dir, "pseudotime_density_plot.png")

png(filename = output_file, width = 1600, height = 1200, res = 200)

filled.contour(dens, 
               xlab = "Actual Pseudotime", 
               ylab = "Predicted Pseudotime", 
               color.palette = colorRampPalette(c("white", "blue", "purple")),
               plot.axes = {
                 axis(1)
                 axis(2)
                 reg_model <- lm(predicted_pseudotime ~ merged_df$pseudotime)
                 abline(reg_model, col = "red", lwd = 2)},
               key.axes = {
                 axis(4)
                 mtext("Point Density", side = 1, line = 1)  
               })

dev.off()
print(paste("Image saved to:", output_file))

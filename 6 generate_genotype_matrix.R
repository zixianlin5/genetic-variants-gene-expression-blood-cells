library(data.table)

setwd("")

input_dir <- file.path(getwd(), "input")
output_dir <- file.path(getwd(), "output")

file_path_gene_snp_list <- file.path(input_dir, "gene_snp_list.rds")
gene_snp_list <- readRDS(file_path_gene_snp_list)

# ---------------------------------------------------------

unique_snps <- unique(unlist(gene_snp_list))

extract_chromosome <- function(snp) {
  #extract the chromosome number from the SNP position
  return(strsplit(snp, ":")[[1]][1])
}

# Group unique_snps by chromosome number
grouped_snps <- split(unique_snps, sapply(unique_snps, extract_chromosome))
keys_list <- names(grouped_snps)

# ---------------------------------------------------------


# Function to process the reverse SNP encoding
flip_encoding <- function(genotype) {
  return((genotype - 1) * -1 + 1)
}

for (CHR in keys_list) {
  cat(paste0("Processing CHR", CHR, "...\n"))
  
  plink_file_name <- paste0("filter_vcf_r08/chr", CHR, "_recoded.raw")
  plink_path <- file.path(input_dir, plink_file_name)
  genotype_data <- fread(plink_path)
  
  snp_columns <- grouped_snps[[as.character(CHR)]]
  
  existing_snp_columns <- snp_columns[snp_columns %in% colnames(genotype_data)]
  missing_snp_columns <- snp_columns[!(snp_columns %in% colnames(genotype_data))]
  
  # Extract existing SNP columns
  snp_data <- genotype_data[, c("IID", existing_snp_columns), with = FALSE]
  
  missing_snps_not_found <- list()
  
  # Preprocess .bim file to create a lookup table
  bim_file_name <- paste0("filter_vcf_r08/chr", CHR, ".bim")
  bim_path <- file.path(input_dir, bim_file_name)
  bim_data <- fread(bim_path, header = FALSE)
  bim_lookup <- setNames(bim_data$V5, bim_data$V2)
  
  # Determine all reverse SNPs needed
  reverse_snp_list <- sapply(missing_snp_columns, function(snp) {
    snp_prefix <- strsplit(snp, "_")[[1]][1]
    ref_allele <- bim_lookup[snp_prefix]
    if (!is.null(ref_allele)) {
      return(paste0(snp_prefix, "_", ref_allele))
    } else {
      return(NA)
    }
  })
  reverse_snp_list <- reverse_snp_list[!is.na(reverse_snp_list)]
  
  # Extract all reverse SNP columns at once
  reverse_snps_data <- genotype_data[, c("IID", reverse_snp_list), with = FALSE, nomatch = 0]
  
  cat("  Processing missing SNP columns...\n")
  # Process missing SNP columns
  if (length(missing_snp_columns) > 0) {
    for (snp in missing_snp_columns) {
      snp_prefix <- strsplit(snp, "_")[[1]][1]
      
      # Get reference allele from lookup table
      ref_allele <- bim_lookup[snp_prefix]
      
      if (!is.null(ref_allele)) {
        reverse_snp <- paste0(snp_prefix, "_", ref_allele)
        
        if (reverse_snp %in% colnames(reverse_snps_data)) {
          temp_data <- reverse_snps_data[, .(IID, temp_snp = get(reverse_snp))]
          
          # Flip 0 and 2 encoding
          temp_data[, temp_snp := flip_encoding(temp_snp)]
          
          setnames(temp_data, "temp_snp", snp)
          snp_data <- merge(snp_data,
                            temp_data,
                            by = "IID",
                            all.x = TRUE)
          
        } else {
          # Store missing SNP not found in genotype_data
          missing_snps_not_found <- c(missing_snps_not_found, snp)
          cat(paste0(
            "  reverse_snp not found in genotype_data: ",
            reverse_snp,
            "\n"
          ))
        }
        
      } else {
        cat(paste0("  snp_prefix not found in bim_lookup: ", snp_prefix, "\n"))
      }
    }
  }
  
  if (length(missing_snps_not_found) > 0) {
    cat(paste0(
      "[!] The following SNPs were not found in genotype_data for CHR",
      CHR,
      ":\n"
    ))
    print(missing_snps_not_found)
  }
  
  # Merge snp_data into genotype_matrix
  if (CHR == 1) {
    genotype_matrix <- snp_data
  } else {
    genotype_matrix <- merge(genotype_matrix, snp_data, by = "IID", all = TRUE)
  }
  
  cat("--------------------------------\n")
  
}

cat("Saving genotype_matrix.csv... \n")

output_file_path <- file.path(input_dir, "genotype_matrix.csv")
fwrite(genotype_matrix, file = output_file_path)

cat("Done!")


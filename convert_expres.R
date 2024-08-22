library(arrow)

setwd("")

input_dir <- file.path(getwd(), "input")
output_dir <- file.path(getwd(), "output")

file_path_adjusted_expres <- file.path(input_dir, "adjusted_expres_df.rds")
adjusted_expres_df <- readRDS(file_path_adjusted_expres)

adjusted_expres_df$rownames <- rownames(adjusted_expres_df)

file_path_feather <- file.path(input_dir, "adjusted_expres_df.feather")
write_feather(adjusted_expres_df, file_path_feather)

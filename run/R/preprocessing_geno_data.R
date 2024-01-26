library(targets)

project_path <- "/proj"
core_paths <- lapply(jsonlite::read_json(sprintf("%s/results/core_paths.json", project_path)), function(x) x[[1]])
data_paths <- list(
  "imputed_matrix" = sprintf("%s/ext_dir/to_vcf_and_imputation/results/kibreed_gmat/filtered/kibreed_gmat.m_m_0.5_mat.qs", project_path),
  "zucht_transition" = sprintf("%s/ext_dir/KIBREED/source_data/BigData_yusheng_v1/zucht_transition.qs", project_path),
  "RD_mat" = sprintf("%s/ext_dir/core_marker_dataset/results/big_mat_RD.qs", project_path)
)

# Define packages needed and default storage format for intermittent data

tar_option_set(packages = c("dplyr", "readr", "ggplot2", 
                            "tidyr", "reshape2", "stringr",
                            "lubridate", "hablar", "readxl",
                            "AGHmatrix",
                            "qs"),
               format = "qs")

tar_source(sprintf("%s/src/R/%s", project_path, "fun_preprocessing_geno_data.R"))

# Create pipeline to process data
preprocessed_phenodata <- tar_read(preprocessing_environ_data, 
                                   store = sprintf("%s/store/preprocessing_environ_data", project_path))[["combined_data_with_env"]]

# Create pipeline to process data

list(
  tar_target(
    name = preprocessing_geno_data,
    command = preprocess_geno_data(existing_data = preprocessed_phenodata, 
                                   paths = data_paths,
                                   log_at = sprintf("%s/%s", core_paths[["logs_R"]], "preprocessing_geno_data"))
  ),
  tar_target(
    name = add_hybrid_data,
    command = generate_hybrid_data(data = preprocessing_geno_data)
  ),
  tar_target(
    name = kinship_data,
    command = generate_kinships(data = add_hybrid_data,
                                write_at = sprintf("%s/%s", core_paths[["results_R"]], "preprocessing_geno_data"))
  ),
  tar_target(
    name = deduplicate_genotypes,
    command = add_deduplication_info(existing_data = preprocessing_geno_data,
                                     data = add_hybrid_data,
                                     paths = data_paths)
  )
)


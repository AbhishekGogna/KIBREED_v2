library(targets)

project_path <- "/proj"
core_paths <- lapply(jsonlite::read_json(sprintf("%s/results/core_paths.json", project_path)), function(x) x[[1]])

# Define packages needed and default storage format for intermittent data

tar_option_set(packages = c("dplyr", "readr", "ggplot2", "tidyr", "reshape2",
                            "hablar", "asreml", 
                            "foreach", "doParallel", "parallel"),
               format = "qs")

tar_source(sprintf("%s/src/R/%s", project_path, "fun_phenotypic_data_kws_processing.R"))

# Create pipeline to process data

list(
  tar_target(
    name = kws_2020_data,
    command = get_pheno_data(path = sprintf("%s/ext_dir/KIBREED/source_data/BigData_pheno_moritz_20_06_2021/", project_path))
  ),
  tar_target(
    name = kws_2021_data,
    command = get_pheno_data(path = sprintf("%s/ext_dir/KIBREED/source_data/BigData_combined_moritz_10_02_2022/211111_I_KWS-2021-QC/data/01 Parse/Read phenotypic data/", project_path))
  ),
  tar_target(
    name = BLUEs_kws_2020_data,
    command = derive_BLUEs(kws_2020_data,
                           log_at = sprintf("%s/%s", core_paths[["logs_R"]], "phenotypic_data_kws_processing"),
                           tag = "2020")
  ),
  tar_target(
    name = BLUEs_kws_2021_data,
    command = derive_BLUEs(kws_2021_data,
                           log_at = sprintf("%s/%s", core_paths[["logs_R"]], "phenotypic_data_kws_processing"),
                           tag = "2021")
  ),
  tar_target(
    name = kws_2020_plots,
    command = produce_plots(raw_data = kws_2020_data, processed_data = BLUEs_kws_2020_data,
                            log_at = sprintf("%s/%s", core_paths[["logs_R"]], "phenotypic_data_kws_processing"),
                            tag = "2020")
  ),
  tar_target(
    name = kws_2021_plots,
    command = produce_plots(raw_data = kws_2021_data, processed_data = BLUEs_kws_2021_data,
                            log_at = sprintf("%s/%s", core_paths[["logs_R"]], "phenotypic_data_kws_processing"),
                            tag = "2021")
  )
)
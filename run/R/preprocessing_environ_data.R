library(targets)

project_path <- "/proj"
core_paths <- lapply(jsonlite::read_json(sprintf("%s/results/core_paths.json", project_path)), function(x) x[[1]])
data_paths <- list(
  "coord_raw" = sprintf("%s/ext_dir/KIBREED/source_data/coord_data_v4.csv", project_path),
  "env_data_sci_adv" = sprintf("%s/ext_dir/KIBREED/source_data/V2_weather_data_kibreed.tabular", project_path),
  "env_data_kws_2020" = sprintf("%s/ext_dir/KIBREED/source_data/kws_locations_weather_data.tabular", project_path),
  "env_data_kws_2021" = sprintf("%s/ext_dir/KIBREED/source_data/weather_2020_2021_val.tabular", project_path),
  "metadata_2020" = sprintf("%s/ext_dir/KIBREED/source_data/BigData_combined_moritz_09_07_2021/01 Parse/kws_2020_2/experiments.csv", project_path),
  "metadata_2021" = sprintf("%s/ext_dir/KIBREED/source_data/BigData_combined_moritz_10_02_2022/211111_I_KWS-2021-QC/data/01 Parse/Read phenotypic data/experiments.csv", project_path)
)

# Define packages needed and default storage format for intermittent data

tar_option_set(packages = c("dplyr", "readr", "ggplot2", 
                            "tidyr", "reshape2", "stringr",
                            "lubridate", "hablar", "readxl"),
               format = "qs")

tar_source(sprintf("%s/src/R/%s", project_path, "fun_preprocessing_environ_data.R"))

# Create pipeline to process data
preprocessed_phenodata <- tar_read(preprocessing_phenodata, 
                                   store = sprintf("%s/store/preprocessing_phenodata", project_path)) # can be done with yaml file?

list(
  tar_target(
    name = preprocessing_environ_data,
    command = preprocess_environ_data(existing_data = preprocessed_phenodata,
                                      paths = data_paths,
                                      log_at = sprintf("%s/%s", core_paths[["logs_R"]], "preprocessing_environ_data"))
  )
)
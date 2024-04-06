library(targets)

project_path <- "/proj"
core_paths <- lapply(jsonlite::read_json(sprintf("%s/results/core_paths.json", project_path)), function(x) x[[1]])

# Define packages needed and default storage format for intermittent data

tar_option_set(packages = c("dplyr", "readr", "ggplot2", 
                            "tidyr", "reshape2", "stringr",
                            "lubridate", "hablar", "readxl",
                            "qs", "feather"),
               format = "qs")

run_name <- "process_cgm_data"

tar_source(sprintf("/%s/src/R/fun_%s.R", project_path, run_name))

# Load exiting data

data <- list(
  "pheno_data" = tar_read(KIBREED_data_full , 
                          store = sprintf("%s/store/KIBREED_data_generation", project_path))[["BLUES_within_env"]]
)
# i actually sent the data - /proj/ext_dir/KIBREED/results_plots/export_for_ZALF/phase_4_yield_and_heading_date.txt, but this had missing values in BLUES_dt column. They removed it and merged the data with management data to get latlongs. Finally they kept only rows uniquely identified by latlong, year and geno

#data_paths <- list("param_metadata" = sprintf("%s/ext_dir/KIBREED/source_data/ZALF_Bonn/16_06_23_data/Parameters_code.csv", project_path),
#                   "param_data" = sprintf("%s/ext_dir/KIBREED/source_data/ZALF_Bonn/16_06_23_data/CombiPar.csv", project_path),
#                   "pheno_data_recieved" = sprintf("%s/ext_dir/KIBREED/source_data/ZALF_Bonn/16_06_23_data/Allresults.csv", project_path))

data_paths <- list("param_metadata" = sprintf("%s/ext_dir/KIBREED/source_data/ZALF_Bonn/17_07_23_data/Parameters_code.csv", project_path),
                   "param_data" = sprintf("%s/ext_dir/KIBREED/source_data/ZALF_Bonn/17_07_23_data/CombiPar.csv", project_path),
                   "pheno_data_used" = sprintf("%s/ext_dir/KIBREED/source_data/ZALF_Bonn/17_07_23_data/GenotypesNewAllupdate.csv", project_path),
                   "pheno_data_recieved" = sprintf("%s/ext_dir/KIBREED/source_data/ZALF_Bonn/17_07_23_data/Allresults.csv", project_path))

# define pipeline
list(
  tar_target(
    name = processed_cgm_data,
    command = process_cgm_data(existing_data = data, 
                               paths = data_paths,
                               write_path_for_R = sprintf("%s/%s", core_paths[["results_R"]], run_name),
                               write_path_for_Py = sprintf("%s/%s", core_paths[["results_Py"]], run_name),
                               log_at = sprintf("%s/%s", core_paths[["logs_R"]], run_name),
                               tmp_at = sprintf("%s/%s", core_paths[["tmp_data_R"]], run_name))
  ),
  tar_target(
    name = BLUES_within_env_cgm,
    command = write_cgm_data(data = processed_cgm_data)
  )
)

#geno <- BLUES_within_env_cgm$BLUEs_within_env_cgm %>% distinct(Env, connect_geno) %>%
#  mutate(present = 1) %>%
#  pivot_wider(id_cols = "connect_geno", names_from = "Env", values_from = "present") %>%
#  rowwise() %>%
#  mutate(freq = sum(c_across(where(is.numeric)), na.rm = TRUE))
#table(cut(geno$freq, breaks = c(0, 4, 8, 16, 32, 64))) / nrow(geno)
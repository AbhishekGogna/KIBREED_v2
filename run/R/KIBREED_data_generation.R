library(targets)

project_path <- "/proj"
core_paths <- lapply(jsonlite::read_json(sprintf("%s/results/core_paths.json", project_path)), function(x) x[[1]])

# Define packages needed and default storage format for intermittent data

tar_option_set(packages = c("dplyr", "readr", "ggplot2", "ggpubr",
                            "tidyr", "reshape2", "stringr",
                            "RColorBrewer", "dendextend", "circlize",
                            "grid", "patchwork", "ggmap",
                            "lubridate", "hablar", "readxl",
                            "foreach", "doParallel", "BGLR",
                            "asreml", "qs", "feather"),
               format = "qs")

tar_source(sprintf("%s/src/R/%s", project_path, "fun_KIBREED_data_generation.R"))

# Load existing data

data <- list(
  "pheno_data" = tar_read(deduplicate_genotypes, 
                         store = sprintf("%s/store/preprocessing_geno_data", project_path))[["data_mod"]],
  
  "geno_data" = tar_read(add_hybrid_data, 
                        store = sprintf("%s/store/preprocessing_geno_data", project_path))[c("GNdata_comb", "GNdata_comb_2")],
  
  "env_data"  = tar_read(preprocessing_environ_data, 
                        store = sprintf("%s/store/preprocessing_environ_data", project_path)),
  
  "meta_geno" = tar_read(add_hybrid_data, 
                        store = sprintf("%s/store/preprocessing_geno_data", project_path))[["names_change_consolidated"]]
  )

hd_data_overview <- qs::qread("/proj/ext_dir/KIBREED/results_plots/HD_pred/overview.qs")
hd_data_overview$res_at <- gsub("~/KIBREED", "/proj/ext_dir/KIBREED", hd_data_overview$res_at)

# Define pipeline

list(
  tar_target(
    name = predicted_HD,
    command = predict_HD_incpl_env(existing_data = data,
                                   log_at = sprintf("%s/%s", core_paths[["logs_R"]], "KIBREED_data_generation"),
                                   tmp_data = sprintf("%s/%s", core_paths[["tmp_data_R"]], "KIBREED_data_generation"))
  ), # 18min
  tar_target(
    name = KIBREED_BLUEs_acr_env,
    command = generate_KIBREED_blues(data = predicted_HD)
  ), # 5 hours
  tar_target(
    name = KIBREED_data_full,
    command = write_KIBREED_data_full(
      existing_data = data,
      data = KIBREED_BLUEs_acr_env,
      write_path_for_R = sprintf("%s/%s", core_paths[["results_R"]], "KIBREED_data_generation"),
      write_path_for_Py = sprintf("%s/%s", core_paths[["results_Py"]], "KIBREED_data_generation")
    ) # 5 minutes
  ),
  tar_target(
    name = overviews,
    command = make_plots_and_tables(
      data = KIBREED_data_full,
      pco_data = tar_read(overviews, store = "/proj/store/feature_importance")[["pco_data"]],
      pco_eig = tar_read(overviews, store = "/proj/store/feature_importance")[["pc_RD"]][["eig"]],
      write_path_for_R = sprintf("%s/%s", core_paths[["results_R"]], "KIBREED_data_generation"),
      write_path_for_Py = sprintf("%s/%s", core_paths[["results_Py"]], "KIBREED_data_generation")
    )
  )
)

## export efraim
#data <- tar_read(KIBREED_data_full)[["BLUES_within_env"]] %>%
#  distinct(Series, Project, Exp, Loc, Year, Site, latlong) %>% 
#  filter(Project %in% c("Hywheat", "Zucht"))
#qsave(data, "/proj/results/R/KIBREED_data_generation/28_02_24_env_info_Hy_Zu.qs")



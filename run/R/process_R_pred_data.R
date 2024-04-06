library(targets)

project_path <- "/proj"
core_paths <- lapply(jsonlite::read_json(sprintf("%s/results/core_paths.json", project_path)), function(x) x[[1]])

# Define packages needed and default storage format for intermittent data

tar_option_set(packages = c("dplyr", "readr", "ggplot2", "ggpubr", 
                            "tidyr", "reshape2", "stringr",
                            "lubridate", "hablar", "readxl",
                            "foreach", "doParallel",
                            "Metrics", "BBmisc",
                            "qs", "feather"),
               format = "qs")

run_name <- "process_R_pred_data"

tar_source(sprintf("/%s/src/R/fun_%s.R", project_path, run_name))

# Load exiting data

data <- list(
  "cv_acr_5f" = tar_read(run_scripts_acr_5f , 
                          store = sprintf("%s/store/generate_prediction_data", project_path))
  , "cv_acr_sce" = tar_read(run_scripts_acr_sce , 
                          store = sprintf("%s/store/generate_prediction_data", project_path))
  #, "cv_acr_str" = tar_read(run_scripts_acr_str , 
  #                        store = sprintf("%s/store/generate_prediction_data", project_path))
  , "cv_wtn_tra" = tar_read(run_scripts_wtn_tra , 
                          store = sprintf("%s/store/generate_prediction_data", project_path))
  , "cv_wtn_LoO" = tar_read(run_scripts_wtn_LoO , 
                          store = sprintf("%s/store/generate_prediction_data", project_path))
  , "cv_wtn_cvL" = tar_read(run_scripts_wtn_cvL , 
                          store = sprintf("%s/store/generate_prediction_data", project_path))
  , "pred_paths_py" = jsonlite::read_json(sprintf("%s/create_slurm_scripts/pred_dirs_paths.json", core_paths[["results_Py"]]))
  , "data_bahareh_at" = sprintf("%s/ext_dir/KIBREED/source_data/ZALF_Bonn/21_11_23_data", project_path)
)

# define pipeline
list(
  tar_target(
    name = pred_file_paths,
    command = get_pred_file_paths(existing_data = data, 
                                  log_at = sprintf("%s/%s", core_paths[["logs_R"]], run_name),
                                  tmp_at = sprintf("%s/%s", core_paths[["tmp_data_R"]], run_name))
  )
  , tar_target(
    name = pred_data,
    command = load_pred_data(data = pred_file_paths,
                             write_at = sprintf("%s/%s", core_paths[["results_R"]], run_name),
                             geno_mapping_from = sprintf("%s/results/R/process_cgm_data/BLUEs_within_env_cgm.qs", project_path))
  )
  #, tar_target(
  #  name = reset_job_meta_data,
  #  command = change_job_meta_data(data = pred_file_paths,
  #                                 new_time = "2-0:0:0",
  #                                 cv_type_vec = c("cv_wtn_tra"),
  #                                 model_name_vec = paste0("M_", c()))
  #)
)

## t-tests
#to_check <- "plot_wtn_with_cgm"
##to_check <- "plot_wtn_no_cgm"
#
#data <- tar_read(pred_data)[[to_check]]$data
#models <- as.vector(unique(data$model_idx))
#cv <- as.vector(unique(data$cv_idx))
#types <- as.vector(unique(data$type))
#check_against <- "M_1"
#
#data_wide <- data %>% 
#  mutate(value = ifelse(value == "NaN", NA, value)) %>%  
#  pivot_wider(id_cols = c(run_type, cv_idx, run_idx, type), 
#              names_from = model_idx, values_from = value) %>%
#  as.data.frame()
#
#cases <- data %>% 
#  distinct(model_idx, cv_idx) %>% 
#  filter(model_idx != check_against) %>%
#  as.data.frame()
#
#out <- NULL
#for (row in 1:nrow(cases)){
#  idx <- as.character(cases[row, "model_idx"])
#  idx_cv <- as.character(cases[row, "cv_idx"])
#  for(typ in types){
#    data_1 <- data_wide %>% filter(cv_idx == idx_cv, type == typ) %>% pull(all_of(check_against))
#    data_2 <- data_wide %>% filter(cv_idx == idx_cv, type == typ) %>% pull(all_of(idx))
#    val <- t.test(data_1, data_2, paired = TRUE)$p.value
#    
#    res <- data.frame("cv" = idx_cv,
#                      "idx" = idx,
#                      "type" = typ,
#                      "p_val" = ifelse(val < 0.05, "sig", "non_sig"))
#    if(is.null(out)){
#      out <- res
#    } else {
#      out <- rbind(out, res)
#    }
#  }
#}
#
#out_wide <- out %>% 
#  pivot_wider(id_cols = c("cv", "idx"), 
#              names_from = "type", 
#              values_from = "p_val") %>%
#  arrange(cv)
#
#write.table(out_wide, 
#            sprintf("/proj/results/R/process_R_pred_data/t_tests_%s_%s.txt", 
#                    to_check, check_against), row.names = FALSE)
#
#pred_data$plot_acr_sce_2$data %>%
#  filter(model_idx %in% c("acr_CNN", "M_3")) %>% 
#  distinct() %>% 
#  pivot_wider(id_cols = c("run_idx", "type", "scenario"), names_from = "model_idx", values_from = "value") %>% 
#  group_by(type, scenario) %>% 
#  summarize(superior_cnn = sum(acr_CNN >= M_3), total = n(), cnn_better = superior_cnn/total)

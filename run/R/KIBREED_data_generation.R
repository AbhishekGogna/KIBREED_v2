library(targets)

project_path <- "/proj"
core_paths <- lapply(jsonlite::read_json(sprintf("%s/results/core_paths.json", project_path)), function(x) x[[1]])

# Define packages needed and default storage format for intermittent data

tar_option_set(packages = c("dplyr", "readr", "ggplot2", "ggpubr",
                            "tidyr", "reshape2", "stringr",
                            "RColorBrewer", "dendextend", "circlize",
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
                                   tmp_data = sprintf("%s/%s", core_paths[["tmp_data_R"]], "KIBREED_data_generation"),
                                   hd_pred_at = hd_data_overview)
  ),
  tar_target(
    name = KIBREED_BLUEs_acr_env,
    command = generate_KIBREED_blues(data = predicted_HD,
                                     asreml_data_at = "/proj/ext_dir/KIBREED/results_plots/env_object_kibreed_v2.qs")
  ),
  tar_target(
    name = KIBREED_data_full,
    command = write_KIBREED_data_full(
      existing_data = data,
      data = KIBREED_BLUEs_acr_env,
      write_path_for_R = sprintf("%s/%s", core_paths[["results_R"]], "KIBREED_data_generation"),
      write_path_for_Py = sprintf("%s/%s", core_paths[["results_Py"]], "KIBREED_data_generation")
    )
  ),
  tar_target(
    name = overviews,
    command = make_plots_and_tables(
      data = KIBREED_data_full,
      write_path_for_R = sprintf("%s/%s", core_paths[["results_R"]], "KIBREED_data_generation"),
      write_path_for_Py = sprintf("%s/%s", core_paths[["results_Py"]], "KIBREED_data_generation")
    )
  )
)

#check <- KIBREED_data_full$BLUES_within_env %>% group_by(Env, Geno_new) %>% 
#  summarize(blues = mean(BLUES_dt, na.rm = TRUE), .groups = "drop") %>% 
#  pivot_wider(id_cols = "Geno_new", names_from = "Env", values_from = "blues") %>%
#  select(-Geno_new)
#
#check_2 <- KIBREED_data_full$BLUES_within_env %>% group_by(Env, Geno_dedup) %>% 
#  summarize(blues = mean(BLUES_dt, na.rm = TRUE), .groups = "drop") %>% 
#  pivot_wider(id_cols = "Geno_dedup", names_from = "Env", values_from = "blues") %>%
#  select(-Geno_dedup)
#
#get_cor_overlap_mat <- function(data) {
#  cols <- colnames(data)
#  mat <- matrix(NA, nrow = length(cols), ncol = length(cols), dimnames = list(cols, cols))
#  cor_coord <- as.data.frame(which(upper.tri(mat, diag = FALSE) == TRUE, arr.ind = TRUE))
#  overlap_coord <- as.data.frame(which(lower.tri(mat, diag = FALSE) == TRUE, arr.ind = TRUE))
#  
#  for (i in 1:nrow(overlap_coord)) {
#    row_id <- overlap_coord[i, 1]
#    col_id <- overlap_coord[i, 2]
#    
#    mat[row_id, col_id] <- length(which(!is.na(data[, row_id]) & !is.na(data[, col_id])))
#  }
#  
#  for (i in 1:nrow(cor_coord)) {
#    row_id <- cor_coord[i, 1]
#    col_id <- cor_coord[i, 2]
#    
#    tryCatch({
#      mat[row_id, col_id] <- round(cor(data[, row_id], data[, col_id], use = "complete.obs"), 2)
#    }, error = function(e) {
#      # Handle the case when there are no complete observations
#      mat[row_id, col_id] <- NA
#    })
#  }
#  
#  return(mat)
#}
#
#mat_1 <- get_cor_overlap_mat(check)
#mat_2 <- get_cor_overlap_mat(check_2)
#
#not_missing_both <- !is.na(mat_1) & !is.na(mat_2)
#only_present_mat_2 <- is.na(mat_1) & !is.na(mat_2)
#not_missing <- not_missing_both | only_present_mat_2
#
#diffs <- as.data.frame(which(not_missing, arr.ind = TRUE), row.names = NA)
#diffs$mat_1 <- NA
#diffs$mat_2 <- NA
#diffs$type <- ifelse(diffs$row < diffs$col, "cor", "overlap")
#
#mat_1_str <- matrix(as.character(mat_1), ncol = ncol(mat_1), dimnames = dimnames(mat_1))
#mat_2_str <- matrix(as.character(mat_2), ncol = ncol(mat_2), dimnames = dimnames(mat_2))
#
#mat_com <- mat_1_str
#  
#for (i in 1:nrow(diffs)){
#  row_id <- diffs[i, 1]
#  col_id <- diffs[i, 2]
#  
#  val_before <- as.numeric(mat_1_str[row_id, col_id])
#  val_after <- as.numeric(mat_2_str[row_id, col_id])
#  
#  diffs$mat_1[i] <- val_before
#  diffs$mat_2[i] <- val_after
#  
#  if(!is.na(val_before)){
#    if(val_before < val_after){
#      mat_com[row_id, col_id] <- paste0(val_before, "_i_", val_after) 
#    } else if (val_before > val_after) {
#      mat_com[row_id, col_id] <- paste0(val_before, "_d_", val_after) 
#    } else if (val_before == val_after) {
#      mat_com[row_id, col_id] <- as.character(val_before)
#    }
#  } else {
#    mat_com[row_id, col_id] <- paste0("NA_", val_after) 
#  }
#}
#
#colnames(mat_com) <- 1:ncol(mat_com)
#
#write.csv2(mat_com, "/proj/tmp_data/dedup.csv", quote = T, row.names = F)
#
#n <- 117
#n_overlaps <- (n*(n-1))/2  
#
#diffs %>% 
#  mutate(type_val = ifelse(mat_2 > mat_1, "inc",
#                           ifelse(mat_2 < mat_1, "dec", 
#                                  ifelse(mat_2 == mat_1, "same", "und")))) %>%
#  group_by(type, type_val) %>% 
#  summarize(val = n()/n_overlaps, .groups = "drop")









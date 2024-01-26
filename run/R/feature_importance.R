library(targets)

project_path <- "/proj"
core_paths <- lapply(jsonlite::read_json(sprintf("%s/results/core_paths.json", project_path)), function(x) x[[1]])

# Define packages needed and default storage format for intermittent data
options(java.parameters = "-Xmx100G")

tar_option_set(packages = c("dplyr", "readr", "ggplot2", "ggmap", 
                            "tidyr", "reshape2", "stringr",
                            "lubridate", "hablar", "readxl",
                            "qs", "feather",
                            "corehunter", "rJava", "BGLR", 
                            "asreml", "RColorBrewer", "dendextend",
                            "ggpubr"),
               format = "qs")

run_name <- "feature_importance"

tar_source(sprintf("/%s/src/R/fun_%s.R", project_path, run_name))

# Load exiting data
ext_parse <- "results/R"
data_paths <- list("g_data" = sprintf("/proj/%s/KIBREED_data_generation/GNdata_comb_add.qs", ext_parse), # based on genetic data
                   "p_wtn" = sprintf("/proj/%s/process_cgm_data/BLUEs_within_env_cgm.qs", ext_parse), # raw p_data
                   "G_a_RM" = sprintf("/proj/%s/preprocessing_geno_data/kin_a.qs", ext_parse),  # based on genetic data
                   "G_d_RM" = sprintf("/proj/%s/preprocessing_geno_data/kin_d.qs", ext_parse),  # based on genetic data
                   "G_aa_RM" = sprintf("/proj/%s/preprocessing_geno_data/kin_aa.qs", ext_parse),  # based on genetic data
                   "ERM_l" = sprintf("/proj/%s/KIBREED_data_generation/ERM_data_linear.qs", ext_parse), # based on environment data
                   "ERM_nl" = sprintf("/proj/%s/KIBREED_data_generation/ERM_data_non_linear.qs", ext_parse), # based on environment data
                   "SRM" = sprintf("/proj/%s/KIBREED_data_generation/SRM.qs", ext_parse), # based on environment data
                   "YRM" = sprintf("/proj/%s/KIBREED_data_generation/YRM.qs", ext_parse), # based on environment data
                   "G_a_S" = sprintf("/proj/%s/process_cgm_data/g_s_mat.qs", ext_parse) # based on CGM output
)


list(
  tar_target(
    name = RD_data,
    command = get_RD_mat(existing_data = data_paths,
                         write_at = sprintf("%s/%s", core_paths[["results_R"]], run_name),
                         log_at = sprintf("%s/%s", core_paths[["logs_R"]], run_name),
                         tmp_at = sprintf("%s/%s", core_paths[["tmp_data_R"]], run_name),
                         key = run_name)
  ), 
  tar_target(
    name = core_set,
    command = get_core_set(existing_data = data_paths,
                           data = RD_data)
  ),
  tar_target(
    name = training_data,
    command = get_training_data(existing_data = data_paths,
                                data = core_set)
  ),
  tar_target(
    name = predicted_data,
    command = predict_missing(existing_data = data_paths,
                              data = training_data,
                              model = "ERM_nl@BRR&G_a@BRR&G_d@BRR&G_aa@BRR&G_a_ERM_nl@RKHS",
                              debug = FALSE)
    
  ),
  tar_target(
    name = residual_vals,
    command = get_residuals(existing_data = data_paths,
                            data = predicted_data)
    
  ),
  tar_target(
    name = env_clusters,
    command = get_env_clusters(existing_data = data_paths,
                               data = residual_vals,
                               ec_path = sprintf("/proj/%s/KIBREED_data_generation/ec_mat.qs", ext_parse))
  ),
  tar_target(
    name = feature_importance,
    command = get_feature_imp_plots(existing_data = data_paths,
                                    data = residual_vals,
                                    scores_at = "/proj/results/Py/feature_importance/feature_imp_scores.csv")
  
  )
)

#data <- predicted_data[["predicted_data"]][["pred"]]%>% filter(is.na(set)) 
#to_color <- data %>% filter(geno %in% c("zf003zm027", "hystarhystar"))
#data %>% 
#  ggplot(aes(x = env, y = pred, group = geno)) +
#  geom_line() +
#  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
#  geom_line(aes(x = env, y = pred, group = geno, color = "red"), data = to_color)

# identify best and worst performers per environment cluster?
#core_dist <- core_set$geno_dist %>% filter(connect_geno %in% core_set$core$sel) %>% 
#  relocate(freq, .after= "connect_geno")
#
#interval_size <- 2
#breaks_in <- seq(0, max(core_dist$freq) + interval_size, interval_size)
#core_dist_series_0 <- core_dist %>% 
#  pivot_longer(!c("connect_geno", "freq"), names_to = "env") %>% 
#  filter(!is.na(value), env != "freq") %>% 
#  mutate(freq_env = freq,
#         freq_env_class_integer = cut(freq_env, breaks = breaks_in, include.lowest = TRUE, labels = FALSE),
#         freq_env_class = paste0("(", (as.numeric(freq_env_class_integer) - 1) * interval_size + 1, ", ", as.numeric(freq_env_class_integer) * interval_size, "]"), #upperbound of lower class +1 to upper bound of this class
#         series = gsub("^([[:alnum:]_]+_\\d_[[:alnum:]]+).*", "\\1", env),
#         value = 1) 
#multi_series <- core_dist_series_0 %>%
#  select(connect_geno, freq_env_class, series) %>%
#  distinct() %>% count(connect_geno) %>% filter(n > 1)
#core_dist_series <- core_dist_series_0 %>%
#  select(connect_geno, freq_env_class, series) %>% 
#  mutate(series = ifelse(connect_geno %in% multi_series$connect_geno, "Exp_multi", series)) %>%
#  distinct() %>%
#  group_by(freq_env_class, series) %>% 
#  summarize(freq_ser = n(), .groups = "drop") %>% 
#  pivot_wider(id_cols = freq_env_class, names_from = series, values_from = freq_ser) %>%
#  rowwise() %>%
#  mutate(class_lower_bound = as.numeric(gsub("\\((\\d+)\\,.*", "\\1", freq_env_class)),
#         total = sum(c_across(starts_with("Exp")), na.rm = TRUE)) %>%
#  ungroup() %>%
#  arrange(class_lower_bound) %>%
#  select(freq_env_class, Exp_2_Zucht, Exp_3_Zucht, Exp_4_Zucht, Exp_6_KWS, Exp_7_GABI, Exp_multi, total)
#
##write.csv2(core_dist_series, "/proj/tmp_data/core_dist_series.csv", row.names = F, quote = T)
#per_series <- core_dist_series %>% select(-freq_env_class, -total) %>% as.matrix()
#per_series_sum <- colSums(per_series, na.rm = T)
#per_series_sum[order(per_series_sum, decreasing = T)]

## how many of the core set are in the extremes
#core <- core_set$core$sel
#blues <- qread(data_paths[["p_wtn"]]) %>%
#  mutate(present_in = ifelse(connect_geno %in% core, "core", "general"))
#
#dist_plot <- ggplot(aes(x = BLUES_dt), data = blues %>% filter(present_in == "general")) +
#  geom_histogram(bins = 50, fill = "#00AFBB", color = "black") +
#  geom_rug(aes(x = BLUES_dt), data = blues %>% filter(present_in == "core"), color = "#E7B800") +
#  labs(x= "yield (quintal per ha)", y = "frequency") +
#  theme_classic()
#  
#ggsave(plot = dist_plot,
#       filename = sprintf("%s/dist_plot.png", "/proj/results/R/feature_importance"), 
#       width = 8.4, height = 8.4, dpi = 600, units = "cm")

# In numbers
#target <- 10
#check <- blues %>% filter(!is.na(BLUES_dt)) %>% group_by(Env) %>%
#  arrange(BLUES_dt) %>%
#  mutate(rank = row_number()) %>%
#  filter(rank %in% 1:target | rank %in% (n() - (target -1)):n()) %>%
#  ungroup() %>%
#  select(-rank) %>%
#  distinct(Geno_new, .keep_all = T)
#
#sum(core %in% check$connect_geno)
  
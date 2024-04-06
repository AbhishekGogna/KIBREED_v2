library(targets)

project_path <- "/proj"
core_paths <- lapply(jsonlite::read_json(sprintf("%s/results/core_paths.json", project_path)), function(x) x[[1]])
data_paths <- list("data_yu" = sprintf("%s/ext_dir/KIBREED/source_data/BigData_yusheng_v1/Bigdata_Blues_within_env_All_nodup_yield.txt", project_path),
                   "transition_table_1" = sprintf("%s/ext_dir/KIBREED/source_data/BigData_yusheng_v1/Translate 25.txt", project_path),
                   "series_1_hd" = sprintf("%s/ext_dir/KIBREED/source_data/Complete_BLUEs_singleEnv_2Step_HeadingDateAsreml 4.txt", project_path),
                   "series_2_hd" = sprintf("%s/ext_dir/KIBREED/source_data/BigData_yusheng_v1/heading date Fac_S1_Yield_plot.csv", project_path),
                   "series_3" = sprintf("%s/ext_dir/KIBREED/source_data/BigData_yusheng_v1/Yield Plant height and heading date Fac_S2_Yield_plot 2 years data.xlsx", project_path),
                   "series_4_top_s1_hd" = sprintf("%s/ext_dir/KIBREED/source_data/Heading date TopS1 2015_2016.xlsx", project_path),
                   "series_4_top_s2_hd" = sprintf("%s/ext_dir/KIBREED/source_data/Heading date TopS2 2017_2018.xlsx", project_path),
                   "series_5_wp1_hd" = sprintf("%s/ext_dir/KIBREED/source_data/HD_bereinigt_corrected_for_block_effects_factorial_Wp2_S1.xlsx", project_path),
                   "series_5_wp2_hd" = sprintf("%s/ext_dir/KIBREED/source_data/Heading date_corrected_for_block_effects_WP2 S2.xlsx", project_path),
                   "series_5_wp3_hd" = sprintf("%s/ext_dir/KIBREED/source_data/Headingdate_corrected_for_block_effects_Wp2_S3.xlsx", project_path),
                   "series_6_old_hd" = sprintf("%s/ext_dir/KIBREED/source_data/BLUEs_PerEnv.txt", project_path),
                   "series_6_2020_yield" = sprintf("%s/ext_dir/KIBREED/results_plots/BLUEs_KWS_2020_acr.csv", project_path),# can be channeled from previous analysis
                   "series_6_2020_hd" = sprintf("%s/ext_dir/KIBREED/results_plots/BLUEs_KWS_2020_acr_hd.csv", project_path),# can be channeled from previous analysis
                   "series_6_2021_yield" = sprintf("%s/ext_dir/KIBREED/results_plots/BLUEs_KWS_2021_acr.csv", project_path),# can be channeled from previous analysis
                   "series_6_2021_hd" = sprintf("%s/ext_dir/KIBREED/results_plots/BLUEs_KWS_2021_acr_hd.csv", project_path),# can be channeled from previous analysis
                   "series_7" = sprintf("%s/ext_dir/GABI/BLUEs/results/BLUES_traits_v2.Rdata", project_path),
                   "transition_for_checks" = sprintf("%s/ext_dir/GABI/source_files/transition_for_checks.txt", project_path),
                   "for_kibreed_coding_to_variety" = sprintf("%s/ext_dir/GABI/source_files/for_kibreed_coding_to_variety.txt", project_path))

# Define packages needed and default storage format for intermittent data

tar_option_set(packages = c("dplyr", "readr", "ggplot2", "tidyr", "reshape2",
                            "hablar", "readxl"),
               format = "qs")

tar_source(sprintf("%s/src/R/%s", project_path, "fun_preprocessing_phenodata.R"))

# Create pipeline to process data

list(
  tar_target(
    name = preprocessing_phenodata,
    command = preprocess_phenodata(paths = data_paths,
                                   log_at = sprintf("%s/%s", core_paths[["logs_R"]], "preprocessing_phenodata"))
  )
)

# Repeatabilities 
## Grain yield
## Exp1
#library(readxl)
#library(dplyr)
#library(hablar)
#
#exp_1_GY <- read_excel("/proj/ext_dir/KIBREED/source_data/Repeatability use classical formula.xlsx", sheet = "Hywheat") %>%
#  rename(loc = LocName, year = Year, rep = Repeatability_all) %>%
#  mutate(series = "Exp_1",
#         trial = NA) %>%
#  select(loc, year, trial, rep, series) %>%
#  convert(chr(loc, trial, series),
#          int(year),
#          num(rep))
#
### Exp6
#exp_6_GY_old <- read_excel("/proj/ext_dir/KIBREED/source_data/Repeatability use classical formula.xlsx", sheet = "KWS_OLD") %>%
#  rename(loc = Location, year = Year, rep = Repeatability) %>%
#  select(loc, year, rep) %>%
#  mutate(series = "Exp_6",
#         trial = NA)%>%
#  select(loc, year, trial, rep, series) %>%
#  convert(chr(loc, trial, series),
#          int(year),
#          num(rep))
#
#exp_6_GY_new <- read.table("/proj/results/R/phenotypic_data_kws_processing/rep_20_21.txt", header = TRUE) %>%
#  filter(grepl("yield", element_name)) %>%
#  filter(!is.na(location)) %>%
#  mutate(loc = location,
#         rep = value,
#         trial = strsplit(element_name, "\\_")[[1]][2],
#         series = "Exp_6")%>%
#  select(loc, year, trial, rep, series) %>%
#  convert(chr(loc, trial, series),
#          int(year),
#          num(rep))
#
### Exp 7
#load(data_paths[["series_7"]])
#exp_7_GY <- BLUES_traits[["YIE"]][["Repeatability"]] %>%
#  mutate(loc = strsplit(Environment, "\\.")[[1]][2],
#         year = strsplit(Environment, "\\.")[[1]][1],
#         rep = Repeatibility,
#         series = "Exp_7",
#         trial = NA) %>%
#  select(loc, year, trial, rep, series) %>%
#  convert(chr(loc, trial, series),
#          int(year),
#          num(rep))
#
#all_rep_GY <- bind_rows(exp_1_GY,
#                     exp_6_GY_old,
#                     exp_6_GY_new,
#                     exp_7_GY) %>% 
#  group_by(series) %>%
#  summarize(min_rep = round(min(rep), 2)
#            , max_rep = round(max(rep), 2)
#            , med_rep = round(median(rep), 2)
#            , .groups = "drop")
#
### Heading date
### Exp 6
#exp_6_HD_new <- read.table("/proj/results/R/phenotypic_data_kws_processing/rep_20_21.txt", header = TRUE) %>%
#  filter(grepl("head", element_name)) %>%
#  filter(!is.na(location)) %>%
#  mutate(loc = location,
#         rep = value,
#         trial = strsplit(element_name, "\\_")[[1]][2],
#         series = "Exp_6")%>%
#  select(loc, year, trial, rep, series) %>%
#  convert(chr(loc, trial, series),
#          int(year),
#          num(rep))
### Exp 7
#load(data_paths[["series_7"]])
#exp_7_HD <- BLUES_traits[["HD"]][["Repeatability"]] %>%
#  mutate(loc = strsplit(Environment, "\\.")[[1]][2],
#         year = strsplit(Environment, "\\.")[[1]][1],
#         rep = Repeatibility,
#         series = "Exp_7",
#         trial = NA) %>%
#  select(loc, year, trial, rep, series) %>%
#  convert(chr(loc, trial, series),
#          int(year),
#          num(rep))
#all_rep_HD <- bind_rows(exp_6_HD_new,
#                        exp_7_HD) %>% 
#  group_by(series) %>%
#  summarize(min_rep = round(min(rep), 2)
#            , max_rep = round(max(rep), 2)
#            , med_rep = round(median(rep), 2)
#            , .groups = "drop")

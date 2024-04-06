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

#convert_to_dataframe <- function(all_rep) {
#  # Initialize empty vectors to store data
#  element_name <- character()
#  location <- character()
#  value <- numeric()
#  # Iterate through the elements of 'all_rep'
#  for (element in names(all_rep)) {
#    if (!any(is.na(all_rep[[element]]))) {
#      # Get the location names
#      locations <- all_rep[[element]][, "Location"]
#      # Get the repeatability values
#      repeatability <- as.numeric(all_rep[[element]][, "repeatability"])
#      # Append data to the vectors
#      element_name <- c(element_name, rep(element, length(locations)))
#      location <- c(location, locations)
#      value <- c(value, repeatability)
#    } else {element_name <- c(element_name, element)
#    location <- c(location, NA)
#    value <- c(value, NA)}
#  }
#  
#  # Create a dataframe
#  dataframe <- data.frame(element_name, location, value)
#  return(dataframe)
#}
#
#tar_load(BLUEs_kws_2020_data)
#tar_load(BLUEs_kws_2021_data)
#all_rep_2020 <- convert_to_dataframe(lapply(BLUEs_kws_2020_data$BLUES_traits_V2, function(x) x[["Repeatability"]]))
#all_rep_2020$year <- "2020"
#all_rep_2021 <- convert_to_dataframe(lapply(BLUEs_kws_2021_data$BLUES_traits_V2, function(x) x[["Repeatability"]]))
#all_rep_2021$year <- "2021"
#all_rep <- rbind(all_rep_2020, all_rep_2021)
#write.table(all_rep, "/proj/results/R/phenotypic_data_kws_processing/rep_20_21.txt", row.names = FALSE)
#all_rep %>% filter(!is.na(value)) %>% 
#  mutate(trait = gsub("BLUEs_\\S\\S\\_(\\S+)", "\\1", element_name)) %>% group_by(trait) %>%
#  summarise(mean_rep = mean(value), min = min(value), max = max(value))

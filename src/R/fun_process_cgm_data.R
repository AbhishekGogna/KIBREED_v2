# core functions
"%!in%" <- Negate("%in%")

prep_as_df_for_feather <- function(data){
  data_df <- cbind("idx" = rownames(data), data)
  data_df <- as.data.frame(data_df)
  return(data_df)
}

mapping_function <- function(.x, ref_df){
  data <- .x
  geno <- unique(data$Genotype)
  #out <- data %>% arrange(SimulationID) %>%
  #  bind_cols(ref_df %>% filter(gens == geno))
  return(print(geno))
}

process_cgm_data <- function(existing_data, paths, 
                             write_path_for_R, write_path_for_Py, 
                             log_at, tmp_at){
  #' given data and paths, process data and return processed data
  
  # Put a log file
  run_instance <- as.character(format(Sys.time(),  format = "%d_%m_%Y_%H_%M"))
  log_at <- sprintf("%s/%s", log_at, run_instance)
  log_file <- sprintf("%s/process_cgm_data.log", log_at)
  
  if(!dir.exists(log_at)){dir.create(log_at, recursive = T)}
  
  # Put a temp_dir
  tmp_data <- sprintf("%s/%s", tmp_at, run_instance)
  
  if(!dir.exists(tmp_data)){dir.create(tmp_data, recursive = T)}
  # Put a write_dir  
  if(!dir.exists(write_path_for_R)){dir.create(write_path_for_R, recursive = T)}
  if(!dir.exists(write_path_for_Py)){dir.create(write_path_for_Py, recursive = T)}
  
  # write log
  cat("Adjusting the KIBREED wtn data to fit the cgm input ------------------",
      file = log_file,
      sep = "\n")
  
  # Sequester data
  pheno_data_sent <- existing_data$pheno_data %>% filter(!is.na(sowing_date)) # this is what i sent
  param_metadata <- read.csv(paths[["param_metadata"]], header = TRUE) %>%
    convert(chr(Number))
  param_metadata$Parameter.name[which(param_metadata$Parameter.name == "DST")] <- paste0("DST", 1:5)
  param_data <- read.csv(paths[["param_data"]], header = TRUE)
  if(!is.null(paths[["pheno_data_used"]]))
    pheno_data_used_0 <- read.csv(paths[["pheno_data_used"]], header = TRUE)
    pheno_data_used <- pheno_data_used_0 %>%
    rowwise() %>% 
    mutate(latlong = paste0(formatC(Latitude, format = "f", digits = 5), ",", formatC(Longitutde, format = "f", digits = 6)),
           Site = gsub("(\\S+)\\d", "\\1", Site, perl = T)) %>%
    ungroup() %>%
    group_by(Genotype) %>% group_modify(~arrange(.x, SimulationID)) %>%
    ungroup() %>%
    select(Genotype, Site, Year, latlong, Yield_dt, Heading_DOY) # this is what they used
  pheno_data_recieved <- read.csv(paths[["pheno_data_recieved"]], header = TRUE) # this is what they produces
  
  # Connect pheno_data_used and pheno_data_received
  data_comb_at <- sprintf("%s/data_comb.qs", write_path_for_R)
  if(!file.exists(data_comb_at)){
    data_comb <- NULL
    uq_gens <- unique(pheno_data_recieved$gens)
    for(gen in uq_gens){
      out <- pheno_data_used %>% filter(Genotype == gen) %>%
        bind_cols(pheno_data_recieved %>% filter(gens == gen) %>% rename(idx_y = X))
      if(is.null(data_comb)){
        data_comb <- out
      } else {
        data_comb <- data_comb %>% bind_rows(out)
      }
      gen_num <- which(uq_gens == gen)
      if(gen_num %% 2000 == 0) cat(sprintf("%s: Done for %s out of %s\n", format(Sys.time(), format = "%Y-%m-%dT%H:%M:%S"), gen_num, length(uq_gens)))
    }
    data_comb <- data_comb %>% mutate(same_HD = ifelse(Heading_DOY == HDobs, TRUE, FALSE), 
                                      same_yie = ifelse(round(Yield_dt/10, 2) == round(Yobs, 2), TRUE, FALSE)) %>%
      select(-site) %>%
      rename(year = Year, site = Site)
    qsave(data_comb, data_comb_at)
  } else {
    data_comb <- qread(data_comb_at)
  }
  
  pheno_data_recieved <- data_comb %>% select(latlong, site, gens, year, Parameterindex)

  # Reshape param data
  param_data_col_names <- data.frame("raw_name" = colnames(param_data),
                                     "link_name" = gsub("X", "", colnames(param_data))) %>%
    mutate(link_name = ifelse(link_name == "", "Parameterindex", link_name)) %>%
    left_join(param_metadata, by = c("link_name" = "Number")) %>%
    mutate(new_name = ifelse(is.na(Parameter.name), link_name, Parameter.name))
  colnames(param_data) <- param_data_col_names$new_name
  
  ## Site renaming
  site_rename <- data.frame("new_site" = c("Grossaitingen", "Nuorvenich", "Suollingen"),
                            "old_site" = c("Großaitingen", "Nörvenich", "Söllingen"))
  non_matching_geno <- unique((pheno_data_recieved$gens[pheno_data_recieved$gens %!in% pheno_data_sent$Geno_new]))
  
  ## Genotype renaming
  geno_rename <- data.frame("new_geno" = non_matching_geno,
                            "Geno_new" = gsub(",", " ", non_matching_geno))
  pheno_data_recieved_reformatted <- pheno_data_recieved %>% left_join(site_rename, by = c("site" = "new_site")) %>%
    left_join(geno_rename, by = c("gens" = "new_geno")) %>%
    mutate(old_site = ifelse(is.na(old_site), site, old_site),
           Geno_new = ifelse(is.na(Geno_new), gens, Geno_new)) %>%
    select(-site, -gens) %>%
    rename(site = old_site)
  
  # recreate zalf input data. i do this to understand what steps they took to 
  # compress the data i sent initially
  
  #pheno_data_multiples <- pheno_data_sent %>% 
  #  count(Site, Year, Geno_new) %>% 
  #  filter(n > 1) %>% select(-n)
#
  #target_env <- pheno_data_sent %>% filter(!is.na(sowing_date)) %>%
  #  count(Env, Altitude, Connect_at, sowing_date, Sowing_density_grain_m2, plot_area_m2, Soil_type, Soil_texture) %>%
  #  select(-n) %>%
  #  mutate(sowing =  paste0(str_sub(Connect_at, end= -5), (as.numeric(str_sub(Connect_at, start= -4)) - 1)),
  #         harvest = Connect_at) %>% select(-Connect_at) %>%
  #  pivot_longer(cols = c("harvest", "sowing"), names_to = "connect_type", values_to = "Connect_at") # For the environments where sowing data are available, i added the next year too since this would allow extraction of climate data for both sowing and harvest year. Wheat, for reference, grows in 10 - 11 months spread over two years.
  
  recreated_pheno_data <- pheno_data_sent %>%filter(!is.na(sowing_date)) %>% 
    filter(!is.na(BLUES_dt)) %>% convert(int(Year)) # 46643 rows. this is the base data for me
  
  # Merge with cgm data
  pheno_data_combined <- recreated_pheno_data %>% 
    left_join(pheno_data_recieved_reformatted, by = c("Site" = "site", "latlong", "Year" = "year", "Geno_new")) %>%
    rename(connect_climate = Connect_at,
           connect_params_idx = Parameterindex,
           connect_geno = connect_geno_data) %>%
    mutate(latlong = gsub(",", "_", latlong),
           connect_param = paste0(latlong, "_", Geno_new),
           idx_cv = row_number()) %>% 
    relocate(tidyselect::contains("connect"), .after = last_col()) %>%
    relocate(idx_cv)# 46643
  
  # recreate input data for zalf
  pheno_data_export_0 <- pheno_data_used_0 %>% 
    mutate(site_mod = gsub("(\\S+)\\d", "\\1", Site, perl = T)) %>%
    left_join(site_rename, by = c("site_mod" = "new_site")) %>%
    left_join(geno_rename, by = c("Genotype" = "new_geno")) %>%
    mutate(old_site = ifelse(is.na(old_site), site_mod, old_site),
           Geno_new = ifelse(is.na(Geno_new), Genotype, Geno_new)) %>%
    select(-site_mod)
  
  pheno_data_export <-  pheno_data_combined %>% 
    select(idx_cv, Env, Site, Year, Geno_new, BLUES_dt, HD, sowing_date, harvest_date) %>% 
    rename(Yield_dt_orig = BLUES_dt, Heading_DOY_orig = HD, old_site = Site) %>%
    rowwise() %>%
    mutate(Yield_dt_orig = round(Yield_dt_orig, 1), 
           Heading_DOY_orig = round(Heading_DOY_orig),
           sowing_date_zalf = format(sowing_date, "%d.%m.%Y"),
           sim_start_zalf = as.Date(sowing_date_zalf, format = "%d.%m.%Y") - 5,
           sim_start_zalf = format(sim_start_zalf, format = "%d.%m.%Y"),
           sim_end_zalf = as.Date(sowing_date_zalf, format = "%d.%m.%Y") + 315,
           sim_end_zalf = format(sim_end_zalf, format = "%d.%m.%Y"),
           harvest_date =  format(harvest_date, format = "%d.%m.%Y")) %>%
    ungroup() %>%
    select(-sowing_date) %>%
    left_join(pheno_data_export_0, by = c("old_site","Year", "Geno_new")) %>%
    select(c("idx_cv", "Env", "harvest_date", "Yield_dt_orig", "Heading_DOY_orig", "sowing_date_zalf", "sim_start_zalf", "sim_end_zalf", colnames(pheno_data_used_0))) %>%
    arrange(idx_cv) %>%
    select(- Yield_dt, -Heading_DOY, -sowing_date, -SimulationStart, -SimulationEnd) %>%
    rename(Yield_dt = Yield_dt_orig,
           Heading_DOY = Heading_DOY_orig,
           sowing_date = sowing_date_zalf,
           SimulationStart = sim_start_zalf,
           SimulationEnd = sim_end_zalf) %>%
    relocate(sowing_date, SimulationStart, SimulationEnd, .after = Genotype) %>%
    relocate(Yield_dt, .after = SimulationID2) %>%
    relocate(Heading_DOY, .after = SimulationEnd)
  
  #check <- pheno_data_export %>% mutate(sowing_sim_diff = as.Date(SimulationEnd, format = "%d.%m.%Y") - as.Date(sowing_date, format = "%d.%m.%Y"),
  #                                      sowing_harvest_diff = as.Date(harvest_date, format = "%d.%m.%Y") - as.Date(sowing_date, format = "%d.%m.%Y")) # better replace it 
  
  # Calculate gxy kinship
  g_x_s_data <- pheno_data_combined %>% distinct(connect_param, connect_params_idx)  %>%
    left_join(param_data, by = c("connect_params_idx" = "Parameterindex")) %>%
    select(-connect_params_idx) %>% 
    relocate(connect_param)
  g_x_s_data_mat <- as.matrix(g_x_s_data[, -1])
  rownames(g_x_s_data_mat) <- g_x_s_data %>% pull(connect_param) %>% as.vector()
  g_x_s_data_mat_scaled <- scale(g_x_s_data_mat, center = T, scale = T)
  #g_s_mat_0 <- g_x_s_data_mat %*% t(g_x_s_data_mat)
  #g_s_kin <- g_s_mat_0 / mean(diag(g_s_mat_0))
  
  # Generate output
  cat("Data prepared",
      file = log_file,
      sep = "\n",
      append = T)
  
  out <- list()
  out[["log_file"]] <- log_file
  out[["tmp_data"]] <- tmp_data
  out[["write_at_R"]] <- write_path_for_R
  out[["write_at_Py"]] <- write_path_for_Py
  out[["pheno_data_wtn"]] <- pheno_data_combined %>% select(-connect_params_idx)
  out[["params_data"]] <- param_data
  out[["g_s_mat"]] <- g_x_s_data_mat_scaled
  out[["pheno_wtn_export"]] <- pheno_data_export
  
  return(out)
}

write_cgm_data <- function(data){
  # Sequester data
  log_file <- data[["log_file"]]
  write_at_R <- data[["write_at_R"]]
  write_at_Py <- data[["write_at_Py"]]
  BLUEs_within_env_cgm <- data[["pheno_data_wtn"]]
  param_data <- data[["params_data"]]
  g_s_mat <- data[["g_s_mat"]]
  pheno_data_export <- data[["pheno_wtn_export"]]
  
  # Produce output
  out <- list()
  out[["log_file"]] <- log_file
  out[["BLUEs_within_env_cgm"]] <- BLUEs_within_env_cgm
  out[["params_data"]] <- param_data
  out[["g_s_mat"]] <- g_s_mat
  out[["pheno_data_export"]] <- pheno_data_export
  
  # Write a section in the log file
  cat("Writing BLUEs within environment but reduced to fit cgm data ------------------",
      file = log_file,
      sep = "\n",
      append = TRUE)
  
  # Write data
  for(data_name in c("BLUEs_within_env_cgm", "params_data", "g_s_mat", "pheno_data_export")){
    qsave(out[[data_name]],
          file = sprintf("%s/%s.qs", write_at_R, data_name))
    write_feather(prep_as_df_for_feather(out[[data_name]]),
                  sprintf("%s/%s.feather", write_at_Py, data_name))
    if(data_name == "pheno_data_export") write.csv2(out[[data_name]], sprintf("%s/%s.csv", write_at_R, data_name), 
                                                   quote = FALSE, row.names = FALSE)
  }
  
  # Produce output
  return(out)
}
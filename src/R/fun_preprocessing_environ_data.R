preprocess_environ_data <- function(paths, log_at, existing_data){
  # Core functions
  "%!in%" <- Negate("%in%")
   
  # Generate output
  out <- list()

  # Put a log file
  run_instance <- as.character(format(Sys.time(),  format = "%d_%m_%Y_%H_%M"))
  log_at <- sprintf("%s/%s", log_at, run_instance)
  
  if(!dir.exists(log_at)){dir.create(log_at, recursive = T)}
  
  # Sequester data
  combined_data <- existing_data
  
  # Combining env data and management data to phenodata-----------------------
  
  # get env data
  coord_raw <- read.csv(paths[["coord_raw"]], header = T) 
  coord_raw[coord_raw == ""] <- NA
  
  coord <- coord_raw %>% as_tibble() %>%
    mutate(Exp = gsub(" ", "", Experiment),
           env = paste0(Code, "_", Exp, "_", Year),
           lat = sprintf("%.5f", lat_degree),
           long = sprintf("%.6f", long_degree),
           latlong = paste0(lat, ",", long),
           year = as.character(Year),
           sowing_date = as.Date(Sowing_date, origin = "1899-12-30"),
           Altitude = ifelse(Altitude == "-" | Altitude == "", NA, Altitude)) %>% 
    select(Series, env, latlong, Altitude, Site, year, sowing_date, sowing_date_raw, 
           Sowing_density_grain_m2, plot_area_m2, Soil_type, Soil_texture)
  
  env_data_sci_adv <- read.delim(paths[["env_data_sci_adv"]], 
                                 header = T) %>% mutate(loc = Loc) %>% select(-Loc)
  
  env_data_kws_2020 <- read.delim(paths[["env_data_kws_2020"]], 
                                  header = T) %>% select(colnames(env_data_sci_adv)) # is a subset of env_sc_adv
  
  env_data_kws_2021 <- read.delim(paths[["env_data_kws_2021"]], 
                                  header = T) %>% dplyr::rename(loc = location) %>% 
    select(colnames(env_data_sci_adv))
  env_data_kws_2021[which(env_data_kws_2021$loc == "Grossaitingen"), "loc"] <- "Großaitingen" # correct name.
  env_data_kws_2021[which(env_data_kws_2021$loc == "Noervenich"), "loc"] <- "Nörvenich"
  env_data_kws_2021[which(env_data_kws_2021$loc == "Wetze_kleine_Horst"), "loc"] <- "Wetze" # I choose this as the value for wetze. Moritz has only for this location two coordinates, i.e. klein-horst and haasburg. These are just ~ 2 km apart so i took only one of them a relevant. 
  options(digits =10)
  env_data_kibreed <- env_data_sci_adv %>% bind_rows(env_data_kws_2020) %>% bind_rows(env_data_kws_2021) %>%
    mutate(lat = sprintf("%.5f", as.numeric(gsub("(\\d+\\.\\d+)\\,\\s?\\d+\\.\\d+", "\\1", latlong, perl = T))),
           long = sprintf("%.6f", as.numeric(gsub("\\d+\\.\\d+\\,\\s?(\\d+\\.\\d+)", "\\1", latlong, perl = T))),
           latlong = paste0(lat, ",", long),
           date = as.Date(date),
           year = format(date, format="%Y")) %>% select(-lat, -long) %>%
    relocate(loc, latlong, year)
  
  #check_miss <- apply(env_data_kibreed, 2, function(x) sum(is.na(x)))
  #check_miss_prop <- check_miss/nrow(env_data_kibreed) #todo - account for missing values
  
  #extra_leap_year_days <- env_data_kibreed %>% distinct(latlong, year, date) %>% 
  #  filter((as.numeric(year) %% 4) == 0) %>%
  #  filter(grepl("-02-29", date))
  
  #env_data_kibreed <- env_data_kibreed %>% anti_join(extra_leap_year_days, by = colnames(extra_leap_year_days)) # 193 days removed
  
  # connect coord data with recieved data from KWS
  # minor edits in coordinates
  env_data_kibreed[which(env_data_kibreed$latlong == "51.83210,11.713056"), "latlong"] <- "51.83208,11.713056" # Bernburg
  env_data_kibreed[which(env_data_kibreed$latlong == "52.73834,9.006702"), "latlong"] <- "52.73832,9.006696" # Asendorf
  env_data_kibreed[which(env_data_kibreed$latlong == "49.54225,10.193251"), "latlong"] <- "49.54222,10.193194" # Aspachhof
  env_data_kibreed[which(env_data_kibreed$latlong == "52.23989,10.170323"), "latlong"] <- "52.23989,10.170306" # Adenstedt
  env_data_kibreed[which(env_data_kibreed$latlong == "51.27528,13.136667"), "latlong"] <-"51.27544,13.136842"  # Nasenberg kws 2021
  env_data_kibreed[which(env_data_kibreed$latlong == "54.40000,9.200000"), "latlong"] <- "50.40000,9.200000" # Wohlde # serios fault. Get it checked and input again.
  
  # currently trials for exp_3 fac_s2 harvested in 2019 lack any env data and this is expected.
  
  # kws data
  # kws 2021 sowing dates
  experiments_2020_2021 <- read.csv(paths[["metadata_2020"]], header = T) %>% 
    bind_rows(read.csv(paths[["metadata_2021"]])) %>%
    rename(Site = `location`) %>%
    mutate(sowing_date = as.Date(sowing_date, origin = "1900-01-01"),
           harvest_date = as.Date(harvest_date, origin = "1900-01-01"),
           year = as.character(year),
           subtrial = gsub("(\\S\\S)\\d\\_.*", "\\1", experiment_name, perl = T),
           Sowing_density_grain_m2 = sowing_density,
           plot_area_m2 = harvest_net_area) %>% 
    group_by(Site, year, Sowing_density_grain_m2) %>% 
    summarize(sowing_date = mean(sowing_date), 
              harvest_date = mean(harvest_date, na.rm = T),
              Sowing_density_grain_m2 = mean(Sowing_density_grain_m2), 
              .groups = "drop") %>% 
    mutate(sowing_date_raw = gsub("-", "@", format(sowing_date, "%d-%m-%Y")),
           harvest_date_raw = gsub("-", "@", format(harvest_date, "%d-%m-%Y")),
           Series = "Exp_6") %>%
    right_join(coord %>% filter(Series == "Exp_6", year %in% c("2020", "2021")) %>% 
                 select(-sowing_date, -sowing_date_raw, -Sowing_density_grain_m2),
               by = c("Series", "Site", "year")) %>%
    relocate(all_of(colnames(coord))) %>%
    convert(chr(Sowing_density_grain_m2))
  
  coord_mod <- coord %>% anti_join(experiments_2020_2021, by = colnames(experiments_2020_2021)[1:6]) %>%
    bind_rows(experiments_2020_2021) %>%
    left_join(env_data_kibreed %>% distinct(loc, year, latlong), by = c("latlong", "year")) %>% 
    arrange(Series, year) # checks if all environments of combined data are present in env data 
  
  # missing values in the last column "loc" indicate that the corresponding env data is absent
  
  #write.table(coord_mod, "~/KIBREED/results_plots/coord_mod.txt", row.names = F)
  
  env_data_kibreed <- env_data_kibreed %>% mutate(loc = gsub(" ", "", loc, perl = T), 
                                                  connect_at = paste0(gsub(",", "_", latlong), "_", loc, "_", year))
  
  
  # combine the two data
  
  #rm(list= ls()[!(ls() %in% c('combined_data', "combined_data_with_env"))])
  combined_data_with_env <- combined_data %>% 
    left_join(coord_mod, by = c("Env" = "env",
                                "Series" = "Series")) %>% select(-year) %>%
    mutate(Env = paste0(Series, "_", Project, "_", Exp, "_", Loc, "_", Year),
           Connect_at = paste0(gsub(",", "_", latlong), "_", Site, "_", Year)) %>%
    convert(fct(Series, Project, Exp, Loc, Year, Env, Type, FemaleNew, MaleNew, Geno_new),
            chr(latlong)) %>% # there are sites where more than one environment has the same year and coordinates hence the same value for connect at eg. 51.99161_11.303554_Hadmersleben_2017 has three environments 
    filter(Connect_at %in% env_data_kibreed$connect_at) # removes those environments where there is no data available
  #NRV_Inbred_2021 is missing too but looks like there is not phenotype data for it.
  
  # get informative years
  informative_years <- combined_data_with_env %>% distinct(Connect_at) %>% 
    mutate(previous_yr = paste0(str_sub(Connect_at, end= -5), (as.numeric(str_sub(Connect_at, start= -4)) - 1)))
  
  unique_env <- unique(c(informative_years %>% pull(Connect_at), informative_years %>% pull(previous_yr)))
  
  check_miss <- env_data_kibreed %>% filter(connect_at %in% unique_env) %>%
    mutate(day = format(date, format="%j")) %>% select(-loc, -latlong, -year, -date) %>% relocate(connect_at) 
  
  check_miss_val <- apply(check_miss, 2, function(x) sum(is.na(x)))/nrow(check_miss)
  
  check_miss_val_rm <- names(check_miss_val[which(check_miss_val > 0.05)]) # 5 of these may be salvaged but discuss this
  
  #check_miss_val_rm <- names(check_miss_val[which(check_miss_val == 0)]) # 5 of these may be salvaged but discuss this
  
  # some previous years overlap with harvest year of a given index value
  
  env_data_kibreed_raw <- env_data_kibreed <- env_data_kibreed %>% select(-all_of(check_miss_val_rm))
  
  cat("Imputing missing values",
      file = sprintf("%s/preprocessing_environ_data.log", log_at),
      sep = "\n")
  
  for(i in 1:nrow(informative_years)){ 
    data_subset_harvest_year <- env_data_kibreed %>% filter(connect_at == informative_years[i, ] %>% pull(Connect_at)) %>%
      mutate(day = format(date, format="%j")) %>% select(-loc, -latlong, -year, -date) %>% relocate(connect_at) %>%
      filter(day <= 243) %>% # assuming harvest before 31 august
      pivot_longer(-c("day", "connect_at")) %>% pivot_wider(names_from = day, values_from = value)
    data_subset_previous_year <- env_data_kibreed %>% filter(connect_at ==  informative_years[i, ] %>% pull(previous_yr)) %>%
      mutate(day = format(date, format="%j")) %>% select(-loc, -latlong, -year, -date) %>% relocate(connect_at) %>%
      filter(day > 306) %>% # assuming sowing after 1 november
      pivot_longer(-c("day", "connect_at")) %>% pivot_wider(names_from = day, values_from = value) 
    if("366" %in% colnames(data_subset_previous_year)){
      data_subset_previous_year <- data_subset_previous_year %>% select(-`366`)
    }
    colnames(data_subset_harvest_year)[3:ncol(data_subset_harvest_year)] <- paste0("h_", colnames(data_subset_harvest_year)[3:ncol(data_subset_harvest_year)])
    colnames(data_subset_previous_year)[3:ncol(data_subset_previous_year)] <- paste0("p_", colnames(data_subset_previous_year)[3:ncol(data_subset_previous_year)])
    
    final_data <- cbind(data_subset_harvest_year[, 1:2],
                        data_subset_previous_year[, 3:ncol(data_subset_previous_year)],
                        data_subset_harvest_year[, 3: ncol(data_subset_harvest_year)])
    
    if(i == 1){
      env_data_kibreed_sub <- final_data
    } else {
      env_data_kibreed_sub <- rbind(env_data_kibreed_sub, final_data)
    }
    if(i %% 10 == 0){cat(paste0("Done for ", i, " row"),
                         file = sprintf("%s/preprocessing_environ_data.log", log_at),
                         append = T,
                         sep = "\n")}
  }
  
  # simple check
  
  unique(combined_data_with_env$Connect_at[which(combined_data_with_env$Connect_at %!in% env_data_kibreed_sub$connect_at)]) # should be zero for full availability
  
  #which(is.na(env_data_kibreed_sub), arr.ind = T) %>% as_tibble() %>% group_by(row) %>% summarise(range = range(col)) %>% mutate(id = row_number()) %>% pivot_wider(names_from = id, values_from = range) %>% rename("col_1" = '1', "col_2" = '2') %>% ungroup() %>% mutate(diff = col_2 - col_1) %>% arrange(desc(diff)) # tells me that in each row the missing values are inblocks of 2 weeks. and no row has two chunks missing.
  
  # check if any row is completely missing
  missing_per_row <- apply(env_data_kibreed_sub[, 3:ncol(env_data_kibreed_sub)], 1, function(x) sum(is.na(x)))/ncol(env_data_kibreed_sub)
  bad_rows <- which(missing_per_row > 0.90)
  
  # impute missing values
  if(length(bad_rows) > 0){
    env_data_kibreed_sub_imp <- env_data_kibreed_sub[-bad_rows, ]  
  } else {
    env_data_kibreed_sub_imp <- env_data_kibreed_sub
  }
  
  
  for(i in 1:nrow(env_data_kibreed_sub_imp)){
    cols <- 3:ncol(env_data_kibreed_sub_imp)
    row_subset <- as.numeric(env_data_kibreed_sub_imp[i, cols])
    ref <- which(!is.na(row_subset))[1]
    impute_val <- mean(row_subset[ref:(ref+7)])
    env_data_kibreed_sub_imp[i, which(is.na(env_data_kibreed_sub_imp[i, ]))] <- impute_val
  }
  
  #sum(is.na(env_data_kibreed_sub))
  #sum(is.na(env_data_kibreed_sub_imp))
  
  # Calculate EC and ERM 
  meta_data <- combined_data_with_env %>% distinct(Connect_at, Site) %>%
    rename(harvest_env = Connect_at) %>%
    mutate(harvest_year = as.numeric(gsub(".*\\_(\\d{4})$", "\\1", harvest_env, perl = T)), 
           sowing_year =  harvest_year - 1,
           sowing_env = sprintf("%s%s", gsub("(.*\\_)\\d{4}$", "\\1", harvest_env, perl = T), sowing_year),
           env = sprintf("Env_%s", row_number())) %>%
    relocate(env, harvest_env)
  
  env_data <- meta_data %>%
    pivot_longer(cols = c("sowing_env", "harvest_env"), names_to = "type", values_to = "connect_at") %>%
    left_join(env_data_kibreed_raw, by = "connect_at") %>%
    mutate(day_of_year = yday(date),
           sowing_leap_year = ifelse(sowing_year %% 4 == 0, TRUE, FALSE),
           harvest_leap_year = ifelse(harvest_year %% 4 == 0, TRUE, FALSE),
           day_from_start_of_sowing_year = ifelse(sowing_year != year & sowing_leap_year == TRUE, day_of_year + 366, 
                                                  ifelse(sowing_year != year & sowing_leap_year == FALSE, day_of_year + 365, day_of_year))) %>%
    relocate(day_of_year, sowing_leap_year, day_from_start_of_sowing_year, .after = "date")
  
  # sanity check for days availability
  #env_data %>% count(env, connect_at) %>% filter(!(n %in% c(365, 366))) # nasenberg has data till 30 of december, so 364 instead of 365, and others in 2021 has data till 30 of november
  #env_data %>% group_by(env, sowing_year) %>% summarize(min = min(day_from_start_of_sowing_year), max = max(day_from_start_of_sowing_year)) # only those env which has leap year in sowing or harvest year have 731 days 
  
  start_leap = yday("2016-10-01")
  start_non_leap = yday("2015-10-01")
  stop_leap = yday("2016-08-31")
  stop_non_leap = yday("2015-08-31") 
  
  env_data_filter <- env_data %>% distinct(env, sowing_year, harvest_year) %>% 
    mutate(start = ifelse(sowing_year %% 4 == 0, start_leap, start_non_leap),
           end = ifelse(sowing_year %% 4 == 0, 366 + stop_non_leap, # harvest year following a leap year will always by non leap
                        ifelse(sowing_year %% 4 != 0 & harvest_year %% 4 == 0, 365 + stop_leap, # harvest year after a non leap year can be leap or non leap
                               ifelse(sowing_year %% 4 != 0 & harvest_year %% 4 != 0, 365 + stop_non_leap, NA))),
           end_fixed = start + 329) %>%
    select(env, start, end, end_fixed)
  
  intervals <- seq(0, 330, 30)
  
  env_data_filtered <- env_data %>% left_join(env_data_filter, by = "env") %>%
    filter(day_from_start_of_sowing_year >= start & day_from_start_of_sowing_year <= end_fixed) %>% # used end fixed for binning days from sowing into 30 days periods
    group_by(env) %>%
    mutate(days_from_sowing  = 1:n(),
           days_from_sowing_intervals = cut(days_from_sowing, breaks = intervals)) %>%
    ungroup() %>%
    select(all_of(c(colnames(env_data_kibreed_raw), "env", "Site", "days_from_sowing", "days_from_sowing_intervals"))) %>%
    relocate(env, Site, connect_at, days_from_sowing, days_from_sowing_intervals, .after = date) %>%
    select(-short_wave_radiation_min..W.m.2.) %>%
    select(-loc) # loc in env data is similar to site in pheno data
  
  # Calculate g_x_e kinship
  
  weather_variables <- colnames(env_data_filtered)[9:ncol(env_data_filtered)]
  
  EC <- env_data_filtered %>% 
    pivot_longer(cols = all_of(weather_variables), names_to = "variable", values_to = "value") %>%
    filter(!is.na(value)) %>% # remove missing
    group_by(env, variable, days_from_sowing_intervals) %>%
    summarize(val_mean = mean(value), .groups = "drop") %>%
    mutate(EC = sprintf("%s_%s", variable, days_from_sowing_intervals)) %>%
    pivot_wider(id_cols = "env", names_from = EC, values_from = val_mean) %>%
    left_join(meta_data %>% select(env, harvest_env), by = "env") %>%
    relocate(harvest_env) %>%
    as.data.frame()
  
  EC_scaled <- scale(EC[, c(-1, -2)], scale = T, center = T)
  
  rownames(EC_scaled) <- EC[, 1]
  
  #EC_miss <- apply(EC_scaled, 2, function(x) sum(is.na(x)))
  
  #EC_scaled_no_miss <- EC_scaled[, which(EC_miss == 0)]
  
  EC_scaled_no_miss <- as.matrix(EC_scaled)
  
  ERM <- tcrossprod(EC_scaled_no_miss)
  
  ERM_data <- list()
  
  ERM_data[["linear"]] <- ERM/(sum(diag(ERM))/nrow(ERM))
  
  guassian_kernel <- function (x, h = NULL) {
    d <- as.matrix(dist(x, upper = TRUE, diag = TRUE))^2
    q <- median(d)
    if (is.null(h)) 
      h <- 1
    return(exp(-h * d/q))
  }
  
  ERM_data[["non_linear"]] <- guassian_kernel(x = EC_scaled_no_miss)
  
  # Calculate site kinship
  EC_s <- env_data_filtered %>% 
    pivot_longer(cols = all_of(weather_variables), names_to = "variable", values_to = "value") %>%
    filter(!is.na(value)) %>% # remove missing
    group_by(Site, variable, days_from_sowing_intervals) %>%
    summarize(val_mean = mean(value), .groups = "drop") %>%
    mutate(EC = sprintf("%s_%s", variable, days_from_sowing_intervals)) %>%
    pivot_wider(id_cols = "Site", names_from = EC, values_from = val_mean) %>%
    as.data.frame()
  
  EC_s_scaled <- scale(EC_s[, -1], scale = T, center = T)
  rownames(EC_s_scaled) <- EC_s[, 1]
  SRM <- tcrossprod(as.matrix(EC_s_scaled))
  SRM <- SRM/(sum(diag(SRM))/nrow(SRM))
  
  # Calculate year kinship
  EC_y <- env_data_filtered %>%
    pivot_longer(cols = all_of(weather_variables), names_to = "variable", values_to = "value") %>%
    filter(!is.na(value)) %>% # remove missing
    group_by(year, variable) %>%
    summarize(val_mean = mean(value), .groups = "drop") %>%
    mutate(EC = sprintf("%s", variable)) %>%
    pivot_wider(id_cols = "year", names_from = EC, values_from = val_mean) %>%
    as.data.frame()
  
  EC_y_scaled <- scale(EC_y[, -1], scale = T, center = T)
  rownames(EC_y_scaled) <- EC_y[, 1]
  YRM <- tcrossprod(as.matrix(EC_y_scaled))
  YRM <- YRM/(sum(diag(YRM))/nrow(YRM))
  
  rm(list = setdiff(ls(), c( "coord_mod", "combined_data_with_env", "env_data_filtered", 
                            "env_data_kibreed_raw", "ERM_data", "EC", "SRM", "YRM", "out", "log_at")))
  
  # sanity check
  # env_data_filtered %>% group_by(env) %>% summarize(sowing_date = min(date), harvest_date = max(date))
  cat("processing completed and file written",
      file = sprintf("%s/preprocessing_environ_data.log", log_at),
      append = T,
      sep = "\n")
  
  # Generate out file -------------------------------------------------------
  out[["coord_mod"]] <- coord_mod
  out[["combined_data_with_env"]] <- combined_data_with_env
  out[["env_data_kibreed_raw"]] <- env_data_kibreed_raw
  out[["ec_mat"]] <- EC
  out[["ERM_data"]] <- ERM_data
  out[["SRM"]] <- SRM
  out[["YRM"]] <- YRM
  
  return(out)
}
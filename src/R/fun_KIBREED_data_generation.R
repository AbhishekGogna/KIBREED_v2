record_raster_grob <- function(width, height, dpi = 300, expr){
  t <- tempfile(fileext = ".png")
  on.exit(unlink(t), add = TRUE)
  png(t, width = width, height = height, res = dpi, units = "cm")
  tryCatch( expr , finally = dev.off() )
  grid::rasterGrob(png::readPNG(t))
}

predict_HD_incpl_env <- function(existing_data, log_at, tmp_data){
  # Core functions
  "%!in%" <- Negate("%in%")
  
  predict_hd <- function(index, A_matrix, D_matrix, pred_info, pheno_data, 
                         names_change_consolidated, write_at, nIter = 12000, 
                         burnIn = 2000, thin = 5){
    
    ## Generate pheno_data
    to_predict <- pred_info %>% filter(order == index)
    run <- to_predict %>% pull(run) %>% unique()
    fold <- to_predict %>% pull(fold) %>% unique()
    hybrids <- to_predict %>% pull(hybrids) %>% unique()
    env <- to_predict %>% pull(Env) %>% unique()
    
    p_data <- pheno_data %>% filter(Env == env) %>%
      select(idx, Geno_new, Geno_dedup, Type, HD) %>%
      left_join(names_change_consolidated %>% select(orig, connect_pheno), by = c("Geno_new" = "connect_pheno")) %>% 
      group_by(idx) %>% slice_head(n = 1) %>% ungroup() # takes the first match if a genotype has multiple genotype hits since this will most likey be the actual genptype and not the deduplicated genotype
    
    missing_prop <- sum(is.na(p_data$HD))/nrow(p_data)
    
    ## Filter markers based on if there are hybrids or inbred in the population
    
    #if(hybrids){
    #  parents <- p_data %>% filter(Type != "Hybrid") %>% pull(orig)
    #  maf <- apply(A_matrix[parents, ], 2, function(x) mean(x)/2)
    #} else{
    #  maf <- apply(A_matrix, 2, function(x) mean(x)/2)
    #}
    #
    #useful_markers <- names(maf[which(maf > 0.05 & maf < 0.95)]) # no need to filter since the input matrix is already filtered
    
    A_matrix_fil <- A_matrix[p_data %>% pull(orig) %>% as.vector(), ]
    D_matrix_fil <- D_matrix[p_data %>% pull(orig) %>% as.vector(), ]
    
    if(!dir.exists(write_at)){
      dir.create(write_at)
    }
    
    if(!dir.exists(paste0(write_at, "/dump"))){
      dir.create(paste0(write_at, "/dump"))
    }
    
    if(!dir.exists(paste0(write_at, "/output"))){
      dir.create(paste0(write_at, "/output"))
    }
    
    # run model
    start.time <- Sys.time()
    
    ETA <- list(A_matrix = list(X = A_matrix_fil, model = 'BRR', saveEffects = TRUE),
                D_matrix = list(X = D_matrix_fil, model = 'BRR', saveEffects = TRUE))
    
    fm_base <- BGLR(y = p_data$HD,
                    ETA = ETA,
                    nIter = nIter,
                    burnIn = burnIn,
                    thin = thin,
                    saveAt = paste0(write_at, "/dump/run_", run, "_fold_", fold, "_"),
                    verbose = TRUE)
    
    end.time <- Sys.time()
    
    pred <- fm_base$yHat
    out <- p_data %>% bind_cols("pred" = pred) %>% mutate(run = run, fold = fold)
    train_set <- out %>% filter(!is.na(HD))
    accuracy <- cor(train_set$HD, train_set$pred)
    
    logs <- data.frame(
      "G_a" = fm_base$ETA[[1]]$varB,
      "G_d" = fm_base$ETA[[2]]$varB,
      "ke" = fm_base$varE,
      "missing_prop" = missing_prop,
      "accuracy" = accuracy,
      "runtime" = difftime(end.time, start.time, units = "min"),
      "log_at" = paste0(write_at, "/output/run_", run, "_fold_", fold, ".qs"),
      "res_at" = paste0(write_at, "/output/run_", run, "_fold_", fold, ".qs")
    )
    
    qsave(out, paste0(write_at, "/output/run_", run, "_fold_", fold, ".qs"))
    qsave(logs, paste0(write_at, "/output/run_", run, "_fold_", fold, ".log"))
    
    return(logs)
  }
  
  # Generate output
  out <- list()
  
  # Put a log file
  run_instance <- as.character(format(Sys.time(),  format = "%d_%m_%Y_%H_%M"))
  log_at <- sprintf("%s/%s", log_at, run_instance)
  
  if(!dir.exists(log_at)){dir.create(log_at, recursive = T)}
  
  # Put a temp_dir
  tmp_data <- sprintf("%s/%s", tmp_data, run_instance)
  
  if(!dir.exists(tmp_data)){dir.create(tmp_data, recursive = T)}
  
  cat("Predicting HD for environments where part of it is missing ------------------",
      file = sprintf("%s/KIBREED_data_generation.log", log_at),
      sep = "\n")
  
  # Sequester data
  data_mod <- existing_data[["pheno_data"]]
  GNdata_comb <- existing_data[["geno_data"]][["GNdata_comb"]]
  GNdata_comb_2 <- existing_data[["geno_data"]][["GNdata_comb_2"]]
  names_change_consolidated <- existing_data[["meta_geno"]]
  
  # generate overviews
  
  missing_overview <- data_mod %>% group_by(Env) %>% 
    summarize(total = n(), miss_prop_HD = sum(is.na(HD))/total, .groups = "drop")
  completely_missing <- missing_overview %>% filter(miss_prop_HD == 1)
  none_missing <- missing_overview %>% filter(miss_prop_HD == 0)
  partially_missing <- missing_overview %>% filter(miss_prop_HD < 1, miss_prop_HD > 0)
  miss_range <- round(range(partially_missing$miss_prop_HD)*100, 2)
  
  cat(sprintf("%s percent to %s percent", miss_range[1], miss_range[2]),
      file = sprintf("%s/KIBREED_data_generation.log", log_at),
      sep = "\n",
      append = T)
  
  # check genomic data availablity
  env_to_impute <- partially_missing %>% pull(Env) %>% as.vector()
  to_impute_data <- data_mod %>% filter(Env %in% env_to_impute) %>% mutate(idx = 1:n())
  to_impute_data_overview <- to_impute_data %>%
    group_by(Env) %>% summarize(total_geno = length(Geno_new),
                                geno_with_gdata = sum(Geno_new %in% names_change_consolidated$connect_pheno))
  
  predict_info <- to_impute_data %>% distinct(Env, Type) %>% group_by(Env) %>% 
    summarize(hybrids = ifelse("Hybrid" %in% Type, TRUE, FALSE), .groups = "drop") %>% 
    distinct(Env, hybrids) %>% group_by(Env) %>% 
    mutate(run_type = "HD_pred", run = cur_group_id(), fold = cur_group_id(), order = cur_group_id()) %>%
    ungroup()
  
  n_cases <- predict_info %>% pull(order) %>% as.vector()
  
  cat(sprintf("Cases to predict %s", length(n_cases)),
      file = sprintf("%s/KIBREED_data_generation.log", log_at),
      sep = "\n",
      append = T)
  
  # Make predictions
  RhpcBLASctl::blas_set_num_threads(1)
  c1 <- makeCluster(length(n_cases), type = "FORK", outfile = sprintf("%s/%s", tmp_data, "cluster.log"))
  registerDoParallel(c1)
  system.time(out_raw <- foreach(parts = n_cases,
                                 .packages = "BGLR",
                                 .export = "%!in%",
                                 .combine = "rbind") %dopar%
                predict_hd(index = parts,
                           A_matrix = GNdata_comb,
                           D_matrix = GNdata_comb_2,
                           pred_info = predict_info,
                           pheno_data = to_impute_data,
                           names_change_consolidated = names_change_consolidated,
                           write_at = tmp_data)) # with a and d
  stopCluster(c1) # 30 minutes

  #qsave(out_raw, "~/KIBREED/results_plots/HD_pred/overview.qs")
  # Read in output
  data <- do.call(rbind, lapply(out_raw[, 8], function(x) qread(x)))
  data_df <- as_tibble(data) %>% arrange(idx) %>% mutate(HD_pred = ifelse(is.na(HD), pred, HD)) %>%
    select(idx, HD_pred) 
  
  data_mod_predicted <- data_mod %>% filter(Env %!in% env_to_impute) %>%
    bind_rows(to_impute_data %>% left_join(data_df, by = "idx") %>% 
                mutate(HD = ifelse(is.na(HD), HD_pred, HD)) %>% select(-idx, -HD_pred))
  missing_overview_post_pred <- data_mod_predicted %>% group_by(Env) %>% 
    summarize(total = n(), miss_prop_HD = sum(is.na(HD))/total, .groups = "drop")
  
  # Generate output
  rm(list = setdiff(ls(), c("data_mod_predicted", "missing_overview_post_pred", "log_at", "out", "tmp_data", "out_raw")))
  
  cat("HD data predicted. Output written",
      file = sprintf("%s/KIBREED_data_generation.log", log_at),
      sep = "\n",
      append = T)
  
  out[["log_at"]] <- log_at
  out[["tmp_data"]] <- tmp_data
  out[["out_raw"]] <- out_raw
  out[["data_mod_predicted"]] <- data_mod_predicted
  out[["missing_overview_post_pred"]] <- missing_overview_post_pred
  
  return(out)
}

generate_KIBREED_blues <- function(data){
  # Core functions
  "%!in%" <- Negate("%in%")
  
  cmd_overlaps_mat <- function(data) {
    data_overlaps_numbers <- lapply(data, function(x) lengths(x))
    df <- do.call(rbind, data_overlaps_numbers)
    colnames(df) <- rownames(df) <- names(data)
    return(df)
  }
  
  cmd_overlaps <- function(file) {
    output <- sapply(file, 
                     function(x) sapply(file, 
                                        function(y) intersect(x, y),
                                        simplify = FALSE, 
                                        USE.NAMES = TRUE), 
                     simplify = FALSE, 
                     USE.NAMES = TRUE)
    output_df <- cmd_overlaps_mat(output)
    return(output_df)
  }
  
  uq_geno <- function(x, blues.3){
    data <- x 
    data_with_blues <- data %>% left_join(blues.3, by = "Geno_dedup") %>% 
      select(Series, Geno_new, BLUEs, Type) %>% 
      distinct(Geno_new, .keep_all = T) 
    
    return(data_with_blues)
  } # extracting BLUEs for geno_new not dedup
  
  # Sequester data
  log_at <- data[["log_at"]]
  tmp_data <- data[["tmp_data"]]
  kibreed_data <- data[["data_mod_predicted"]]
  #kibreed_data$Env <- as.factor(kibreed_data$Env)
  #kibreed_data$Geno_dedup <- as.factor(kibreed_data$Geno_dedup)
  kibreed_data$Type <- ifelse(kibreed_data$Type != "Hybrid", "Line", "Hybrid")
  # adjust data types
  kibreed_data <- kibreed_data %>% 
    convert(fct(Env, Geno_dedup, Male_dedup, Female_dedup, Type, Series)) %>%
    mutate(BLUES_T = BLUES_dt/10) %>%
    arrange(Env) # this means that the order of genotypes will be different win the final data
  
  missing_overview_post_pred <- data[["missing_overview_post_pred"]]

  asreml.options(maxit = 50,             
                 workspace = "25gb",       
                 pworkspace = "25gb",
                 trace=T)
  env <- asreml(fixed = BLUES_dt ~ Geno_dedup, 
                random = ~ Env,
                data = kibreed_data) #  takes 4.5 hours
  
  env_BLUP_T <- asreml(fixed = BLUES_T ~ Type + Series, 
                        random = ~ Env + 
                          ~at(Type, "Line"):Geno_dedup + 
                          ~at(Type, "Hybrid"):Male_dedup +
                          ~at(Type, "Hybrid"):Female_dedup +
                          ~at(Type, "Hybrid"):Geno_dedup +
                          ~at(Type, "Line"):Geno_dedup:Env + 
                          ~at(Type, "Hybrid"):Male_dedup:Env +
                          ~at(Type, "Hybrid"):Female_dedup:Env +
                          ~at(Type, "Hybrid"):Geno_dedup:Env,
                        residual = ~dsum(~units | Env),
                        data = kibreed_data) #  takes 5 minutes
  
  env_BLUP_dt <- asreml(fixed = BLUES_dt ~ Type + Series, 
                     random = ~ Env + 
                       ~at(Type, "Line"):Geno_dedup + 
                       ~at(Type, "Hybrid"):Male_dedup +
                       ~at(Type, "Hybrid"):Female_dedup +
                       ~at(Type, "Hybrid"):Geno_dedup +
                       ~at(Type, "Line"):Geno_dedup:Env + 
                       ~at(Type, "Hybrid"):Male_dedup:Env +
                       ~at(Type, "Hybrid"):Female_dedup:Env +
                       ~at(Type, "Hybrid"):Geno_dedup:Env,
                     residual = ~dsum(~units | Env),
                     data = kibreed_data) #  takes 0.5 hours
  
  var <- summary(env_BLUP_dt)$varcomp
  varcomp <- var$component
  names(varcomp) <- rownames(var)
  err <- sum(varcomp[grep("!R", rownames(var))])/length(grep("!R", rownames(var)))
  
  test_lines <- kibreed_data %>% filter(Type == "Line") %>% distinct(Env, Geno_dedup) %>% 
    count(Geno_dedup) %>% summarize(mean_env = mean(n)) %>% as.numeric()
  test_hybrid <- kibreed_data %>% filter(Type == "Hybrid") %>% distinct(Env, Geno_dedup) %>% 
    count(Geno_dedup) %>% summarize(mean_env = mean(n)) %>% as.numeric()
  
  line_var <- varcomp["at(Type, Line):Geno_dedup"]
  herit_lines <- as.numeric(line_var/ (line_var + err/test_lines))
  
  hybrid_var <- varcomp["at(Type, Hybrid):Male_dedup"] + varcomp["at(Type, Hybrid):Female_dedup"] + varcomp["at(Type, Hybrid):Geno_dedup"]
  herit_hybrids <- as.numeric(hybrid_var / (hybrid_var + err/test_hybrid))
  
  blues <- data.frame(env$coefficients$fixed)
  blues.1 <- data.frame(rownames(blues),blues)
  
  blues.2 <- data.frame (blues.1[,1],
                         blues.1[,2] + blues.1[which(blues.1[,1]== "(Intercept)"),2]) 
  blues.3 <- blues.2[which(blues.2[,1]!= "(Intercept)"),] 
  
  blues.3$blues.1...1. <- gsub("Geno_dedup_(\\S+)", "\\1", blues.3$blues.1...1., perl = T)
  colnames(blues.3) <- c("Geno_dedup", "BLUEs")
  #blues.3$unique_idx <- 1:nrow(blues.3) # check_1check_1 in exp_7 has name premiopremio. effetevely what is happening is that during left_join the premio is getting the value of premio as present in Geno_dedup. Should i correct it?
  
  BLUES_for_series <- kibreed_data %>% mutate(Series = paste0(Series, "_",Project),
                                              group_fct = Series) %>% 
    group_by(group_fct) %>% 
    group_map(~uq_geno(.x, blues.3))
  
  names(BLUES_for_series) <- unlist(lapply(BLUES_for_series, function(x) unique(x$Series)))
  
  Geno_connectivity <- cmd_overlaps(lapply(BLUES_for_series, function(x) x$Geno_new))
  
  BLUES_for_series <- do.call(rbind, BLUES_for_series)
  BLUES_for_series_0 <- BLUES_for_series %>% distinct(Geno_new) %>% mutate(unique_idx = row_number())
  BLUES_for_series <- BLUES_for_series %>% mutate(idx_with_series = row_number()) %>% 
    left_join(BLUES_for_series_0, by = "Geno_new")
  
  # generate output
  out <- list()
  out[["log_at"]] <- log_at
  out[["BLUE_model"]] <- env
  out[["BLUP_model"]] <- env_BLUP_dt
  out[["BLUP_model_T"]] <- env_BLUP_T
  out[["BLUES_acr_env"]] <- as.data.frame(BLUES_for_series)
  out[["BLUES_within_env"]] <- kibreed_data
  out[["missing_overview_post_pred"]] <- missing_overview_post_pred
  
  return(out)
}

write_KIBREED_data_full <- function(existing_data,
                                    data,
                                    write_path_for_R,
                                    write_path_for_Py){
  
  generate_geno_connection <- function(pheno_data, connect_data){
    ## define connection with marker_data
    geno_data <- c()
    for(i in pheno_data$Geno_new){
      geno_data_name <- connect_data$orig[which(connect_data$connect_pheno == i)][1]
      geno_data <- c(geno_data, geno_data_name)
    }
    return(geno_data)
  }
  
  # core function
  prep_as_df_for_feather <- function(data){
    data_df <- cbind("idx" = rownames(data), data)
    data_df <- as.data.frame(data_df)
    return(data_df)
  }
  
  # Sequester data
  log_at <- data[["log_at"]]
  
  # Produce data
  pheno_acr_recurrent_genotypes <- data[["BLUES_acr_env"]] %>% count(Geno_new) # add a filtering step here if needed
  pheno_acr <- data[["BLUES_acr_env"]] %>% 
    filter(Geno_new %in% pheno_acr_recurrent_genotypes$Geno_new) 
  pheno_acr$connect_geno_data <- generate_geno_connection(pheno_acr, existing_data[["meta_geno"]])
  
  pheno_wtn_recurrent_genotypes <- data[["BLUES_within_env"]] %>% count(Geno_new) # add a filtering step here if needed 
  pheno_wtn <-   data[["BLUES_within_env"]] %>% 
    filter(Geno_new %in% pheno_wtn_recurrent_genotypes$Geno_new)
  pheno_wtn$connect_geno_data <- generate_geno_connection(pheno_wtn, existing_data[["meta_geno"]])
  
  ## Notes from revision of data. 
  # 30.03.24 -  I found that for gabi wheat genomic duplication was not done correctly. i corrected it and did all the steps again till here. 
  # The results were checked for identity to past results 
  # BLUES_within_env - only few values in HD changed. the difference is numeric accuracy less than one so i do not change anything related to this data
  ### check <- KIBREED_data_full$BLUES_within_env # from the function output with updated results for GABI wheat dedup
  ### check_2 <- qread("/proj/ext_dir/qg_10_backups/zfs-auto-snap_daily-2024-03-29-0226/AGR-QG/Gogna/KIBREED_v2/results/R/KIBREED_data_generation/BLUES_within_env.qs") # results i saved from the past
  ### identical(check[, c(1:7, 11, 13:25, 29)], check_2[, c(1:7, 11, 13:25, 29)]) # return TRUE
  # BLUES_acr_env - This changed the most and results will be generated again. 
  # check <- KIBREED_data_full$BLUES_acr_env
  # check_2 <- qread("/proj/ext_dir/qg_10_backups/zfs-auto-snap_daily-2024-03-29-0226/AGR-QG/Gogna/KIBREED_v2/results/R/KIBREED_data_generation/BLUES_acr_env.qs")
  # identical(check[, c(1:2, 4:7)], check_2[, c(1:2, 4:7)])
  # correlation between BLUEs - 0.97
  # all others - No change
  
  pheno_wtn_old <- qread("/proj/results/R/KIBREED_data_generation/old_res/BLUES_within_env.qs")
  
  # Generate output
  out <- list()
  out[["BLUES_within_env"]] <- pheno_wtn_old
  out[["BLUES_acr_env"]] <- pheno_acr
  out[["climate_data"]] <- existing_data[["env_data"]][["env_data_kibreed_raw"]]
  out[["GNdata_comb_add"]] <- existing_data[["geno_data"]][["GNdata_comb"]]
  out[["GNdata_comb_dom"]] <- existing_data[["geno_data"]][["GNdata_comb_2"]]
  out[["ERM_data_linear"]] <- existing_data[["env_data"]][["ERM_data"]][["linear"]]
  out[["ERM_data_non_linear"]] <- existing_data[["env_data"]][["ERM_data"]][["non_linear"]]
  out[["SRM"]] <- existing_data[["env_data"]][["SRM"]]
  out[["YRM"]] <- existing_data[["env_data"]][["YRM"]]
  out[["ec_mat"]] <- existing_data[["env_data"]][["ec_mat"]]

  if(!dir.exists(write_path_for_R)){dir.create(write_path_for_R, recursive = T)}
  if(!dir.exists(write_path_for_Py)){dir.create(write_path_for_Py, recursive = T)}
  
  for (i in names(out)){
    data_name <- i
    save_path_r <- sprintf("%s/%s.qs", write_path_for_R, data_name)
    save_path_py <- sprintf("%s/%s.feather", write_path_for_Py, data_name)
    if (!file.exists(save_path_r)) qs::qsave(out[[i]], save_path_r)
    if (!file.exists(save_path_py)) feather::write_feather(prep_as_df_for_feather(out[[i]]), save_path_py)
    print(sprintf("done for %s", i))
  }

  
  # Put a log file
  if(!dir.exists(log_at)){dir.create(log_at, recursive = T)}
  
  cat("Generating outputs for further work ------------------",
      file = sprintf("%s/KIBREED_data_generation.log", log_at),
      sep = "\n", 
      append = TRUE)
  
  return(out)
}

make_plots_and_tables <- function(data,
                                  pco_data,
                                  pco_eig,
                                  write_path_for_R,
                                  write_path_for_Py){
  # Sequester data -------------------------------------------------------------
  log_at <- data[["log_at"]]
  pheno_wtn <- data[["BLUES_within_env"]]
  ERM_linear <- data[["ERM_data_linear"]]
  ec_mat <- data[["ec_mat"]] %>%
    mutate(Site = gsub("^\\S+\\_\\S+\\_(\\S+)\\_\\d+", "\\1", harvest_env),
           Year = gsub("^\\S+\\_\\S+\\_\\S+\\_(\\d+)", "\\1", harvest_env),
           abb_1 = gsub("(\\S{4}).*", "\\1", Site), 
           abb_2 = gsub("\\d{2}(\\d+)", "\\1", Year),
           abb = paste0(abb_1, abb_2)) %>%
    group_by(abb) %>%
    mutate(ids_mod = ifelse(row_number() != 1, paste0(abb, "_", row_number()), abb)) %>%
    ungroup() %>%
    select(-abb_1, -abb_2, -abb)
  data_acr <- data[["BLUES_acr_env"]] %>% 
    mutate(Type = ifelse(Type != "Hybrid", "Lines", "Hybrid"),
           Series = gsub("(Exp_\\d).*", "\\1", Series, perl = TRUE))
  
  colors_years <- ec_mat %>% distinct(Year)
  n_years <- nrow(colors_years)
  #coul <- brewer.pal(n = 12, name = "Paired")[1:n_years]
  coul <- c("#6EB5FF", "#1E5DAB", "#7FFF7F", "#2E8B57", "#FF6347", "#CD5C5C", "#BA55D3", "#800080", "#FFD700", "#DAA520", "#FF6961")
  colors_years$color <- coul
  colors_years <- colors_years[order(as.numeric(colors_years$Year)),] 
  
  colors_years$Year_fct <- factor(colors_years$Year, levels = colors_years$Year)
  colors_years$color_fct <- factor(colors_years$color, levels = colors_years$color)
  
  ec_mat <- ec_mat %>% left_join(colors_years %>% select(Year, color), by = "Year")
  
  # Check for existing directories ---------------------------------------------
  if(!dir.exists(write_path_for_R)){dir.create(write_path_for_R, recursive = T)}
  if(!dir.exists(write_path_for_Py)){dir.create(write_path_for_Py, recursive = T)}
  
  # Make plots -----------------------------------------------------------------
  overview <- pheno_wtn %>% distinct(Connect_at, Site, Loc, Year) %>%
    mutate(lat = as.numeric(gsub("(\\S+)_\\S+_\\S+_\\S+", "\\1", Connect_at)),
           long = as.numeric(gsub("\\S+_(\\S+)_\\S+_\\S+", "\\1", Connect_at)),
           loc = gsub("\\S+_\\S+_(\\S+)_\\S+", "\\1", Connect_at),
           year = Year) %>%
    left_join(colors_years, by = "Year")
  ec_mat <- ec_mat %>% filter(harvest_env %in% overview$Connect_at)
  
  ## Site Map ------------------------------------------------------------------
  rectangle <- c(left = 1, bottom = 47, right = 18, top = 55)
  
  #bw_map_2 <- get_stadiamap(rectangle, zoom = 6, maptype = "stamen_toner_lite", force = TRUE)
  bw_map_2 <- qread("/proj/results/R/KIBREED_data_generation/stadia_map.qs")
  cluster_map <- ggmap(bw_map_2) +
    geom_point(aes(x = long, y = lat, color = color),
               data = overview, size = 2, shape = 13) +
    labs(x = "Longitude", y = "Latitude") +
    guides(color = guide_legend(title = "Year")) +
    scale_colour_identity(guide = "legend",
                          labels = levels(colors_years$Year_fct),
                          breaks = levels(colors_years$color_fct))+
    theme_classic()+
    theme(legend.spacing.y = unit(0.25, "cm"),
          legend.box.margin = margin(0,0,0,0, "pt"),
          legend.key.size = unit(0.4, "cm"),
          plot.margin = margin(5,5,0,1, "pt"),
          legend.position = "right")
  
  ggsave(cluster_map, filename = sprintf("%s/%s", write_path_for_R, "cluster_map_physical.png"),
         width = 11, height = 8.4, units = "cm", dpi = 600, bg = "white")
  
  ## Site PCoA -----------------------------------------------------------------
  cmd <- cmdscale(dist(ERM_linear), eig = T, k = 3)
  comps <- 100*(cmd$eig[1:3]/sum(cmd$eig)) 
  plot_data_0 <- tibble("x" = cmd$points[, 1],
                       "y" = cmd$points[, 2],
                       "y2" = cmd$points[, 3],
                       "Connect_at" = rownames(cmd$points)) %>% 
      right_join(overview, by = "Connect_at") %>%
      select(Connect_at, x, y, y2,Site, Year, color) %>%
      convert(fct(Site, Year)) %>%
      filter(!is.na(x) & !is.na(y))
  
  #outliers <- c("54.38339_10.468389_Schmoel_2015", "54.38339_10.468389_Schmoel_2013", "54.40000_9.850000_Harzhof_2012", "54.40000_9.850000_Harzhof_2013")
  
  #outliers_pc <- plot_data_0 %>% filter(Connect_at %in% outliers)
  
  plot_data_1 <- plot_data_0  %>%
    ggplot(aes(x, y, color = color)) +
    geom_point(size = 2, shape = 13) +
    labs(x = paste0("PCo1 ", as.character(round(comps[1], 2)), "%"), 
         y = paste0("PCo2 ", as.character(round(comps[2], 2)), "%")) +
    guides(color = guide_legend(title = "Year")) +
    scale_colour_identity(guide = "legend",
                          labels = levels(colors_years$Year_fct),
                          breaks = levels(colors_years$color_fct)) +
    theme_classic(base_size = 10) +
    theme(legend.spacing.y = unit(0, "pt"),
          legend.box.margin = margin(0,0,0,0, "pt"),
          legend.key.size = unit(0.5, "cm"),
          #legend.position = "none",
          plot.margin = margin(0,0,0,0, "pt")) +
    coord_fixed(xlim = c(-5, 5), ylim = c(-5, 5))
  
  plot_data_2 <- plot_data_0  %>%
    ggplot(aes(x, y2, color = color)) +
    geom_point(size = 2, shape = 13) +
    labs(x = paste0("PCo1 ", as.character(round(comps[1], 2)), "%"), 
         y = paste0("PCo3 ", as.character(round(comps[3], 2)), "%")) +
    guides(color = guide_legend(title = "Year")) +
    scale_colour_identity(guide = "legend",
                          labels = levels(colors_years$Year_fct),
                          breaks = levels(colors_years$color_fct))+
    theme_classic(base_size = 10)+
    theme(legend.spacing.y = unit(0, "pt"),
          legend.box.margin = margin(0,0,0,0, "pt"),
          legend.key.size = unit(0.5, "cm"),
          #legend.position = "none",
          plot.margin = margin(0,0,0,0, "pt"))  +
    coord_fixed(ylim = c(-5, 5), xlim = c(-5, 5))
  
  plot_data <- ggarrange(plot_data_1, plot_data_2, legend = "right",
                         ncol = 2, nrow = 1, align = "hv", common.legend = TRUE)
  
  ggsave(plot_data_1, filename = sprintf("%s/%s", write_path_for_R, "env_pca_plot.png"),
         width = 8.4, height = 7.4, units = "cm", dpi = 600, bg = "white")
  
  ## Site phylogenetic analysis ------------------------------------------------
  ec_mat_plot <- ec_mat %>% select(-harvest_env, -env, -color, -Year, -Site, -ids_mod) %>% 
    as.data.frame() 
  rownames(ec_mat_plot) <- ec_mat$ids_mod
  ec_mat_plot_scaled <- scale(ec_mat_plot, center = T, scale = T)
  
  tree_data <- ec_mat_plot_scaled %>% 
    dist(method = "euclidean") %>% 
    hclust() %>% as.dendrogram() %>%
    set("branches_lwd", c(2,1,2)) %>%
    set("branches_lty", c(1,1,3,1,1,2)) %>%
    set("labels_cex", c(.9,1.2))
  label_order = data.frame("ids_mod" = tree_data %>% labels) 
  colors_to_use <- label_order %>% 
    left_join(ec_mat %>% distinct(harvest_env, ids_mod, color), 
              by = "ids_mod") %>%
    #mutate(color = ifelse(harvest_env %in% outliers, "#000000", color)) %>%
    pull(color) %>% as.vector()
  labels_colors(tree_data) <- colors_to_use
  
  #ggd1 <- as.ggdend(tree_data)

  #ggplot(ggd1, labels = FALSE) + scale_y_reverse(expand = c(0.2, 0)) + coord_polar(theta="x")
  
  cic_dend <-  record_raster_grob(8.4, 8.4, 600, expr = {
    par(mfrow=c(1,1), mar = c(2,2.25,2.1,2), bg="white", xpd = NA, cex = 0.7) #(b,l,t,r)
    circlize_dendrogram(tree_data, labels = TRUE)
  })
  
  ggsave(wrap_elements(cic_dend), filename = sprintf("%s/%s", write_path_for_R, "env_dend.png"),
         width = 8.4, height = 8.4, units = "cm", dpi = 600, bg = "white")
  
  ## Genetic divergence Fst plot -----------------------------------------------
  parent_path <- "/proj/ext_dir/KIBREED/dump/fst/"
  fst_val <- "weighted"
  files_0 <- list.files(parent_path)
  files <- grep("kibreed_gmat.m_m_0.5.Beagle.fst\\S+\\.log", files_0, value = T, perl = T)
  mat_out <- matrix(NA, 7, 7, dimnames = list(paste0("Exp_", 1:7), paste0("Exp_", 1:7)))
  for(i in files){
    data <- readLines(paste0(parent_path, i))
    coord_row <- as.numeric(gsub(".*fst\\.(\\d+)\\.(\\d+)\\.log", "\\1", i, perl = T))
    coord_col <- as.numeric(gsub(".*fst\\.(\\d+)\\.(\\d+)\\.log", "\\2", i, perl = T))
    mat_out[coord_row, coord_col] <- as.numeric(gsub(paste0(".*", fst_val," Fst estimate\\: (-?\\d+)"), 
                                                     "\\1", 
                                                     grep(paste0(fst_val," Fst"), data, value = T), 
                                                     perl = T))
  }
  
  mat_out[upper.tri(mat_out)] <- t(mat_out)[upper.tri(t(mat_out), diag = F)]
  diag(mat_out) <- 0
  
  # with dendextend
  dend_ob <- as.dendrogram(hclust(as.dist(mat_out))) %>%
    set("branches_lwd", 2) %>%
    set("labels_cex", 0.5) %>%
    hang.dendrogram()
  all_col <- brewer.pal(n = 8, name = "Dark2")
  dend_label_order = data.frame("ids_mod" = c(dend_ob %>% labels, "multi"),
                                color = all_col)
  labels_colors(dend_ob) <- dend_label_order$color[1:7]
  fst_plot <- record_raster_grob(8.4, 4, 300, expr = {
    par(mfrow=c(1,1), mar = c(1.5,0,1,2), bg="white", xpd = NA, cex = 1.5, tcl = NA, mgp = c(0, 0, 0))
    dend_ob %>% 
      color_branches(col = dend_label_order$color[1:7]) %>%
      plot(cex.axis=0.5, horiz = T)})
  
  ggsave(wrap_elements(fst_plot), filename = sprintf("%s/%s", write_path_for_R, "fst_dend.png"),
         width = 8.4, height = 8.4, units = "cm", dpi = 600, bg = "white")
  
  ## Genetic PCoA plot ---------------------------------------------------------
  dend_label_ordered <- dend_label_order %>% arrange(ids_mod)
  
  pco_data_col <- pco_data %>% left_join(dend_label_order, by = c("Series" = "ids_mod"))
  pc_eig <- round(100*(pco_eig[1:2]/sum(pco_eig)), 2)
  
  pco_data_1 <- pco_data_col %>% filter(!in_core) %>%
    filter(Series != "Exp_6")
  
  pco_data_2 <- pco_data_col %>% filter(!in_core) %>%
    filter(Series == "Exp_6")
  
  pco_data_3 <- pco_data %>% filter(in_core)
  
  pco_plot <- pco_data_2 %>%
    ggplot(aes(x = V1, y = V2, color = color))+
    geom_point(size = 0.3, alpha = 0.5) +
    geom_point(aes(x = V1, y = V2, color = color), size = 0.3, alpha = 0.5,
               data = pco_data_1) +
    geom_point(aes(x = V1, y = V2), color = "red", size = 0.3, shape = 4,
               data = pco_data_3) +
    labs(x = paste0("PCo1 ", pc_eig[1], " %"), 
         y = paste0("PCo2 ", pc_eig[2], " %")) +
    scale_colour_identity(guide = "legend",
                          labels = dend_label_ordered$ids_mod,
                          breaks = dend_label_ordered$color) +
    coord_fixed(xlim = c(-0.2, 0.3), ylim = c(-0.2, 0.3)) +
    theme_classic(base_size = 10) +
    theme(legend.position = "bottom",
          legend.spacing.x = unit(0, "pt"),
          legend.spacing.y = unit(0, "pt"),
          legend.box.margin = margin(t = -10, r = 0, b = 0, l = 0, "pt"),
          legend.key.size = unit(0.5, "cm"),
          plot.margin = margin(0,0,0,0, "pt")) +
    guides(colour = guide_legend(override.aes = list(size=2),
                                 title = "Series"))
  
  ggsave(pco_plot, filename = sprintf("%s/%s", write_path_for_R, "RD_PCo.png"),
         width = 8.4, height = 8.4, units = "cm", dpi = 600, bg = "white")
  
#  layout <- "
#AABBB
#AADDD
#CCDDD
#CCEEE
#"

  #(plot_data_2 + plot_layout(tag_level = 'new')) 
#  joint_plot <- cluster_map + 
#    (plot_data_1 + coord_fixed()) +
#    wrap_elements(cic_dend) + 
#    #wrap_elements(panel = textGrob('Here are some text')) +
#    (pco_plot + coord_fixed()) +
#    fst_plot +
#    plot_annotation(tag_levels = "a") + 
#    plot_layout(design = layout)
  
  joint_plot <- cluster_map + labs(tag = 'a') +
    #plot_data_1 + 
    pco_plot + labs(tag = 'c') +
    wrap_elements(cic_dend) + labs(tag = 'b') + 
    #wrap_elements(panel = textGrob('Here are some text')) +
    fst_plot + labs(tag = 'd') +
    #plot_annotation(tag_levels = "a") + 
    plot_layout(nrow = 2, ncol = 2, heights = c(1, 1), widths = c(1, 1))
  
  ggsave(joint_plot, filename = sprintf("%s/%s", write_path_for_R, "env_phylo_plot.png"),
         width = 16.8, height = 12, units = "cm", dpi = 600, bg = "white")
  
  # supplementary figure -------------------------------------------------------
  data_acr_means <- data_acr %>% 
    aggregate(BLUEs ~ Series + Type, ., function(i) round(mean(i), 3)) %>%
    mutate(value = ifelse(BLUEs > 1, round(BLUEs, 2), NA))
  
  supp_plot <- data_acr %>% ggplot(aes(x = Type, y = BLUEs, color = Type)) + 
    geom_boxplot() + 
    labs(x = "", y = "Grain yield (Quintals per hectare)") +
    facet_wrap(~ Series) + 
    geom_text(aes_string(label = "value", y = 30), 
              data = data_acr_means, position = position_dodge(1), show.legend = FALSE) + 
    coord_cartesian(ylim = c(20, 110)) + 
    theme_classic()
  
  ggsave(supp_plot, filename = sprintf("%s/%s", write_path_for_R, "BLUEs_plot.png"),
         width = 16.8, height = 16.8, units = "cm", dpi = 600, bg = "white")
  
  # produce output
  out <- list()
  out[["pca_plot"]] <- plot_data
  out[["phylo_plot"]] <- cic_dend
  out[["combined_plot"]] <- joint_plot
  
  return(out)
}


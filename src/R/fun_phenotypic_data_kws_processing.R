# this script has functions to process KWS field trial data

get_pheno_data <- function(path){
  
  # get file data
  
  traits <- c("yield", "heading_date", "plant_height")
  
  files <- grep(".csv", list.files(path), value = T)
  
  pheno_data <- list()
  
  for (i in files){pheno_data[[gsub("(kws\\.)?(\\S+)\\.csv", "\\2", i, perl = T)]] <- read_csv(paste0(path,
                                                                                                      i, 
                                                                                                      col_types = cols(.default = "c")),
                                                                                               guess_max = 8000)}
  #names(pheno_data)
  
  pheno_data$scoring <- pheno_data$scoring %>% mutate(trial = ifelse(trial == "AS1", "AP2", trial))
  
  pheno_data$scoring$trial <- gsub("(\\w\\w)\\d", 
                                   "\\1", 
                                   pheno_data$scoring$trial, 
                                   perl = T) #adjust the bracket to get results for "A" and "B"
  #pheno_data$scoring$location <- paste0(pheno_data$scoring$location, ".", pheno_data$scoring$field)
  #pheno_data$scoring$location <- gsub(" ", ".", pheno_data$scoring$location)
  
  #Index
  
  experiments_names <- unique(pheno_data$scoring$experiment_name)
  
  trials <- unique(pheno_data$scoring$trial)
  
  para_indices <- NULL
  
  for (i in traits){
    for (j in trials){
      feature <- cbind(i, j)
      para_indices <- rbind(para_indices, 
                            feature)
    }
  }
  
  # generate output 
  out <- list()
  out[["pheno_data"]] <- pheno_data
  out[["para_indices"]] <- para_indices
  
  return(out)
}

derive_BLUEs <- function(data, log_at, tag){
  # store output
  out <- list()
  
  # put a log file
  run_instance <- as.character(format(Sys.time(),  format = "%d_%m_%Y_%H_%M"))
  log_at <- sprintf("%s/%s/%s", log_at, tag, run_instance)
  
  if(!dir.exists(log_at)){dir.create(log_at, recursive = T)}

  # core functions
  detect_outlier<-function(arg1,arg2, only_outliers = T)
  {
    source("https://code.bioconductor.org/browse/multtest/raw/RELEASE_3_16/R/mt.basic.R")
    studresid.data <- arg1$resid/sd(arg1$resid, na.rm=TRUE)
    resi <- cbind(residuals(arg1, type="response"))
    medi <- median(resi, na.rm=TRUE)
    MAD<-median((abs(resi-medi)), na.rm=TRUE)
    re_MAD<-MAD*1.4828 
    res_MAD<-resi/re_MAD 
    rawp.BHStud <- 2 * (1 - pnorm(abs(res_MAD)))  
    rawp.BHStud.all <- cbind(arg2, studresid.data, rawp.BHStud)
    test.BHStud<-mt.rawp2adjp(rawp.BHStud,proc=c("Holm"))  # get this function from the source link above
    adjp <- cbind(test.BHStud[[1]][,1])
    bholm <- cbind(test.BHStud[[1]][,2])
    index <- cbind(test.BHStud[[2]]) 
    out_flag <- ifelse(bholm<0.05, "OUTLIER ", ".")  
    BHStud_test <- cbind(adjp,bholm,index,out_flag) 
    BHStud_test2 <- BHStud_test[order(index),]  
    nam <- c("rawp","bholm","index","out_flag")
    colnames(BHStud_test2) <- nam  
    total.m2_data <- cbind(rawp.BHStud.all, BHStud_test2) 
    outliers_BH <- total.m2_data[which(total.m2_data$out_flag!="."),]  
    if (only_outliers == T){
      return(outliers_BH)
    } else {
      return(total.m2_data)
    }
  }
  
  process_outliers <- function(asreml_obj, data){
    
    outliers_un <- detect_outlier(asreml_obj, data) # 12 values
    
    outliers <- detect_outlier(asreml_obj, data, only_outliers = F)
    
    outliers_omit <- outliers_un %>% select(experiment_name, location, replication, genotype, row, col)
    
    data_w_outliers <- data %>% anti_join(outliers_omit, by = colnames(outliers_omit))
    
    return(data_w_outliers)
  }
  
  post_outlier_detection <- function(arg1, arg2, design_el)
  {
    a <- arg1$coefficients$random
    row_sub <- apply(a, 1, function(row) all(row !=0 ))
    a2<-data.frame(a[row_sub,])
    colnames(a2) <- "effect"
    
    arg2$row <- factor(arg2$row)
    arg2$design_effekt_row <- NA
    arg2$design_effekt_col <- NA
    
    for(j in 1:nrow(a2)){
      sp_vec<-strsplit(rownames(a2)[j], "\\_|\\:| ")
      sp_vec
      if("row"%in%sp_vec[[1]] == T){
        loc <- sp_vec[[1]][2]
        exp <- paste0(sp_vec[[1]][5], "_", sp_vec[[1]][6], "_", sp_vec[[1]][7])
        row <- as.numeric(sp_vec[[1]][9])
        out <- cbind(a2$effect[j], loc, exp, row)
        temp<-intersect(intersect(which(arg2$location==loc),
                                  which(arg2$experiment_name==exp)),
                        which(arg2$row==row))
        arg2$design_effekt_row[temp] <- as.numeric(a2[j,1])
      } else if ("col"%in%sp_vec[[1]] == T) {
        loc <- sp_vec[[1]][2]
        exp <- paste0(sp_vec[[1]][5], "_", sp_vec[[1]][6], "_", sp_vec[[1]][7])
        col <- as.numeric(sp_vec[[1]][9])
        out <- cbind(a2$effect[j], loc, exp, col)
        temp<-intersect(intersect(which(arg2$location==loc),
                                  which(arg2$experiment_name==exp)),
                        which(arg2$col==col))
        arg2$design_effekt_col[temp] <- as.numeric(a2[j,1])
      }
    }
    
    arg2 <- arg2 %>% mutate(trait_corrected = trait - (design_effekt_row + design_effekt_col))
    
    return(arg2)
  }
  
  BLUEs_V2 <- function(data, index, heritability = T, 
                       repeatab = T, BLUE = T, resi = T, miss_count = 20, log_at)
  {
    `%!in%` = Negate(`%in%`)
    # Get data ----------------------------------------------------------------
    
    trait <- para_indices[index, "i"]
    trial <- para_indices[index, "j"]
    
    cat(paste0("Calculating BLUEs for trait = ", 
               trait, " and trial = ", trial), 
        file = sprintf("%s/%s_%s.log", log_at, trait, trial),
        sep = "\n")
    
    if (trial %in% c("A0", "AP", "A")){
      data_sub <- data[which(data$trial ==  trial) , c("experiment_name", "location", "replication", "block" ,"genotype", trait)] %>%
        convert(fct(experiment_name, location, replication, block, genotype))
    } else if (trial %in% c("BD", "B")){
      data_sub <- data[which(data$trial ==  trial) , c("experiment_name", "location", "replication", "block" ,"genotype", "row", "col", trait)] %>%
        filter(replication == 1) %>%
        convert(fct(experiment_name, location, replication, block, row, col, genotype)) %>%
        select(-block)
    }
    colnames(data_sub)[length(colnames(data_sub))] <- "trait"
    #trait %>% group_by(id) %>% summarise(n = n()) # shows values per environment
    
    
    # Generate information for repetitions ------------------------------------
    
    reps_per_loc <- data_sub %>% group_by(experiment_name, location) %>% 
      group_map(~ length(unique(.x$replication))) %>% 
      unlist()
    
    missing_per_loc <- data_sub %>% group_by(experiment_name, location) %>% 
      group_map(~ sum(is.na(.x$trait))) %>% 
      unlist()
    
    check_list<-cbind(data_sub %>% group_by(experiment_name, location) %>% summarise(n = n(), .groups = "drop"), 
                      "reps" = reps_per_loc, 
                      "missing" = missing_per_loc)
    
    cat(paste0("Total location and experiment name combinations = ", dim(check_list)[1]), 
        file = sprintf("%s/%s_%s.log", log_at, trait, trial), 
        sep = "\n", 
        append = T)
    
    omit <- check_list %>% 
      filter(missing > miss_count)
    
    if (dim(omit)[1] > 0) {
      cat(paste0("Those with missing values more than ", miss_count, " = ", dim(omit)[1]), 
          file = sprintf("%s/%s_%s.log", log_at, trait, trial), 
          sep = "\n", 
          append = T)
      
      #suppressWarnings(write.table(omit , file = sprintf("%s/%s_%s.log", log_at, trait, trial),
      #                             append = T))
    } else {
      cat(paste0("No trials with missing values more than ", miss_count), 
          file = sprintf("%s/%s_%s.log", log_at, trait, trial), 
          sep = "\n", 
          append = T)}
    
    check_list_2 <- check_list %>% 
      filter(missing < miss_count)
    
    loc.rep <- unique(check_list_2$location[which(check_list_2$reps %in% c(2, 3))]) %>%
      droplevels() %>% 
      as.character()
    
    loc.unrep <- unique(check_list_2$location[which(check_list_2$reps %in% c(1))]) %>% 
      droplevels() %>% 
      as.character()
    
    if (length(loc.rep) > 0){
      cat(paste0("only ", length(loc.rep), " locations included instead of ", length(unique(check_list$location))), 
          file = sprintf("%s/%s_%s.log", log_at, trait, trial), 
          sep = "\n", 
          append = T)
    } else {
      cat(paste0("ALL locations are unreplicated. So all ", length(loc.unrep), " locations were included"), 
          file = sprintf("%s/%s_%s.log", log_at, trait, trial), 
          sep = "\n", 
          append = T)
    }
    
    #Remove locations with missing values
    exclude <- omit %>%
      select(experiment_name, location)
    
    data_sub <- data_sub %>% anti_join(exclude, by = colnames(exclude))
    
    # Step 1------------------------------------------------------
    #Outlier Correction and correction for design effects per trial per location
    asreml.options(maxit = 50,             
                   workspace = "512mb",       
                   pworkspace = "512mb",
                   trace=F,
                   extra = 10)
    
    cat("ASREML options set", 
        file = sprintf("%s/%s_%s.log", log_at, trait, trial), 
        sep = "\n", 
        append = T)
    # For replicated environments
    if (length(loc.rep) >0){
      
      dta.rep <- data_sub
      
      env.1 <- asreml(fixed = trait ~ genotype, 
                      random = ~ location +
                        location:genotype +
                        location:experiment_name +
                        location:experiment_name:replication + 
                        location:experiment_name:replication:block,
                      data = dta.rep
      )
      
      outliers_rep <- detect_outlier(env.1, data_sub)
      
      outliers <- detect_outlier(env.1, data_sub, only_outliers = F)
      
      outliers_omit <- outliers_rep %>% select(experiment_name, location, replication, block, genotype)
      
      dta.rep_w_outliers <- dta.rep %>% anti_join(outliers_omit, by = colnames(outliers_omit))
      
      #replicated environments are passed to next step
      
    }  else if (length(loc.unrep) > 0){
      # For unreplicated data 
      dta.unrep <- data_sub
      
      env.1 <- asreml(fixed = trait ~ genotype,
                      random = ~ location + 
                        location:genotype +
                        location:experiment_name +
                        location:experiment_name:row +
                        location:experiment_name:col,
                      data = dta.unrep)
      
      outliers_unrep <- detect_outlier(env.1, dta.unrep) 
      
      outliers <- detect_outlier(env.1, dta.unrep, only_outliers = F)
      
      outliers_omit <- outliers_unrep %>% select(experiment_name, location, replication, genotype, row, col)
      
      dta.unrep_w_outliers <- dta.unrep %>% anti_join(outliers_omit, by = colnames(outliers_omit)) # 12 values
      
      #correct design effects
      env.2 <- asreml(fixed = trait ~ genotype,
                      random = ~ location + 
                        location:genotype +
                        location:experiment_name +
                        location:experiment_name:row +
                        location:experiment_name:col,
                      data = dta.unrep_w_outliers)
      
      dta.w_design_effects <- post_outlier_detection(env.2, 
                                                     dta.unrep_w_outliers)
      
    }
    
    cat("Step 1 complete", 
        file = sprintf("%s/%s_%s.log", log_at, trait, trial), 
        sep = "\n", 
        append = T)
    
    # Step 2------------------------------------------------------
    # BLUEs within locations
    if (length(loc.rep) > 0){
      #Replicated trials
      repeatability <- rep(NA, length(loc.rep))
      residual <- rep(NA,length(loc.rep))
      repeatability_df <- NULL
      
      for (i in 1:length(loc.rep))
      {
        print(i)
        field.red <- dta.rep_w_outliers[which(dta.rep_w_outliers$location == loc.rep[i]),] %>% 
          droplevels()
        
        env.1 <- asreml(fixed = trait ~  1 + experiment_name, 
                        random = ~ genotype + 
                          experiment_name:replication + 
                          experiment_name:replication:block,
                        data = field.red
        )
        
        residual[i] <- summary(env.1)$varcomp["units!R", "component"]
        
        repeatability[i] <- summary(env.1)$varcomp["genotype","component"]/(summary(env.1)$varcomp["genotype","component"] + 
                                                                              summary(env.1)$varcomp["units!R","component"])
        
        rep <- cbind(loc.rep[i], repeatability[i])
        
        repeatability_df <- rbind(repeatability_df, rep)
        
        env.2 <- asreml(fixed = trait ~  experiment_name + 
                          genotype, 
                        random = ~ experiment_name:replication + 
                          experiment_name:replication:block,
                        data = field.red
        )
        
        Geno <- predict.asreml(env.2, classify = "genotype",)$pvals[,1:2]
        
        colnames(Geno)[2] <- loc.rep[i]
        
        if (i == 1) {BLUES <- Geno} else {BLUES <- merge(BLUES,Geno, all = T)} 
        
        cat(paste0("REP; BLUEs done for ", loc.rep[i]), file = sprintf("%s/%s_%s.log", log_at, trait, trial), 
            append = T, 
            sep = "\n")
      }
      
      colnames(repeatability_df) <- c("Location", "repeatability")
      
      BLUES[,1] <- as.character(BLUES[,1])
      
      dta <-  pivot_longer(BLUES, !genotype, names_to = "loc", values_to = "trait") %>% 
        convert(fct(loc, genotype))
    } else if (length(loc.unrep) > 0) {
      
      for (i in 1:length(loc.unrep)){
        
        field.red <- dta.w_design_effects[which(dta.w_design_effects$location == loc.unrep[i]),] %>% 
          droplevels()
        
        env.1 <- asreml(fixed = trait_corrected ~ 1,
                        random = ~ genotype +
                          experiment_name,
                        data = field.red)
        
        env.2 <- asreml(fixed = trait_corrected ~ genotype,
                        random = ~ experiment_name,
                        data = field.red)
        Geno <- predict.asreml(env.2, classify = "genotype",)$pvals[,1:2]
        
        colnames(Geno)[2] <- loc.unrep[i]
        
        if (i == 1) {BLUES <- Geno} else {BLUES <- merge(BLUES,Geno, all = T)} 
        
      }
      BLUES[,1] <- as.character(BLUES[,1])
      
      dta <-  pivot_longer(BLUES, !genotype, names_to = "loc", values_to = "trait") %>% 
        convert(fct(loc, genotype))
    }
    
    # Step 3------------------------------------------------------
    # Blues across locations
    
    env.3 <- asreml(fixed = trait ~ 1 , 
                    random = ~ genotype + loc,
                    data = dta)
    
    env.4 <- asreml(fixed = trait ~ genotype, 
                    random = ~ loc,
                    data = dta)
    
    Blues.acr <- predict.asreml(env.4, classify = "genotype",)$pvals[, 2]
    
    com_BLUEs <- cbind(BLUES, Blues.acr)
    
    # Calculate heritabilities
    
    sigma.g <- summary(env.3)$varcomp["genotype","component"]
    sigma.e <- summary(env.3)$varcomp["units!R","component"]
    
    if (length(loc.rep) > 0){
      sigma.g.e <- sigma.e - (mean(residual)/2)
      herit <- sigma.g/(sigma.g + sigma.g.e/length(loc.rep) + sigma.e/(2*length(loc.rep))) 
    } else if (length(loc.unrep) > 0) {
      herit <- sigma.g/(sigma.g + sigma.e/(length(loc.unrep))) 
    }
    
    cat(paste0("BLUEs done for all locations"), file = sprintf("%s/%s_%s.log", log_at, trait, trial), 
        append = T, 
        sep = "\n")
    
    cat(paste0("heritability = ", herit), file = sprintf("%s/%s_%s.log", log_at, trait, trial), 
        append = T, 
        sep = "\n")
    
    # Producing outputs -------------------------------------------------------
    
    output <- list()
    
    if (repeatab == T) {
      if (length(loc.rep) > 0){
        output["Repeatability"] <- list(repeatability_df)
      } else {
        repeatability <- NA
        output["Repeatability"] <- repeatability
      }
    } 
    
    if (heritability == T) {
      output[["Heritability"]] <- herit
    }
    
    if (BLUE == T) {
      output[paste0("BLUEs_", trial, "_", trait)] <- list(com_BLUEs)
    }
    
    if (resi == T){
      output[paste0("resi_step_2_", trial, "_", trait)] <- list(cbind(env.4$mf[, c(1:4)], env.4$residuals))
    }
    
    output[["outlier_step_1"]] <- outliers
    
    return(output)
    
  } # keeps in the trial effect. each elementof the out
  
  BLUEs_acr_trials <- function (data, repetibility = T, 
                                BLUE = T, trait = "yield")
  {
    
    # Reshape data ------------------------------------------------------------
    BLUEs_long <- list()
    for (i in names(data)){
      if (length(grep(trait, i)) != 0){
        subset_data <- data[[i]][[i]]
        subset_data$trial <- i
        BLUEs_long[[i]] <- subset_data%>% select(-Blues.acr)
      }
    }
    
    BLUEs_long_comb <- BLUEs_long[[1]] %>%
      bind_rows(BLUEs_long[[2]]) %>%
      bind_rows(BLUEs_long[[3]]) %>%
      relocate(trial) %>%
      convert(fct(genotype, trial))
    
    names_checks <- names(which(table(BLUEs_long_comb$genotype)>1))
    locations <- BLUEs_long_comb %>% select(-trial, -genotype) %>% colnames()
    
    # Fit model ---------------------------------------------------------------
    asreml.options(maxit = 50,             
                   workspace = "512mb",       
                   pworkspace = "512mb",
                   trace=F,
                   extra = 10)
    
    repe <- rep(NA, length(locations))
    names(repe) <- locations
    
    for(i in 1:length(locations)){
      loc <- locations[i]
      field.red <- BLUEs_long_comb %>% select(trial, genotype, all_of(loc))
      n_trials <- as.character(unique(field.red[!is.na(field.red[, loc]), ] %>% pull(trial)))
      if(length(n_trials) == 1){
        #print(paste0("Only one trial for ", loc, ". No correction therefore done for trial effect."))
        Blues_3s <- field.red %>% select(-trial) %>% distinct(genotype, .keep_all = T)
        repe[loc] <- NA
      } else if(length(n_trials) > 1) {
        field.red <- field.red %>% rename(value = all_of(loc))
        env.1 <- asreml(fixed = value ~ 1,
                        random = ~ genotype + trial,
                        data = field.red)
        env.2 <- asreml(fixed = value ~ genotype,
                        random = ~ trial,
                        data = field.red)
        Blues_3s <- predict.asreml(env.2, classify = "genotype",)$pvals[, c(1, 2)]
        colnames(Blues_3s) <- c("genotype", loc)
        sigma.g <- summary(env.1)$varcomp["genotype","component"]
        sigma.e <- summary(env.1)$varcomp["units!R","component"]
        
        repe[loc] <- sigma.g/(sigma.g + sigma.e)
      }
      if(i == 1){BLUES_out <- Blues_3s} else {BLUES_out <- merge(BLUES_out, Blues_3s, by = "genotype")}
    }
    
    #Produce outputs
    
    output <- list()
    
    if (repetibility == T) {
      output[["repetibility"]] <- repe
    }
    
    if (BLUE == T) {
      output[paste0("BLUEs_", trait)] <- list(BLUES_out)  
    }
    
    return(output)
  }
  
  # get data
  pheno_data <- data[["pheno_data"]]
  para_indices <- data[["para_indices"]]
  
  # BLUEs estimation
  
  BLUES_traits_V2 <- list()
  
  #for (i in 1:nrow(para_indices)){
  #  data_out <- NULL
  #  data_out <- BLUEs_V2(pheno_data$scoring, index = i, log_at = log_at)
  #  BLUES_traits_V2[[names(data_out[3])]]  <- data_out
  #}
  
  cl <- makeCluster(dim(para_indices)[1])
  
  registerDoParallel(cl)
  
  system.time(BLUES_traits_V2 <- foreach(parts = 1:dim(para_indices)[1],
                                         .packages = c("tidyverse", 
                                                       "hablar", 
                                                       "asreml")) %dopar% 
                BLUEs_V2(pheno_data$scoring, 
                         index = parts, 
                         log_at = log_at))
  
  stopImplicitCluster()
  
  names <- c()
  
  for (i in 1:length(BLUES_traits_V2)){
    names <- c(names, names(BLUES_traits_V2[[i]][3]))
  }
  
  names(BLUES_traits_V2) <- names

  out[["BLUES_traits_V2"]] <- BLUES_traits_V2
  
  # Calculate BLUEs across sub trials
  out[["BLUEs_acr"]] <- BLUEs_acr_trials(BLUES_traits_V2, trait = "yield")
  
  out[["BLUEs_acr_hd"]] <- BLUEs_acr_trials(BLUES_traits_V2, trait = "heading_date")
  
  out[["data_long"]] <- out[["BLUEs_acr"]]$BLUEs_yield %>% 
    pivot_longer(!genotype,
                 names_to = "location",
                 values_to = "value")
  
  out[["data_long_heading_date"]] <- out[["BLUEs_acr_hd"]]$BLUEs_heading_date %>%
    pivot_longer(!genotype,
                 names_to = "location",
                 values_to = "value")
  return(out)
}

produce_plots <- function(raw_data, processed_data, tag, 
                          results_dir = core_paths[["results_R"]],
                          width = 16.8, height = 16.8, units = "cm", dpi = 600,
                          log_at){
  # store output
  out <- list()
  
  # put a log file
  run_instance <- as.character(format(Sys.time(),  format = "%d_%m_%Y_%H_%M"))
  log_at <- sprintf("%s/%s/%s", log_at, tag, run_instance)
  
  if(!dir.exists(log_at)){dir.create(log_at, recursive = T)}
  
  # core functions
  BLUEs_acr_trials_v1 <- function (data, heritability = T, 
                                   BLUE = T, trait = "yield")
  {
    
    # Reshape data ------------------------------------------------------------
    BLUEs_long <- NULL
    for (i in names(data)){
      if (length(grep(trait, i)) != 0){
        subset_data <- data[[i]][[i]][, c("genotype","Blues.acr")]
        subset_data$trial <- i
        BLUEs_long <- rbind(BLUEs_long, subset_data)
      }
    }
    
    BLUEs_long <- BLUEs_long %>% convert(fct(genotype, trial))
    
    names_checks <- names(which(table(BLUEs_long$genotype)>1))
    
    # Fit model ---------------------------------------------------------------
    asreml.options(maxit = 50,             
                   workspace = "512mb",       
                   pworkspace = "512mb",
                   trace=F)
    env.1 <- asreml(fixed = Blues.acr ~ 1,
                    random = ~ genotype + trial,
                    data = BLUEs_long)
    
    env.2 <- asreml(fixed = Blues.acr ~ genotype,
                    random = ~ trial,
                    data = BLUEs_long)
    
    Blues_3s <- predict.asreml(env.2, classify = "genotype",)$pvals[, c(1, 2)]
    
    sigma.g <- summary(env.1)$varcomp["genotype","component"]
    sigma.e <- summary(env.1)$varcomp["units!R","component"]
    
    Herit <- sigma.g/(sigma.g + sigma.e)
    
    #Produce outputs
    
    Blues_3s_out <- BLUEs_long %>% left_join(Blues_3s, by = "genotype")
    
    output <- list()
    
    if (heritability == T) {
      output[["Heritability"]] <- Herit
      
    }
    
    if (BLUE == T) {
      output[paste0("BLUEs_", trait)] <- list(Blues_3s_out)  
    }
    
    return(output)
  }
  
  repeatability_plot <- function(data, trait, rows = 2)
  {
    # remove traits with no repeatability
    index <- which(is.na(sapply(data, "[[", 1)) == T)
    list <- data[-index]
    
    #Herit data
    herit_data <- data.frame("Trial_trait" = gsub("\\S{6}(.*)", "\\1", names(list), perl = T),
                             "heritability" = sapply(list, "[[", 2), row.names = 1:length(list),
                             "orig_env" = names(list))
    
    # creating labels
    herit_data$labels <- paste0("(h^2 = ", round(herit_data$heritability, 2), ")")
    
    # Create data for plot
    output <- NULL
    for (i in names(list)){
      data <- list[[i]]$Repeatability
      data <- cbind(data, "env" = rep(i, dim(data)[1]))
      output <- rbind(output, data)
    }
    
    output <- as.data.frame(output) %>% 
      left_join(herit_data, by = c("env" = "orig_env")) %>%
      select(-env) 
    
    # Add missing environments to traits
    matrix <- matrix(0, length(unique(output$Trial_trait)), 9)
    colnames(matrix) <- unique(output$Location)
    rownames(matrix) <- unique(output$Trial_trait)
    for (i in rownames(matrix)){
      data <- output %>% filter(Trial_trait == i)
      matrix[i, which(colnames(matrix) %in% data$Location)] <- data$repeatability
    }
    
    output_long <- melt(matrix) 
    colnames(output_long) <- c("Trial_trait", "Location", "Repeatibility")
    output_long <- output_long %>% convert(fct("Trial_trait", "Location"))
    
    # Adding labels for heritabilities
    output_plot <- output_long %>% left_join(herit_data, by = c("Trial_trait")) %>% 
      select(-orig_env) %>% convert(fct(labels, Trial_trait))
    
    output_plot$Repeatibility <- as.numeric(output_plot$Repeatibility)
    
    #output_plot$Trial_trait <- paste0(output_plot$Trial_trait, " ", output_plot$labels)
    
    #Produce plot
    plot <- ggplot(data = output_plot, aes(x = Location, y = Repeatibility, fill = Location))+
      geom_bar(stat="identity", width=0.5)+
      facet_wrap(~Trial_trait, scales = "free_x", nrow = rows)+
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            panel.spacing = unit(0.2, "lines"),
            text = element_text(size = 9))+
      coord_cartesian(xlim = c(0, 9), ylim = c(0, 1))
    return(plot)
  }
  
  
  BLUEs_corrplot <- function(data_list, trial_loc)
  {
    output <- NULL
    for (i in names(data_list)){
      data <- data_list[[i]][[3]]
      data_cor <- data[, -c(1, dim(data)[2])]
      plot_cor <- cor(data_cor, use="complete.obs")
      plot_cor[upper.tri(plot_cor, diag = T)] <- NA
      
      matrix_loc <- matrix(NA, length(trial_loc), length(trial_loc))
      
      all_names <- trial_loc
      
      names <- c(rownames(plot_cor), setdiff(all_names, rownames(plot_cor)))
      
      colnames(matrix_loc) <- rownames(matrix_loc) <- names
      
      for (j in rownames(plot_cor)){
        for (k in colnames(plot_cor)){
          matrix_loc[which(rownames(matrix_loc) == j),
                     which(colnames(matrix_loc) == k)] <- plot_cor[j, k]
        }
      }
      
      check <- melt(matrix_loc)
      check <- check %>% convert(fct(Var1, Var2))
      
      out <- data.frame(check, i)
      output <- rbind(output, out)
    }
    
    a <- ggplot(output %>% na.omit(), aes(x = Var1, y = Var2, fill = value)) +
      geom_tile() +
      geom_text(aes(label=round(value, 2))) +
      scale_fill_gradientn(colours=c("red", "white","darkgreen"),
                           limits=c(0,1), 
                           na.value = "white")+
      facet_wrap(. ~i, scales='fixed', nrow = 3) +
      theme_bw() +
      theme(axis.line        = element_blank(),
            axis.text.x      = element_text(angle = 45, hjust=1),
            axis.text.y      = element_text(hjust=1),
            axis.ticks       = element_blank(),
            strip.background = element_rect(fill = 'white'), 
            aspect.ratio = 1,
            axis.text = element_text(size = 11))+
      labs(x = "", y = "", fill = "Pearson\nCorrelation" )
    
    return(a)
  }
  
  
  residual_plot <- function(data_list)
  {
    output <- NULL
    for (i in names(data_list)) {
      trait_trial = i
      data <- cbind(unname(data_list[[i]][["outlier_step_1"]][, c("genotype", "location", "index", "e", "out_flag")]), trait_trial)
      colnames(data) <- c("geno", "loc", "index", "resi", "flag","triat_trail")
      data <- data %>% arrange(loc)
      data$index_2 <- 1:nrow(data)
      output <- rbind(output, data)
    }
    
    highlight <- output %>% filter(flag == "OUTLIER ")
    
    output$triat_trail <- as.factor(output$triat_trail)
    #str(output)
    
    a <- ggplot(aes(x = index_2, y = resi, color = loc), data = output)+
      geom_point()+
      geom_point(data = highlight,
                 aes(x = index_2, y = resi),
                 color = "red")+
      facet_wrap(~triat_trail, scales = "free")
    
    return(a)
  }
  
  
  cor_plot <- function(data_list)
  {
    traits <- c("yield", "heading_date", "plant_height")
    output <- NULL
    for (i in traits){
      data <- BLUEs_acr_trials_v1(data_list, heritability = F, trait = i)[[1]]
      output <- rbind(output, data)
    }
    
    output$trait <- gsub("\\S{9}(.*)", "\\1", output$trial, perl = T)
    output$trial <- gsub("\\S{6}(\\S{2}).*", "\\1", output$trial, perl = T)
    
    a <- ggplot(data = output %>% na.omit(), aes(x = predicted.value, y = Blues.acr, color = trial))+
      geom_point()+
      facet_wrap(~trait, scales = "free")+
      stat_cor(method = "pearson")
    
    b <- ggplot(output %>% na.omit(), aes(x = predicted.value, fill = trial)) + 
      geom_density(alpha = 0.5) + 
      facet_wrap(~trait, scales = "free")
    
    #c <- ggplot(output %>% na.omit(), aes(x = predicted.value)) + geom_density(alpha = 0.5) + facet_wrap(~trait, scales = "free")
    cor_1 <- output %>% na.omit() %>% filter(grepl("A0", trial)) %>% 
      pivot_wider(!c(trial, Blues.acr), names_from= trait, values_from = predicted.value) %>% 
      ggpairs(columns = c("yield", "heading_date", "plant_height"), title = "A0") 
    #ggsave(filename = "~/KIBREED/results_plots/3_step_BLUEs_A0_cor.png", width = 16.8, height = 16.8, units = "cm", dpi = 600)
    
    cor_2 <- output %>% na.omit() %>% filter(grepl("AP", trial)) %>% 
      pivot_wider(!c(trial, Blues.acr), names_from= trait, values_from = predicted.value) %>% 
      ggpairs(columns = c("yield", "heading_date", "plant_height"), title = "AP") 
    #ggsave(filename = "~/KIBREED/results_plots/3_step_BLUEs_AP_cor.png", width = 16.8, height = 16.8, units = "cm", dpi = 600)
    
    cor_3 <- output %>% na.omit() %>% filter(grepl("BD", trial)) %>% 
      pivot_wider(!c(trial, Blues.acr), names_from= trait, values_from = predicted.value) %>% 
      ggpairs(columns = c("yield", "heading_date", "plant_height"), title = "BD") 
    #ggsave(filename = "~/KIBREED/results_plots/3_step_BLUEs_BD_cor.png", width = 16.8, height = 16.8, units = "cm", dpi = 600)
    
    output <- list()
    
    output[["cor_plot"]] <- ggarrange(b, a, nrow = 3, ncol = 1)
    output[["pairwise_cor_A0"]] <- cor_1
    output[["pairwise_cor_AP"]] <- cor_2
    output[["pairwise_cor_BD"]] <- cor_3
    
    return(output)
  }
  
  # to store outputs
  out <- list()
  
  # get data
  BLUES_traits_V2 <- processed_data[["BLUES_traits_V2"]]
  pheno_data <- raw_data[["pheno_data"]]
  
  out[["rep_plot"]] <- repeatability_plot(BLUES_traits_V2)
  cat(paste0("plot 1 done"), 
      file = sprintf("%s/plots.log", log_at), 
      sep = "\n")
  
  
  out[["BLUEs_cor"]] <- BLUEs_corrplot(BLUES_traits_V2, trial_loc = unique(pheno_data$scoring$location))
  cat(paste0("plot 2 done"), 
      file = sprintf("%s/plots.log", log_at), 
      sep = "\n", 
      append = T)
  
  out[["BLUEs_resi"]] <- residual_plot(BLUES_traits_V2) 
  cat(paste0("plot 3 done"), 
      file = sprintf("%s/plots.log", log_at), 
      sep = "\n", 
      append = T)
  
  #out[["BLUEs_cor_2_3"]] <- cor_plot(BLUES_traits_V2)[["cor_plot"]]
  #cat(paste0("plot 4 done"), 
  #    file = sprintf("%s/plots.log", log_at), 
  #    sep = "\n", 
  #    append = T)
  
  # save_plots
  save_at <- sprintf("%s/phenotypic_data_kws_processing", results_dir)
  if(!dir.exists(save_at)){dir.create(save_at, recursive = T)}
  
  for(i in names(out)){
    file_name <- sprintf("%s/%s_%s.png", save_at, tag, i)
    ggsave(out[[i]], filename = file_name, width = width, height = height, units = units, dpi = dpi)
  }
  return(out)
}


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
#      } else {element_name <- c(element_name, element)
#      location <- c(location, NA)
#      value <- c(value, NA)}
#    }
#     
#  # Create a dataframe
#  dataframe <- data.frame(element_name, location, value)
#  return(dataframe)
#}
#all_rep_2020 <- convert_to_dataframe(lapply(BLUEs_kws_2020_data$BLUES_traits_V2, function(x) x[["Repeatability"]]))
#all_rep_2020$year <- "2020"
#all_rep_2021 <- convert_to_dataframe(lapply(BLUEs_kws_2021_data$BLUES_traits_V2, function(x) x[["Repeatability"]]))
#all_rep_2021$year <- "2021"
#all_rep <- rbind(all_rep_2020, all_rep_2021) %>% filter(!grepl("height", element_name))
#all_rep %>% filter(!is.na(value)) %>% 
#  mutate(trait = gsub("BLUEs_\\S\\S\\_(\\S+)", "\\1", element_name)) %>% group_by(trait) %>%
#  summarise(mean_rep = mean(value), min = min(value), max = max(value))
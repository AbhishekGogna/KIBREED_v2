# Core functions
"%!in%" <- Negate("%in%")

preprocess_geno_data <- function(existing_data,
                                 paths,
                                 log_at){

  # Generate output
  out <- list()
  
  # Put a log file
  run_instance <- as.character(format(Sys.time(),  format = "%d_%m_%Y_%H_%M"))
  log_at <- sprintf("%s/%s", log_at, run_instance)
  
  if(!dir.exists(log_at)){dir.create(log_at, recursive = T)}
  
  cat("Merging pheno and geno data ------------------",
      file = sprintf("%s/preprocessing_geno_data.log", log_at),
      sep = "\n")
  
  # Sequester data
  combined_data_with_env <- existing_data
  
  # Combine the data with genodata ------------------------------------------
  # since i do not include kws 2012 genodata i remove this
  
  kws_2012 <- combined_data_with_env %>% group_by(Series, Project, Exp, Year) %>% 
    count() %>%
    filter(Series == "Exp_6", 
           Project == "KWS",
           Exp == "Inbred",
           Year == "2012")
  
  combined_data_with_env_fil <- combined_data_with_env %>% anti_join(kws_2012, by = colnames(kws_2012)[1:4])
  
  
  to_check <- combined_data_with_env_fil %>% filter(Type != "Hybrid") %>% distinct(Geno_new) %>% pull(Geno_new) %>% as.character()
  
  # load genomat 
  
  imp_mat <- qread(paths[["imputed_matrix"]]) # my imputed matrix
  
  # filter based on MAF
  # Get informative markers
  mono_info <- apply(imp_mat, 2, function(x)length(unique(x)))
  
  to_keep <- names(mono_info[which(mono_info > 1)])
  
  imp_mat <- imp_mat[, to_keep]
  
  # filter maf
  maf <- apply(imp_mat, 2, mean)/2
  maf_to_keep <- names(maf[which(maf>=0.01)])
  
  imp_mat <- imp_mat[, maf_to_keep]
  
  # generate D matrix
  
  execute <- TRUE
  
  if (execute) {# Generate dominance matrix for parents 
    
    imp_mat_changed_coding <- imp_mat -1
    
    imp_mat_dom <- matrix(NA, dim(imp_mat_changed_coding)[1], dim(imp_mat_changed_coding)[2], 
                          dimnames = list(rownames(imp_mat_changed_coding), 
                                          colnames(imp_mat_changed_coding)))
    
    #options(warn = 0)
    for(i in colnames(imp_mat_changed_coding)){
      marker <- imp_mat_changed_coding[, i]
      counts <- table(marker)
      if("-1" %in% names(counts)){
        p11 <- as.numeric(counts["-1"])/nrow(imp_mat_changed_coding)
      } else {p11 = 0}
      
      if("0" %in% names(counts)){
        p12 <- as.numeric(counts["0"])/nrow(imp_mat_changed_coding)
      } else {p12 = 0}
      
      if("1" %in% names(counts)){
        p22 <- as.numeric(counts["1"])/nrow(imp_mat_changed_coding)
      } else {p22 = 0}
      
      common_den <- (p11 + p22 - (p11 - p22)^2)
      
      placeholder_0 <- -(2*(p12*p22)/common_den)
      placeholder_1 <- (4*(p11*p22)/common_den)
      placeholder_2 <- -(2*(p11*p12)/common_den)
      
      marker[which(marker == -1)] <- as.numeric(placeholder_0)
      marker[which(marker == 0)] <- as.numeric(placeholder_1)
      marker[which(marker == 1)] <- as.numeric(placeholder_2)
      
      imp_mat_dom[, i] <- marker
    }
    
    #qsave(imp_mat_dom, "~/KIBREED/results_plots/kibreed_gmat.m_m_0.5_mat_dom.qs")
  } #ref - 
  
  cat("A and D matrices loaded",
      file = sprintf("%s/preprocessing_geno_data.log", log_at),
      sep = "\n",
      append = T)
  
  #imp_mat_dom <- qread("~/KIBREED/results_plots/kibreed_gmat.m_m_0.5_mat_dom.qs")
  
  # ceate and add genodata for composites and mixes
  
  # same genotypes are named different in the genodata since all this data is not coming from one source. We can therefore get the sync up by incorporating the transition data
  
  zucht_trans <- qread(paths[["zucht_transition"]]) # transtion data for zucht project
  
  additional_geno <- data.frame(add_geno = as.vector(do.call(c, zucht_trans[c("mixes", "composite")])),
                                to_combine = c("m808_m810", 
                                               "Apache_Rocky_KWSFerrum", 
                                               "Piko_Henrik_Hermann", 
                                               "malemixearly_Mironovskaja66",
                                               "SUR.99820.SUR_SUR.186",
                                               "sur99820sur186_X910102"),
                                type = c("com", "com", "com", "hyb", "com", "hyb"))
  
  additional_geno <- additional_geno[order(additional_geno$type),]
  
  # for additive matrix
  out_mat <- matrix(NA, nrow(additional_geno),dim(imp_mat)[2], dimnames = list(additional_geno$add_geno, 
                                                                               colnames(imp_mat)))
  
  for (i in 1:6){
    type_geno <- additional_geno[i, "type"]
    name_geno <- additional_geno[i, "add_geno"]
    to_split <- additional_geno[i, "to_combine"]
    if(type_geno == "com"){
      genotype <- imp_mat[str_split(to_split, "_")[[1]], ]
      out_mat[name_geno, ] <- colSums(genotype)/nrow(genotype)
    } else if (type_geno == "hyb"){
      splits <- str_split(to_split, "_")[[1]]
      p1 <- out_mat[splits[1], ]
      p2 <- imp_mat[splits[2], ]
      out_mat[name_geno, ] <- (p1+p2)/2
    }
  }
  
  imp_mat <- rbind(imp_mat, out_mat)
  
  # for dominance matrix
  out_mat_dom <- matrix(NA, nrow(additional_geno),dim(imp_mat_dom)[2], dimnames = list(additional_geno$add_geno, 
                                                                                       colnames(imp_mat_dom)))
  
  for (i in 1:6){
    type_geno <- additional_geno[i, "type"]
    name_geno <- additional_geno[i, "add_geno"]
    to_split <- additional_geno[i, "to_combine"]
    if(type_geno == "com"){
      genotype <- imp_mat_dom[str_split(to_split, "_")[[1]], ]
      out_mat_dom[name_geno, ] <- colSums(genotype)/nrow(genotype)
    } else if (type_geno == "hyb"){
      splits <- str_split(to_split, "_")[[1]]
      p1 <- out_mat_dom[splits[1], ]
      p2 <- imp_mat_dom[splits[2], ]
      out_mat_dom[name_geno, ] <- (p1+p2)/2
    }
  }
  
  imp_mat_dom <- rbind(imp_mat_dom, out_mat_dom)
  
  cat("Genodata for mixtures and synthetics added",
      file = sprintf("%s/preprocessing_geno_data.log", log_at),
      sep = "\n",
      append = T)
  
  # source script at ~/KIBREED/source_data/BigData_yusheng_v1/
  
  hywheat_trans <- data.frame("in_geno" = c("Norin", "Oakley", "Humber", "Scout", "Horizon", "Podium", "Colonia", "m808"),
                              "name_ph_yu" = c("f110", "f006", "f007", "f009", "f010", "f011", "f107", "m006"))
  
  geno_connection <- rbind(zucht_trans$conv_table_zu_final, hywheat_trans)
  #names_change <- rbind(names_change, data.frame("orig" = additional_geno$add_geno,
  #                                               "orig_lower" = NA,
  #                                               "in_pat_data" = NA,
  #                                               "name_ph_yu" = additional_geno$add_geno,
  #                                               "new_name" = additional_geno$add_geno,
  #                                               "name_uni" = additional_geno$add_geno))
  
  
  # create a data frame to hold the transitions 
  names_change <- data.frame("orig" = rownames(imp_mat),
                             "orig_lower" = tolower(rownames(imp_mat)),
                             "orig_dup" = paste0(tolower(rownames(imp_mat)), tolower(rownames(imp_mat))))
  
  to_check_abs <- to_check[which(to_check %!in% names_change$orig_dup)]
  
  # where does these come from?
  
  origin_of_misfits <- combined_data_with_env_fil %>% filter(Type != "Hybrid") %>% 
    filter(Geno_new %in% to_check_abs) %>% 
    distinct(Geno_new, .keep_all = T) %>%
    dplyr::count(Series, Project, Exp, Loc, Year) %>%
    arrange(Year)
  
  
  names_change_consolidated <- names_change %>%
    left_join(geno_connection, by = c("orig" = "in_geno")) %>% # gives a column called names_ph_zu which just corrects for zuchtwert project naming of genotypes
    mutate(connect_pheno = ifelse(is.na(name_ph_yu), orig_dup, paste0(name_ph_yu, name_ph_yu)),
           connect_pheno_undup  = ifelse(is.na(name_ph_yu), orig_lower, name_ph_yu)) # connect_pheno and connect_pheno_undup connects matrix names to Geno_new column
  
  # create a mapping file for fst calculation
  population_str <- names_change_consolidated %>% select(orig, connect_pheno) %>% distinct() %>%
    left_join(combined_data_with_env %>% select(Series, Geno_new) %>% distinct(), by = c("connect_pheno" = "Geno_new")) %>%
    select(orig, Series) %>% distinct()
  
  #qsave(population_str, "~/KIBREED/results_plots/population_str.qs")
  # check if problem is little better now
  
  to_check_abs_2 <- to_check[which(to_check %!in% names_change_consolidated$connect_pheno)] # it improves a lot
  
  # other ways to look at what is actually missing here
  
  origin_of_misfits_2 <- combined_data_with_env_fil %>% filter(Type != "Hybrid") %>% 
    filter(Geno_new %in% to_check_abs_2) %>% 
    distinct(Geno_new, .keep_all = T) %>%
    count(Series, Project, Exp, Loc, Year) %>%
    arrange(Year) # majority lies in kws 2013
  
  # remove further environments based on missing data
  
  geno_absent <- combined_data_with_env_fil %>% 
    filter(Type != "Hybrid") %>% 
    filter(Geno_new %!in% names_change_consolidated$connect_pheno) %>%
    group_by(Project, Env) %>% 
    dplyr::count() %>% 
    dplyr::rename(absent = n)
  
  geno_stats <- combined_data_with_env_fil %>% 
    filter(Type != "Hybrid") %>% 
    filter(Geno_new %in% names_change_consolidated$connect_pheno) %>%
    group_by(Series, Project, Env) %>% 
    dplyr::count() %>% 
    dplyr::rename (present = n) %>% 
    left_join(geno_absent, by = c("Project", "Env")) %>%
    mutate(total = sum(absent, present, na.rm = T),
           missing = ifelse(!is.na(absent), (absent/total), 0)) %>%
    arrange(desc(missing)) %>%
    filter(missing >= 0.1) # tells you which envt is missing how many genotypes from the geno data
  
  # lets remove the phdata rows that does not have Gdata
  
  to_remove_lines <- combined_data_with_env %>% 
    filter(Type != "Hybrid") %>%
    filter(Geno_new %!in% names_change_consolidated$connect_pheno)
  
  to_remove_hybrids <- combined_data_with_env %>% 
    filter(Type == "Hybrid") %>%
    filter(paste0(MaleNew, MaleNew) %!in% names_change_consolidated$connect_pheno | paste0(FemaleNew, FemaleNew)  %!in% names_change_consolidated$connect_pheno)
  
  to_remove <- to_remove_lines %>% bind_rows(to_remove_hybrids)
  
  combined_data_with_env_connected <- combined_data_with_env_fil %>% anti_join(to_remove, by = colnames(to_remove))
  #nrow(combined_data_with_env) == (nrow(to_remove) + nrow(combined_data_with_env_connected))
  
  # check status again 
  removed_status <- combined_data_with_env_fil %>% 
    group_by(Env) %>% count() %>% rename(total = n) %>%
    left_join(combined_data_with_env_connected %>% 
                group_by(Env) %>% count() %>% rename(total_new = n), "Env") %>%
    mutate(removed = total - total_new)
  
  rm(list = setdiff(ls(), c("combined_data_with_env_connected", "names_change_consolidated", 
                            "population_str", "imp_mat", "imp_mat_dom", "log_at", "out")))
  
  cat("Geno data curated. Output written",
      file = sprintf("%s/preprocessing_geno_data.log", log_at),
      sep = "\n",
      append = T)
  
  # generate output
  out[["log_at"]] <- log_at
  out[["population_str"]] <- population_str
  out[["combined_data_with_env_connected"]] <- combined_data_with_env_connected
  out[["names_change_consolidated"]] <- names_change_consolidated
  out[["imp_mat"]] <- imp_mat
  out[["imp_mat_dom"]] <- imp_mat_dom
  
  return(out)
}

generate_hybrid_data <- function(data){

  # Generate output
  out <- list()
  
  # Sequester data
  log_at <- data[["log_at"]]
  combined_data_with_env_connected <- data[["combined_data_with_env_connected"]]
  names_change_consolidated <- data[["names_change_consolidated"]]
  imp_mat <- data[["imp_mat"]]
  imp_mat_dom <- data[["imp_mat_dom"]]
  
  # Put a log file
  
  if(!dir.exists(log_at)){dir.create(log_at, recursive = T)}
  
  cat("Generate hybrid data -------------------",
      file = sprintf("%s/preprocessing_geno_data.log", log_at),
      sep = "\n",
      append = T)
  
  # Generate hybrid data ----------------------------------------------------
  
  # i extract the genotypes of hybrids based on the assumption that geno dedup offers a good approximation to the Geno_new, after using geno dedup i link the respective dedup parent names to the geno new hybrid names since i will eventually use geno new for predictions. But this may be discussed further. 
  
  # generate gneotypes for hybrids
  
  hybrid_data <- combined_data_with_env_connected %>% filter(Type == "Hybrid") %>% select(FemaleNew, MaleNew, Geno_new) %>% distinct()
  
  hybrids_with_parent_data <- hybrid_data %>% filter((MaleNew %in% names_change_consolidated$connect_pheno_undup) & 
                                                       (FemaleNew %in% names_change_consolidated$connect_pheno_undup)) # match 100%
  
  hybrids_with_parent_data <- as.data.frame(hybrids_with_parent_data) %>% as.data.frame()
  
  mle <- c()
  fle <- c()
  
  for(i in 1:nrow(hybrids_with_parent_data)){
    mle <- c(mle, names_change_consolidated$orig[which(names_change_consolidated$connect_pheno_undup == hybrids_with_parent_data[i, "MaleNew"])[1]])
    fle <- c(fle, names_change_consolidated$orig[which(names_change_consolidated$connect_pheno_undup == hybrids_with_parent_data[i, "FemaleNew"])[1]])
  }  
  
  hybrids_with_parent_data$name_gdta_f <- fle
  hybrids_with_parent_data$name_gdta_m <- mle
  
  # for additive matrix
  male_mat <- imp_mat[hybrids_with_parent_data$name_gdta_m, colnames(imp_mat)]
  female_mat <- imp_mat[hybrids_with_parent_data$name_gdta_f, colnames(imp_mat)] # if the genotypes do not exist in GNData then try working with dedup data instead of Geno_new data
  hyb_mat <- (male_mat + female_mat)/2
  rownames(hyb_mat) <- hybrids_with_parent_data$Geno_new
  # append to the imp_mat
  GNdata_comb <- rbind(imp_mat, hyb_mat)
  
  # for dominance matrix
  male_mat_d <- imp_mat_dom[hybrids_with_parent_data$name_gdta_m, colnames(imp_mat_dom)]
  female_mat_d <- imp_mat_dom[hybrids_with_parent_data$name_gdta_f, colnames(imp_mat_dom)] # if the genotypes do not exist in GNData then try working with dedup data instead of Geno_new data
  hyb_mat_d <- (male_mat_d + female_mat_d)/2
  rownames(hyb_mat_d) <- hybrids_with_parent_data$Geno_new
  # append to the imp_mat_dom
  GNdata_comb_2 <- rbind(imp_mat_dom, hyb_mat_d)
  
  cat("Respective hybrid data added",
      file = sprintf("%s/preprocessing_geno_data.log", log_at),
      sep = "\n",
      append = T)
  
  # extend the names_change_consolidated data frame
  names_change_consolidated  <- names_change_consolidated %>% bind_rows(data.frame(
    orig = hybrids_with_parent_data$Geno_new,
    orig_lower = NA,
    orig_dup = NA,
    name_ph_yu = NA,
    connect_pheno = hybrids_with_parent_data$Geno_new,
    connect_pheno_undup = NA
  ))
  
  rm(list = setdiff(ls(), c("GNdata_comb", "GNdata_comb_2", "names_change_consolidated", "out", "log_at")))
  
  cat("Output written",
      file = sprintf("%s/preprocessing_geno_data.log", log_at),
      sep = "\n",
      append = T)
  
  # Generate output
  out[["log_at"]] <- log_at
  out[["GNdata_comb"]] <- GNdata_comb
  out[["GNdata_comb_2"]] <- GNdata_comb_2
  out[["names_change_consolidated"]] <- names_change_consolidated
  return(out)
}

generate_kinships <- function(data, write_at){
  # Sequester data
  log_at <- data[["log_at"]]
  g_data <- data[["GNdata_comb"]]
  g_data_dom <- data[["GNdata_comb_2"]]
  
  # Put a log file
  if(!file.exists(sprintf("%s/preprocessing_geno_data.log", log_at))){file.create(sprintf("%s/preprocessing_geno_data.log", log_at), recursive = T)}
  
  # Put a results folder
  if(!dir.exists(write_at)) dir.create(write_at, recursive = TRUE)
  
  # Compute and save G_a_mat
  cat("Kinship G_a -------------------",
      file = sprintf("%s/preprocessing_geno_data.log", log_at),
      sep = "\n",
      append = T)
  save_at_a <- sprintf("%s/kin_a.qs", write_at)
  g_a_kin <- Gmatrix(as.matrix(g_data), method = "VanRaden", integer = FALSE)
  qsave(g_a_kin, save_at_a)
  
  # Compute and save G_aa_mat
  cat("Kinship G_aa -------------------",
      file = sprintf("%s/preprocessing_geno_data.log", log_at),
      sep = "\n",
      append = T)
  save_at_aa <- sprintf("%s/kin_aa.qs", write_at) 
  g_aa_kin <- g_a_kin * g_a_kin
  qsave(g_aa_kin, save_at_aa)
  
  # Compute and save G_d_mat
  cat("Kinship G_d -------------------",
      file = sprintf("%s/preprocessing_geno_data.log", log_at),
      sep = "\n",
      append = T)
  save_at_d <- sprintf("%s/kin_d.qs", write_at)
  d_mat_0 <- g_data_dom %*% t(g_data_dom)
  g_d_kin <- d_mat_0 / mean(diag(d_mat_0))
  qsave(g_d_kin, save_at_d)
  
  file_list <- c(save_at_a, save_at_aa, save_at_d)
  write_check <- all(file.exists(unlist(file_list)))
  
  # Generate output
  cat("Kinship calculation complete -------------------",
      file = sprintf("%s/preprocessing_geno_data.log", log_at),
      sep = "\n",
      append = T)
  out <- list()
  out[["write_check"]] <- write_check
  out[["kin_paths"]] <- file_list
  
  return(out)
}

add_deduplication_info <- function(existing_data,
                                   paths,
                                   data){
  
  # Generate output
  out <- list()
  
  # Sequester data
  log_at <- data[["log_at"]]
  combined_data_with_env_connected <- existing_data[["combined_data_with_env_connected"]]
  names_change_consolidated <- existing_data[["names_change_consolidated"]]
  GNdata_comb <- data[["GNdata_comb"]]
  GNdata_comb_2 <- data[["GNdata_comb_2"]]
  
  # Put a log file
  if(!dir.exists(log_at)){dir.create(log_at, recursive = T)}
  
  cat("Deduplication -------------------",
      file = sprintf("%s/preprocessing_geno_data.log", log_at),
      sep = "\n",
      append = T)
  
  # Deduplication -----------------------------------------------------------
  
  RD <- qread(paths[["RD_mat"]])
  
  DistM_snp<-1-as.matrix(RD)
  DistM_snp[1:5,1:5]
  Geno<-rownames(RD)
  DupSave<-NULL
  Level<-0.97 # 
  
  for(i in 1:dim(RD)[1]){
    Temp<-array(0,30)
    TT<-which(DistM_snp[i,]>Level)
    if(length(TT)>1){
      if(TT[1]>=i){
        Temp[1:length(c(TT))]<-as.character(Geno[TT])
        DupSave<-rbind(DupSave,Temp)
      }
    }
  }
  
  cat(paste0("Total duplicates found = ", dim(DupSave)[1]),
      file = sprintf("%s/preprocessing_geno_data.log", log_at),
      sep = "\n",
      append = T)
  
  
  
  SaveRG<-matrix(0,dim(DupSave)[1],2)
  colnames(SaveRG) <- c("max_RD_bw_any_pair", "members in the family")
  
  rename_df <- NULL
  for(i in 1:dim(DupSave)[1]){
    Dup <- setdiff(DupSave[i,],0) # members in the family
    PosD <- which(is.element(Geno,Dup)) # locations of the members in the RD matrix
    
    if(length(PosD)<2) {
      cat(print(i),
          file = sprintf("%s/preprocessing_geno_data.log", log_at),
          sep = "\n",
          append = T)
    } # basically those where there is just one member in the family
    
    SaveRG[i,1] <- max(RD[PosD,PosD]) # maximum pairwise distance amongst the family members
    SaveRG[i,2] <- length(PosD) # totat menbers in the family
    
    all_pairs <- RD[Dup, Dup]
    all_pairs[lower.tri(all_pairs, diag = T)] <- NA
    all_pairs_df <- na.omit(melt(all_pairs))
    rownames(all_pairs_df) <- NULL
    
    out_res <- cbind(rep(i, nrow(all_pairs_df)), 
                     rep(DupSave[i, 1], nrow(all_pairs_df)), 
                     all_pairs_df)
    
    colnames(out_res) <- c("family", "source", "G_1", "G_2", "similarity")
    rename_df <- rbind(rename_df, out_res)  
    
    if(SaveRG[i,1]>0.03){
      write.table(round(RD[PosD,PosD],3),
                  file = sprintf("%s/preprocessing_geno_data.log", log_at),
                  append = T)
    } # prints those families where at least one pairwise distance is more than 0.03
  }
  #plot(SaveRG[,1])
  
  ## these are few groups need to check manually
  SaveRG[which(SaveRG[,1]>0.03),]
  
  families_to_look_manually <- rename_df %>% filter(similarity > 0.03) %>% pull(family) %>% as.character() # have to decide something for these 
  
  rename_df_no_conflict <- rename_df %>% filter(family %!in% families_to_look_manually) %>%
    distinct(G_2, source) %>%
    left_join(names_change_consolidated[ ,c("orig", "connect_pheno_undup")], by = c("source" = "orig")) %>%
    left_join(names_change_consolidated[ ,c("orig", "connect_pheno_undup")], by = c("G_2" = "orig"), suffix = c(".s", ".g")) %>%
    distinct(connect_pheno_undup.s, connect_pheno_undup.g)
  
  data_mod <- combined_data_with_env_connected  %>%
    left_join(rename_df_no_conflict, by = c("FemaleNew" = "connect_pheno_undup.g"))  %>%
    left_join(rename_df_no_conflict, by = c("MaleNew" = "connect_pheno_undup.g"), suffix = c("Fem", "Mal")) %>% 
    distinct(Env, Geno_new, BLUES_dt, .keep_all = T) %>%
    mutate(Female_dedup = ifelse(is.na(connect_pheno_undup.sFem), FemaleNew, connect_pheno_undup.sFem),
           Male_dedup = ifelse(is.na(connect_pheno_undup.sMal), MaleNew, connect_pheno_undup.sMal),
           Geno_dedup = paste0(Female_dedup, Male_dedup)) %>% select(-connect_pheno_undup.sFem, -connect_pheno_undup.sMal) 
  
  rm(list = setdiff(ls(), c("data_mod", "log_at", "out")))
  
  # Generate overview
  before_dedup <- data_mod %>% distinct(Series, Geno_new) %>%
    mutate(present = 1) %>%
    pivot_wider(id_cols = "Geno_new", names_from = "Series", values_from = "present") %>%
    select(-Geno_new)
  
  after_dedup <-  data_mod %>% distinct(Series, Geno_dedup) %>%
    mutate(present = 1) %>%
    pivot_wider(id_cols = "Geno_dedup", names_from = "Series", values_from = "present") %>%
    select(-Geno_dedup)
  series <- colnames(before_dedup)
  dedup_res <- matrix(NA, nrow = length(series),
                      ncol = length(series),
                      dimnames = list(series, series))
  diag(dedup_res) <- 0
  coord_before <- which(upper.tri(dedup_res, diag = F), arr.ind = T)
  for(row in 1:nrow(coord_before)){
    row_id <- as.integer(coord_before[row, "row"])
    col_id <- as.integer(coord_before[row, "col"])
    
    data <- before_dedup[, c(rownames(dedup_res)[row_id], 
                             colnames(dedup_res)[col_id])]
    data$sum <- rowSums(data, na.rm = T)
    dedup_res[row_id, col_id] <- length(which(data$sum > 1))
  }
  
  coord_after <- which(lower.tri(dedup_res, diag = F), arr.ind = T)
  for(row in 1:nrow(coord_after)){
    row_id <- as.integer(coord_after[row, "row"])
    col_id <- as.integer(coord_after[row, "col"])
    
    data <- after_dedup[, c(rownames(dedup_res)[row_id],
                            colnames(dedup_res)[col_id])]
    data$sum <- rowSums(data, na.rm = T)
    dedup_res[row_id, col_id] <- length(which(data$sum > 1))
  }
  
  # Generate output
  out[["data_mod"]] <- data_mod
  out[["dedup_res"]] <- dedup_res
  
  cat("Output written",
      file = sprintf("%s/preprocessing_geno_data.log", log_at),
      sep = "\n",
      append = T)
  
  return(out)
}
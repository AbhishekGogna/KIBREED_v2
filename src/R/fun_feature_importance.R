# Core functions
"%!in%" <- Negate("%in%")

create_dir_file <- function(paths, file = TRUE) {
  for (path in paths) {
    if (file) {
      if (!file.exists(path)) {
        tryCatch({
          file.create(path, recursive = TRUE)
        }, error = function(e) {
          message("Error creating file '", path, "': ", e$message)
        })
      }
    } else {
      if (!dir.exists(path)) {
        tryCatch({
          dir.create(path, recursive = TRUE)
        }, error = function(e) {
          message("Error creating directory '", path, "': ", e$message)
        })
      }
    }
  }
  return(paths)
}

cmd_RD_func<-function(SNP, print_status = T){
  ## calculate the matrix of Roger's distance
  ## SNP must as 0,1,2 
  ## Transfer element class from integer/character to numeric
  ## Data may have missing values, use the na.rm function
  nmarker <- ncol(SNP)
  No.Geno<-dim(SNP)[1]             # Auto input number of genotypes
  Geno<-rownames(SNP)
  RD <- matrix(0,No.Geno,No.Geno)
  if(print_status == T){
    print(paste0("Start at = ", Sys.time()))
  }
  system.time(
    for(i in 1:No.Geno)
    {
      if(is.element("2",names(table(SNP[,1:10])))=="FALSE"){
        print("SNP marker may not fit 0,1,2")
        break
      }
      if(length(SNP[i,])!=nmarker){
        print("input matrix is wrong")
        break}
      for (j in 1:No.Geno)
      {
        RD[i,j] <- mean(abs(SNP[i,]-SNP[j,])/2,na.rm=TRUE)
        if(length(SNP[i,])!=nmarker){
          print("input matrix is wrong")
          break}
      }
      if(print_status == T){
        if(i%%1000==0){print(paste("row",i,"finished",Sys.time()))}
      }
    }
  )
  if(print_status == T){
    print(paste0("End at = ", Sys.time()))
  }
  rownames(RD)<-Geno
  colnames(RD)<-Geno
  return(RD)
}

get_product <- function(incidence_mat, kin_mat, mat_names, save_at){
  kin_mat_ordered <- kin_mat[colnames(incidence_mat), colnames(incidence_mat)]
  t0 <- Sys.time()
  mat_eigen <- eigen(kin_mat_ordered)
  mat_eigen$vectors<-mat_eigen$vectors[,mat_eigen$values>1e-8]
  mat_eigen$values<-mat_eigen$values[mat_eigen$values>1e-8]
  mat_PC <- sweep(mat_eigen$vectors,MARGIN=2,STATS=sqrt(mat_eigen$values),FUN='*')
  # mat_PC<-t(apply(mat_eigen$vectors, 1, function(x) x*sqrt(mat_eigen$values)))
  out_mat <- incidence_mat %*% mat_PC
  if (!is.null(mat_names)){
    rownames(out_mat) <- mat_names
  }
  colnames(out_mat) <- paste0("var_", 1:ncol(out_mat))
  t1 <- Sys.time()
  period <- round(as.numeric(t1 - t0, units = "mins"), 2)
  ## generate output
  out <- list()
  out[["mat"]] <- out_mat
  out[["time_min"]] <- period
  
  if (!is.null(save_at)){
    if(!file.exists(save_at)) {
      out_df <- prep_as_df_for_feather(out[["mat"]])
      write_feather(out_df, save_at)
    }
  }
  
  return(out)
}

get_had_prod <- function(mat_1, mat_2, save_at){
  t0 <- Sys.time()
  mat_1_sym <- mat_1 %*% t(mat_1)
  mat_2_sym <- mat_2 %*% t(mat_2)
  mat_had <- mat_1_sym * mat_2_sym
  t1 <- Sys.time()
  period <- round(as.numeric(t1 - t0, units = "mins"), 2)
  ## generate output
  out <- list()
  out[["mat"]] <- mat_had
  out[["time_min"]] <- period
  
  if (!is.null(save_at)){
    if(!file.exists(save_at)) {
      out_df <- prep_as_df_for_feather(out[["mat"]])
      write_feather(out_df, save_at)
    }
  }
  
  return(out)
}

get_incidence_matrix <- function(data, column_name) {
  
  if (length(column_name) > 1){
    subset_0 <- paste(data[, column_name[1]], data[, column_name[2]], sep = ":")
    subset <- unique(subset_0)
    column_name <- paste0(column_name, collapse = ":")
  } else {
    subset <- NULL
  }
  
  formula <- sprintf("~ -1 + %s", column_name)
  t0 <- Sys.time()
  matrix <- model.matrix(as.formula(formula), data)
  
  
  if (is.null(subset)){
    colnames(matrix) <- gsub(column_name, '', colnames(matrix), perl = TRUE)  
  } else {
    part_1 <- gsub(strsplit(column_name, ":")[[1]][1], '', colnames(matrix), perl = TRUE)
    colnames(matrix) <- gsub(strsplit(column_name, ":")[[1]][2], '', part_1, perl = TRUE)
    matrix <- matrix[, subset]
  }
  
  t1 <- Sys.time()
  period <- round(as.numeric(t1 - t0, units = "mins"), 2)
  ## generate output
  out <- list()
  out[["mat"]] <- matrix
  out[["time_min"]] <- period
  return(out)
}

get_BRR_pred <- function(pheno_data, col_name, kin_path, var_path, 
                         tmp_data_dir, dump_dir, 
                         cv_id, nIter, burnIn){
  
  # create tmp dirs
  dump_at_temp <- paste0(dump_dir, "/par_pred")
  if(!dir.exists(dump_at_temp)) dir.create(dump_at_temp, recursive = T)
  
  # create tmp_dir
  tmp_data <- paste0(tmp_data_dir, "/par_pred")
  if(!dir.exists(tmp_data)) dir.create(tmp_data)
  
  pheno_data_0 <- pheno_data %>%
    mutate(idx = row_number())
  
  g_data <- qread(kin_path)
  param_data <- qread(var_path)
  param_data_df <- as.data.frame(param_data, row.names = NA)
  param_data_df[, col_name] <- rownames(param_data)
  
  pheno_data_mod <- pheno_data_0 %>% mutate(obs_set = ifelse(is.na(blues_var), "test", "train")) %>%
    arrange(desc(obs_set)) %>% left_join(param_data_df, by = col_name)
  
  envs <- as.character(unique(pheno_data_mod$env))
  meta_info <- NULL
  param_data_out <- NULL
  
  for (env in envs){
    ifrm(pheno_data_0_site, env = environment())
    pheno_data_0_site <- pheno_data_0 %>% filter(env == !!as.character(env))
    param_data_site <- pheno_data_mod %>% filter(env == !!as.character(env)) %>%
      mutate(across(all_of(colnames(param_data)), ~ifelse(!is.na(blues_var), ., NA), .names = "to_pred_{.col}"))
    geno_site <- param_data_site %>% pull(connect_geno) %>% as.character()
    geno_data_site <- g_data[geno_site, geno_site]
    for (param in colnames(param_data)){
      cols_to_select <- c("idx", colnames(pheno_data), "obs_set", param, paste0("to_pred_", param))
      ifrm(param_data_site_param, env = environment())
      param_data_site_param <- param_data_site %>% select(all_of(cols_to_select)) %>%
        rename("param_raw" = all_of(param),
               "param" = all_of(paste0("to_pred_", param)))
      model_fit <- BGLR(y = param_data_site_param[, "param"],
                        ETA = list(G=list(K=geno_data_site,model='RKHS')),
                        nIter = nIter,
                        burnIn = burnIn,
                        thin = 5,
                        saveAt = sprintf("%s/%s_%s_%s_", dump_at_temp, env, param, cv_id),
                        verbose = FALSE)
      param_data_site_param$pred_param <- model_fit$yHat
      
      meta <- data.frame("env" = env,
                         "param" = param,
                         "cv_id" = cv_id,
                         "var_G" = model_fit$ETA$G$varU,
                         "var_E" = model_fit$varE,
                         "rep" = model_fit$ETA$G$varU/(model_fit$ETA$G$varU+model_fit$varE),
                         "train_set_cor_avg" = param_data_site_param %>% 
                           filter(obs_set == "train") %>% group_by(env) %>% 
                           summarize(corr_v = cor(param_raw, pred_param)) %>% ungroup() %>% 
                           pull(corr_v) %>% as.numeric() %>% mean(),
                         "test_set_cor_avg" =  param_data_site_param %>% 
                           filter(obs_set == "test") %>% group_by(env) %>% 
                           summarize(corr_v = cor(param_raw, pred_param)) %>% ungroup() %>% 
                           pull(corr_v) %>% as.numeric() %>% mean())
      
      meta_info <- rbind(meta_info, meta)
      
      param_data_site_param <- param_data_site_param %>% mutate(param = ifelse(is.na(param), pred_param, param)) %>%
        select(-param_raw, -pred_param, -obs_set) %>%
        relocate(param, .after = last_col())
      colnames(param_data_site_param)[ncol(param_data_site_param)] <- param
      pheno_data_0_site <- pheno_data_0_site %>% left_join(param_data_site_param, by = c("idx", colnames(pheno_data))) # binds columns
    }
    if(is.null(param_data_out)) {
      param_data_out <- pheno_data_0_site
    } else {
      param_data_out <- param_data_out %>% bind_rows(pheno_data_0_site) # binds rows
    }
  }
  
  # remove temp_data
  system(sprintf("rm -rf %s", dump_at_temp))
  
  # derive_output
  param_data_out_mat_0 <- param_data_out %>% select(-all_of(c(colnames(pheno_data), "idx")))
  param_data_out_mat <- as.matrix(param_data_out_mat_0)
  rownames(param_data_out_mat) <- param_data_out[, col_name]
  out <- list("mat" = param_data_out_mat)
  
  # write_data
  write_feather(param_data_out, sprintf("%s/%s.feather", tmp_data, cv_id))
  write.table(meta_info, sprintf("%s/meta_%s.csv", tmp_data, cv_id))
  
  return(out)
}

get_mat <- function(pheno_data, col_name, BRR_mat){
  
  # calulate kinship
  param_data_df <- as.data.frame(BRR_mat, row.names = NA)
  param_data_df[, col_name] <- rownames(BRR_mat)
  
  pheno_data_mod <- pheno_data %>% 
    mutate(obs_set = ifelse(is.na(blues_var), "test", "train")) %>%
    arrange(desc(obs_set)) %>% left_join(param_data_df, by = col_name) %>% 
    filter(obs_set == "train")
  param_data_kin <- pheno_data_mod %>% 
    select(all_of(colnames(BRR_mat)), "latlong", "geno") %>% 
    group_by(latlong) %>% 
    summarize(across(all_of(colnames(BRR_mat)), ~ mean(as.numeric(.x), na.rm = TRUE)))
  
  param_data_kin_scaled <- scale(param_data_kin[, -1], scale = T, center = T)
  rownames(param_data_kin_scaled) <- param_data_kin %>% pull(latlong) %>% as.vector()
  param_data_kin_scaled_miss <- as.matrix(param_data_kin_scaled)
  
  RM <- tcrossprod(param_data_kin_scaled_miss)
  RM_scaled <- RM/(sum(diag(RM))/nrow(RM))
  
  # calculate incidence
  inc_mat <- get_incidence_matrix(pheno_data, "latlong")[["mat"]]
  
  # generate output
  if(any(colnames(RM_scaled) %in% colnames(inc_mat))) {
    mat_out <- get_product(inc_mat, RM_scaled, NULL, NULL)
  } else {
    mat_out <- NULL
  }
  return(mat_out)
}

create_and_get_product <- function(col_name, pheno_data, kin_mat_path) {
  incidence_mat <- get_incidence_matrix(data = pheno_data,
                                        column_name = col_name)
  kin_mat <- qread(kin_mat_path)
  mat_out <- get_product(incidence_mat[["mat"]], kin_mat, NULL, NULL)
  return(mat_out)
}

wrapper_functions <- list(
  wrap_inc = function(pheno_data, matrix_params) {
    col_name <- matrix_params[["col_name"]]
    out_0 <- get_incidence_matrix(pheno_data, col_name)
    out <- list(X = out_0[["mat"]])
    return(out)
  },
  
  wrap_BRR_kin = function(pheno_data, matrix_params, paths) {
    col_name <- matrix_params[["col_name"]]
    kin_path <- paths[[matrix_params[["path"]]]]
    incidence_mat <- get_incidence_matrix(data = pheno_data,
                                          column_name = col_name)
    kin_mat <- qread(kin_path)
    out_0 <- get_product(incidence_mat[["mat"]], kin_mat, NULL, NULL)
    out <- list(X = out_0[["mat"]])
    return(out)
  },
  
  wrap_RKHS = function(pheno_data, matrix_params, paths) {
    
    matrix_params <- matrix_params[1:2]
    
    matrix_params_1 <- sapply(matrix_params, function(x) x[[1]], USE.NAMES = TRUE, simplify = FALSE)
    matrix_params_2 <- sapply(matrix_params, function(x) x[[2]], USE.NAMES = TRUE, simplify = FALSE)
    
    matrix1 <- create_and_get_product(matrix_params_1[["col_name"]], pheno_data, paths[[matrix_params_1[["path"]]]])
    matrix2 <- create_and_get_product(matrix_params_2[["col_name"]], pheno_data, paths[[matrix_params_2[["path"]]]])
    
    out_0 <- get_had_prod(matrix1[["mat"]], matrix2[["mat"]], NULL)
    out <- list(K = out_0[["mat"]])
    return(out)
  },
  
  wrap_BRR = function(pheno_data, matrix_params, paths) {
    col_name <- matrix_params[["col_name"]]
    BRR_mat <- qread(paths[[matrix_params[["path"]]]])
    out <- list(X = BRR_mat[pheno_data[, col_name], ])
    return(out)
  },
  
  wrap_BRR_2 = function(pheno_data, matrix_params, paths) {
    col_name <- matrix_params[["col_name"]]
    kin_path <- paths[[matrix_params[["path"]][1]]]
    var_path <- paths[[matrix_params[["path"]][2]]]
    tmp_data_dir <- matrix_params[["tmp_data_dir"]]
    dump_dir <- matrix_params[["dump_dir"]]
    cv_id  <- matrix_params[["cv_id"]]
    out_0 <- get_BRR_pred(pheno_data, col_name, kin_path, var_path, 
                          tmp_data_dir, dump_dir, 
                          cv_id, nIter = 15000, burnIn = 2000)
    out <- list(X = out_0[["mat"]])
    return(out)
  },
  
  wrap_BRR_3 = function(pheno_data, matrix_params, paths) {
    col_name <- matrix_params[["col_name"]]
    BRR_mat <- qread(paths[[matrix_params[["path"]]]])
    out_0 <- get_mat(pheno_data, col_name, BRR_mat)
    out <- list(X = out_0[["mat"]])
    return(out)
  },
  
  wrap_RKHS_1 = function(pheno_data, matrix_params, paths) {
    
    matrix_params <- matrix_params[1:2]
    
    matrix_params_1 <- sapply(matrix_params, function(x) x[[1]], USE.NAMES = TRUE, simplify = FALSE)
    matrix_params_2 <- sapply(matrix_params, function(x) x[[2]], USE.NAMES = TRUE, simplify = FALSE)
    
    matrix1 <- create_and_get_product(matrix_params_1[["col_name"]], pheno_data, paths[[matrix_params_1[["path"]]]])
    
    col_name <- matrix_params_2[["col_name"]]
    BRR_mat <- qread(paths[[matrix_params_2[["path"]]]])
    matrix2 <- get_mat(pheno_data, col_name, BRR_mat) 
    
    out_0 <- get_had_prod(matrix1[["mat"]], matrix2[["mat"]], NULL)
    out <- list(K = out_0[["mat"]])
    return(out)
  }
)

return_elements <- function(element, pheno_data, paths, 
                            tmp_data_dir, dump_dir, cv_id) {
  matrix_mapping <- list(
    "E_i@BRR"        = list(col_name = "connect_climate", 
                            type = "BRR", wrapper = "wrap_inc"),
    "G_i@BRR"        = list(col_name = "connect_geno", 
                            type = "BRR", wrapper = "wrap_inc"),
    "S_i@BRR"        = list(col_name = "site", 
                            type = "BRR", wrapper = "wrap_inc"),
    "Y_i@BRR"        = list(col_name = "year", 
                            type = "BRR", wrapper = "wrap_inc"),
    "G_a_S_i@BRR"    = list(col_name = c("connect_geno", "latlong"), 
                            type = "BRR", wrapper = "wrap_inc"),
    "S@BRR"          = list(col_name = "site", path = "SRM",
                            type = "BRR", wrapper = "wrap_BRR_kin"),
    "Y@BRR"          = list(col_name = "year", path = "YRM",
                            type = "BRR", wrapper = "wrap_BRR_kin"),
    "G_a@BRR"        = list(col_name = "connect_geno", path = "G_a_RM",
                            type = "BRR", wrapper = "wrap_BRR_kin"),
    "G_d@BRR"        = list(col_name = "connect_geno", path = "G_d_RM",
                            type = "BRR", wrapper = "wrap_BRR_kin"),
    "G_aa@BRR"       = list(col_name = "connect_geno", path = "G_aa_RM",
                            type = "BRR", wrapper = "wrap_BRR_kin"),
    "ERM_l@BRR"      = list(col_name = "connect_climate", path = "ERM_l",
                            type = "BRR", wrapper = "wrap_BRR_kin"),
    "ERM_nl@BRR"     = list(col_name = "connect_climate", path = "ERM_nl",
                            type = "BRR", wrapper = "wrap_BRR_kin"),
    "G_a_ERM_l@RKHS" = list(col_name = c("connect_geno", "connect_climate"),
                            path = c("G_a_RM", "ERM_l"), 
                            type = "RKHS", wrapper = "wrap_RKHS"),
    "G_d_ERM_l@RKHS" = list(col_name = c("connect_geno", "connect_climate"),
                            path = c("G_d_RM", "ERM_l"),
                            type = "RKHS", wrapper = "wrap_RKHS"),
    "G_a_ERM_nl@RKHS"= list(col_name = c("connect_geno", "connect_climate"),
                            path = c("G_a_RM", "ERM_nl"),
                            type = "RKHS", wrapper = "wrap_RKHS"),
    "G_d_ERM_nl@RKHS"= list(col_name = c("connect_geno", "connect_climate"),
                            path = c("G_d_RM", "ERM_nl"),
                            type = "RKHS", wrapper = "wrap_RKHS"),
    "G_a_Y@RKHS"      = list(col_name = c("connect_geno", "year"), 
                             path = c("G_a_RM", "YRM"),  
                             type = "RKHS", wrapper = "wrap_RKHS"),
    "G_a_S@BRR"      = list(col_name = "connect_param", path = "G_a_S",
                            type = "BRR", wrapper = "wrap_BRR"),
    "G_a_S_p@BRR"    = list(col_name = "connect_param", 
                            path = c("G_a_RM", "G_a_S"),
                            type = "BRR", wrapper = "wrap_BRR_2"),
    "S_al@BRR"       = list(col_name = "connect_param", path = "G_a_S", 
                            type = "BRR", wrapper = "wrap_BRR_3"),
    "G_a_S_al@RKHS"  = list(col_name = c("connect_geno", "connect_param"), 
                            path = c("G_a_RM", "G_a_S"),
                            type = "RKHS", wrapper = "wrap_RKHS_1")
  )
  
  matrix_params <- matrix_mapping[[element]] # get needed params
  
  wrapper <- matrix_params[["wrapper"]]
  wrapper_function <- wrapper_functions[[wrapper]]
  if (element == "G_a_S_p@BRR"){
    matrix_params[["tmp_data_dir"]] <- tmp_data_dir 
    matrix_params[["dump_dir"]] <- dump_dir 
    matrix_params[["cv_id"]] <- cv_id
  }
  
  if (wrapper == "wrap_inc"){
    mat_out <- wrapper_function(pheno_data, matrix_params)
  } else {
    mat_out <- wrapper_function(pheno_data, matrix_params, paths)
  }
  
  output_0 <- list(model = matrix_params[["type"]],
                   saveEffects = TRUE,
                   name = as.character(element))
  
  output <- c(output_0, mat_out)
  
  print(paste("ETA element for column =", element, "defined"))
  
  return(output)
}

calculate_pred <- function(paths, pheno_data, index, log_at, tmp_data){
  # define paths
  create_dir_file(paste0(tmp_data, "/tmp"), file = FALSE)
  create_dir_file(paste0(tmp_data, "/dump"), file = FALSE)
  
  tmp_data_dir <- paste0(tmp_data, "/tmp/", index, "_")
  dump_dir <- paste0(tmp_data, "/dump/", index, "_")
  
  # define matrices
  elements <- str_split(index, "&")[[1]]
  
  ETA <- lapply(elements, return_elements,
                pheno_data = pheno_data, 
                paths = paths,
                tmp_data_dir = tmp_data_dir,
                dump_dir = dump_dir,
                cv_id = "vars_pred")
  names(ETA) <- sapply(ETA, function(x) x[["name"]])
  
  fm_base <- BGLR(y = pheno_data$blues,
                  ETA = ETA,
                  nIter = 15000,
                  burnIn = 2000,
                  thin = 5,
                  saveAt = dump_dir,
                  verbose = FALSE)
  
  pheno_data$pred <- as.numeric(fm_base$yHat)
  
  ## vars and sd
  out_var_sd <- NULL
  out_var_sd_names <- NULL
  for (i in elements){
    if (i %in% names(fm_base[["ETA"]])){
      fit <- fm_base[["ETA"]][[i]]
      interest_names <- grep("var\\S", names(fit), value = T, perl = T)
      sd_id <- grep("SD.", interest_names)
      sd <- fit[[interest_names[sd_id]]]
      var <- fit[[interest_names[-sd_id]]]
      out_var_sd <- c(out_var_sd, var, sd)
      out_var_sd_names <- c(out_var_sd_names, paste0(i, c("_var", "_sd")))
    } else {
      var <- NA
      sd <- NA
      out_var_sd <- c(out_var_sd, var, sd)
      out_var_sd_names <- c(out_var_sd_names, paste0(i, c("_var", "_sd")))
    }
  }
  out_var_sd <- c(out_var_sd, fm_base$varE)
  out_var_sd_names <- c(out_var_sd_names, "err_var")
  
  combined_vars_sd <- as.data.frame(rbind(out_var_sd), row.names = run_name)
  colnames(combined_vars_sd) <- out_var_sd_names
  
  # Write a section in the log file
  cat("Done with predicting data",
      file = log_at,
      sep = "\n",
      append = TRUE)
  
  # produce output
  out <- list()
  out[["vars"]] <- combined_vars_sd 
  out[["pred"]] <- pheno_data
  return(out)
}

create_cluster_plot <- function(mds.df, wss_values, label, write_at, n_clust, labs, col_vector) {
  kmclusters <- wss_values[[n_clust]][["clusters"]]
  mds.df$groups <- kmclusters
  mds.df$env <- rownames(mds.df)
  
  clust_plot <- ggscatter(mds.df,
                          x = "V1", 
                          y = "V2", 
                          xlab = sprintf("PC1 = %s%%", as.character(round(labs[1], 1))),
                          ylab = sprintf("PC2 = %s%%", as.character(round(labs[2], 1))),
                          label = label,
                          color = "groups", 
                          #palette = col_vector[1:n_clust], 
                          size = 1, 
                          ellipse = TRUE, 
                          ellipse.type = "convex", 
                          repel = TRUE,
                          ggtheme = theme_gray(),
                          title = sprintf("Clusters = %s", n_clust))
  
  ggsave(plot = clust_plot, filename = sprintf("%s/clust_plot_%s_%s.png", write_at, n_clust, ifelse(is.null(label), "no_label", label)), 
         width = 16.8, height = 16.8, units = "cm", dpi = 600)
  return(clust_plot)
}

# Task functions
get_RD_mat <- function(existing_data,
                       write_at,
                       log_at,
                       tmp_at,
                       key){
  # Put a log file
  run_instance <- as.character(format(Sys.time(),  format = "%d_%m_%Y_%H_%M"))
  log_at <- sprintf("%s/%s", log_at, run_instance)
  log_file <- sprintf("%s/%s.log", log_at, key)
  create_dir_file(log_at, file = FALSE)
  
  # Put a temp_dir
  tmp_data <- sprintf("%s/%s", tmp_at, run_instance)
  
  create_dir_file(tmp_data, file = FALSE)
  
  cat(sprintf("Logs for task: %s ------------------", key),
      file = log_file,
      sep = "\n")
  
  # Put a results dir
  create_dir_file(write_at, file = FALSE)
  
  # Load data
  g_data <- qread(existing_data[["g_data"]])
  
  # Write a section in the log file
  cat("GN data loaded",
      file = log_file,
      sep = "\n",
      append = TRUE)
  
  # Produce output
  save_at <- sprintf("%s/%s/RD_data.qs", core_paths[["results_R"]], run_name)
  if(!file.exists(save_at)){
    RD <- cmd_RD_func(g_data)
    qsave(RD_data, save_at)
  } else {
    RD <- qread(save_at)
  }
  
  # Write a section in the log file
  cat("RD data generated",
      file = log_file,
      sep = "\n",
      append = TRUE)
  
  # Check availability
  check <- sum(rownames(RD) == row.names(g_data)) == nrow(g_data)
  
  if(check) {
    RD_out <- RD
  } else {
    RD_out <- NULL
  }

  # return output
  out <- list()
  out[["log_file"]] <- log_file
  out[["tmp_data"]] <- tmp_data
  out[["write_at"]] <- write_at
  out[["RD_out"]] <- RD_out
  return(out)
}

get_core_set <- function(existing_data, data) {
  # Sequester data
  log_file <- data[["log_file"]]
  tmp_data <- data[["tmp_data"]]
  write_at <- data[["write_at"]]
  RD_out <- data[["RD_out"]]
  p_wtn <- qread(existing_data[["p_wtn"]])
  
  # Write a section in the log file
  cat("Starting with core sampling",
      file = log_file,
      sep = "\n",
      append = TRUE)
  
  # Produce output: Subset RD_mat for genotypes of interest
  geno <- p_wtn %>% distinct(Env, connect_geno) %>%
    mutate(present = 1) %>%
    pivot_wider(id_cols = "connect_geno", names_from = "Env", values_from = "present") %>%
    rowwise() %>%
    mutate(freq = sum(c_across(where(is.numeric)), na.rm = TRUE))
  
  overview <- table(cut(geno$freq, breaks = seq(0, 120, 10)))
  
  RD_sub <- RD_out[geno$connect_geno, geno$connect_geno]
  
  core_obj <- list(
    objective("EN", "PD", weight = 1)
    #, objective("AN", "GD", weight = 1)
  ) # multiple objectives (custom weight)
  
  set.seed(12345)
  core <- sampleCore(data = distances(RD_sub),
                     obj = core_obj, # default. can include more than one, in which case normalization is done for combined results
                     size = 500,
                     #always.selected = integer(0), #geno vector to be always be selected in the core collection
                     #never.selected = integer(0), #geno vector to be never be selected in the core collection
                     mode = "default",
                     normalize = TRUE, # For single-objective configurations, this argument is ignored.
                     time = NA,
                     impr.time = NA, # Defaults to ten or two seconds, mode = default or fast
                     steps = NA,
                     impr.steps = NA,
                     indices = FALSE,
                     verbose = FALSE)
  
  # Write a section in the log file
  cat("Done with core sampling",
      file = log_file,
      sep = "\n",
      append = TRUE)
  
  # Generate output
  out <- list()
  out[["log_file"]] <- log_file
  out[["tmp_data"]] <- tmp_data
  out[["write_at"]] <- write_at
  out[["core"]] <- core
  out[["geno_dist"]] <- geno
  out[["geno_dist_intervals"]] <- as.data.frame(overview)
  
  return(out)
}

get_training_data <- function(existing_data, data) {
  # Sequester data
  log_file <- data[["log_file"]]
  tmp_data <- data[["tmp_data"]]
  write_at <- data[["write_at"]]
  core_set <- data[["core"]]
  p_wtn <- qread(existing_data[["p_wtn"]]) %>%
    mutate(Type = as.character(Type)) %>%    # Convert 'Type' column to character type
    mutate(Type = ifelse(Type != "Hybrid", "Non_hybrid", Type)) %>%  # Replace non-"Hybrid" values with "Non_hybrid"
    mutate(Type = as.factor(Type)) %>% # Convert 'Type' column to a factor
    rename(env = Env,
           geno = Geno_new,
           type = Type,
           site = Site,
           year = Year,
           blues = BLUES_dt) %>%
    select(env, type, 
           site, latlong, year, geno, 
           connect_climate, connect_geno,
           blues) %>%
    convert(fct(env, site, latlong, year, geno),
            chr(type, starts_with("connect")),
            num(blues))
  
  # Write a section in the log file
  cat("Starting with defining a training data",
      file = log_file,
      sep = "\n",
      append = TRUE)
  
  # Produce output
  core <- p_wtn %>% 
    filter(connect_geno %in% core_set$sel) %>%
    mutate(set = "core")
  rest <- p_wtn %>% 
    filter(connect_geno %!in% core_set$sel) %>%
    mutate(set = "rest")
  
  geno_connection <- p_wtn %>%
    distinct(geno, connect_geno, type)

  core_reformatted <- core %>%
    select(-connect_geno) %>%
    group_by(env, site, latlong, year, connect_climate, geno) %>%
    summarise(blues = mean(blues), .groups = "drop") %>%
    pivot_wider(id_cols = c("env",
                            "site", "latlong", 
                            "year", "connect_climate"),
                names_from = "geno",
                values_from = "blues") %>%
    pivot_longer(!(c("env",
                     "site", "latlong", 
                     "year", "connect_climate")),
                 names_to = "geno",
                 values_to = "blues") %>%
    left_join(geno_connection, by = "geno") %>%
    select(all_of(colnames(p_wtn)))
  
  training_data <- rest %>% bind_rows(core_reformatted)
  #training_data <- core_reformatted
  
  miss_prop <- round(sum(is.na(training_data$blues))/nrow(training_data), 2)
  print(sprintf("Miss prop = %s", miss_prop))
  cat(sprintf("Miss prop = %s", miss_prop),
      file = log_file,
      sep = "\n",
      append = TRUE)
  # Write a section in the log file
  cat("Done with defining a training data",
      file = log_file,
      sep = "\n",
      append = TRUE)
  
  # Generate output
  out <- list()
  out[["log_file"]] <- log_file
  out[["tmp_data"]] <- tmp_data
  out[["write_at"]] <- write_at
  out[["training_data"]] <- training_data
  return(out)
}

predict_missing <- function(existing_data, data, model, debug){
  # Sequester data
  log_file <- data[["log_file"]]
  tmp_data <- data[["tmp_data"]]
  write_at <- data[["write_at"]]
  training_data <- data[["training_data"]]
  
  if (debug) {
    all_geno <- unique(training_data$geno)[1:10]
    
    training_data <- training_data %>%
      filter(geno %in% all_geno)
  }
  
  # Write a section in the log file
  cat("Start with predicting data",
      file = log_file,
      sep = "\n",
      append = TRUE)
  
  # Produce output
  predicted_data <- calculate_pred(paths = existing_data,
                                   pheno_data = training_data,
                                   index = model,
                                   log_at = sprintf("%s/predict_missing.log", tmp_data),
                                   tmp_data = tmp_data)
  
  train_pred <- predicted_data$pred %>% filter(set == "core")

  # Write a section in the log file
  cat(sprintf("Done with predicting data. obs/pred = %s/%s. train acc = %s", 
              train_pred %>% filter(!is.na(blues)) %>% nrow(),
              train_pred %>% filter(is.na(blues)) %>% nrow(),
              train_pred %>% filter(!is.na(blues)) %>% 
                summarize(corr = cor(blues, pred)) %>%
                pull(corr) %>%
                as.numeric()),
      file = log_file,
      sep = "\n",
      append = TRUE)

  # Generate output
  out <- list()
  out[["log_file"]] <- log_file
  out[["tmp_data"]] <- tmp_data
  out[["write_at"]] <- write_at
  out[["predicted_data"]] <- predicted_data

  return(out)
}

get_residuals <- function(existing_data, data){
  # Sequester data
  log_file <- data[["log_file"]]
  tmp_data <- data[["tmp_data"]]
  write_at <- data[["write_at"]]
  predicted_data <- data[["predicted_data"]]$pred %>%
    convert(fct(geno)) %>% 
    filter(is.na(set)) %>%
    select(-set) # modify this part in training set creation where set == "core" is lost
  
  # Write a section in the log file
  cat("Start with deriving residuals",
      file = log_file,
      sep = "\n",
      append = TRUE)
  
  asr_fit <- asreml(fixed = pred ~ 1,
                    random = ~ geno + env,
                    data = predicted_data)
  
  asr_fir_sum <- summary(asr_fit)$varcomp
  
  fit <- asr_fir_sum["geno", "component"]/(asr_fir_sum["geno", "component"]+asr_fir_sum["units!R", "component"])
  
  predicted_data$resid <- asr_fit$residuals
  
  # Write a section in the log file
  cat(sprintf("Done with deriving residuals. h2 = %s", fit),
      file = log_file,
      sep = "\n",
      append = TRUE)
  
  # Generate output
  out <- list()
  out[["log_file"]] <- log_file
  out[["tmp_data"]] <- tmp_data
  out[["write_at"]] <- write_at
  out[["residual_vals"]] <- predicted_data
  return(out)
}

get_env_clusters <- function(existing_data, data, ec_path){
  # Sequester data
  log_file <- data[["log_file"]]
  tmp_data <- data[["tmp_data"]]
  write_at <- data[["write_at"]]
  residual_vals <- data[["residual_vals"]]
  ec_mat <- qread(ec_path)
  
  # Write a section in the log file
  cat("Start with deriving env_clusters",
      file = log_file,
      sep = "\n",
      append = TRUE)
  
  # reformat to exg
  residual_vals_wide <- residual_vals %>%
    distinct(env, connect_climate, geno, year, resid) %>%
    pivot_wider(id_cols = c("env", "connect_climate", "year"), names_from = "geno", values_from = "resid") %>%
    select(-env) %>%
    distinct(connect_climate, .keep_all = T) %>%
    mutate(Site = gsub("^\\S+\\_\\S+\\_(\\S+)\\_\\d+", "\\1", connect_climate),
           Year = gsub("^\\S+\\_\\S+\\_\\S+\\_(\\d+)", "\\1", connect_climate),
           abb_1 = gsub("(\\S{4}).*", "\\1", Site), 
           abb_2 = gsub("\\d{2}(\\d+)", "\\1", Year),
           abb = paste0(abb_1, abb_2)) %>%
    group_by(abb) %>%
    mutate(ids_mod = ifelse(row_number() != 1, paste0(abb, "_", row_number()), abb)) %>%
    ungroup() %>%
    select(-Site, -Year, -abb_1, -abb_2, -abb) %>%
    relocate(ids_mod)
  
  colors_years <- residual_vals_wide %>% distinct(year)
  n_years <- nrow(colors_years)
  coul <- brewer.pal(n = 12, name = "Paired")[1:n_years]
  colors_years$color <- coul
  
  residual_vals_wide <- residual_vals_wide %>%
    left_join(colors_years, by = "year") %>%
    relocate(color, .after = "year")
  
  # derive relationship matrix
  mat <- as.matrix(residual_vals_wide[, 5:ncol(residual_vals_wide)])
  row.names(mat) <- residual_vals_wide$ids_mod
  ERM_resid <- tcrossprod(mat)
  ERM_resid_scaled <- ERM_resid/(sum(diag(ERM_resid))/nrow(ERM_resid))
  
  # perform clustering and produce plots
  tree_data <- mat %>%
    dist(method = "euclidean") %>% 
    hclust() %>% as.dendrogram() %>%
    set("branches_lwd", c(3,2,3)) %>%
    set("branches_lty", c(1,1,3,1,1,2)) %>%
    set("labels_cex", c(.9,1.2))
  label_order = data.frame("ids_mod" = tree_data %>% labels) 
  colors_to_use <- label_order %>% 
    left_join(residual_vals_wide %>% distinct(ids_mod, color), 
              by = "ids_mod") %>%
    pull(color) %>% as.vector()
  labels_colors(tree_data) <- colors_to_use
  
  pdf(sprintf("%s/%s", write_at, "env_geno_resid_phylo_plot.pdf"), width = 9.5, height = 9.5)
  par(mar = c(2,2,2,2), bg=NA)
  (cic_plot <- circlize_dendrogram(tree_data))
  dev.off()
  
  # cluster analysis
  options("ggrepel.max.overlaps" = Inf)
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  
  dist.object <- dist(mat, method="euclidean", diag=T, upper=T)
  dist.matrix <- as.matrix(dist.object)
  mds <- cmdscale(dist.matrix,eig=TRUE, k=2)
  mds.df <- as.data.frame(mds$points)
  labs <- 100*(mds$eig/sum(mds$eig))[1:2]
  
  k.values <- 1:20
  wss <- function(df, k) {
    k_fit <- kmeans(df, k, nstart = 10 )
    out <- list()
    out[["wss"]] <- k_fit$tot.withinss
    out[["n_clust"]] <- k
    out[["min_size"]] <- min(k_fit$size) 
    out[["clusters"]] <- as.factor(k_fit$cluster)
    return(out)
  }
  set.seed(2345)
  wss_values <- lapply(k.values, wss, df = mds.df)
  names(wss_values) <- k.values
  
  n_clust_qual <- sapply(wss_values, function(x) x[["min_size"]])
  n_clust_0 <- n_clust_qual[which(n_clust_qual >= 3)]
  n_clust <- as.numeric(names(n_clust_0[length(n_clust_0)]))
  
  kmean_n_clust <- data.frame("clusters" = k.values,
                              "wss" = as.numeric(sapply(wss_values, function(x) x[["wss"]]))) %>%
    ggplot(aes(x = clusters, y = wss)) +
    geom_line() +
    geom_point() +
    geom_vline(xintercept = n_clust) +
    theme_classic()
  
  ggsave(plot = kmean_n_clust, filename = sprintf("%s/kmean_n_clust.png", write_at), 
         width = 8.4, height = 8.4, units = "cm", dpi = 600)
  
  clust_plot_0 <- create_cluster_plot(mds.df = mds.df, label = "env",
                                      wss_values = wss_values,
                                      write_at = write_at, 
                                      n_clust = 10, 
                                      labs = labs, 
                                      col_vector = col_vector)
  
  clust_plot_1 <- create_cluster_plot(mds.df = mds.df, 
                                      wss_values = wss_values,
                                      label = NULL, 
                                      write_at = write_at, 
                                      n_clust = 10, 
                                      labs = labs, 
                                      col_vector = col_vector)
  
  clust_plot_2 <- create_cluster_plot(mds.df = mds.df, label = "env",
                                      wss_values = wss_values,
                                      write_at = write_at, 
                                      n_clust = 8, 
                                      labs = labs, 
                                      col_vector = col_vector)
  
  clust_plot_3 <- create_cluster_plot(mds.df = mds.df, label = "env",
                                      wss_values = wss_values,
                                      write_at = write_at, 
                                      n_clust = 12, 
                                      labs = labs, 
                                      col_vector = col_vector)
  
  #kmclusters <- wss_values[[n_clust]][["clusters"]]
  #mds.df$groups <- kmclusters
  #mds.df$env <- rownames(mds.df)
  #  
  #clust_plot_0 <- ggscatter(mds.df,
  #                        x = "V1", 
  #                        y = "V2", 
  #                        xlab = sprintf("PC1 = %s", as.character(round(labs[1], 1))),
  #                        ylab = sprintf("PC2 = %s", as.character(round(labs[2], 1))),
  #                        label = "env", 
  #                        color = "groups", 
  #                        palette = col_vector[1:n_clust], 
  #                        size = 1, 
  #                        ellipse = TRUE, 
  #                        ellipse.type = "convex", 
  #                        repel = TRUE,
  #                        ggtheme = theme_gray(),
  #                        title = sprintf("Clusters = %s", n_clust))
  #
  #clust_plot_1 <- ggscatter(mds.df,
  #                          x = "V1", 
  #                          y = "V2", 
  #                          xlab = sprintf("PC1 = %s", as.character(round(labs[1], 1))),
  #                          ylab = sprintf("PC2 = %s", as.character(round(labs[2], 1))),
  #                          label = NULL, 
  #                          color = "groups", 
  #                          palette = col_vector[1:n_clust], 
  #                          size = 1, 
  #                          ellipse = TRUE, 
  #                          ellipse.type = "convex", 
  #                          repel = TRUE,
  #                          title = sprintf("Clusters = %s", n_clust))
  ##clust_plot <- clust_plot_0 + theme(legend.position = "none")
  #
  #ggsave(plot = clust_plot_0, filename = sprintf("%s/clust_plot_env.png", write_at), 
  #       width = 16.8, height = 16.8, units = "cm", dpi = 600)
  #ggsave(plot = clust_plot_1, filename = sprintf("%s/clust_plot.png", write_at), 
  #       width = 16.8, height = 16.8, units = "cm", dpi = 600)
  
  # Write a section in the log file
  cat("Done with deriving env_clusters",
      file = log_file,
      sep = "\n",
      append = TRUE)
  
  # add groups to original mat
  mds.df_0 <- as.data.frame(clust_plot_1$data) %>% 
    tibble::rownames_to_column(var = "ids_mod") %>%
    left_join(residual_vals_wide %>% select(ids_mod, connect_climate, year, color), by = "ids_mod") %>%
    select(ids_mod, connect_climate, year, color, groups) %>%
    mutate(lat = as.numeric(gsub("(\\S+)_\\S+_\\S+_\\S+", "\\1", connect_climate)),
           long = as.numeric(gsub("\\S+_(\\S+)_\\S+_\\S+", "\\1", connect_climate)),
           loc = gsub("\\S+_\\S+_(\\S+)_\\S+", "\\1", connect_climate),
           year = as.numeric(gsub("\\S+_\\S+_\\S+_(\\S+)", "\\1", connect_climate)))
  
  # visualize on a map
  #bw_map <- get_googlemap(center = c(9.9019, 49.8430), zoom = 4, scale = 2,
  #                        style = "feature:road|visibility:off&style=element:labels|visibility:off&style=feature#:administrative|visibility:on")
  
  #load("/proj/ext_dir/KIBREED/source_data/map.Rdata")
  
  rectangle <- c(left = 5, bottom = 45, right = 20, top = 55)
  
  #bw_map <- get_stadiamap(rectangle, zoom = 6)
  #qsave(bw_map, "/proj/results/R/feature_importance/stadia_map.qs")
  bw_map <- qread("/proj/results/R/feature_importance/stadia_map.qs")
  
  cluster_map <- ggmap(bw_map) +
    geom_point(aes(x = long, y = lat, color = as.factor(year)),
               data = mds.df_0, size = 1) +
    #scale_color_discrete(type = mds.df_0$color)+
    labs(x = "Latitude", y = "Longitude") +
    facet_wrap(~groups, nrow = 2, ncol = 5) +
    guides(color = guide_legend(title = "Year")) +
    scale_colour_manual(values=coul)
  
  ggsave(plot = cluster_map,
         filename = sprintf("%s/clust_map.png", write_at), 
         width = 16.4, height = 8.4, dpi = 600, units = "cm") 
  
  final_data <- ec_mat %>%
    left_join(mds.df_0, by = c("harvest_env" = "connect_climate")) %>%
    filter(!is.na(groups)) %>%
    relocate(all_of(c("ids_mod", "year", "color", "groups")), .after = "env")
  
  write.csv(final_data, file = sprintf("%s/clust_data.csv", write_at), row.names = F)
  
  # Generate output
  out <- list()
  out[["log_file"]] <- log_file
  out[["tmp_data"]] <- tmp_data
  out[["write_at"]] <- write_at
  out[["input_mat"]] <- mat
  out[["dend_plot"]] <- cic_plot
  out[["kmean_n_clust"]] <- kmean_n_clust
  out[["kmean_raw"]] <- wss_values
  out[["clust_plot"]] <- clust_plot_0
  out[["final_data"]] <- final_data
  
  return(out)
}

get_feature_imp_plots <- function(existing_data, data, scores_at){
  # Sequester data
  log_file <- data[["log_file"]]
  tmp_data <- data[["tmp_data"]]
  write_at <- data[["write_at"]]
  feature_imp_0 <- read.csv(scores_at, header = T)
  
  # Write a section in the log file
  cat("Start with feature importance scores",
      file = log_file,
      sep = "\n",
      append = TRUE)
  
  # Produce output
  period_to_month <- data.frame(period = c("0-30", "30-60", "60-90", 
                                           "90-120", "120-150", "150-180", 
                                           "180-210", "210-240", "240-270", 
                                           "270-300", "300-330"),
                                month = c(month.name[10:12], month.name[1:8]))
  feature_imp <- feature_imp_0 %>%
    mutate(variable = gsub("(\\S+)\\_\\(\\S+", "\\1", old, perl = T),
           period = gsub("(\\S+)\\_\\((\\S+)\\]", "\\2", old, perl = T),
           period = gsub(",", "-", period)) %>%
    left_join(period_to_month, by = "period") %>%
    #filter(period != "BLUEs_raw") %>%
    convert(fct(variable)) %>%
    arrange(desc(gain_scaled))
  feature_imp$period <- factor(feature_imp$period, levels = period_to_month$period)
  feature_imp$month <- factor(feature_imp$month, levels = period_to_month$month)
  
  ## mean feature importances
  feature_imp_mean <- feature_imp %>%
    group_by(variable) %>%
    summarize(mean_gain_scaled = mean(gain_scaled), .groups = "drop") %>%
    arrange(desc(mean_gain_scaled))
  feature_imp_mean$variable <- factor(feature_imp_mean$variable, levels = feature_imp_mean$variable)
  
  mean_imp_plot <- ggplot(data=feature_imp_mean, aes(x=variable, y=mean_gain_scaled)) +
    geom_bar(stat="identity") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  ggsave(mean_imp_plot, filename = sprintf("%s/%s", write_at, "mean_imp_plot.png"),
         width = 16.8, height = 16.8, units = "cm", dpi = 600, bg = "white")
  
  ## feature importance per growth stage
  feature_imp$variable <- factor(feature_imp$variable, levels = feature_imp_mean$variable)
  imp_stage_plot <- ggplot(feature_imp, aes(x = variable, y = month, fill = gain_scaled)) +
    geom_tile() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    scale_fill_binned(limits=c(0, 0.04), breaks=seq(0,0.04,by=0.005), type = "viridis")
  
  ggsave(imp_stage_plot, filename = sprintf("%s/%s", write_at, "imp_stage_plot.png"),
         width = 16.8, height = 16.8, units = "cm", dpi = 600, bg = "white")
  
  # Write a section in the log file
  cat("Done with feature importance scores",
      file = log_file,
      sep = "\n",
      append = TRUE)
  
  # Generate output
  out <- list()
  out[["log_file"]] <- log_file
  out[["tmp_data"]] <- tmp_data
  out[["write_at"]] <- write_at
  out[["mean_imp_plot"]] <- mean_imp_plot
  out[["imp_stage_plot"]] <- imp_stage_plot
}
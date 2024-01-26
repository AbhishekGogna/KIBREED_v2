#!/usr/bin/env Rscript

source("/proj/renv/activate.R")

# R version 4.0.1
args = commandArgs(trailingOnly=TRUE)

# load arguments
cv_type = args[[1]]
cv_id = args[[2]] # index of the run info file
model_name = args[[3]]
model_spec = args[[4]]
model_spec = args[[4]]
in_cc = as.logical(args[[5]])
options(scipen = 999)

# functions
ifrm <- function(obj, env = globalenv()) {
  obj <- deparse(substitute(obj))
  if(exists(obj, envir = env)) {
    rm(list = obj, envir = env)
  }
}

prep_as_df_for_feather <- function(data){
  data_df <- cbind("idx" = rownames(data), data)
  data_df <- as.data.frame(data_df)
  return(data_df)
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
  formula <- sprintf("~ -1 + %s", column_name)
  t0 <- Sys.time()
  matrix <- model.matrix(as.formula(formula), data)
  colnames(matrix) <- gsub(column_name, '', colnames(matrix), perl = TRUE)
  t1 <- Sys.time()
  period <- round(as.numeric(t1 - t0, units = "mins"), 2)
  ## generate output
  out <- list()
  out[["mat"]] <- matrix
  out[["time_min"]] <- period
  return(out)
}

create_and_get_product <- function(matrix_params, pheno_data) {
  col_name <- matrix_params$col_name
  kin_mat_path <- matrix_params$kin_mat_path
  incidence_mat <- get_incidence_matrix(data = pheno_data,
                                        column_name = col_name)
  kin_mat <- qread(kin_mat_path)
  mat_out <- get_product(incidence_mat[["mat"]], kin_mat, NULL, NULL)
  return(mat_out)
}

create_and_get_had_prod <- function(matrix_params, pheno_data) {
  if("type" %in% names(matrix_params)) matrix_params <- matrix_params[which(names(matrix_params) != "type")]
      
  matrix_params_1 <- sapply(matrix_params, function(x) x[[1]], USE.NAMES = TRUE, simplify = FALSE)
  matrix_params_2 <- sapply(matrix_params, function(x) x[[2]], USE.NAMES = TRUE, simplify = FALSE)
  matrix1 <- create_and_get_product(matrix_params_1, pheno_data)
  matrix2 <- create_and_get_product(matrix_params_2, pheno_data)
  return(get_had_prod(matrix1[["mat"]], matrix2[["mat"]], NULL))
}

get_BRR_pred <- function(pheno_data, matrix_params, tmp_data_dir, dump_dir, 
                         cv_id, nIter, burnIn){
  
  # create tmp dirs
  dump_at_temp <- paste0(dump_dir, "/par_pred")
  if(!dir.exists(dump_at_temp)) dir.create(dump_at_temp, recursive = T)
  
  # create tmp_dir
  tmp_data <- paste0(tmp_data_dir, "/par_pred")
  if(!dir.exists(tmp_data)) dir.create(tmp_data)
  
  pheno_data_0 <- pheno_data %>%
    mutate(idx = row_number())
  
  col_name <- "connect_param"
  g_data <- qread(matrix_params$kin_mat_path)
  param_data <- qread(matrix_params$BRR_mat_path)
  param_data_df <- as.data.frame(param_data, row.names = NA)
  param_data_df[, col_name] <- rownames(param_data)
  
  pheno_data_mod <- pheno_data_0 %>% mutate(obs_set = ifelse(is.na(blues), "test", "train")) %>%
    arrange(desc(obs_set)) %>% left_join(param_data_df, by = col_name)
  
  envs <- as.character(unique(pheno_data_mod$env))
  meta_info <- NULL
  param_data_out <- NULL
  
  for (env in envs){
    ifrm(pheno_data_0_site, env = environment())
    pheno_data_0_site <- pheno_data_0 %>% filter(env == !!as.character(env))
    param_data_site <- pheno_data_mod %>% filter(env == !!as.character(env)) %>%
      mutate(across(all_of(colnames(param_data)), ~ifelse(!is.na(blues), ., NA), .names = "to_pred_{.col}"))
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

return_elements <- function(element, pheno_data, paths, tmp_data_dir, dump_dir, cv_id) {
  mat_name <- str_split(element, "@")[[1]][1]
  model_type <- str_split(element, "@")[[1]][2]
  
  matrix_mapping <- list(
    E_i        = list(col_name = "connect_climate", type = "inc"),             # Incidence_mat
    G_i        = list(col_name = "connect_geno", type = "inc"),                # Incidence_mat
    S_i        = list(col_name = "site", type = "inc"),                        # Incidence_mat
    Y_i        = list(col_name = "year", type = "inc"),                        # Incidence_mat
    S          = list(col_name = "site",
                      kin_mat_path = paths[["SRM"]]),                          # BRR
    S_al        = list(col_name = "latlong",
                      BRR_mat_path = paths[["G_a_S"]]),                        # BRR
    Y          = list(col_name = "year", 
                      kin_mat_path = paths[["YRM"]]),                          # BRR
    G_a        = list(col_name = "connect_geno", 
                      kin_mat_path = paths[["G_a_RM"]]),                       # BRR
    G_d        = list(col_name = "connect_geno", 
                      kin_mat_path = paths[["G_d_RM"]]),                       # BRR
    G_aa       = list(col_name = "connect_geno", 
                      kin_mat_path = paths[["G_aa_RM"]]),                      # BRR
    ERM_l      = list(col_name = "connect_climate", 
                      kin_mat_path = paths[["ERM_l"]]),                        # BRR
    ERM_nl     = list(col_name = "connect_climate", 
                      kin_mat_path = paths[["ERM_nl"]]),                       # BRR
    G_a_ERM_l  = list(col_name = c("connect_geno", "connect_climate"), 
                      kin_mat_path = c(paths[["G_a_RM"]], paths[["ERM_l"]])),  # RKHS
    G_d_ERM_l  = list(col_name = c("connect_geno", "connect_climate"), 
                      kin_mat_path = c(paths[["G_d_RM"]], paths[["ERM_l"]])),  # RKHS
    G_a_ERM_nl = list(col_name = c("connect_geno", "connect_climate"), 
                      kin_mat_path = c(paths[["G_a_RM"]], paths[["ERM_nl"]])), # RKHS
    G_d_ERM_nl = list(col_name = c("connect_geno", "connect_climate"), 
                      kin_mat_path = c(paths[["G_d_RM"]], paths[["ERM_nl"]])), # RKHS
    G_a_S      = list(col_name = "connect_param", 
                      BRR_mat_path = paths[["G_a_S"]]),                        # BRR
    G_a_S_i    = list(col_name = c("connect_geno", "site"),
                      kin_mat_path = c(paths[["G_a_RM"]], paths[["SRM"]])),    # RKHS
    G_a_S_p    = list(col_name = NA,
                      kin_mat_path = paths[["G_a_RM"]], 
                      BRR_mat_path = paths[["G_a_S"]]),                        # BRR
    G_a_Y      = list(col_name = c("connect_geno", "year"), 
                      kin_mat_path = c(paths[["G_a_RM"]], paths[["YRM"]]))     # BRR
  )
  
  matrix_params <- matrix_mapping[[mat_name]] # get needed paths
  if (!"type" %in% names(matrix_params)) matrix_params[["type"]] <- model_type # add model_types
  
  # define needed values for conditional implementation
  BRR <- matrix_params$BRR_mat_path
  kin <- matrix_params$kin_mat_path
  
  mat_out_final <- list(model = model_type, saveEffects = TRUE, name = mat_name)
  
  if(matrix_params$type == "inc"){
    mat_out <- get_incidence_matrix(data = pheno_data,column_name = matrix_params$col_name)
    mat_out_final[["X"]] <- mat_out[["mat"]]
  } else if (matrix_params$type == "BRR"){
    if (!is.null(BRR) && is.null(kin)) {
      BRR_mat <- qread(BRR)
      mat_out <- list("mat" = BRR_mat[pheno_data[, matrix_params$col_name], ])
      mat_out_final[["X"]] <- mat_out[["mat"]]
      
    } else if (!is.null(BRR) && !is.null(kin)) {
      mat_out <- get_BRR_pred(pheno_data = pheno_data, matrix_params = matrix_params, 
                              tmp_data_dir = tmp_data_dir, dump_dir = dump_dir, 
                              cv_id = cv_id, nIter = 15000, burnIn = 2000)
      mat_out_final[["X"]] <- mat_out[["mat"]]
    } else if (is.null(BRR) && !is.null(kin)) {
      mat_out <- create_and_get_product(matrix_params, pheno_data)
      mat_out_final[["X"]] <- mat_out[["mat"]]
    }
  } else if (matrix_params$type == "RKHS"){
    mat_out <- create_and_get_had_prod(matrix_params, pheno_data)
    if (dim(mat_out[["mat"]])[1] != dim(mat_out[["mat"]])[2]) {
      print("Model type is set to RKHS but the matrix produced is not square")
      mat_out_final <- NULL
    } else {
      mat_out_final[["K"]] <- mat_out[["mat"]]
    }
  }
  
  log_info(paste("ETA element for column =", element, "defined"))
  
  return(mat_out_final)
}

# to_debug
if (in_cc){
  ext_parse <- "results/R"
  write_at <- sprintf("/proj/results/R/generate_prediction_data/%s/run_data/%s", cv_type, cv_id)
} else {
  ext_parse <- "ext_dir"
  write_at <- "/proj"
}

# generate paths for files
log_at <- sprintf("%s/logs/%s.log", write_at, model_name)
result_at <- sprintf("%s/preds/%s.qs", write_at, model_name)
dump_dir <- sprintf("%s/tmp_data/%s_dump", write_at, model_name)
tmp_data_dir <- sprintf("%s/tmp_data/%s_tmp", write_at, model_name)

for (dir in c(dump_dir, tmp_data_dir)) if(!dir.exists(dir)){dir.create(dir, recursive = T)}
dump_at <- paste0(dump_dir, "/", model_name, "_")

# generate environment
from_cran <- c("dplyr",
               "tidyr",
               "hablar",
               "stringr",
               "BGLR",
               "logger",
               "qs",
               "jsonlite",
               "AGHmatrix",
               "feather") # done installation from rserver

if(sum(unlist(lapply(from_cran, function(x) suppressPackageStartupMessages(require(x, character.only = TRUE))))) == length(from_cran)) {
  cat("All required packages are present and are loaded. Version check was not done.", file = log_at, sep = "\n")
} else {
  cat("Some packages were not loaded or are not installed. Please install and load packages manually", file = log_at, sep = "\n")
}

# register the log file for further logs
log_appender(appender_file(log_at))

# load necessary data
log_info(sprintf("ext parse = %s", ext_parse))

kibreed_data_blues_within <- qread(sprintf("/proj/%s/process_cgm_data/BLUEs_within_env_cgm.qs", ext_parse))
run_data <- read_json(sprintf("/proj/%s/generate_prediction_data/%s/%s.json", ext_parse, cv_type, cv_type))

paths <- list(
  "G_a_RM" = sprintf("/proj/%s/preprocessing_geno_data/kin_a.qs", ext_parse),
  "G_d_RM" = sprintf("/proj/%s/preprocessing_geno_data/kin_d.qs", ext_parse),
  "G_aa_RM" = sprintf("/proj/%s/preprocessing_geno_data/kin_aa.qs", ext_parse),
  "ERM_l" = sprintf("/proj/%s/KIBREED_data_generation/ERM_data_linear.qs", ext_parse),
  "ERM_nl" = sprintf("/proj/%s/KIBREED_data_generation/ERM_data_non_linear.qs", ext_parse),
  "SRM" = sprintf("/proj/%s/KIBREED_data_generation/SRM.qs", ext_parse),
  "YRM" = sprintf("/proj/%s/KIBREED_data_generation/YRM.qs", ext_parse),
  "G_a_S" = sprintf("/proj/%s/process_cgm_data/g_s_mat.qs", ext_parse)
)
log_info("All data loaded")

# generate run data
specific_run <- run_data[[cv_id]]
run_name <- paste0(cv_type, "_", model_name, "_", cv_id)
test_index <- do.call(c, specific_run$test)
train_index <- do.call(c, c(specific_run$train, specific_run$val))
all_data <- sort(c(test_index, train_index))
log_info("Run data defined")

## set missing  
phenoGE_with_NA_subset <- kibreed_data_blues_within %>%
  mutate(obs = BLUES_dt) %>%               # Create a new column 'obs' with 'BLUES_dt' values
  mutate(BLUES_dt = ifelse(idx_cv %in% test_index, NA, BLUES_dt)) %>%  # Replace 'BLUES_dt' values with NA where 'test_index' is TRUE
  filter(idx_cv %in% all_data)  %>%   # Subset rows based on 'all_data' condition
  mutate(Type = as.character(Type)) %>%    # Convert 'Type' column to character type
  mutate(Type = ifelse(Type != "Hybrid", "Non_hybrid", Type)) %>%  # Replace non-"Hybrid" values with "Non_hybrid"
  mutate(Type = as.factor(Type), # Convert 'Type' column to a factor
         run_name = run_name) %>%
  rename(env = Env,
         geno = Geno_new,
         type = Type,
         site = Site,
         year = Year,
         blues = BLUES_dt) %>%
  select(run_name, env, type, 
         site, latlong, year, geno, 
         connect_climate, connect_param, connect_geno,
         blues, obs) %>%
  convert(fct(env, site, latlong, year, geno),
          chr(type, run_name, starts_with("connect")),
          num(blues, obs))
## for debugging
##phenoGE_with_NA_subset <- phenoGE_with_NA_subset[sample(1:nrow(phenoGE_with_NA_subset), 1000), ]

#execute = TRUE # switch for quick debugging
#if (execute){
t0 = Sys.time()

if(!file.exists(result_at)){
  # define matrices
  elements <- str_split(model_spec, "&")[[1]]
  
  mat_names <- paste0(phenoGE_with_NA_subset$env, ":", phenoGE_with_NA_subset$geno)
  matrices <- unlist(lapply(str_split(elements, "@"), function(x) x[[1]]))
  
  ETA <- lapply(elements, return_elements,
                pheno_data = phenoGE_with_NA_subset, 
                paths = paths,
                tmp_data_dir = tmp_data_dir,
                dump_dir = dump_dir,
                cv_id = cv_id)
  names(ETA) <- sapply(ETA, function(x) x[["name"]])
  
  t1 = Sys.time()
  log_info(sprintf("ETAs defined for model %s i.e. %s . The whole process starting with defining of E_i Took %s min",
                   model_name, 
                   model_spec,
                   round(as.numeric(t1 - t0, units = "mins"), 2)))
  
  # check if any element has missing values 
  check_na <- sapply(ETA, function(x) ifelse(sum(is.na(x[[grep("X|K", names(x), value = T)]])) == 0, TRUE, FALSE))
  
  if(any(check_na)){
    # run model
    #RhpcBLASctl::blas_set_num_threads(1)
    fm_base <- BGLR(y = phenoGE_with_NA_subset$blues,
                    ETA = ETA,
                    nIter = 15000,
                    burnIn = 2000,
                    thin = 5,
                    saveAt = dump_at,
                    verbose = FALSE)
    t2 <- Sys.time()
    log_info(sprintf("%s took %s minutes.\nPredictions finished.",
                     model_name, 
                     round(as.numeric(t2 - t1, units = "mins"), 2)))
    
    # extract outputs
    
    ## Yield
    phenoGE_with_NA_subset$pred <- as.numeric(fm_base$yHat)
    
    ## vars and sd
    out_var_sd <- NULL
    out_var_sd_names <- NULL
    for (i in matrices){
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
    out_var_sd_names <- c(out_var_sd_names, "Var.E")
    
    combined_vars_sd <- as.data.frame(rbind(out_var_sd), row.names = run_name)
    colnames(combined_vars_sd) <- out_var_sd_names
    
    # save results
    output <- list()
    output[["vars"]] <- combined_vars_sd
    output[["preds"]] <- phenoGE_with_NA_subset
    
    qsave(output, result_at)
    log_info("Results written")
  } else {
    log_info("Predictors have missing values, execution skipped")
    log_info(paste(names(check_na), check_na, sep = ": ", collapse = ", "))
  }
} else {
  log_info("Result file was present, execution skipped")
}

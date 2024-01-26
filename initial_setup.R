# make a project in rstudio with name eg KIBREED_gen_3
project_dir <- "/proj" # if you have a different path to the project directory. Change the string

# initialize dependency management with renv
renv::init(project = project_dir, bare = TRUE,
           settings = list(use.cache = FALSE,
                           external.libraries = c("/proj/rst/R/x86_64-pc-linux-gnu-library/4.0", 
                                                  "/usr/local/lib/R/library")),
           force = TRUE)

# Modify rprofile 
cat('if(getwd() == "/proj") source("renv/activate.R")\n', 
    file = sprintf("%s/%s", project_dir, ".Rprofile"))

# Packages ----------------------------------------------------------------
packs <- c("remotes", # for version specific install
           "tidyverse:1.3.2", "hablar", "readxl", "reshape2", # for data wrangling
           "targets", "tarchetypes", "visNetwork", "usethis", # for reproducible analysis
           "rmarkdown", "knitr", "pbkrtest:0.5.1", "ggpubr:0.5.0", # for reporting
           "foreach", "doParallel", "RhpcBLASctl", # for parallel processes
           "qs", "jsonlite", "feather", # for saving and writing files
           "BGLR", "corehunter:3.2.2") # for statistical analysis

## remotes is the key package
if(!sapply("remotes", require, character.only = TRUE)){
  install.packages("remotes")
}

## extract packs info and install them
pack_list <- strsplit(packs, ":")
pack_names <- do.call(c, lapply(pack_list, function(x) x[[1]]))
pack_versions <- unlist(lapply(pack_list, function(x) if(!is.na(x[2])) x[2] else "NULL"))
names(pack_versions) <- pack_names

success <- suppressWarnings(sapply(pack_names, require, character.only = TRUE))
for(pack in names(success[which(success == FALSE)])){
  vrn = pack_versions[pack]
  if (vrn == "NULL") {vrn <- NULL}
  remotes::install_version(pack, version = vrn,  Ncpus = 30) # modify available cpu's if you have low computing power
}

# to install asreml you need to put the tar file ar /proj/renv/cellar formatted as asreml_version.tar.gz
cellar_at <- sprintf("%s/renv/cellar", project_dir) # make a cellar for user defined packages
asreml_tar_name <-grep("asreml", list.files("/proj/renv/cellar/"), value = TRUE)[1]
if(identical(asreml_tar_name, character(0))){
  print("Asreml tar file not found. Please check.")
} else {
  install.packages(sprintf("%s/%s", cellar_at, asreml_tar_name), repos = NULL, type = "source")
}
#multtest # use source code from bioconductor

# OpenBLAS ----------------------------------------------------------------
# You need to install OpenBLAS on a host from within the container. 
# To do this start shell in the rserver and execute the following
#git clone -b v0.3.23 https://github.com/xianyi/OpenBLAS
#cd OpenBLAS
#ln -s /usr/lib/x86_64-linux-gnu/libmpfr.so.6 /qg-10/data/AGR-QG/Gogna/computing_containers/lib_symlinks/libmpfr.so.4
#export LD_LIBRARY_PATH="/qg-10/data/AGR-QG/Gogna/computing_containers/lib_symlinks:$LD_LIBRARY_PATH"
# make DYNAMIC_ARCH=1
# make install PREFIX="inst/$(hostname)"
# now use 
#sessionInfo()
# to get the default location of the BLAS library in R. Then go outside and bind 
# the absolute path of installed blas library to the path given by sessionInfo.
# eg "${OpenBLAS_lib}:/usr/local/lib/R/lib/libRblas.so" where 
# OpenBLAS_lib="/qg-10/data/AGR-QG/Gogna/computing_containers/OpenBLAS/inst/qg-10.ipk-gatersleben.de/lib/libopenblas.so"
# you can benchmark OpenBLAS with https://mac.r-project.org/benchmarks/R-benchmark-25.R # With OpenBLAS was completed in 8 secs on qg-10
# Set target management ----------------------------------------------------

targets::use_targets(open = FALSE) # will initialize targets package

# Define core paths
my_dirs <- list()

core_dirs <- c("src", "run", "results", "tmp_data", "logs")
raw_data <-  sprintf("%s/raw_data", project_dir)
store_path <- sprintf("%s/store", project_dir)

if(!dir.exists(raw_data)){dir.create(raw_data, recursive = T)}
if(!dir.exists(store_path)){dir.create(store_path, recursive = T)}

for (i in core_dirs){
  for_R <- sprintf("%s/%s/%s", project_dir, i, "R")
  for_Py <- sprintf("%s/%s/%s", project_dir, i, "Py")
  my_dirs[[paste0(i, "_R")]] <- for_R
  my_dirs[[paste0(i, "_Py")]] <- for_Py
  if(!dir.exists(for_R)){dir.create(for_R, recursive = T)}
  if(!dir.exists(for_Py)){dir.create(for_Py, recursive = T)}
}

# store dir data
jsonlite::write_json(my_dirs, sprintf("%s/results/%s", project_dir, "core_paths.json"))

# For R: Modify targets.yaml -----------------------------------------------------
scripts <- c("phenotypic_data_kws_processing",
             "preprocessing_phenodata",
             "preprocessing_environ_data",
             "preprocessing_geno_data",
             "KIBREED_data_generation",
             "process_cgm_data",
             "generate_prediction_data",
             "process_R_pred_data",
             "get_vars",
             "feature_importance"
             ) # define script names here
# add common lines to scripts

common_chunk_script <- 'library(targets)

project_path <- "/proj"
core_paths <- lapply(jsonlite::read_json(sprintf("%s/results/core_paths.json", project_path)), function(x) x[[1]])

# Define packages needed and default storage format for intermittent data

tar_option_set(packages = c("dplyr", "readr", "ggplot2", 
                            "tidyr", "reshape2", "stringr",
                            "lubridate", "hablar", "readxl",
                            "qs", "feather"),
               format = "qs")

run_name <- "%s"

tar_source(sprintf("/%s/src/R/fun_%s.R", project_path, run_name))
'

for (i in scripts){
  script_path <- sprintf("%s/%s.R", my_dirs[["run_R"]], i)
  function_path <- sprintf("%s/fun_%s.R", my_dirs[["src_R"]], i)
  if(!file.exists(script_path)){
    file.create(script_path)
    cat(sprintf(common_chunk_script, "%s", i, "%s", "%s"), file = script_path)
  }
  if(!file.exists(function_path)){
    file.create(function_path)
  }
  if(!file.exists(function_path)){
    file.create(function_path)
  }
  
  targets::tar_config_set(script = script_path,
                          store = sprintf("%s/%s", store_path, i), 
                          project = sprintf("%s", i),
                          inherits = NULL,
                          reporter_make = NULL)
}


# For python -------------------------------------------------------------------

cat("##for DL
tensorflow==2.8.0
tensorboard==2.8.0
pyarrow==5.0.0
matplotlib==3.5.1
pandas==1.4
scikit-learn==1.0.2
patsy==0.5.2
protobuf==3.20.*
keras-tuner==1.1.3
ipykernel==6.22.0
## for ML
xgboost==2.0.0
##for doit
graphviz==0.20
doit==0.36.0
pygraphviz==1.9
import_deps==0.2.0
", file = sprintf("%s/requirements.txt", project_dir))
# Run these from the terminal at /proj
#python3 -m venv /proj/py_env
#source /proj/py_env/bin/activate
#pip3 --no-cache-dir install -r /proj/requirements.txt
#echo "source /proj/py_env/bin/activate" > /proj/.bash_profile
#rm /proj/requirements.txt
# you can then create an enironment with 
#python -m venv /path/to/new/virtual/environment
# This creates a directory at /path/to/new/virtual/environment. 
# You can add source /path/to/new/virtual/environment/bin/activate to .bash_profile to enable it for bash
# To add this environment to jupyter run 
#python -m ipykernel install --user --name=py_env
# Then you can choose this environment from dropdown at top right of the ipynb notebook

# Modify nature of dependency discovery ----------------------------------------------------
cat(list.files(project_dir, recursive = FALSE, all.files = TRUE, no.. = TRUE), sep = "\n",
    file = sprintf("%s/%s", project_dir, ".gitignore")) # remove this file, src and run since these contain scripts you write.

# Write lock file ---------------------------------------------------------
## check status
renv::status()

## update lockfile
renv::snapshot(project = project_dir, force = TRUE)

# Produce R library list --------------------------------------------------
scripts <- c("phenotypic_data_kws_processing",
             "preprocessing_phenodata",
             "preprocessing_environ_data",
             "preprocessing_geno_data",
             "KIBREED_data_generation",
             "process_cgm_data",
             "generate_prediction_data",
             "process_R_pred_data"
) # define script names here

extract_package_names <- function(script_text) {
  # Extract package names within double quotes using regular expression
  package_names <- unique(regmatches(script_text, gregexpr("\"[^\"]+\"", script_text))[[1]])
  
  # Remove double quotes from package names
  package_names <- unique(gsub('"', '', package_names))
  
  all_packs <- as.data.frame(installed.packages())
  
  # Create a dataframe with the package names
  df <- all_packs[which(all_packs$Package %in% package_names), c("Package", "Version", "LibPath")]
  
  return(df)
}

used_packs <- lapply(scripts, function(x) {
  script_text <- readLines(paste0("/proj/run/R/", x, ".R"))
  extracted_packages <- extract_package_names(paste(script_text, collapse = " "))
  return(extracted_packages)
})
names(used_packs) <- scripts
used_packs_df <- do.call(rbind, used_packs)
used_packs_df$script <- rownames(used_packs_df)
used_packs_df$script <- gsub("(\\S+)\\.\\S+", "\\1", used_packs_df$script, perl = T)
rownames(used_packs_df) <- NULL
used_packs_df %>% distinct(Package, .keep_all = T) %>% select(-LibPath, -script)

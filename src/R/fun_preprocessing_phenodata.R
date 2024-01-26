preprocess_phenodata <- function(paths, log_at){
  # put a log file
  run_instance <- as.character(format(Sys.time(),  format = "%d_%m_%Y_%H_%M"))
  log_at <- sprintf("%s/%s", log_at, run_instance)
  
  if(!dir.exists(log_at)){dir.create(log_at, recursive = T)}
  
  # Sci. Adv. data
  
  data_yu <- read.table(paths[["data_yu"]]) %>%
    mutate(Series = ifelse(Exp == "HBD", "Exp_1",
                           ifelse(Exp == "Fac1_S1", "Exp_2",
                                  ifelse(Exp == "Fac1_S2", "Exp_3",
                                         ifelse(grepl("Top", Exp), "Exp_4",
                                                ifelse(grepl("Ap2", Exp), "Exp_5",
                                                       ifelse(Exp == "Inbred", "Exp_6", NA))))))) %>%
    mutate(Env = paste0(Loc, "_", Exp, "_", Year))
  
  cat("Sci. adv data loaded", 
      file = sprintf("%s/preprocess_phenodata.log", log_at), 
      sep = "\n")
  
  # Justifications for column removal
  
  execute_justification <- FALSE
  if(execute_justification){
    data_yu_fil_1 <- data_yu %>% filter(Genotypen_uniform == Geno_new)
    data_yu_fil_1_out <- data_yu %>% filter(Genotypen_uniform != Geno_new)
    
    #nrow(data_yu) == (nrow(data_yu_fil_1_out) + nrow(data_yu_fil_1))
    
    transition_table_1 <- read.table(paths[["transition_table_1"]]) %>%
      select(Name_in_Fac1S2, Name_in_Fac1S2.1) %>%
      mutate(Name_in_Fac1S2 = tolower(Name_in_Fac1S2),
             Name_in_Fac1S2.1 = tolower(Name_in_Fac1S2.1))
    
    data_yu_fil_1_out_MF <- data_yu_fil_1_out %>% filter(paste0(Male, "&&", Female) == Genotypen_uniform) %>% # genotypen uniform formed by joining male then female
      left_join(transition_table_1, by = c("Male" = "Name_in_Fac1S2.1")) %>%
      left_join(transition_table_1, by = c("Female" = "Name_in_Fac1S2.1"), suffix = c(".M", ".F")) %>% # no values for females. Only male names were changed
      mutate(mod = paste0(Name_in_Fac1S2.M, "&&", Female)) %>%
      filter(Geno_new == mod) # after merging the mod is equal to geno new so this 1200 are sorted meaning geno new also has for these male first then female.
    
    data_yu_fil_1_out_FM <- data_yu_fil_1_out %>% filter(paste0(Female, "&&", Male) == Genotypen_uniform) %>% # genotypen uniform formed by joining female then male
      #mutate(mod = paste0(Male, "&&", Female)) %>%
      #filter(Geno_new != mod) %>% # majority is resolved 
      mutate(Female = gsub("_", "", Female)) %>%
      #mutate(mod = paste0(Male, "&&", Female)) %>%
      #filter(Geno_new != mod) %>% # remainig are resolved are correcting female names for special characters
      mutate(Female = ifelse(Female %in% c("turkis", "tã½²kis"), "tuerkis", Female),
             mod = paste0(Male, "&&", Female)) %>%
      filter(Geno_new != mod) # remaining are sorted by correcting for tuerkis
    
    data_yu_fil_1_out_lines <- data_yu_fil_1_out %>% filter(Female == Genotypen_uniform & Male == Genotypen_uniform) # genotypen uniform formed by putting male or female name. Basically correct males/female for transition table and then correct special characters and tuerkis variants. 
    
    nrow(data_yu_fil_1_out) == (nrow(data_yu_fil_1_out_MF) + nrow(data_yu_fil_1_out_FM) + nrow(data_yu_fil_1_out_lines))
    
    # basically this suggests that Male, Female, Genotypen_uniform, and Geno_new are useless!
    
    data_yu_mod <- data_yu %>% 
      select(-Order, -Male, -Female, -Genotypen_uniform,  -Geno_new, -Group_detail, -BLUEs, -ChangeN, -Check, -FemaleNew_uni, -MaleNew_uni, -Geno_uni) %>% # Group_detail is basically giving the same information as Type except for the few cases 
      relocate(Series) %>% mutate(Geno_new = paste0(FemaleNew, MaleNew)) %>%
      select(Series, Project, Exp, Loc, Year, Env, Type, FemaleNew, MaleNew, Geno_new, BLUES_dt)
    
    #data_yu %>% distinct(Group_detail, Type)
    #LineLSA Female
    #LineIPK Female
    #LineIPK   Male
    #LineLSA   Male
  }
  
  # Series wise data processing
  ## Series 1
  series_1 <- data_yu %>% filter(Series == "Exp_1")
  
  old <- series_1 %>% group_by(Genotypen_uniform, Female, Male, Type, Geno_new, FemaleNew, MaleNew, FemaleNew_uni, MaleNew_uni, Geno_uni) %>% 
    count() %>% ungroup() %>% select(Geno_new, Type) %>% 
    mutate(Geno = ifelse(Type == "Hybrid",
                         paste0(gsub("(.*)&&(.*)", "\\2", Geno_new, perl = T), 
                                gsub("(.*)&&(.*)", "\\1", Geno_new, perl = T)),
                         Geno_new))
  
  series_1_hd <- read.table(paths[["series_1_hd"]], header = T) %>% 
    mutate(Loc = gsub("([A-Za-z]?)([0-9]?)", "\\1",  Location,  perl = T),
           Year = as.integer(gsub("([A-Za-z]?)([0-9]?)", "\\2",  Location,  perl = T)),
           Env = paste0(Loc, "_HBD_", Year),
           Geno = tolower(Genotype),
           HD = Blues) %>%
    left_join(old, by = "Geno") %>%
    select(-Location, -Genotype, -status, -TraitName, -std.error, -Blues) %>%
    mutate(Geno = ifelse(is.na(Type), substr(Geno, 3, nchar(Geno)), Geno)) %>%
    select(-Type, -Geno_new) %>%
    left_join(old, by = "Geno") %>%
    select(-Geno) %>%
    relocate(Env, Loc, Year, Geno_new, Type, HD)
  
  series_1_processed <- series_1 %>%
    left_join(series_1_hd, by = c("Env", "Loc", "Year", "Geno_new", "Type")) %>%
    select(-Order, -Male, -Female, -Genotypen_uniform,  -Geno_new, -Group_detail, -BLUEs, -ChangeN, -Check, -FemaleNew_uni, -MaleNew_uni, -Geno_uni) %>% # Group_detail is basically giving the same information as Type except for the few cases 
    relocate(Series) %>% mutate(Geno_new = paste0(FemaleNew, MaleNew)) %>%
    select(Series, Project, Exp, Loc, Year, Env, Type, FemaleNew, MaleNew, Geno_new, BLUES_dt, HD)
  
  rm(series_1, old, series_1_hd)
  cat("Series 1 processing done", 
      file = sprintf("%s/preprocess_phenodata.log", log_at), 
      sep = "\n", 
      append = T)
  ## Series 2
  series_2 <- data_yu %>% filter(Series == "Exp_2") %>% 
    select(-Order, -Male, -Female, -Genotypen_uniform,  -Geno_new, -Group_detail, -BLUEs, -ChangeN, -Check, -FemaleNew_uni, -MaleNew_uni, -Geno_uni) %>% # Group_detail is basically giving the same information as Type except for the few cases 
    relocate(Series) %>% mutate(Geno_new = paste0(FemaleNew, MaleNew)) %>%
    mutate(yield = BLUES_dt) %>%
    select(Series, Project, Exp, Loc, Year, Env, Type, FemaleNew, MaleNew, Geno_new, yield)
  
  pheno_data_HD_Exp_2 <- read.csv(paths[["series_2_hd"]], header = T) %>% select(-HD) %>%
    dplyr::rename(HD_raw = HD_design_corrected) %>%
    mutate(Genotype = gsub(" ", "", tolower(Genotype)),
           Env = gsub(" ", "", paste0(Loc, "_Fac1_S1_", Year)),
           Type = ifelse(grepl("_x_", Genotype), "Hybrid",
                         ifelse(grepl("zm", Genotype), "Male",
                                ifelse(grepl("zf", Genotype), "Female", "Check"))),
           Male = gsub("(.*)\\_x\\_(.*)", "\\1", Genotype, perl = T),
           Female = gsub("(.*)\\_x\\_(.*)", "\\2", Genotype, perl = T),
           Geno_new = ifelse(Type == "Hybrid", paste0(Female, Male), paste0(Genotype, Genotype))) %>%
    select(Env, Loc, Year, Geno_new, Type, HD_raw)
  
  series_2_processed <- series_2 %>% 
    group_by(Series, Project, Exp, Loc, Year, Env, Type, FemaleNew, MaleNew, Geno_new) %>% 
    summarize(BLUES_dt = mean(yield), .groups = "drop") %>%
    left_join(pheno_data_HD_Exp_2 %>% 
                group_by(Env, Loc, Year, Type, Geno_new) %>% 
                summarize(HD = mean(HD_raw), .groups = "drop")%>%
                select(-Env), 
              by = c("Loc", "Year", "Type", "Geno_new"))
  
  rm(series_2, pheno_data_HD_Exp_2)
  cat("Series 2 processing done", 
      file = sprintf("%s/preprocess_phenodata.log", log_at), 
      sep = "\n", 
      append = T)
  
  ## Series 3
  series_3 <- data_yu %>% filter(Series == "Exp_3")
  
  data_yu_yield_Exp_3_new_year_1_and_2 <- read_xlsx(paths[["series_3"]], sheet = 1, col_types = "text") %>% select(-`...1`) # todo: BLUEs_dt of data_yu maps to Yield_design_corrected2. This column also has NA coded as text so correct it.
  
  data_yu_HD_Exp_3_new_year_1_and_2 <- read_xlsx(paths[["series_3"]], sheet = 3, col_types = "text")  %>% select(-`...1`, - HD_BBCH) # todo: HD corrected of data_yu maps to Yield_design_corrected2. This column also has NA coded as text so correct it.
  
  data_yu_Exp_3_new_year_1_and_2 <- data_yu_yield_Exp_3_new_year_1_and_2 %>% 
    left_join(data_yu_HD_Exp_3_new_year_1_and_2, by = all_of(colnames(data_yu_yield_Exp_3_new_year_1_and_2)[c(1:15, 17)]))
  
  old <- series_3 %>% group_by(Genotypen_uniform, Female, Male, Type, Geno_new, FemaleNew, MaleNew, FemaleNew_uni, MaleNew_uni, Geno_uni) %>% 
    count() %>% ungroup() %>% select(Geno_new, Type) 
  
  series_3_processed <- data_yu_Exp_3_new_year_1_and_2 %>% 
    as_tibble() %>%
    mutate(Genotype_new = tolower(gsub("_x_", "&&", Genotype))) %>% 
    left_join(old, by = c("Genotype_new" = "Geno_new")) %>%
    mutate(Order = as.integer(Order),
           Project = "Zucht",
           Exp = "Fac1_S2", #since versuch has levels = "FactorialSerie_2_Jahr_1_LP",  "FactorialSerie_2_Jahr_2_LP" and "Factorial_Serie_2_Jahr_2_LP"
           Loc = Location, # 4 new locations "Bauer"    "Streng"   "Strube"   "Nordsaat"
           Year = as.integer(Year),
           Env = paste0(Loc, "_", Exp, "_", Year),
           Group_detail = ifelse(Type =="Check", "CHK",
                                 ifelse(Type == "Female", "FEM",
                                        ifelse(Type == "Male", "MALE",
                                               ifelse(Type == "Hybrid" , "HYB", NA)))),
           Genotypen_uniform = Genotype_new,
           Female = gsub("(\\S+)\\&\\&(\\S+)", "\\1", Genotype_new, perl = T),
           Male = gsub("(\\S+)\\&\\&(\\S+)", "\\2", Genotype_new, perl = T),
           Type = ifelse(is.na(Type) & grepl("&&", Genotype_new), "Hybrid", Type),
           Yield_design_corrected2 = ifelse(Yield_design_corrected2 == "NA", NA, Yield_design_corrected2),
           Headingdate_design_corrected2 = ifelse(Headingdate_design_corrected2 == "NA", NA, Headingdate_design_corrected2),
           BLUES_dt = as.numeric(Yield_design_corrected2),
           BLUEs = BLUES_dt/10,
           HD = as.numeric(Headingdate_design_corrected2),
           Geno_new = Genotype_new,
           FemaleNew = Female,
           MaleNew = Male,
           ChangeN = as.integer(0),
           Check = TRUE,
           FemaleNew_uni = Female,
           MaleNew_uni = Male,
           Geno_uni = gsub("(\\S+)\\&\\&(\\S+)", "\\1\\2", Genotype_new, perl = T),
           Series = "Exp_3") %>%
    select(all_of(colnames(series_3)), HD) %>%
    select(-Order, -Male, -Female, -Genotypen_uniform,  -Geno_new, -Group_detail, -BLUEs, -ChangeN, -Check, -FemaleNew_uni, -MaleNew_uni, -Geno_uni) %>% # Group_detail is basically giving the same information as Type except for the few cases 
    relocate(Series) %>% mutate(Geno_new = paste0(FemaleNew, MaleNew)) %>%
    select(Series, Project, Exp, Loc, Year, Env, Type, FemaleNew, MaleNew, Geno_new, BLUES_dt, HD) %>%
    group_by(Series, Project, Exp, Loc, Year, Env, Type, FemaleNew, MaleNew, Geno_new) %>%  # multiple counts of checks per environments so i will take the mean value here
    summarise(BLUES_dt = mean(BLUES_dt), HD = mean(HD), .groups = "drop")
  
  rm(series_3, data_yu_yield_Exp_3_new_year_1_and_2, data_yu_HD_Exp_3_new_year_1_and_2, data_yu_Exp_3_new_year_1_and_2, old)
  cat("Series 3 processing done", 
      file = sprintf("%s/preprocess_phenodata.log", log_at), 
      sep = "\n", 
      append = T)
  
  ## Series 4
  series_4 <- data_yu %>% filter(Series == "Exp_4")
  
  old <- series_4 %>% group_by(Genotypen_uniform, Female, Male, Type, Geno_new, FemaleNew, MaleNew, FemaleNew_uni, MaleNew_uni, Geno_uni) %>% 
    count() %>% ungroup() %>% select(Geno_new, Type) # to check genotype overlap
  
  transition_table_1 <- read.table(paths[["transition_table_1"]]) %>%
    mutate(Male_new = Name_in_Fac1S2,
           old_name = Name_in_Fac1S2.1) %>%
    select(old_name, Male_new) # since the name change was only done for series 2
  
  series_4_top_s1_hd <- read_xlsx(paths[["series_4_top_s1_hd"]]) %>% 
    select(Loc, Year, Genotype, Type, trait_corrected) %>% 
    mutate(Project = "Zucht", Exp = "TopS1",  HD = as.numeric(ifelse(trait_corrected != "NA", trait_corrected, NA))) %>% select(-trait_corrected) %>%
    relocate(Project, Exp, Loc, Year, Genotype, Type)
  
  series_4_top_s2_hd <- read_xlsx(paths[["series_4_top_s2_hd"]]) %>%
    select(Location, Year, Genotype, Type, trait_corrected) %>%
    mutate(Project = "Zucht", Exp = "TopS2", HD = as.numeric(ifelse(trait_corrected != "NA", trait_corrected, NA)), Loc = Location) %>% select(-trait_corrected, -Location) %>%
    relocate(Project, Exp, Loc, Year, Genotype, Type)
  
  series_4_hd <- series_4_top_s1_hd %>% bind_rows(series_4_top_s2_hd) %>% 
    mutate(Male = gsub("(\\S+)\\_x\\_(\\S+)", "\\1", Genotype, perl = T),
           Female = gsub("(\\S+)\\_x\\_(\\S+)", "\\2", Genotype, perl = T)) %>%
    left_join(transition_table_1, by = c("Male" = "old_name")) %>%
    mutate(Genotype_new = ifelse(!is.na(Male_new) & Type == "HYB", paste0(Male_new, "_x_", Female), 
                                 ifelse(!is.na(Male_new) & Type == "MALE", Male_new, Genotype)),
           Genotype = Genotype_new, 
           Male = ifelse(is.na(Male_new), Male, Male_new)) %>%
    select(-Genotype_new, -Male_new, -Male, -Female) %>%
    mutate(Genotype_new = tolower(gsub("_x_", "&&", Genotype))) %>%
    select(-Genotype) %>%
    left_join(old, by = c("Genotype_new" = "Geno_new")) %>%
    mutate(Group_detail = Type.x,
           Type = Type.y,
           Geno_new = Genotype_new) %>%
    select(-Type.x, -Type.y, -Genotype_new)
  
  series_4_processed <-  series_4 %>%
    left_join(series_4_hd, by = c("Project", "Exp", "Loc", "Year", "Group_detail", "Type", "Geno_new")) %>%
    select(all_of(colnames(series_4)), HD) %>%
    select(-Order, -Male, -Female, -Genotypen_uniform,  -Geno_new, -Group_detail, -BLUEs, -ChangeN, -Check, -FemaleNew_uni, -MaleNew_uni, -Geno_uni) %>% # Group_detail is basically giving the same information as Type except for the few cases 
    relocate(Series) %>% mutate(Geno_new = paste0(FemaleNew, MaleNew)) %>%
    select(Series, Project, Exp, Loc, Year, Env, Type, FemaleNew, MaleNew, Geno_new, BLUES_dt, HD)
  
  rm(series_4, transition_table_1, series_4_top_s1_hd, series_4_top_s2_hd, series_4_hd, old)
  cat("Series 4 processing done", 
      file = sprintf("%s/preprocess_phenodata.log", log_at), 
      sep = "\n", 
      append = T)
  
  ## Series 5
  unique_columns <- colnames(data_yu)[grep("BLUE|Order", colnames(data_yu), invert = T)]
  
  series_5 <- data_yu %>% filter(Series == "Exp_5") %>% 
    group_by_at(unique_columns) %>%
    summarize(BLUES_dt = mean(BLUES_dt), BLUEs = mean(BLUEs), .groups = "drop") %>%
    mutate(Geno_new = gsub("franz_piko", "piko&&franz", Geno_new),
           MaleNew = gsub("franz_piko", "piko", MaleNew),
           FemaleNew = gsub("franz_piko", "franz", FemaleNew))
  
  old <- series_5 %>% group_by(Exp, Genotypen_uniform, Female, Male, Group_detail, Type, Geno_new, FemaleNew, MaleNew, FemaleNew_uni, MaleNew_uni, Geno_uni) %>% 
    count() %>% ungroup() %>% select(Exp, Geno_new, Group_detail, Type) %>% distinct(Exp, Geno_new, Group_detail, Type)
  
  HD_s1 <- as.matrix(read_xlsx(paths[["series_5_wp1_hd"]], col_types = "text"))
  HD_s1[which(HD_s1 == "NA")] <- NA
  HD_s1 <- as_tibble(HD_s1) %>% mutate(Loc = Ort, 
                                       Year = as.integer(Year), 
                                       Env = paste0(Loc, "_", Year),
                                       Genotype_ML = gsub("Turkis", "Tuerkis", Genotype_ML),
                                       Genotype_ML = gsub("_kurz", "", Genotype_ML),
                                       Genotype_ML = gsub("_lang", "", Genotype_ML),
                                       Genotype = tolower(ifelse(grepl("_x_", Genotype_ML),
                                                                 paste0(gsub("(\\S+)\\_x\\_(\\S+)", "\\2", Genotype_ML, perl = T),
                                                                        "&&",
                                                                        gsub("(\\S+)\\_x\\_(\\S+)", "\\1", Genotype_ML, perl = T)),
                                                                 Genotype_ML)),
                                       Type = Type2,
                                       HD = as.numeric(trait_corrected),
                                       Exp = "Ap2_S1") %>%
    group_by(Exp, Env, Loc, Year, Genotype, Type) %>%
    summarize(HD = mean(HD), .groups = "drop")
  
  HD_s2 <- as.matrix(read_xlsx(paths[["series_5_wp2_hd"]], col_types = "text"))
  HD_s2[which(HD_s2 == "NA")] <- NA
  HD_s2 <- as_tibble(HD_s2) %>% mutate(Loc = Ort, 
                                       Year = as.integer(Jahr), 
                                       Env = paste0(Loc, "_", Year), 
                                       Type = Status_Final,
                                       Genotyp_neu = gsub("IPKF1_Franz_Piko", "IPKF1_Franz_x_Piko", Genotyp_neu),
                                       Genotyp_neu_edit = ifelse(grepl("Hybrid|Check", Type), gsub("^(\\w{1,6}\\_)?(.*)", "\\2", Genotyp_neu, perl = T), Genotyp_neu),
                                       Genotyp_neu_edit = gsub("Türkis", "tuerkis", Genotyp_neu_edit),
                                       p1 = ifelse(grepl("Hybrid|Check", Type, perl = T) & grepl("_x_", Genotyp_neu_edit), 
                                                   gsub("(.*)\\_x\\_(.*)", "\\1", Genotyp_neu_edit, perl = T),
                                                   gsub("(.*)(\\_TRI.*)", "\\1", Genotyp_neu_edit, perl = T)),
                                       p2 = ifelse(grepl("Hybrid|Check", Type, perl = T) & grepl("_x_", Genotyp_neu_edit), 
                                                   gsub("(.*)\\_x\\_(.*)", "\\2", Genotyp_neu_edit, perl = T),
                                                   gsub("(.*\\_)(TRI.*)", "\\2", Genotyp_neu_edit, perl = T)),
                                       p1 = tolower(gsub("[^[:alnum:]]", "", p1, perl = T)),
                                       p2 = tolower(gsub("[^[:alnum:]]", "", p2, perl = T)),
                                       Genotype = ifelse(grepl("Hybrid", Type) | grepl("_x_", Genotyp_neu_edit), paste0(p2, "&&", p1), p1),
                                       Genotype = gsub("mvberes1", "mvberes(1)", Genotype),
                                       HD = as.numeric(trait_corrected),
                                       Exp = "Ap2_S2") %>%
    select(Exp, Env, Loc, Year, Genotype, Type, HD)
  
  HD_s3 <- as.matrix(read_xlsx(paths[["series_5_wp3_hd"]], col_types = "text")) 
  HD_s3[which(HD_s3 == "NA")] <- NA
  HD_s3 <- as_tibble(HD_s3) %>% mutate(Loc = Ort, 
                                       Year = as.integer("2018"), 
                                       Env = paste0(Loc, "_", Year),
                                       Genotype_F1 = ifelse(grepl("F1", Genotype), gsub(".*F1(MT|TC)?(.*)", "\\2", Genotype, perl = T), Genotype), 
                                       Genotype_old = Genotype,
                                       Genotype = tolower(ifelse(grepl("_x_", Genotype_F1),
                                                                 paste0(gsub("(\\S+)\\_x\\_(\\S+)", "\\2", Genotype_F1, perl = T),
                                                                        "&&",
                                                                        gsub("(\\S+)\\_x\\_(\\S+)", "\\1", Genotype_F1, perl = T)),
                                                                 Genotype_F1)),
                                       Type = Status_Final,
                                       HD = as.numeric(trait_corrected),
                                       Exp = "Ap2_S3") %>%
    select(Exp, Env, Loc, Year, Genotype, Type, HD, Genotype_F1, Genotype_old)
  
  HD_comb <- HD_s1 %>% bind_rows(HD_s2) %>% bind_rows(HD_s3) %>% 
    group_by(Exp, Loc, Year, Genotype, Type) %>% summarize(HD = mean(HD), .groups = "drop") %>%
    left_join(old, c("Exp", "Type" = "Group_detail", "Genotype" = "Geno_new")) %>%
    select(-Type.y)
  
  # Sanity check
  length(HD_comb$Genotype %in% old$Geno_new) == nrow(HD_comb) # all the names in combined HD data are in the series_5
  
  series_5_processed <- series_5 %>% 
    left_join(HD_comb, by = c("Exp", "Loc", "Year", "Geno_new" = "Genotype", "Group_detail" = "Type")) %>%
    select(-Male, -Female, -Genotypen_uniform,  -Geno_new, -Group_detail, -BLUEs, -ChangeN, -Check, -FemaleNew_uni, -MaleNew_uni, -Geno_uni) %>% # Group_detail is basically giving the same information as Type except for the few cases 
    relocate(Series) %>% mutate(Geno_new = paste0(FemaleNew, MaleNew)) %>%
    select(Series, Project, Exp, Loc, Year, Env, Type, FemaleNew, MaleNew, Geno_new, BLUES_dt, HD)
  
  rm(series_5, HD_comb, HD_s1, HD_s2, HD_s3, old, unique_columns)
  cat("Series 5 processing done", 
      file = sprintf("%s/preprocess_phenodata.log", log_at), 
      sep = "\n", 
      append = T)
  
  ## Series 6
  series_6 <- data_yu %>% filter(Series == "Exp_6")
  
  old <- series_6 %>% group_by(Exp, Genotypen_uniform, Female, Male, Group_detail, Type, Geno_new, FemaleNew, MaleNew, FemaleNew_uni, MaleNew_uni, Geno_uni) %>% 
    count() %>% ungroup() %>% select(Exp, Geno_new, Group_detail, Type) %>% distinct(Geno_new, Type) # only lines in this one
  
  series_6_HD <- read.table(paths[["series_6_old_hd"]], header = T) %>% 
    mutate(Geno_new = tolower(Entry)) %>%
    left_join(old, by = "Geno_new") %>%
    mutate(Geno_new = ifelse(is.na(Type), paste0("x", Geno_new), Geno_new)) %>%
    left_join(old, by = "Geno_new") %>%
    mutate(HD = BLUE,
           Year = as.integer(paste0("20", gsub("(\\d{2})\\_(\\S+)", "\\1", Environment, perl = T))),
           Loc =  gsub("(\\d{2})\\_(\\S+)", "\\2", Environment, perl = T)) %>%
    select(-Type.x, -Entry, -Environment, -BLUE) %>%
    mutate(Type = Type.y) %>%
    select(-Type.y) %>%
    relocate(Loc, Year, Geno_new, Type, HD)
  
  series_6_pre <- series_6 %>% 
    left_join(series_6_HD, by = c("Loc", "Year", "Geno_new", "Type")) %>%
    select(-Male, -Female, -Genotypen_uniform,  -Geno_new, -Group_detail, -BLUEs, -ChangeN, -Check, -FemaleNew_uni, -MaleNew_uni, -Geno_uni) %>% # Group_detail is basically giving the same information as Type except for the few cases 
    relocate(Series) %>% mutate(Geno_new = paste0(FemaleNew, MaleNew)) %>%
    select(Series, Project, Exp, Loc, Year, Env, Type, FemaleNew, MaleNew, Geno_new, BLUES_dt, HD)
  
  kws_formatted <- NULL
  
  for(i in c("KWS_2020", "KWS_2021")){
    if(i == "KWS_2020"){
      kws_pheno_yield <- read.csv(paths[["series_6_2020_yield"]], header = T) %>% filter(!is.na(value)) %>% mutate(YIE = value) %>% select(-value)
      kws_pheno_hd <- read.csv(paths[["series_6_2020_hd"]], header = T) %>% filter(!is.na(value)) %>% mutate(HD = value) %>% select(-value)
    } else if (i == "KWS_2021"){
      kws_pheno_yield <- read.csv(paths[["series_6_2021_yield"]], header = T) %>% filter(!is.na(value)) %>% mutate(YIE = value) %>% select(-value)
      kws_pheno_hd <- read.csv(paths[["series_6_2021_hd"]], header = T) %>% filter(!is.na(value)) %>% mutate(HD = value) %>% select(-value)
    }
    
    
    kws_pheno <- kws_pheno_yield %>% left_join(kws_pheno_hd, by = c("location", "genotype"))
    #kws_pheno$subtrial <- unlist(lapply(str_split(kws_pheno$data_type, "_"), function(x) x[[2]])) # was there do process output from Blues...acr functions
    
    #check <- lapply(str_split(kws_pheno$data_type, "_"), function(x) {
    #  sub_data <- x[c(-1, -2)]
    #  return(sub_data)})
    
    #kws_pheno$trait <- unlist(lapply(check, function(x) {
    #  input <- x
    #  if (length(input) == 1){
    #    out <- input
    #  } else if (length(input) > 1){
    #    out <- paste0(input, collapse = "_")
    #  }
    #  return(out)
    #}))
    
    year <- gsub("\\S+(\\d{4})", "\\1", i, perl = T)
    
    kws_pheno_formatted <- kws_pheno %>% 
      #select(-data_type) %>%
      #pivot_wider(id_cols = c(location, subtrial, genotype), names_from = trait, values_from = value) %>%
      mutate(Series = "Exp_6",
             Project = "KWS_Bg",
             Exp = "Inbred",
             Year = as.integer(year),
             FemaleNew = tolower(genotype),
             MaleNew = tolower(genotype),
             Geno_new = paste0(FemaleNew, MaleNew),
             Type = "Lines",
             BLUES_dt = as.numeric(YIE),
             HD = as.numeric(HD),
             Code = ifelse(grepl("Großaitingen", location), "GAI",
                           ifelse(grepl("Kondratowice", location), "KON",
                                  ifelse(grepl("Nörvenich", location), "NRV",
                                         ifelse(grepl("Nasenberg", location), "NAS",
                                                ifelse(grepl("Bernburg", location), "BER",
                                                       ifelse(grepl("Seligenstadt", location), "SEL",
                                                              ifelse(grepl("Wetze", location), "WET",
                                                                     ifelse(grepl("Wohlde", location), "WOH",
                                                                            ifelse(grepl("Schmoel", location), "SML", NA))))))))),
             Loc = Code,
             Env = paste0(Loc, "_", Exp, "_", Year)) %>%
      select(Series, Project, Exp, Loc, Year, Env, Type, FemaleNew, MaleNew, Geno_new, BLUES_dt, HD)
    
    kws_formatted <- rbind(kws_formatted, kws_pheno_formatted)
  }
  
  series_6_processed <- series_6_pre %>%
    bind_rows(kws_formatted)
  
  rm(series_6, kws_formatted, kws_pheno, kws_pheno_formatted, series_6_HD, series_6_pre, old, kws_pheno_hd, kws_pheno_yield)
  cat("Series 6 processing done", 
      file = sprintf("%s/preprocess_phenodata.log", log_at), 
      sep = "\n", 
      append = T)
  
  ## Series 7
  load(paths[["series_7"]])
  
  traits <- c("YIE", "HD")
  
  series_7_processed <- NULL
  
  for(trait in traits){
    gabi_blues <- BLUES_traits[[trait]][["BLUEs"]]
    
    cols_within_loc <- grep(trait, colnames(gabi_blues), value = T)
    checks_data_trial_1 <- read.table(paths[["transition_for_checks"]], header = F, col.names = c("Coding", "name", "correct_name"))
    
    recoding <- read.table(paths[["for_kibreed_coding_to_variety"]], 
                           header = F, col.names = c("IDT", "genodata")) 
    
    for (i in 1:nrow(checks_data_trial_1)){
      gabi_blues$Coding[which(gabi_blues$Coding == checks_data_trial_1$Coding[i])] <- checks_data_trial_1$name[i]
    }
    
    if(trait == "YIE"){trait_name <- "BLUES_dt"} else {trait_name = "HD"}
    
    gabi_blues_within_env <- gabi_blues[, c("Coding", cols_within_loc)] %>%
      left_join(recoding, by = c("Coding" = "IDT")) %>%
      pivot_longer(!c("Coding", "genodata"), names_to = "id", values_to = "value") %>%
      mutate(Series = "Exp_7",
             Project = "GABI",
             Exp = "Inbred",
             Loc = gsub("(\\S+)\\.(\\S+)\\.(\\S+)", "\\2", id, perl = T),
             Year = as.integer(gsub("(\\S+)\\.(\\S+)\\.(\\S+)", "\\1", id, perl = T)),
             Env = paste0(Loc, "_", Exp, "_", Year),
             Type = "Lines",
             FemaleNew = tolower(Coding),
             MaleNew = tolower(Coding),
             Geno_new = paste0(tolower(genodata), tolower(genodata)),
             "{trait_name}" := value) %>%
      select(Series, Project, Exp, Loc, Year, Env, Type, FemaleNew, MaleNew, Geno_new, all_of(trait_name))
    
    if (is.null(series_7_processed)){
      series_7_processed <- gabi_blues_within_env
    } else {
      series_7_processed <- series_7_processed %>% left_join(gabi_blues_within_env, by = all_of(colnames(series_7_processed)[1:10]))
    }
  }
  
  rm(BLUES_traits, traits, recoding, checks_data_trial_1, gabi_blues, gabi_blues_within_env)
  cat("Series 7 processing done", 
      file = sprintf("%s/preprocess_phenodata.log", log_at), 
      sep = "\n", 
      append = T)
  
  ## merge data 
  combined_data <- series_1_processed %>%
    bind_rows(series_2_processed) %>%
    bind_rows(series_3_processed) %>%
    bind_rows(series_4_processed) %>%
    bind_rows(series_5_processed) %>%
    bind_rows(series_6_processed) %>%
    bind_rows(series_7_processed) 
  
  # Correct BBCH heading date for KWS 2012 to 2015
  
  series_6 <- combined_data %>% filter(Series == "Exp_6")
  
  HD_BBCH <- series_6 %>% filter(Project == "KWS")
  
  non_HD_BBCH <- combined_data %>% filter(Series != "Exp_6") %>% 
    bind_rows(series_6 %>% filter(Project == "KWS_Bg"))
  
  HD_BBCH_geno <- unique(HD_BBCH$Geno_new)
  
  HD_BBCH_geno_overlap <- HD_BBCH_geno[HD_BBCH_geno %in% non_HD_BBCH$Geno_new] # 40 checks
  
  non_HD_BBCH_check_data <- non_HD_BBCH %>% filter(Geno_new %in% HD_BBCH_geno_overlap) %>% 
    group_by(Geno_new) %>% summarize(HD_mean_jan_1 = mean(HD, na.rm = T), .groups = "drop")
  
  HD_BBCH_check_data <- HD_BBCH %>% filter(Geno_new %in% HD_BBCH_geno_overlap) %>% 
    group_by(Geno_new) %>% summarize(HD_mean_BBCH = mean(HD, na.rm = T), .groups = "drop")
  
  check_data <- non_HD_BBCH_check_data %>% left_join(HD_BBCH_check_data, by = "Geno_new")
  
  l_model <- lm( HD_mean_jan_1 ~ HD_mean_BBCH, data = check_data)
  
  #summary(l_model)
  
  HD_BBCH$HD_pred <- predict(l_model, newdata = HD_BBCH %>% mutate(HD_mean_BBCH = HD))
  
  corrected_data <- non_HD_BBCH %>% bind_rows(HD_BBCH %>% select(-HD) %>% mutate(HD = HD_pred) %>% select(-HD_pred))
  
  rm(list = setdiff(ls(), c("%!in%", "corrected_data", "log_at")))
  cat("BBCH stage corrected for series 6 and output saved.", 
      file = sprintf("%s/preprocess_phenodata.log", log_at), 
      sep = "\n", 
      append = T)
  
  combined_data <- corrected_data
  
  return(combined_data)
  
}

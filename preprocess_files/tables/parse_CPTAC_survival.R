# Yige Wu @ WashU 2019 Jan
## parse CPTAC to a data frame with showing survial info for Cox model


if (cancer_tmp == "CCRCC") {
  # input clinical file -----------------------------------------------------
  clinical_case_file_tmp <- fread(input = "./Ding_Lab/Projects_Current/CPTAC/PGDAC/CPTAC_Biospecimens_Clinical_Data/CPTAC3_Clinical_Data/CCRCC/20190606/CCRCC_June2019_case.tsv", data.table = F)
  clinical_discovery_set <- read_xlsx(path = "./Ding_Lab/Projects_Current/CPTAC/CPTACIII/CCRCC_AWG/CCRCC_shared_data/manuscripts/ccRCC Manuscript/CPTAC3-ccRCC-SupplementaryTables_Final/Table S1.xlsx", sheet = "ccrcc_clinical_characteristics")
  
  # make the data frame for cox model ---------------------------------------
  partIDs_tmp <- unique(clinical_case_file_tmp$Case_ID)
  
  df_cox_tmp <- NULL
  for (partID_tmp in partIDs_tmp) {
    # partID_tmp <- "C3L-01106"
    clinical_case_file_partID_tmp <- clinical_case_file_tmp %>%
      filter(Case_ID == partID_tmp)
    last_contact_day_tmp <- as.numeric(unlist(strsplit(x = clinical_case_file_partID_tmp$`follow-up/Path_Diag_to_Last_Contact_Day`, split = "\\|")))
    last_contact_day_tmp <- unique(last_contact_day_tmp[length(last_contact_day_tmp)])
    
    event_day_tmp <- as.numeric(unlist(strsplit(x = clinical_case_file_partID_tmp$`follow-up/Path_Diag_to_Death_days`, split = "\\|")))
    event_day_tmp <- unique(event_day_tmp[!is.na(event_day_tmp)])
    
    if (length(last_contact_day_tmp) == 0) {
      next()
    } else {
      if (length(event_day_tmp) > 0) { ## status 2: event happened
        time_tmp <- event_day_tmp
        df_tmp <- data.frame(partID = partID_tmp, time = time_tmp, status = 2)
        
      } else { ## status 1: cencorted
        time_tmp <- last_contact_day_tmp
        df_tmp <- data.frame(partID = partID_tmp, time = time_tmp, status = 1)
      }
      df_cox_tmp <- rbind(df_tmp, df_cox_tmp)
    }
    
  }
  
  
  # incorporate age ---------------------------------------------------------
  ages_tmp <- data.frame(partID = clinical_case_file_tmp$Case_ID, age_string = clinical_case_file_tmp$`consent/Age`)
  tmp <- str_split_fixed(string = ages_tmp$age_string, pattern = ">=", 2)
  tmp[tmp[,1] == "", 1] <- tmp[tmp[,1] == "", 2]
  ages_tmp$age <- as.numeric(tmp[,1])
  df_cox_tmp <- merge(df_cox_tmp, ages_tmp[, c("partID", "age")], all.x = T)
  
  # incorporate stage -------------------------------------------------------
  stages_tmp <- data.frame(partID = clinical_case_file_tmp$Case_ID)
  stages_tmp$stage_numeric <- NA
  stages_tmp$stage_numeric[clinical_case_file_tmp$`baseline/Tumor_Stage_Pathological` == "Stage I"] <- 1
  stages_tmp$stage_numeric[clinical_case_file_tmp$`baseline/Tumor_Stage_Pathological` == "Stage II"] <- 2
  stages_tmp$stage_numeric[clinical_case_file_tmp$`baseline/Tumor_Stage_Pathological` == "Stage III"] <- 3
  stages_tmp$stage_numeric[clinical_case_file_tmp$`baseline/Tumor_Stage_Pathological` == "Stage IV"] <- 4
  df_cox_tmp <- merge(df_cox_tmp, stages_tmp, all.x = T)
  
  # incorporate grade -------------------------------------------------------
  grades_tmp <- data.frame(partID = clinical_discovery_set$Case_ID)
  grades_tmp$grade_numeric <- NA
  grades_tmp$grade_numeric[clinical_discovery_set$Grade == "G1"] <- 1
  grades_tmp$grade_numeric[clinical_discovery_set$Grade == "G2"] <- 2
  grades_tmp$grade_numeric[clinical_discovery_set$Grade == "G3"] <- 3
  grades_tmp$grade_numeric[clinical_discovery_set$Grade == "G4"] <- 4
  df_cox_tmp <- merge(df_cox_tmp, grades_tmp, all.x = T)
  
  # incoporate metastasis ---------------------------------------------------
  mets_tmp <- data.frame(partID = clinical_case_file_tmp$Case_ID, metastasis = NA)
  mets_tmp$metastasis[clinical_case_file_tmp$`baseline/Clin_Stage_Dist_Mets_cM` == "cM1" | clinical_case_file_tmp$`baseline/Path_Stage_Dist_Mets_pM` == "pM1"] <- 1
  mets_tmp$metastasis[clinical_case_file_tmp$`baseline/Clin_Stage_Dist_Mets_cM` == "cM0" | clinical_case_file_tmp$`baseline/Path_Stage_Dist_Mets_pM` == "pM0" | clinical_case_file_tmp$`baseline/Path_Stage_Dist_Mets_pM` == "No pathologic evidence of distant metastasis"] <- 0
  df_cox_tmp <- merge(df_cox_tmp, mets_tmp, all.x = T)
  
}


# Yige Wu @ WashU 2019 Jan
## parse CPTAC to a data frame with showing survial info for Cox model


if (cancer_tmp == "CCRCC") {
  # input clinical file -----------------------------------------------------
  clinical_case_file_tmp <- fread(input = "./Ding_Lab/Projects_Current/CPTAC/PGDAC/CPTAC_Biospecimens_Clinical_Data/CPTAC3_Clinical_Data/CCRCC/20190606/CCRCC_June2019_case.tsv", data.table = F)
  
  # make the data frame for cox model ---------------------------------------
  partIDs_tmp <- unique(clinical_case_file_tmp$Case_ID)
  
  df_cox_tmp <- NULL
  for (partID_tmp in partIDs_tmp) {
    # partID_tmp <- "C3L-01106"
    clinical_case_file_partID_tmp <- clinical_case_file_tmp %>%
      filter(Case_ID == partID_tmp)
    last_contact_day_tmp <- as.numeric(unlist(strsplit(x = clinical_case_file_partID_tmp$`follow-up/Path_Diag_to_Last_Contact_Day`, split = "\\|")))
    last_contact_day_tmp <- last_contact_day_tmp[length(last_contact_day_tmp)]
    
    event_day_tmp <- as.numeric(unlist(strsplit(x = clinical_case_file_partID_tmp$`follow-up/Path_Diag_to_new_Tumor_days`, split = "\\|")))
    event_day_tmp <- unique(event_day_tmp[!is.na(event_day_tmp)])
    
    if (length(event_day_tmp) > 0) { ## status 2: event happened
      time_tmp <- event_day_tmp
      df_tmp <- data.frame(partID = partID_tmp, time = time_tmp, status = 2)
      df_cox_tmp <- rbind(df_tmp, df_cox_tmp)
      
    } else if (length(last_contact_day_tmp) == 1){ ## status 1: censorted
      time_tmp <- last_contact_day_tmp
      df_tmp <- data.frame(partID = partID_tmp, time = time_tmp, status = 1)
      df_cox_tmp <- rbind(df_tmp, df_cox_tmp)
    }
  }
  
}


# Yige Wu @ WashU Mar 2018 
## merge summary reports from batch1-2 Summary reports from CDAP

# source -------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/preprocess_files/preprocess_files.R')


# inputs ------------------------------------------------------------------
summaryreport_suffixs <- c("Proteome.tmt10.tsv", "Phosphoproteome.phosphosite.tmt10.tsv")
cdap_dirs <- vector("list")
## input directory/file names for each dataset
cdap_dirs[["CCRCC"]] <- vector("list")
cdap_dirs[["CCRCC"]][["masterDirName"]] <- "6_CPTAC3_Clear_Cell_Renal_Cell_Carcinoma"
for (fs in summaryreport_suffixs) {
  cdap_dirs[["CCRCC"]][[fs]] <- paste0("CPTAC3_Clear_Cell_Renal_Cell_Carcinoma_", fs)
}

# process per cancer type -------------------------------------------------
for (cancer in c("CCRCC")) {
  num_batchs <- 1:2
  masterpath <- paste(inputD, cdap_dirs[[cancer]][["masterDirName"]], "/", sep = "")
  ppoutpath <- paste0(inputD, cancer)
  dir.create(ppoutpath)
  mergeoutpath <- paste0(ppoutpath, "/", cancer, "_batch", head(num_batchs, n = 1), "-", tail(num_batchs, n = 1), "_SummaryReports/")
  dir.create(mergeoutpath)
  
  for (fs in summaryreport_suffixs) {
    exp_table <- fread(input = paste0(masterpath, cancer, "_batch", 1, "_SummaryReports/", cdap_dirs[[cancer]][[fs]]),
                   data.table = F)
    for (num_batch in 2:tail(num_batchs, n = 1)) {
      ## input
      exp_table_tmp <- fread(input = paste0(masterpath, cancer, "_batch", num_batch, "_SummaryReports/", cdap_dirs[[cancer]][[fs]]))
      
      ## merge
      exp_table <- merge(exp_table, exp_table_tmp, by = intersect(colnames(exp_table), colnames(exp_table_tmp)), all = T)
    }
    write.table(exp_table, row.names = F, quote=F, sep = '\t', file=paste(mergeoutpath, cdap_dirs[[cancer]][[fs]],sep=""))
  }
}

# Yige Wu @ WashU 2019 Jun
## calculate the association between kinase-substrate score with days to new tumor


# source ------------------------------------------------------------------
library("survival")
# library("survminer")


# set variables -----------------------------------------------------------
enzyme_type <- "kinase"
reg_nonNA <- 20

# input regression result table -------------------------------------------
# regression <- fread(input = paste0(ppnD, "regression/tables/", "change_regression_nonNA/", 
#                                    "regression", "_", "cptac2p_cptac3", "_", "tumor", "_reg_nonNA", reg_nonNA, ".txt"), data.table = F)
regression <- fread(input = paste0(ppnD, "regression/tables/annotate_regression_with_mut_impact/", 
                                   "regression_cptac2p_cptac3_tumor_reg_nonNA20_mut_impact_cancer_specificity_annotated.txt"),
                    data.table = F)
## already annotated with the supporting evidence for each pair

## divide into kinase and phosphatase related
regression <- regression[regression$enzyme_type == enzyme_type,]
## get the significant pairs
regression <- markSigSiteCan(regression = regression, sig_thres = ifelse(enzyme_type == "kinase" ,fdr_pk, fdr_pp), enzyme_type = enzyme_type)
regression$regulated <- (regression$fdr_sig & regression$coef_sig)

for (cancer_tmp in "CCRCC") {

# parse survival data -----------------------------------------------------
  source("./cptac2p_analysis/preprocess_files/tables/parse_CPTAC_new_tumor_days.R")
  

# input kinase-substrate pair score ---------------------------------------

  esscore_tab_fn <- paste0("./cptac2p/analysis_results/phospho_network/es_score/table/calculate_es_pair_score/kinase/esscore_tab", "_", cancer_tmp, "_", enzyme_type, "_reg_nonNA", reg_nonNA, ".txt")
  esscore_tab <- fread(input = esscore_tab_fn, data.table = F)  
  partIDs_exp <- colnames(esscore_tab); partIDs_exp <- partIDs_exp[!(partIDs_exp %in% c("GENE", "SUB_GENE", "SUB_MOD_RSD", "SELF", "pair", "site_id"))]
  esscore_tab <- esscore_tab[!is.na(esscore_tab$pair) & esscore_tab$pair %in% regression$pair[regression$Cancer == cancer_tmp & regression$regulated == T],]
  esscore_tab <- esscore_tab[!duplicated(esscore_tab$pair),]
  
  cox_stat_summary <- NULL
  for (pair_tmp in esscore_tab$pair) {
    esscore_tab_tmp <- esscore_tab[esscore_tab$pair == pair_tmp,]
    esscore_tab_tmp.m <- melt(esscore_tab_tmp[, c("pair", partIDs_exp)], id.vars = "pair")
    esscore_tab_tmp.m <- esscore_tab_tmp.m %>%
      mutate(esscore = value) %>%
      select(variable, esscore)
    
    # merge es score with survival data table ---------------------------------
    df_cox <- merge(df_cox_tmp, esscore_tab_tmp.m, by.x = c("partID"), by.y = c("variable"))
    df_cox <- df_cox %>%
      filter(!is.na(esscore))
    
    if (nrow(df_cox) > 10) {
      # Run cox model -----------------------------------------------------------
      res.cox <- coxph(Surv(time, status) ~ esscore, data = df_cox)
      res.cox$score
      cox_stat_tmp <- summary(res.cox)
      cox_stat_df <- data.frame(pair = pair_tmp, 
                                esscore_z = cox_stat_tmp$coefficients["esscore","z"], 
                                esscore_coef = cox_stat_tmp$coefficients["esscore","coef"],
                                esscore_HR = cox_stat_tmp$conf.int["esscore", "exp(coef)"],
                                esscore_HR_bottom = cox_stat_tmp$conf.int["esscore", "lower .95"],
                                esscore_HR_upper = cox_stat_tmp$conf.int["esscore", "upper .95"],
                                esscore_logrank_pvalue = cox_stat_tmp$logtest["pvalue"])
      cox_stat_summary <- rbind(cox_stat_df, cox_stat_summary)
    }
  }
  cox_stat_summary <- merge(cox_stat_summary, esscore_tab[, c("pair", "SELF", "GENE")], all.x = T)
  cox_stat_summary$cancer <- cancer_tmp
  cox_stat_summary$esscore_logrank_FDR <- FDR_by_id_columns(p_vector = cox_stat_summary$esscore_logrank_pvalue, id_columns = c("SELF", "cancer"), df = cox_stat_summary)
  
  
  fn <- paste0(makeOutDir(resultD = resultD),
               "cox_new_tumor_days_with_esscore_", cancer_tmp, "_reg_nonNA", reg_nonNA, ".txt")
  write.table(x = cox_stat_summary, file = fn, 
              row.names = F, quote = F, sep = "\t")
}


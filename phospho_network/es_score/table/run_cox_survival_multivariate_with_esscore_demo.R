# Yige Wu @ WashU 2019 Jan
## calculate the association between kinase-substrate score with survival


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

tmp <- regression %>%
  filter(regulated == T) %>%
  filter(GENE == "MET") %>%
  filter(Cancer == "CCRCC")

for (cancer_tmp in "CCRCC") {
  
  # parse survival data -----------------------------------------------------
  source("./cptac2p_analysis/preprocess_files/tables/parse_CPTAC_survival.R")
  
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
      filter(!is.na(esscore)) %>%
      filter(!is.na(age)) %>%
      filter(!is.na(stage_numeric)) %>%
      filter(!is.na(grade_numeric)) %>%
      filter(!is.na(metastasis))
    
    if (nrow(df_cox) > 10) {
      # Run cox model -----------------------------------------------------------
      res.cox <- coxph(Surv(time, status) ~ esscore + age + stage_numeric + grade_numeric + metastasis, data = df_cox)
      cox_stat_tmp <- summary(res.cox)
      
      cox_stat_df <- data.frame(pair = pair_tmp, 
                                esscore_coef = cox_stat_tmp$coefficients["esscore","coef"],
                                esscore_HR = cox_stat_tmp$conf.int["esscore", "exp(coef)"],
                                esscore_HR_bottom = cox_stat_tmp$conf.int["esscore", "lower .95"],
                                esscore_HR_upper = cox_stat_tmp$conf.int["esscore", "upper .95"],
                                esscore_logrank_pvalue = cox_stat_tmp$coefficients["esscore", "Pr(>|z|)"],
                                age_coef = cox_stat_tmp$coefficients["age","coef"],
                                age_HR = cox_stat_tmp$conf.int["age", "exp(coef)"],
                                age_HR_bottom = cox_stat_tmp$conf.int["age", "lower .95"],
                                age_HR_upper = cox_stat_tmp$conf.int["age", "upper .95"],
                                age_logrank_pvalue = cox_stat_tmp$coefficients["age", "Pr(>|z|)"],
                                stage_numeric_coef = cox_stat_tmp$coefficients["stage_numeric","coef"],
                                stage_numeric_HR = cox_stat_tmp$conf.int["stage_numeric", "exp(coef)"],
                                stage_numeric_HR_bottom = cox_stat_tmp$conf.int["stage_numeric", "lower .95"],
                                stage_numeric_HR_upper = cox_stat_tmp$conf.int["stage_numeric", "upper .95"],
                                stage_numeric_logrank_pvalue = cox_stat_tmp$coefficients["stage_numeric", "Pr(>|z|)"],
                                grade_numeric_coef = cox_stat_tmp$coefficients["grade_numeric","coef"],
                                grade_numeric_HR = cox_stat_tmp$conf.int["grade_numeric", "exp(coef)"],
                                grade_numeric_HR_bottom = cox_stat_tmp$conf.int["grade_numeric", "lower .95"],
                                grade_numeric_HR_upper = cox_stat_tmp$conf.int["grade_numeric", "upper .95"],
                                grade_numeric_logrank_pvalue = cox_stat_tmp$coefficients["grade_numeric", "Pr(>|z|)"],
                                metastasis_coef = cox_stat_tmp$coefficients["metastasis","coef"],
                                metastasis_HR = cox_stat_tmp$conf.int["metastasis", "exp(coef)"],
                                metastasis_HR_bottom = cox_stat_tmp$conf.int["metastasis", "lower .95"],
                                metastasis_HR_upper = cox_stat_tmp$conf.int["metastasis", "upper .95"],
                                metastasis_logrank_pvalue = cox_stat_tmp$coefficients["metastasis", "Pr(>|z|)"],
                                sample_size = nrow(df_cox))
      cox_stat_summary <- rbind(cox_stat_df, cox_stat_summary)
    }
  }
  cox_stat_summary <- merge(cox_stat_summary, esscore_tab[, c("pair", "SELF", "GENE")], all.x = T)
  cox_stat_summary$cancer <- cancer_tmp
  cox_stat_summary$esscore_logrank_FDR <- FDR_by_id_columns(p_vector = cox_stat_summary$esscore_logrank_pvalue, id_columns = c("SELF", "cancer"), df = cox_stat_summary)
  
  fn <- paste0(makeOutDir(resultD = resultD),
               "survival_cox_multivariate_with_esscore_", cancer_tmp, "_reg_nonNA", reg_nonNA, ".txt")
  write.table(x = cox_stat_summary, file = fn, 
              row.names = F, quote = F, sep = "\t")
}


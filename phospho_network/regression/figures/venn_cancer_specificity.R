# Yige Wu @ WashU 2019 Feb
# venn diagram showing whether regulated enzyme-substrate-phosphosite pairs are also the sites previously reported

# source ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_plotting.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
library(eulerr)

# set variables -----------------------------------------------------------
reg_nonNA <- 20
# fdr_thres <- c(0.05, 0.1); names(fdr_thres) <- c("kinase", "phosphatase")
reg_sig <- c(0.05, 0.05); names(reg_sig) <- c("kinase", "phosphatase")
cancersprocess <- c("BRCA", "OV", "CO", "UCEC", "CCRCC")

# inputs regression result ------------------------------------------------------------------
regression_sup <- fread(input = paste0(ppnD, "regression/tables/annotate_regression_with_mut_impact/", 
                                       "regression_cptac2p_cptac3_tumor_reg_nonNA20_mut_impact_cancer_specificity_annotated.txt"),
                        data.table = F)
regression_sup <- regression_sup %>%
  filter(pair_pro %in% omnipath_tab$pair_pro[!(omnipath_tab$Source %in% c("NetKIN", "PhosphoNetworks", "MIMP"))] | pair_pro %in% psp_tab$pair_pro)

regression_sup <- adjust_regression_by_nonNA(regression = regression_sup, reg_nonNA = 20, reg_sig = reg_sig)
regression_sup <- annotate_ks_source(regression = regression_sup)

for (cancer_tmp in cancersprocess) {
  # for (SELF_tmp in c("cis", "trans")) {
  for (SELF_tmp in c("trans")) {
    for (enzyme_type_tmp in c("kinase")) {
      
      regression <- regression_sup %>%
        filter(enzyme_type == enzyme_type_tmp) %>%
        filter(SELF == SELF_tmp) %>%
        filter(Cancer == cancer_tmp)
      
      dat <- data.frame(pair = unique(regression$pair[regression$regulated]))
      dat[, paste0("regulated_in_", cancer_tmp)] <- T
      dat$regulated_in_other_cohort <- (dat$pair %in% regression_sup$pair[regression_sup$regulated & regression_sup$Cancer != cancer_tmp])
      dat$detected_in_other_cohort <- (dat$pair %in% regression_sup$pair[regression_sup$Cancer != cancer_tmp])
      dat[, paste0("uniquely_regulated_in_", cancer_tmp)] <- (dat$pair %in% regression$pair[regression$regulated_uniq & regression$Cancer == cancer_tmp]) & !(dat$pair %in% regression_sup$pair[regression_sup$regulated & regression_sup$Cancer != cancer_tmp])
      
      fn = paste(makeOutDir(resultD = resultD), cancer_tmp, "_", SELF_tmp, "_",  enzyme_type_tmp, '_substrate_site_level_pairs.pdf',sep ="")
      grid.newpage()
      pdf(fn, height = 12, width = 20, useDingbats = FALSE)
      fit <- euler(combinations = dat[, c(paste0("regulated_in_", cancer_tmp), "regulated_in_other_cohort", "detected_in_other_cohort", paste0("uniquely_regulated_in_", cancer_tmp)), ], input = "disjoint", shape = 'circle')
      p <-plot(fit, quantities = list(fontsize = 80), legend = list(fontsize = 50))
      grid.draw(p)
      dev.off()
    }
  }
}
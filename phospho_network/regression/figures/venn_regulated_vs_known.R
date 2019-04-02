# Yige Wu @ WashU 2018 Apr
# venn diagram showing whether regulated enzyme-substrate-phosphosite pairs are also the sites previously reported

# source ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_plotting.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
library(eulerr)


# set variables -----------------------------------------------------------
reg_nonNA <- 20
# fdr_thres <- c(0.05, 0.1); names(fdr_thres) <- c("kinase", "phosphatase")
reg_sig <- c(0.05, 0.05); names(reg_sig) <- c("kinase", "phosphatase")
SELF = "trans"

# inputs regression result ------------------------------------------------------------------

# regression <- fread(input = paste0(ppnD, "regression/tables/", "change_regression_nonNA/", 
#                                    "regression", "_", "cptac2p_cptac3", "_", "tumor", "_reg_nonNA", reg_nonNA, ".txt"), data.table = F)
# ## filtering
# regression <- regression[regression$enzyme_type_tmp == enzyme_type_tmp & regression$SELF == SELF,]
# regression <- markSigKS(regression = regression, sig_thres = fdr_thres[enzyme_type_tmp], enzyme_type_tmp = enzyme_type_tmp)

regression_sup <- fread(input = paste0(ppnD, "regression/tables/annotate_regression_with_mut_impact/", 
                                       "regression_cptac2p_cptac3_tumor_reg_nonNA20_mut_impact_cancer_specificity_annotated.txt"),
                        data.table = F)
regression_sup <- regression_sup %>%
  filter(pair_pro %in% omnipath_tab$pair_pro[!(omnipath_tab$Source %in% c("NetKIN", "PhosphoNetworks", "MIMP"))] | pair_pro %in% psp_tab$pair_pro)

regression_sup <- adjust_regression_by_nonNA(regression = regression_sup, reg_nonNA = 20, reg_sig = reg_sig)
regression_sup <- annotate_ks_source(regression = regression_sup)

for (SELF_tmp in c("cis", "trans")) {
  # first plot site-level comparison with known ones ------------------------
  for (enzyme_type_tmp in c("kinase")) {
    regression <- regression_sup %>%
      filter(enzyme_type == enzyme_type_tmp) %>%
      filter(SELF == SELF_tmp)
    regression %>%
      filter(pair == "ERBB2:ERBB2:Y1109")
    
    dat <- data.frame(pair = unique(regression$pair))
    dat$covariation <- (dat$pair %in% regression$pair[regression$regulated])
    dat$experimental <- (dat$pair %in% regression$pair[!is.na(regression$Source) & !(regression$Source %in% c("PhosphoNetworks", "MIMP", "NetKIN"))])
    dat$PhosphoNetworks <- (dat$pair %in% regression$pair[!is.na(regression$Source) & grepl(pattern = "PhosphoNetworks", x = regression$Source)])
    dat$MIMP <- (dat$pair %in% regression$pair[!is.na(regression$Source) & grepl(pattern = "MIMP", x = regression$Source)])
    dat$NetKIN <- (dat$pair %in% regression$pair[!is.na(regression$Source) & grepl(pattern = "NetKIN", x = regression$Source)])
    
    for (sequence_prediction in c("PhosphoNetworks", "MIMP", "NetKIN")) {
      fn = paste(makeOutDir(resultD = resultD), SELF_tmp, "_", enzyme_type_tmp, '_substrate_site_level_pairs_covariation_vs_experimental_vs', sequence_prediction, '.pdf',sep ="")
      grid.newpage()
      pdf(fn, height = 12, width = 20, useDingbats = FALSE)
      fit <- euler(combinations = dat[, c("covariation", "experimental", sequence_prediction)], input = "disjoint", shape = 'circle')
      p <-plot(fit, quantities = list(fontsize = 30), legend = list(fontsize = 30))
      grid.draw(p)
      dev.off()
    }
    
    fn = paste(makeOutDir(resultD = resultD), SELF_tmp, "_",  enzyme_type_tmp, '_substrate_site_level_pairs_covariation_vs_experimental.pdf',sep ="")
    grid.newpage()
    pdf(fn, height = 12, width = 20, useDingbats = FALSE)
    fit <- euler(combinations = dat[, c("covariation", "experimental")], input = "disjoint", shape = 'circle')
    p <-plot(fit, quantities = list(fontsize = 80), legend = F)
    grid.draw(p)
    dev.off()
    stop()
    
    dat <- data.frame(pair_pro = unique(regression$pair_pro))
    dat$covariation <- (dat$pair_pro %in% regression$pair_pro[regression$regulated])
    dat$experimental <- T
    dat$PhosphoNetworks <- (dat$pair_pro %in% omnipath_tab$pair_pro[!is.na(omnipath_tab$Source) & grepl(pattern = "PhosphoNetworks", x = omnipath_tab$Source)])
    dat$MIMP <- (dat$pair_pro %in% omnipath_tab$pair_pro[!is.na(omnipath_tab$Source) & grepl(pattern = "MIMP", x = omnipath_tab$Source)])
    dat$NetKIN <- (dat$pair_pro %in% omnipath_tab$pair_pro[!is.na(omnipath_tab$Source) & grepl(pattern = "NetKIN", x = omnipath_tab$Source)])
    
    for (sequence_prediction in c("PhosphoNetworks", "MIMP", "NetKIN")) {
      fn = paste(makeOutDir(resultD = resultD), SELF_tmp, "_", enzyme_type_tmp, '_substrate_protein_level_pairs_covariation_vs_experimental_vs', sequence_prediction, '.pdf',sep ="")
      grid.newpage()
      pdf(fn, height = 12, width = 20, useDingbats = FALSE)
      fit <- euler(combinations = dat[, c("covariation", "experimental", sequence_prediction)], input = "disjoint", shape = 'circle')
      p <-plot(fit, quantities = list(fontsize = 30), legend = list(fontsize = 30))
      grid.draw(p)
      dev.off()
    }
    
    fn = paste(makeOutDir(resultD = resultD), SELF_tmp, "_",  enzyme_type_tmp, '_substrate_protein_level_pairs_covariation_vs_experimental.pdf',sep ="")
    grid.newpage()
    pdf(fn, height = 12, width = 20, useDingbats = FALSE)
    fit <- euler(combinations = dat[, c("covariation", "experimental")], input = "disjoint", shape = 'circle')
    p <-plot(fit, quantities = list(fontsize = 80), legend = list(fontsize = 80))
    grid.draw(p)
    dev.off()
  }
}


# enrichement of direct targets within regulated sites --------------------
fisher.test(table(dat[, c("covariation", "experimental")]))



# regression <- fread(input = paste0(ppnD, "regression/tables/", "change_regression_nonNA/", 
#                                    "regression", "_", "cptac2p_cptac3", "_", "tumor", "_reg_nonNA", reg_nonNA, ".txt"), data.table = F)
# ## filtering
# regression <- regression[regression$enzyme_type_tmp == enzyme_type_tmp & regression$SELF == SELF,]
# regression <- markSigKS(regression = regression, sig_thres = fdr_thres[enzyme_type_tmp], enzyme_type_tmp = enzyme_type_tmp)

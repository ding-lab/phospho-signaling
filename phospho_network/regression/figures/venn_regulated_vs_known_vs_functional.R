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



# input regulatory sites --------------------------------------------------
regulatory_sites <- fread(input = "./Ding_Lab/Projects_Current/CPTAC/pan3can_shared_data/Phospho_databases/PhosphositePlus/Mar_04_2019/Regulatory_sites", data.table = F)
regulatory_sites <- regulatory_sites %>%
  mutate(SUB_MOD_RSD = str_split_fixed(string = MOD_RSD, pattern = "-", n = 2)[,1])
regulatory_sites2merge <- regulatory_sites %>%
  select(GENE, SUB_MOD_RSD, ON_FUNCTION, ON_PROCESS, ON_PROT_INTERACT, ON_OTHER_INTERACT) %>%
  unique
regression_sup <- merge(regression_sup, regulatory_sites2merge, by.x = c("SUBSTRATE", "SUB_MOD_RSD"), by.y = c("GENE", "SUB_MOD_RSD"), all.x = T)

# plot showing only regulated sites ---------------------
# for (SELF_tmp in c("cis", "trans")) {
for (SELF_tmp in c("trans")) {
  for (enzyme_type_tmp in c("kinase")) {
    ## first plot site-level comparison with known ones
    regression <- regression_sup %>%
      filter(enzyme_type == enzyme_type_tmp) %>%
      filter(SELF == SELF_tmp) %>%
      mutate(is.functional = (!is.na(ON_FUNCTION) | !is.na(ON_PROCESS) | !is.na(ON_PROT_INTERACT) | !is.na(ON_OTHER_INTERACT)))
    
    dat <- data.frame(pair = unique(regression$pair[regression$regulated]))
    dat$associated <- (dat$pair %in% regression$pair[regression$regulated])
    # dat$known <- (dat$pair %in% regression$pair[!is.na(regression$Source) & !(regression$Source %in% c("PhosphoNetworks", "MIMP", "NetKIN"))])
    dat$known <- (dat$pair %in% regression$pair[regression$is.direct])
    dat$functional <- (dat$pair %in% regression$pair[regression$is.functional])
    
    fn = paste(makeOutDir(resultD = resultD), SELF_tmp, "_",  enzyme_type_tmp, '_substrate_site_level_pairs_associated_vs_known_regulated_only.pdf',sep ="")
    grid.newpage()
    pdf(fn, height = 12, width = 20, useDingbats = FALSE)
    fit <- euler(combinations = dat[, c("associated", "known", "functional")], input = "disjoint", shape = 'circle')
    p <-plot(fit, quantities = list(fontsize = 50), legend = list(fontsize = 50))
    grid.draw(p)
    dev.off()
    
    dat <- data.frame(pair_pro = unique(regression$pair_pro))
    dat$associated <- (dat$pair_pro %in% regression$pair_pro[regression$regulated])
    dat$known <- T
    
    fn = paste(makeOutDir(resultD = resultD), SELF_tmp, "_",  enzyme_type_tmp, '_substrate_protein_level_pairs_associated_vs_known_regulated_only.pdf',sep ="")
    grid.newpage()
    pdf(fn, height = 12, width = 20, useDingbats = FALSE)
    fit <- euler(combinations = dat[, c("associated", "known")], input = "disjoint", shape = 'circle')
    p <-plot(fit, quantities = list(fontsize = 80), legend = list(fontsize = 50))
    grid.draw(p)
    dev.off()
  }
}



# plot showing all sites including not regulated ones ---------------------
# for (SELF_tmp in c("cis", "trans")) {
for (SELF_tmp in c("trans")) {
  for (enzyme_type_tmp in c("kinase")) {
    ## first plot site-level comparison with known ones
    regression <- regression_sup %>%
      filter(enzyme_type == enzyme_type_tmp) %>%
      filter(SELF == SELF_tmp) %>%
      mutate(is.functional = (!is.na(ON_FUNCTION) | !is.na(ON_PROCESS) | !is.na(ON_PROT_INTERACT) | !is.na(ON_OTHER_INTERACT)))
    regression %>%
      filter(pair == "ERBB2:ERBB2:Y1109")
    
    dat <- data.frame(pair = unique(regression$pair))
    dat$associated <- (dat$pair %in% regression$pair[regression$regulated])
    # dat$known <- (dat$pair %in% regression$pair[!is.na(regression$Source) & !(regression$Source %in% c("PhosphoNetworks", "MIMP", "NetKIN"))])
    dat$known <- (dat$pair %in% regression$pair[regression$is.direct])
    dat$functional <- (dat$pair %in% regression$pair[regression$is.functional])
    
    fn = paste(makeOutDir(resultD = resultD), SELF_tmp, "_",  enzyme_type_tmp, '_substrate_site_level_pairs_associated_vs_known.pdf',sep ="")
    grid.newpage()
    pdf(fn, height = 12, width = 20, useDingbats = FALSE)
    fit <- euler(combinations = dat[, c("associated", "known", "functional")], input = "disjoint", shape = 'circle')
    p <-plot(fit, quantities = list(fontsize = 50), legend = list(fontsize = 50))
    grid.draw(p)
    dev.off()
    
    dat <- data.frame(pair_pro = unique(regression$pair_pro))
    dat$associated <- (dat$pair_pro %in% regression$pair_pro[regression$regulated])
    dat$known <- T
    
    fn = paste(makeOutDir(resultD = resultD), SELF_tmp, "_",  enzyme_type_tmp, '_substrate_protein_level_pairs_associated_vs_known.pdf',sep ="")
    grid.newpage()
    pdf(fn, height = 12, width = 20, useDingbats = FALSE)
    fit <- euler(combinations = dat[, c("associated", "known")], input = "disjoint", shape = 'circle')
    p <-plot(fit, quantities = list(fontsize = 80), legend = list(fontsize = 50))
    grid.draw(p)
    dev.off()
  }
}

# Yige Wu @ WashU Mar 2019
## write out regulated substrates for overrepresetation gene set enrichment analysis

# source ------------------------------------------------------------------
baseD = "/Users/yigewu/Box\ Sync/"
setwd(baseD)
source("./cptac2p_analysis/phospho_network/phospho_network_shared.R")

# input regression --------------------------------------------------------
regression <- fread(input = paste0(ppnD, "regression/tables/annotate_regression_with_mut_impact/", 
                                   "regression_cptac2p_cptac3_tumor_reg_nonNA20_mut_impact_cancer_specificity_annotated.txt"),
                    data.table = F)
reg_sig <- c(0.05, 0.05); names(reg_sig) <- c("kinase", "phosphatase")
regression %>% nrow()
regression <- regression %>%
  filter(pair_pro %in% omnipath_tab$pair_pro[!(omnipath_tab$Source %in% c("NetKIN", "PhosphoNetworks", "MIMP"))] | pair_pro %in% psp_tab$pair_pro)
regression %>% nrow()

regression <- adjust_regression_by_nonNA(regression = regression, reg_nonNA = 20, reg_sig = reg_sig)
regression <- annotate_ks_source(regression = regression)
table(regression$Cancer)


# write table -------------------------------------------------------------
for (cancer_tmp in unique(regression$Cancer)) {
  regulated_substrates <- regression %>%
    filter(Cancer == cancer_tmp) %>%
    filter(SELF == "trans") %>%
    filter(regulated == T) %>%
    select(SUB_GENE) %>%
    unique()
  write.table(x = regulated_substrates$SUB_GENE, file = paste0(makeOutDir(resultD = resultD), cancer_tmp, "_regulated_substrates.txt"), quote = F, row.names = F, col.names = F)
  
  
  detected_substrates <- regression %>%
    filter(Cancer == cancer_tmp) %>%
    filter(SELF == "trans") %>%
    select(SUB_GENE) %>%
    unique()
  write.table(x = detected_substrates$SUB_GENE, file = paste0(makeOutDir(resultD = resultD), cancer_tmp, "_detected_substrates.txt"), quote = F, row.names = F, col.names = F)
}



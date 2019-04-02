# Yige Wu @ WashU Mar 2019
## choose validation pairs

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

regression %>%
  filter(regulated_uniq == T) %>%
  filter(GENE == "AKT1")

regression %>%
  filter(GENE == "AKT1") %>%
  filter(SUB_GENE == "GSK3B") %>%
  filter(SUB_MOD_RSD == "S9")

# input druggability ------------------------------------------------------
esscore_tab_outlier_drug <- fread(input = "./cptac2p/analysis_results/phospho_network/druggability/figures/grid_percent_patient_with_druggable_outlier/esscore_tab_outlier_drug_genes.txt", data.table = F)
esscore_tab_outlier_drug <- annotate_ks_source(regression = esscore_tab_outlier_drug)
esscore_tab_outlier_drug <- esscore_tab_outlier_drug %>%
  mutate(pair_cancer = paste0(pair, ":", cancer)) %>%
  mutate(regulated = (pair_cancer %in% regression$pair_cancer[regression$regulated == T])) %>%
  mutate(regulated_uniq = (pair_cancer %in% regression$pair_cancer[regression$regulated_uniq == T]))


# select pairs overlapping draggability and cancer-type-specificity -------
esscore_tab_outlier_drug_uniq <- esscore_tab_outlier_drug %>%
  filter(regulated_uniq == T)



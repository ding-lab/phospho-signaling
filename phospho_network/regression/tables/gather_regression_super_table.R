# Yige Wu @ WashU 2019 Jan
## gather regression results from multiple cancer types (multiple runs) to generate a super table

# source ------------------------------------------------------------------
baseD = "/Users/yigewu/Box\ Sync/"
wd <- getwd()
if (wd != baseD) {
  setwd(baseD)
}
source('cptac2p_analysis/phospho_network/phospho_network_shared.R')

# bussiness ---------------------------------------------------------------
## gather BRCA, OV, CO
sup_tab <- NULL
for (enzyme_type in c("kinase", "phosphatase")) {
  tab_tmp <- fread(input = paste0("./cptac2p/analysis_results/phospho_network/regression/tables/regression_omnipath&newworkin&depod&signor_protein_level/", 
                                  enzyme_type, "_substrate_regression_", "cptac2p_3can", "_tumor.txt"), data.table = F)
  sup_tab <- rbind(sup_tab, tab_tmp)
}


## gather UCEC and CCRCC
for (cancer in c("UCEC", "CCRCC")) {
  tab_tmp <- fread(input = paste0("./cptac2p/analysis_results/phospho_network/regression/tables/regression_omnipath&newworkin&depod&signor&manual_protein_level/", 
                                  "regression", "_", "cptac3", "_", cancer, "_", "tumor", "_", "PGDAC", "_", ifelse(cancer == "CCRCC", "MD_MAD", "median_polishing"), ".txt"), data.table = F)
  sup_tab <- rbind(sup_tab, tab_tmp[, colnames(sup_tab)])
}

## adjust columns
sup_tab$GENE <- as.vector(sup_tab$KINASE)
sup_tab$SUB_GENE <- as.vector(sup_tab$SUBSTRATE)
sup_tab$pair_pro <- paste0(sup_tab$GENE, ":", sup_tab$SUB_GENE)
sup_tab$pair <- paste0(sup_tab$pair_pro, ":", sup_tab$SUB_MOD_RSD)

# clean up kinase/phosphatase annotation -----------------------------------------------------------
## clean up kinase/phosphatase annotation
sup_tab$enzyme_type <- "kinase"
sup_tab$enzyme_type[!(sup_tab$GENE %in% kinases)] <- "phosphatase"

write.table(x = sup_tab, file = paste0(makeOutDir(resultD = resultD), 
                                       "regression", "_", "cptac2p_cptac3", "_", "tumor", ".txt"), 
            row.names = F, quote = F, sep = "\t")
# Yige Wu @ WashU 2019 Jan
## gather regression results from multiple cancer types (multiple runs) to generate a super table

# source ------------------------------------------------------------------
baseD = "/Users/yigewu/Box\ Sync/"
setwd(baseD)

source('cptac2p_analysis/phospho_network/phospho_network_shared.R')

# variable ----------------------------------------------------------------
data2process <- matrix(data = c("BRCA", "CDAP", "tumor", "scaled", "cptac2p",
                                "CO", "CDAP", "tumor", "scaled", "cptac2p",
                                "UCEC", "PGDAC", "tumor", "median_polishing", "cptac3",
                                "OV", "CDAP", "tumor", "scaled", "cptac2p",
                                "CCRCC", "PGDAC", "tumor", "MD_MAD", "cptac3"), ncol = 5, byrow = T)

reg_sig <- c(0.05, 0.05); names(reg_sig) <- c("kinase", "phosphatase")

# gather simulated protein pairs ---------------------------------------------------------------
sup_tab <- NULL
for (i in 1:nrow(data2process)) {
  cancer <- data2process[i,1]
  pipeline_type <- data2process[i,2]
  sample_type <- data2process[i,3]
  norm_type <- data2process[i,4]
  cptac_phase <- data2process[i,5]
  
  tn = paste0("./cptac2p/analysis_results/phospho_network/regression/tables/regression_simulated_PPI/iteration1/", cancer, "/", "regression_", cptac_phase , "_", cancer, "_", sample_type, "_", pipeline_type, "_", norm_type, "_iteration1.txt")
  tab_tmp <- fread(input = tn, data.table = F)
  tab_tmp <- tab_tmp %>%
    filter(!is.na(KINASE))
  sup_tab <- rbind(sup_tab, tab_tmp)
}
sup_tab <- adjust_regression_by_nonNA(regression = sup_tab, reg_nonNA = 20, reg_sig = reg_sig)
sup_tab %>%
  head()
## adjust columns
sup_tab$GENE <- as.vector(sup_tab$KINASE)
sup_tab$SUB_GENE <- as.vector(sup_tab$SUBSTRATE)
sup_tab$pair_pro <- paste0(sup_tab$GENE, ":", sup_tab$SUB_GENE)
sup_tab$pair <- paste0(sup_tab$pair_pro, ":", sup_tab$SUB_MOD_RSD)
sup_tab %>%
  filter(FDR_pho_kin < 0.05) %>%
  head()
sup_tab %>%
  nrow()
write.table(x = sup_tab, file = paste0(makeOutDir(resultD = resultD), 
                                       "regression_simulated_PPI", "_", "cptac2p_cptac3", "_", "tumor", ".txt"), 
            row.names = F, quote = F, sep = "\t")


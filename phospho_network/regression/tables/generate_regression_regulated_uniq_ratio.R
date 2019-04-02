# Yige Wu @ WashU 2018 Jan
# draw a grid showing the distribution of correlated kinase-substrate pairs are distributed across oncogenic pathways

# source ------------------------------------------------------------------
baseD = "/Users/yigewu/Box\ Sync/"
wd <- getwd()
if (wd != baseD) {
  setwd(baseD)
}
source("./cptac2p_analysis/phospho_network/phospho_network_shared.R")
source("./cptac2p_analysis/phospho_network/phospho_network_plotting.R")
library(pheatmap)

# set variables -----------------------------------------------------------
reg_nonNA <- 20
fdr_thres <- c(0.05, 0.05); names(fdr_thres) <- c("kinase", "phosphatase")
SELF = "trans"
num_top <- 25
enzyme_type <- "kinase"
# enzyme_type <- "phosphatase"
enzyme_types2process <- c("kinase", "phosphatase")
size <- 83
cancers2process <- c("BRCA", "OV", "CO", "UCEC", "CCRCC")

for (enzyme_type in enzyme_types2process) {
  # source("./cptac2p_analysis/phospho_network/regression/tables/generate_regression_regulated_uniq_marked.R")
  file_path_tmp <- paste0("./cptac2p/analysis_results/phospho_network/regression/figures/grid_regression_kinase_top_ratio_cance_specific/", "regression_", enzyme_type, "_substrate_size", size, "_detected_in_", paste0(cancers2process, collapse = "_"),".txt")
  regression <- fread(input = file_path_tmp, data.table = F)
  
  for (col2evaluate in c("SUB_GENE", "GENE")) {
    tab_pairs <- data.frame(table(regression[, c(col2evaluate, "Cancer", "regulated_uniq")]))
    tab_pairs %>%
      tail()
    
    tab_pairs.wratio <- merge(tab_pairs[tab_pairs$regulated_uniq == "TRUE",],
                              tab_pairs[tab_pairs$regulated_uniq == "FALSE", c(col2evaluate, "Cancer", "Freq")],
                              by = c(col2evaluate, "Cancer"), suffixes = c(".regulated_uniqT", ".regulated_uniqF"), all.y = T)
    tab_pairs.wratio %>%
      tail()
    tab_pairs.wratio$ratio.regulated_uniqT <- tab_pairs.wratio$Freq.regulated_uniqT/(tab_pairs.wratio$Freq.regulated_uniqT +  tab_pairs.wratio$Freq.regulated_uniqF)
    tab_pairs.wratio %>%
      head()
    stop()
    for (cancer_tmp in cancers2process) {
      tab4GSEA <- tab_pairs.wratio %>%
        filter(Cancer == cancer_tmp) %>%
        filter(ratio.regulated_uniqT > 0)
      tab4GSEA <- tab4GSEA[, c(col2evaluate, "ratio.regulated_uniqT")]
      write.table(x = tab4GSEA, file = paste0(makeOutDir(resultD = resultD), cancer_tmp, "_", enzyme_type, "_", col2evaluate, "_regulated_uniq_ratio2detected.txt"), 
                  quote = F, row.names = F, col.names = F, sep = "\t")
    }
  }
}
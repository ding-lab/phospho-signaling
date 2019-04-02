# Yige Wu @ WashU 2018 Jan
# generate tables marking the kinase-substrate pairs regulated uniquely across cancer types

# source ------------------------------------------------------------------
baseD = "/Users/yigewu/Box\ Sync/"
wd <- getwd()
if (wd != baseD) {
  setwd(baseD)
}
source("./cptac2p_analysis/phospho_network/phospho_network_shared.R")

# set variables -----------------------------------------------------------
reg_nonNA <- 20
SELF = "trans"
num_top <- 25
size <- 83

for (fdr_thres in c(0.05)) {
  reg_sig <- rep(0.05, 2); names(reg_sig) <- c("kinase", "phosphatase")
  summary_tab_tmp <- NULL
  for (cancer in cancers2process) {
    ## input the regression table for only pairs detectable in all cancers
    file_path_tmp <- paste0("./cptac2p/analysis_results/phospho_network/regression/figures/plot_downsize_affect_regression/", "regression_", cancer, "_size", size, "_detected_in_", paste0(cancers2process, collapse = "_"),".txt")
    tab_tmp <- fread(input = file_path_tmp, data.table = F)
    tab_tmp <- unique(tab_tmp)
    tab_tmp <- tab_tmp[order(tab_tmp$P_pho_kin, tab_tmp$pair, decreasing = T),]
    tab_tmp$GENE <- tab_tmp$KINASE
    tab_tmp$SUB_GENE <- tab_tmp$SUBSTRATE
    # tab_test <- tab_tmp[tab_tmp$pair %in% tab_tmp$pair[duplicated(tab_tmp$pair)],]
    tab_tmp <- tab_tmp[!duplicated(tab_tmp$pair),]
    tab_tmp$enzyme_type <- NULL
    tab_tmp <- annotate_enzyme_type(regression = tab_tmp, kinases = kinases, phosphatases = phosphatases)
    tab_tmp %>% nrow()
    tab_tmp <- tab_tmp %>%
      filter(enzyme_type != "") %>%
      filter(Size >= reg_nonNA)
    tab_tmp %>% nrow()
    
    ## clean up
    tab_tmp <- adjust_regression_by_nonNA(regression = tab_tmp, reg_nonNA = reg_nonNA, reg_sig = reg_sig)
    rownames(tab_tmp) <- tab_tmp$pair
    
    ## mark whether each cancer type is significantly correlated for each pair
    if (is.null(summary_tab_tmp)) {
      summary_tab_tmp <- tab_tmp[,c("pair", "SELF", "GENE", "SUB_GENE", "SUB_MOD_RSD")]
      rownames(summary_tab_tmp) <- summary_tab_tmp$pair
    }
    summary_tab_tmp[, paste0("regulated_", cancer)] <- as.numeric(tab_tmp[as.vector(summary_tab_tmp$pair), "regulated"])
  }
  
  regression <- NULL
  for (cancer in cancers2process) {
    ## input the regression table for only pairs detectable in all cancers
    file_path_tmp <- paste0("./cptac2p/analysis_results/phospho_network/regression/figures/plot_downsize_affect_regression/", "regression_", cancer, "_size", size, "_detected_in_", paste0(cancers2process, collapse = "_"),".txt")
    tab_tmp <- fread(input = file_path_tmp, data.table = F)
    tab_tmp <- tab_tmp[order(tab_tmp$P_pho_kin, tab_tmp$pair, decreasing = T),]
    tab_tmp$GENE <- tab_tmp$KINASE
    tab_tmp$SUB_GENE <- tab_tmp$SUBSTRATE
    tab_tmp %>% head()
    tab_tmp <- tab_tmp[!duplicated(tab_tmp$pair),]
    tab_tmp <- adjust_regression_by_nonNA(regression = tab_tmp, reg_nonNA = reg_nonNA, reg_sig = reg_sig)
    rownames(tab_tmp) <- tab_tmp$pair
    
    tab_tmp$regulated_uniq <- ifelse(tab_tmp$pair %in% summary_tab_tmp$pair[rowSums(summary_tab_tmp[, paste0("regulated_", cancers2process)]) == 1 & summary_tab_tmp[, paste0("regulated_", cancer)] == 1], T, F)
    
    if (is.null(regression)) {
      regression <- tab_tmp
      ordered_cols <- colnames(regression)
    } else {
      regression <- rbind(regression, tab_tmp[,ordered_cols])
      
    }
  }
  fn <- paste0(makeOutDir(resultD = resultD), "regression_size", size, "_FDR", fdr_thres, "_detected_in_", paste0(cancers2process, collapse = "_"),".txt")
  write.table(x = regression, file = fn, quote = F, row.names = F, sep = "\t")
}

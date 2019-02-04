# Yige Wu @ WashU 2019 Feb
# look at the intersections of regulated kinase-substrate pairs in BRCA, OV and CO datasets


# input -------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
source(paste0(baseD, "cptac2p_analysis/pan3can_aes.R")) # aes for general purposes
library(UpSetR)



# Plot site level cancer specificity --------------------------------------
size <- 83
tab2p <- NULL
for (cancer in cancers2process) {
  ## input the regression table for only pairs detectable in all cancers
  file_path_tmp <- paste0("./cptac2p/analysis_results/phospho_network/regression/figures/plot_downsize_affect_regression/", "regression_", cancer, "_size", size, "_detected_in_", paste0(cancers2process, collapse = "_"),".txt")
  tab_tmp <- fread(input = file_path_tmp, data.table = F)
  print(nrow(tab_tmp))
  tab_tmp <- tab_tmp[order(tab_tmp$FDR_pho_kin, tab_tmp$pair, decreasing = F),]
  tab_tmp$GENE <- tab_tmp$KINASE
  tab_tmp$SUB_GENE <- tab_tmp$SUBSTRATE
  tab_tmp %>% head()
  tab_tmp <- tab_tmp[!duplicated(tab_tmp$pair),]
  rownames(tab_tmp) <- tab_tmp$pair
  print(nrow(tab_tmp))
  
  ## mark whether each cancer type is significantly correlated for each pair
  if (is.null(tab2p)) {
    tab2p <- tab_tmp[,c("pair", "SELF", "GENE", "SUB_GENE", "SUB_MOD_RSD")]
    rownames(tab2p) <- tab2p$pair
  }
  tab2p[, paste0("regulated_", cancer)] <- as.numeric(tab_tmp[as.vector(tab2p$pair), "regulated"])
}
tab2p %>% head()
stop()
regulated_cols <- colnames(tab2p)[!(colnames(tab2p) %in% c("pair", "SELF", "GENE", "SUB_GENE", "SUB_MOD_RSD"))]
for (SELF in c("cis", "trans")) {
  fn = paste(makeOutDir(resultD = resultD), SELF, '_site_level', '_regulated_pair_overlap.pdf',sep ="")
  pdf(fn, height = 3, width = 6, onefile = F)
  upset(tab2p[tab2p$SELF == SELF,], sets.bar.color = "#56B4E9", sets = regulated_cols,
        empty.intersections = NULL, order.by = "freq")
  dev.off()
  
}


# test SMG/cancer driver gene enrichment ----------------------------------
for (cancer in cancers2process) {
  fisher_tab <- table(tab2p$SUB_GENE %in% driver_genes$Gene, (tab2p[, paste0("regulated_", cancer)]) == 1)
  print(fisher_tab)
  
  if (all(dim(fisher_tab) == c(2,2))) {
    print(fisher.test(fisher_tab, alternative = "greater"))
  }
}

for (cancer in cancers2process) {
  fisher_tab <- table(tab2p$GENE %in% driver_genes$Gene, (tab2p[, paste0("regulated_", cancer)]) == 1)
  print(fisher_tab)
  
  if (all(dim(fisher_tab) == c(2,2))) {
    print(fisher.test(fisher_tab, alternative = "greater"))
  }
}

for (cancer in cancers2process) {
  fisher_tab <- table(tab2p$GENE %in% driver_genes$Gene[driver_genes$Cancer == "KIRC"], (tab2p[, paste0("regulated_", cancer)]) == 1)
  print(fisher_tab)
  
  if (all(dim(fisher_tab) == c(2,2))) {
    print(fisher.test(fisher_tab, alternative = "greater"))
  }
}

## no enrichment of SMGs in regulated pairs
for (cancer in cancers2process) {
  fisher_tab <- table(tab2p$SUB_GENE %in% SMGs[[cancer]], (rowSums(tab2p[,regulated_cols]) == 1 & tab2p[, paste0("regulated_", cancer)]))
  if (all(dim(fisher_tab) == c(2,2))) {
    # print(fisher.test(fisher_tab))
    
  }
  print(tab2p$pair[tab2p$SUB_GENE %in% SMGs[[cancer]] & tab2p[, paste0("regulated_", cancer)]])
  
}

for (cancer in cancers2process) {
  fisher_tab <- table(tab2p$GENE %in% SMGs[[cancer]], (rowSums(tab2p[,regulated_cols]) == 1 & tab2p[, paste0("regulated_", cancer)]))
  if (all(dim(fisher_tab) == c(2,2))) {
    print(fisher_tab)
    print(fisher.test(fisher_tab))
    
  }
}

# Plot site level cancer specificity --------------------------------------
size <- 83
tab2p <- NULL
for (cancer in cancers2process) {
  ## input the regression table for only pairs detectable in all cancers
  file_path_tmp <- paste0("./cptac2p/analysis_results/phospho_network/regression/figures/plot_downsize_affect_regression/", "regression_", cancer, "_size", size, "_detected_in_", paste0(cancers2process, collapse = "_"),".txt")
  tab_tmp <- fread(input = file_path_tmp, data.table = F)
  print(nrow(tab_tmp))
  tab_tmp <- tab_tmp[order(tab_tmp$FDR_pho_kin, tab_tmp$pair, decreasing = F),]
  tab_tmp$GENE <- tab_tmp$KINASE
  tab_tmp$SUB_GENE <- tab_tmp$SUBSTRATE
  tab_tmp$pair_pro <- paste0(tab_tmp$GENE, ":", tab_tmp$SUB_GENE)
    
  tab_tmp %>% head()
  tab_tmp <- tab_tmp[!duplicated(tab_tmp$pair),]
  rownames(tab_tmp) <- tab_tmp$pair
  print(nrow(tab_tmp))
  
  ## mark whether each cancer type is significantly correlated for each pair
  if (is.null(tab2p)) {
    tab2p <- unique(tab_tmp[,c("pair_pro", "SELF", "GENE", "SUB_GENE")])
  }
  tab2p[, paste0("regulated_", cancer)] <- as.numeric(tab2p$pair_pro %in% tab_tmp$pair_pro[tab_tmp$regulated])
}
tab2p %>% head()
regulated_cols <- colnames(tab2p)[!(colnames(tab2p) %in% c("pair_pro", "SELF", "GENE", "SUB_GENE"))]
for (SELF in c("trans", "cis")) {
  fn = paste(makeOutDir(resultD = resultD), SELF, '_protein_level', '_regulated_pair_overlap.pdf',sep ="")
  pdf(fn, height = 3, width = 6, onefile = F)
  upset(tab2p[tab2p$SELF == SELF,], sets.bar.color = "#56B4E9", sets = regulated_cols,
        empty.intersections = NULL, order.by = "freq")
  dev.off()
  
}

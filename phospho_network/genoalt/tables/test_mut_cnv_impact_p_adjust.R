# Yige Wu @ WashU 2018 Aug
## take the wilcox test p-values from compari the substrate phosphorylation level of kinase-mutated tumors, kinase-amplified tumors, kinase-deleted tumors against other tumors for each kinase-substrate pair
## adjust to FDR 

# source ------------------------------------------------------------------
source("./cptac2p_analysis/phospho_network/phospho_network_shared.R")
source('/Users/yigewu/Box Sync/cptac2p_analysis/dependencies/tables/load_TCGA_pathways.R')

# set variables -----------------------------------------------------------
sig_thres <- 0.3
fdr_thres <- 0.3

# inputs ------------------------------------------------------------------
mut_cnv_cans <- fread(input = paste0(ppnD, "genoalt/figures/barplot_mut_cnv_impact/mut_cnv_sig_cans.txt"), data.table = F)

## annotate mutated genes that are oncogenes/tsgs
driver_gene_type.en <- vector(mode = "character", length = nrow(mut_cnv_cans))
driver_gene_type.en[mut_cnv_cans$GENE %in% oncogenes] <- "oncogene"
driver_gene_type.en[mut_cnv_cans$GENE %in% tsgs] <- "tsg"
mut_cnv_cans$driver_gene_type.en <- driver_gene_type.en 

# adjust p-value based by limiting to SMGs of each cancer thypes ----------------------------------------------------------------------
mut_cnv_cans_p_adjusted <- NULL
mut_cnv_cans_smgs_p_adjusted <- NULL

for (cancer in cancers_sort) {
  # cancer <- "BRCA"
  tab_can <- mut_cnv_cans
  print(nrow(tab_can))
  tab_can <- tab_can[tab_can$cancer == cancer,]
  print(nrow(tab_can))
  tab_can_smgs <- tab_can[tab_can$GENE %in% SMGs[[cancer]],]
  print(nrow(tab_can_smgs))
  
  for (self in c("cis", "trans")) {
    # self <- "cis"
    tab_can_smgs_self <- tab_can_smgs[tab_can_smgs$SELF == self,]
    tab_can_self <- tab_can[tab_can$SELF == self,]
    
    print(nrow(tab_can_smgs_self))
    
    for (genoalt_type in unique(tab_can_self$genoalt_type)) {
      # genoalt_type <- "mutation"
      tab_can_smgs_self_genoalt <- tab_can_smgs_self[tab_can_smgs_self$genoalt_type == genoalt_type,]
      tab_can_self_genoalt <- tab_can_self[tab_can_self$genoalt_type == genoalt_type,]
      
      print(nrow(tab_can_smgs_self_genoalt))
      tab_can_smgs_self_genoalt$fdr <- p.adjust(tab_can_smgs_self_genoalt$p, method = "fdr")
      tab_can_self_genoalt$fdr <- p.adjust(tab_can_self_genoalt$p, method = "fdr")
      
      ## annotate the substrates
      sub_genes2pathways <- map2TCGApathwaways(gene_list = unique(tab_can_self_genoalt$SUB_GENE), pathway_list = tcga_pathways_pluskegg_and_pathway)
      sub_genes2pathways <- data.frame(SUB_GENE = rep(x = names(sub_genes2pathways), sapply(X = sub_genes2pathways, FUN = function(x) length(x))), 
                                       SUB_GENE.path = unlist(sub_genes2pathways, use.names = F))
      tab_can_self_genoalt <- merge(tab_can_self_genoalt, sub_genes2pathways, all.x = T)
      tmp <- as.vector(tab_can_self_genoalt$SUB_GENE.path)
      tmp[is.na(tmp)] <- "NA"
      tab_can_self_genoalt$SUB_GENE.path <- tmp
      
      mut_cnv_cans_smgs_p_adjusted <- rbind(mut_cnv_cans_smgs_p_adjusted, tab_can_smgs_self_genoalt)
      mut_cnv_cans_p_adjusted <- rbind(mut_cnv_cans_p_adjusted, tab_can_self_genoalt)
    }
  }
}
write.table(x = mut_cnv_cans_p_adjusted, file = paste0(makeOutDir(resultD = resultD), "mut_cnv_cans_p_adjusted.csv"), quote = F, row.names = F, sep = ",")
write.table(x = mut_cnv_cans_smgs_p_adjusted, file = paste0(makeOutDir(resultD = resultD), "mut_cnv_cans_smgs_p_adjusted.csv"), quote = F, row.names = F, sep = ",")

mut_cans_p_adjusted_sigP <- mut_cnv_cans_p_adjusted[mut_cnv_cans_p_adjusted$p < sig_thres & mut_cnv_cans_p_adjusted$genoalt_type == "mutation",]
nrow(mut_cans_p_adjusted_sigP)
write.table(x = mut_cans_p_adjusted_sigP, file = paste0(makeOutDir(resultD = resultD), "mut_cans_p_adjusted_sigP", sig_thres, ".csv"), quote = F, row.names = F, sep = ",")

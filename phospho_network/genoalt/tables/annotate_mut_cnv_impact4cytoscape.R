# Yige Wu @ WashU 2018 July
## overview figure for mutational impact of kinases on substrates

# source ------------------------------------------------------------------
source("./Box Sync/cptac2p_analysis/phospho_network/phospho_network_plotting.R")
source("./Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R")
source("./Box Sync/cptac2p_analysis/dependencies/tables/load_TCGA_pathways.R")


# set variables -----------------------------------------------------------
show_p_thres <- 0.2
sig_p_thres <- 0.05

# inputs ------------------------------------------------------------------
mut_cnv_cans <- fread(input = paste0(ppnD, "genoalt/tables/test_mut_cna_impact_cptac2p/mut_cnv_sig_cans.txt"), data.table = F)
mut_cnv_cans$log10_P <- -log10(mut_cnv_cans$p)
## annotate the substrates
sub_genes2pathways <- map2TCGApathwaways(gene_list = unique(c(as.vector(mut_cnv_cans$SUB_GENE), as.vector(mut_cnv_cans$GENE))), pathway_list = tcga_pathways_pluskegg_and_pathway)
sub_genes2pathways <- data.frame(SUB_GENE = rep(x = names(sub_genes2pathways), sapply(X = sub_genes2pathways, FUN = function(x) length(x))), 
                                 SUB_GENE.path = unlist(sub_genes2pathways, use.names = F))

## input enzyme-substrate table
ptms_site_pairs_sup <- read_csv(paste0(ppnD, "compile_enzyme_substrate/tables/compile_omnipath_networkin_depod_signor_manual/Omnipath_NetworKin_DEPOD_SignorNotSiteMapped_manualAdded.csv"))
ptms_site_pairs_sup <- ptms_site_pairs_sup[ptms_site_pairs_sup$Source != "NetKIN" | (ptms_site_pairs_sup$Source == "NetKIN" & ptms_site_pairs_sup$networkin_score >= 5),]
ptms_site_pairs_sup <- data.frame(ptms_site_pairs_sup)


# write network file for MAPK pathway ---------------------------------------------
for (path_name in c("MAPK")) {
  path_genes <- KEGG[[keggs2add[pathway_name]]]
  path_genes <- c("ARAF", path_genes)
  network <- mut_cnv_cans[mut_cnv_cans$SUB_GENE %in% path_genes & mut_cnv_cans$GENE %in% path_genes,]
  network <- network[(network$p < show_p_thres & network$SELF == "trans") | (network$p < sig_p_thres & network$SELF == "cis"),]
  network <- unique(network)
  network <- network[network$GENE %in% c("AKT1", "ARAF", "BRAF", "MAP2K1", "MAP2K2", "MAP2K4", "MAP3K1", "MAP3K5", "MAP3K8", "MAPK1", "MAPK3"),]
  network <- network[network$SUB_GENE %in% c("AKT1", "ARAF", "BRAF", "MAP2K1", "MAP2K2", "MAP2K4", "MAP3K1", "MAP3K5", "MAP3K8", "MAPK1", "MAPK3"),]
  network$pair_pro <- paste0(network$GENE, ":", network$SUB_GENE)
  network <- merge(network, ptms_site_pairs_sup[, c("pair_pro", "enzyme_type")], all.x = T)
  network <- network[!(is.na(network$enzyme_type) & network$GENE == ""),]
  network <- network[order(network$p, network$pair_pro),]
  network <- network[!(duplicated(network$pair_pro)),]
  network$site <- paste0(network$SUB_GENE, "_", network$SUB_MOD_RSD)
  network$target_arrow_shape <- 1
  network$target_arrow_shape[!is.na(network$sig_BRCA) & (network$uniq_BRCA | network$uniq_OV | network$uniq_CO | network$shared3can)] <- 2
  network$sig_Pvalue <- (network$p < sig_p_thres)
  
  ## make another data table for showing phosphoprotein
  network2add <- network[!duplicated(network$site),]
  network2add$GENE <- as.vector(network2add$site)
  network2add$site <- as.vector(network2add$SUB_GENE)
  network2add$target_arrow_shape <- 0
  network2add$log10_P <- 1
  network2add$meddiff <- 0
  
  network4cyto <- rbind(network, network2add)
  network4cyto <- network4cyto[, c("GENE", "site", "log10_P", "meddiff", "target_arrow_shape", "cancer", "SELF", "sig_Pvalue")]
  write.table(x = network4cyto, file = paste0(makeOutDir(resultD = resultD), path_name, "_4cytoscape_network.txt"),
              quote = F, row.names = F, col.names = T, sep = "\t")
  
  node_tab_kin <- data.frame(shared_name = unique(network$GENE), text = unique(network$GENE))
  node_tab_kin$type <- ifelse(node_tab_kin$shared_name %in% ptms_site_pairs_sup$GENE[ptms_site_pairs_sup$enzyme_type == "phosphatase"], "phosphatase", "kinase")
  node_tab_site <- data.frame(shared_name = network$site[!duplicated(network$site)], text = network$SUB_MOD_RSD[!duplicated(network$site)])
  node_tab_site$type <- "phosphosite"
  node_tab_sub <- data.frame(shared_name = unique(network$SUB_GENE), text = unique(network$SUB_GENE))
  node_tab_sub$type <- "substrate"
  node_tab_sub <- node_tab_sub[!(node_tab_sub$shared_name %in% node_tab_kin$shared_name),]
  node_tab <- rbind(node_tab_kin, node_tab_site, node_tab_sub)
  write.table(x = node_tab, file = paste0(makeOutDir(resultD = resultD), path_name, "_4cytoscape_node.txt"), 
              quote = F, row.names = F, col.names = T, sep = "\t")
}

# write network file for PI3K pathway ---------------------------------------------
for (path_name in c("PI3K")) {
  path_genes <- sub_genes2pathways$SUB_GENE[sub_genes2pathways$SUB_GENE.path == path_name]
  network <- mut_cnv_cans[mut_cnv_cans$SUB_GENE %in% path_genes & mut_cnv_cans$GENE %in% path_genes,]
  network <- network[(network$p < show_p_thres & network$SELF == "trans") | (network$p < sig_p_thres & network$SELF == "cis"),]
  network <- unique(network)
  network <- network[network$GENE %in% c("AKT1", "BAD", "GSK3B", "FOXO3", "PIK3R1", "PIK3CA", "PTEN"),]
  network <- network[network$SUB_GENE %in% c("AKT1", "BAD", "GSK3B", "FOXO3", "PIK3R1", "PIK3CA", "PTEN"),]
  network$pair_pro <- paste0(network$GENE, ":", network$SUB_GENE)
  network <- merge(network, ptms_site_pairs_sup[, c("pair_pro", "enzyme_type")], all.x = T)
  network <- network[!(is.na(network$enzyme_type) & network$GENE == ""),]
  network <- network[order(network$p, network$pair_pro),]
  network <- network[!(duplicated(network$pair_pro)),]
  network$site <- paste0(network$SUB_GENE, "_", network$SUB_MOD_RSD)
  network$target_arrow_shape <- 1
  network$target_arrow_shape[!is.na(network$sig_BRCA) & (network$uniq_BRCA | network$uniq_OV | network$uniq_CO | network$shared3can)] <- 2
  network$sig_Pvalue <- (network$p < sig_p_thres)
  
  ## make another data table for showing phosphoprotein
  network2add <- network[!duplicated(network$site),]
  network2add$GENE <- as.vector(network2add$site)
  network2add$site <- as.vector(network2add$SUB_GENE)
  network2add$target_arrow_shape <- 0
  network2add$log10_P <- 1
  network2add$meddiff <- 0
  
  network4cyto <- rbind(network, network2add)
  network4cyto <- network4cyto[, c("GENE", "site", "log10_P", "meddiff", "target_arrow_shape", "cancer", "SELF", "sig_Pvalue")]
  write.table(x = network4cyto, file = paste0(makeOutDir(resultD = resultD), path_name, "_4cytoscape_network.txt"),
              quote = F, row.names = F, col.names = T, sep = "\t")
  
  node_tab_kin <- data.frame(shared_name = unique(network$GENE), text = unique(network$GENE))
  node_tab_kin$type <- ifelse(node_tab_kin$shared_name %in% ptms_site_pairs_sup$GENE[ptms_site_pairs_sup$enzyme_type == "phosphatase"], "phosphatase", "kinase")
  node_tab_site <- data.frame(shared_name = network$site[!duplicated(network$site)], text = network$SUB_MOD_RSD[!duplicated(network$site)])
  node_tab_site$type <- "phosphosite"
  node_tab_sub <- data.frame(shared_name = unique(network$SUB_GENE), text = unique(network$SUB_GENE))
  node_tab_sub$type <- "substrate"
  node_tab_sub <- node_tab_sub[!(node_tab_sub$shared_name %in% node_tab_kin$shared_name),]
  node_tab <- rbind(node_tab_kin, node_tab_site, node_tab_sub)
  write.table(x = node_tab, file = paste0(makeOutDir(resultD = resultD), path_name, "_4cytoscape_node.txt"), 
              quote = F, row.names = F, col.names = T, sep = "\t")
}

# write network file for WNT pathway ---------------------------------------------
for (path_name in c("WNT")) {
  path_genes <- as.vector(sub_genes2pathways$SUB_GENE[sub_genes2pathways$SUB_GENE.path == path_name])
  path_genes <- c(path_genes, "AXIN1", "LRP6", "LRP5", "FZD10", "FAM123B", "AXIN2", "TCF7L2", "FBXW7", "SOX9", "TCF7", "MYC", "ARID1A", "DKK1", "DKK2", "DKK3", "DKK4")
  network <- mut_cnv_cans[mut_cnv_cans$SUB_GENE %in% path_genes & mut_cnv_cans$GENE %in% path_genes,]
  network <- network[(network$p < show_p_thres & network$SELF == "trans") | (network$p < sig_p_thres & network$SELF == "cis"),]
  network <- unique(network)
  network <- network[network$GENE %in% c("APC", "CTNNB1"),]
  network <- network[network$SUB_GENE %in% c("APC", "CTNNB1"),]
  network$pair_pro <- paste0(network$GENE, ":", network$SUB_GENE)
  network <- merge(network, ptms_site_pairs_sup[, c("pair_pro", "enzyme_type")], all.x = T)
  network <- network[!(is.na(network$enzyme_type) & network$GENE == ""),]
  network <- network[order(network$p, network$pair_pro),]
  network <- network[!(duplicated(network$pair_pro)),]
  network$site <- paste0(network$SUB_GENE, "_", network$SUB_MOD_RSD)
  network$target_arrow_shape <- 1
  network$target_arrow_shape[!is.na(network$sig_BRCA) & (network$uniq_BRCA | network$uniq_OV | network$uniq_CO | network$shared3can)] <- 2
  network$sig_Pvalue <- (network$p < sig_p_thres)
  
  ## make another data table for showing phosphoprotein
  network2add <- network[!duplicated(network$site),]
  network2add$GENE <- as.vector(network2add$site)
  network2add$site <- as.vector(network2add$SUB_GENE)
  network2add$target_arrow_shape <- 0
  network2add$log10_P <- 1
  network2add$meddiff <- 0
  
  network4cyto <- rbind(network, network2add)
  network4cyto <- network4cyto[, c("GENE", "site", "log10_P", "meddiff", "target_arrow_shape", "cancer", "SELF", "sig_Pvalue")]
  write.table(x = network4cyto, file = paste0(makeOutDir(resultD = resultD), path_name, "_4cytoscape_network.txt"),
              quote = F, row.names = F, col.names = T, sep = "\t")
  
  node_tab_kin <- data.frame(shared_name = unique(network$GENE), text = unique(network$GENE))
  node_tab_kin$type <- ifelse(node_tab_kin$shared_name %in% ptms_site_pairs_sup$GENE[ptms_site_pairs_sup$enzyme_type == "phosphatase"], "phosphatase", "kinase")
  node_tab_site <- data.frame(shared_name = network$site[!duplicated(network$site)], text = network$SUB_MOD_RSD[!duplicated(network$site)])
  node_tab_site$type <- "phosphosite"
  node_tab_sub <- data.frame(shared_name = unique(network$SUB_GENE), text = unique(network$SUB_GENE))
  node_tab_sub$type <- "substrate"
  node_tab_sub <- node_tab_sub[!(node_tab_sub$shared_name %in% node_tab_kin$shared_name),]
  node_tab <- rbind(node_tab_kin, node_tab_site, node_tab_sub)
  write.table(x = node_tab, file = paste0(makeOutDir(resultD = resultD), path_name, "_4cytoscape_node.txt"), 
              quote = F, row.names = F, col.names = T, sep = "\t")
}

# write network file for PRKDC ---------------------------------------------
enzyme <- "PRKDC"
network <- mut_cnv_cans[mut_cnv_cans$genoalt_type == "mutation" & mut_cnv_cans$GENE == "PRKDC" & mut_cnv_cans$SUB_GENE %in% c("TP53BP1", "VIM"),]
network <- unique(network)
network$site <- paste0(network$SUB_GENE, "_", network$SUB_MOD_RSD)
network <- network[network$site %in% network$site[network$p < 0.05],]

## make another data table for showing phosphoprotein
network2add <- network[!duplicated(network$site) & network$SELF == "trans",]
network2add$GENE <- as.vector(network2add$site)
network2add$site <- as.vector(network2add$SUB_GENE)
network4cyto <- rbind(network, network2add)
network4cyto$log10_p <- -log10(network4cyto$p)
write.table(x = network4cyto, file = paste0(makeOutDir(resultD = resultD), "mut_cnv_", enzyme, "_4cytoscape_network.txt"),
            quote = F, row.names = F, col.names = T, sep = "\t")

node_tab_kin <- data.frame(shared_name = unique(network$GENE), text = unique(network$GENE))
node_tab_kin$type <- "kinase"
node_tab_site <- data.frame(shared_name = network$site[!duplicated(network$site)], text = network$SUB_MOD_RSD[!duplicated(network$site)])
node_tab_site$type <- "phosphosite"
node_tab_sub <- data.frame(shared_name = unique(network$SUB_GENE), text = unique(network$SUB_GENE))
node_tab_sub$type <- "substrate"
node_tab <- rbind(node_tab_kin, node_tab_site, node_tab_sub)
write.table(x = node_tab, file = paste0(makeOutDir(resultD = resultD), "mut_cnv_", enzyme, "_4cytoscape_node.txt"), 
            quote = F, row.names = F, col.names = T, sep = "\t")
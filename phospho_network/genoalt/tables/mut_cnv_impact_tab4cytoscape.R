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
mut_cnv_cans <- fread(input = paste0(ppnD, "genoalt/tables/test_mut_cnv_impact_p_adjust/mut_cnv_cans_p_adjusted.csv"), data.table = F)
## annotate the substrates
sub_genes2pathways <- map2TCGApathwaways(gene_list = unique(mut_cnv_cans$SUB_GENE), pathway_list = tcga_pathways_pluskegg_and_pathway)
sub_genes2pathways <- data.frame(SUB_GENE = rep(x = names(sub_genes2pathways), sapply(X = sub_genes2pathways, FUN = function(x) length(x))), 
                                 SUB_GENE.path = unlist(sub_genes2pathways, use.names = F))
sub_genes2pathways <- sub_genes2pathways[!(sub_genes2pathways$SUB_GENE == "APC" & sub_genes2pathways$SUB_GENE.path != "WNT") & !(sub_genes2pathways$SUB_GENE == "GSK3B" & sub_genes2pathways$SUB_GENE.path != "PI3K") & !(sub_genes2pathways$SUB_GENE == "IRS1" & sub_genes2pathways$SUB_GENE.path != "RTK RAS"),]
mut_cnv_cans$SUB_GENE.path <- NULL
mut_cnv_cans <- merge(mut_cnv_cans, sub_genes2pathways, all.x = T)

# write network file for BRCA ---------------------------------------------
for (cancer in c("BRCA")) {
  network <- mut_cnv_cans[mut_cnv_cans$cancer == cancer & mut_cnv_cans$genoalt_type == "mutation" & mut_cnv_cans$GENE == "AKT1" & mut_cnv_cans$SUB_GENE %in% c("GSK3B", "BAD", "FOXO3", "ITGB4", "NFKB1", "SRSF1", "ACIN1", "VIM"),]
  network <- unique(network)
  network$site <- paste0(network$SUB_GENE, "_", network$SUB_MOD_RSD)
  network <- network[network$p < show_p_thres,]

  ## make another data table for showing phosphoprotein
  network2add <- network[!duplicated(network$site),]
  network2add$GENE <- as.vector(network2add$site)
  network2add$site <- as.vector(network2add$SUB_GENE)
  network4cyto <- rbind(network, network2add)
  network4cyto$log10_p <- -log10(network4cyto$p)
  write.table(x = network4cyto, file = paste0(makeOutDir(resultD = resultD), "mut_cnv_", cancer, "_4cytoscape_network.txt"),
              quote = F, row.names = F, col.names = T, sep = "\t")
  
  node_tab_kin <- data.frame(shared_name = unique(network$GENE), text = unique(network$GENE))
  node_tab_kin$type <- "kinase"
  node_tab_site <- data.frame(shared_name = network$site[!duplicated(network$site)], text = network$SUB_MOD_RSD[!duplicated(network$site)])
  node_tab_site$type <- "phosphosite"
  node_tab_sub <- data.frame(shared_name = unique(network$SUB_GENE), text = unique(network$SUB_GENE))
  node_tab_sub$type <- "substrate"
  node_tab <- rbind(node_tab_kin, node_tab_site, node_tab_sub)
  write.table(x = node_tab, file = paste0(makeOutDir(resultD = resultD), "mut_cnv_", cancer, "_4cytoscape_node.txt"), 
              quote = F, row.names = F, col.names = T, sep = "\t")
}

# write network file for CO ---------------------------------------------
for (cancer in c("CO")) {
  network <- mut_cnv_cans[mut_cnv_cans$cancer == cancer & mut_cnv_cans$genoalt_type == "mutation" & mut_cnv_cans$GENE == "APC" & mut_cnv_cans$SUB_GENE %in% c("APC", "CTNNB1"),]
  network <- unique(network)
  network$site <- paste0(network$SUB_GENE, "_", network$SUB_MOD_RSD)
  network <- network[network$p < show_p_thres,]
  
  ## make another data table for showing phosphoprotein
  network2add <- network[!duplicated(network$site) & network$SELF == "trans",]
  network2add$GENE <- as.vector(network2add$site)
  network2add$site <- as.vector(network2add$SUB_GENE)
  network4cyto <- rbind(network, network2add)
  network4cyto$log10_p <- -log10(network4cyto$p)
  write.table(x = network4cyto, file = paste0(makeOutDir(resultD = resultD), "mut_cnv_", cancer, "_4cytoscape_network.txt"),
              quote = F, row.names = F, col.names = T, sep = "\t")
  
  node_tab_kin <- data.frame(shared_name = unique(network$GENE), text = unique(network$GENE))
  node_tab_kin$type <- "kinase"
  node_tab_site <- data.frame(shared_name = network$site[!duplicated(network$site)], text = network$SUB_MOD_RSD[!duplicated(network$site)])
  node_tab_site$type <- "phosphosite"
  node_tab_sub <- data.frame(shared_name = unique(network$SUB_GENE), text = unique(network$SUB_GENE))
  node_tab_sub$type <- "substrate"
  node_tab <- rbind(node_tab_kin, node_tab_site, node_tab_sub)
  write.table(x = node_tab, file = paste0(makeOutDir(resultD = resultD), "mut_cnv_", cancer, "_4cytoscape_node.txt"), 
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
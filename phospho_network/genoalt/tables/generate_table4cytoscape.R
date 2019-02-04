# Yige Wu @ WashU 2019 Jan
## genearte table to plot in cytoscape

if (!exists("tab4cytoscape")) {
  stop("input data frame for mutational impact on proteome first!")
}


# reformat the whole data from to cytoscape -------------------------------
network <- tab4cytoscape
network <- network[network$p < sig_p_thres,]
network <- unique(network)
network$pair_pro <- paste0(network$GENE, ":", network$SUB_GENE)
# network <- merge(network, ptms_site_pairs_sup[, c("pair_pro", "enzyme_type")], all.x = T)
network <- network[order(network$p, network$pair_pro),]
network <- network[!(duplicated(network$pair_pro)),]
network$site <- paste0(network$SUB_GENE, "_", network$SUB_MOD_RSD)
network$target_arrow_shape <- 1
# network$target_arrow_shape[!is.na(network$sig_BRCA) & (network$uniq_BRCA | network$uniq_OV | network$uniq_CO | network$shared3can)] <- 2
network$sig_Pvalue <- (network$p < sig_p_thres)

## make another data table for showing phosphoprotein
network2add <- network[!duplicated(network$site),]
network2add$GENE <- as.vector(network2add$site)
network2add$site <- as.vector(network2add$SUB_GENE)
network2add$target_arrow_shape <- 0
network2add$log10_Pvalue <- 1
network2add$meddiff <- 0

network4cyto <- rbind(network, network2add)
network4cyto <- network4cyto[, c("GENE", "site", "log10_Pvalue", "meddiff", "target_arrow_shape", "cancer", "SELF", "sig_Pvalue")]
write.table(x = network4cyto, file = paste0(fn, "_4cytoscape_network.txt"),
            quote = F, row.names = F, col.names = T, sep = "\t")

node_tab_kin <- data.frame(shared_name = unique(network$GENE), text = unique(network$GENE), group = unique(network$GENE))
# node_tab_kin$type <- ifelse(node_tab_kin$shared_name %in% ptms_site_pairs_sup$GENE[ptms_site_pairs_sup$enzyme_type == "phosphatase"], "phosphatase", "kinase")
node_tab_kin$type <- "kinase"
node_tab_site <- data.frame(shared_name = network$site[!duplicated(network$site)], text = network$SUB_MOD_RSD[!duplicated(network$site)], group = network$SUB_GENE[!duplicated(network$site)])
node_tab_site$type <- "phosphosite"
node_tab_sub <- data.frame(shared_name = unique(network$SUB_GENE), text = unique(network$SUB_GENE), group = unique(network$SUB_GENE))
node_tab_sub$type <- "substrate"
node_tab_sub <- node_tab_sub[!(node_tab_sub$shared_name %in% node_tab_kin$shared_name),]
node_tab <- rbind(node_tab_kin, node_tab_site, node_tab_sub)
write.table(x = node_tab, file = paste0(fn, "_4cytoscape_node.txt"), 
            quote = F, row.names = F, col.names = T, sep = "\t")


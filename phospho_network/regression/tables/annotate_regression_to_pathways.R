# Yige Wu @ WashU 2018 Jul
# barplot for the top kinases with high and low kinase-substrate pairs

# source ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_plotting.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/dependencies/tables/load_TCGA_pathways.R')

library(dplyr)
library(ggrepel)
library(colorspace)
# variables ---------------------------------------------------------------
# reg_nonNA2test <- seq(from = 5, to = 30, by = 5)
reg_nonNAs <- c(15, 25); names(reg_nonNAs) <- c("phosphatase", "kinase")
fdr_thres <- c(0.2, 0.05); names(fdr_thres) <- c("phosphatase", "kinase")


# Input phosphatase -------------------------------------------------------
enzyme_type <- "phosphatase"
reg_nonNA <- reg_nonNAs[enzyme_type]
tab_pp <- fread(input = paste0(ppnD, "regression/tables/", "change_regression_nonNA/", 
                               enzyme_type, "_substrate_regression_cptac2p_3can_tumor",
                               "_reg_nonNA", reg_nonNA, ".txt"), data.table = F)
tab_pp <- markSigKS(regression = tab_pp, sig_thres = fdr_thres[enzyme_type], enzyme_type = enzyme_type)
tab_pp$regulated <- (tab_pp$fdr_sig & tab_pp$coef_sig)
tab_pp_reg <- tab_pp[tab_pp$regulated,] 

tab_pp_common <- fread(input = paste0(ppnD, "regression/tables/", "change_regression_nonNA/", 
                               enzyme_type, "_substrate_regression_cptac2p_3can_tumor",
                               "_reg_nonNA", reg_nonNA, "_commonsite_nonNA", reg_nonNA, ".txt"), data.table = F)
tab_pp_common <- markSigSiteCan(regression = tab_pp_common, sig_thres = fdr_thres[enzyme_type], enzyme_type = enzyme_type)

# Input kinase -------------------------------------------------------
enzyme_type <- "kinase"
reg_nonNA <- reg_nonNAs[enzyme_type]
tab_pk <- fread(input = paste0(ppnD, "regression/tables/", "change_regression_nonNA/", 
                               enzyme_type, "_substrate_regression_cptac2p_3can_tumor",
                               "_reg_nonNA", reg_nonNA, ".txt"), data.table = F)
tab_pk <- markSigKS(regression = tab_pk, sig_thres = fdr_thres[enzyme_type], enzyme_type = enzyme_type)
tab_pk$regulated <- (tab_pk$fdr_sig & tab_pk$coef_sig)
tab_pk_reg <- tab_pk[tab_pk$regulated,] 

tab_pk_common <- fread(input = paste0(ppnD, "regression/tables/", "change_regression_nonNA/", 
                                      enzyme_type, "_substrate_regression_cptac2p_3can_tumor",
                                      "_reg_nonNA", reg_nonNA, "_commonsite_nonNA", reg_nonNA, ".txt"), data.table = F)
tab_pk_common <- markSigSiteCan(regression = tab_pk_common, sig_thres = fdr_thres[enzyme_type], enzyme_type = enzyme_type)

# Annotate substrates -----------------------------------------------------
sup_cans_tab_en_reg <- rbind(tab_pp_reg, tab_pk_reg)
tab_common <- rbind(tab_pp_common, tab_pk_common)

sub_genes2pathways <- map2TCGApathwaways(gene_list = unique(c(as.vector(sup_cans_tab_en_reg$SUB_GENE), as.vector(sup_cans_tab_en_reg$GENE))), pathway_list = tcga_pathways_pluskegg_and_pathway)
sub_genes2pathways <- data.frame(SUB_GENE = rep(x = names(sub_genes2pathways), sapply(X = sub_genes2pathways, FUN = function(x) length(x))), 
                                 SUB_GENE.path = unlist(sub_genes2pathways, use.names = F))
sub_genes2pathways <- sub_genes2pathways[!(sub_genes2pathways$SUB_GENE == "CDKN1B" & sub_genes2pathways$SUB_GENE.path != "Cell Cycle") & !(sub_genes2pathways$SUB_GENE == "TP53" & sub_genes2pathways$SUB_GENE.path != "TP53") & !(sub_genes2pathways$SUB_GENE == "APC" & sub_genes2pathways$SUB_GENE.path != "WNT") & !(sub_genes2pathways$SUB_GENE == "GSK3B" & sub_genes2pathways$SUB_GENE.path != "PI3K") & !(sub_genes2pathways$SUB_GENE == "IRS1" & sub_genes2pathways$SUB_GENE.path != "RTK RAS") & !(sub_genes2pathways$SUB_GENE == "AKT1" & sub_genes2pathways$SUB_GENE.path != "PI3K") & !(sub_genes2pathways$SUB_GENE == "PTEN" & sub_genes2pathways$SUB_GENE.path != "PI3K"),]
sub_genes2pathways <- rbind(data.frame(SUB_GENE = c("MAP3K1", "MAP2K4", "CDK12"), SUB_GENE.path = c("p38", "p38", "DNA Damage")), sub_genes2pathways)
sup_cans_tab_en_reg <- merge(sup_cans_tab_en_reg, sub_genes2pathways, all.x = T)
sup_cans_tab_en_reg <- merge(sup_cans_tab_en_reg, tab_common[, c("pair", "Cancer", 
                                                                    "sig_BRCA", "sig_OV", "sig_CO", 
                                                                    "shared3can", "uniq_BRCA", "uniq_OV", "uniq_CO")], by = c("pair", "Cancer"), all.x = T)

# plot --------------------------------------------------------------------
tab_gene_cis2p <- NULL
tab_gene_as_en2p <- NULL
tab_gene_as_sub2p <- NULL

# for (gene in c("AKT1", "ERBB2", "MAP2K4", "RB1", "CTNNB1", "PTEN")) {
for (gene in unique(c(unlist(SMGs), as.vector(driver_genes$Gene[driver_genes$Cancer %in% c("BRCA", "OV", "COADREAD")])))) {
    
  tab_gene_cis <- sup_cans_tab_en_reg[sup_cans_tab_en_reg$KINASE == gene & sup_cans_tab_en_reg$SELF == "cis",]
  tab_gene_cis_anno <- tab_gene_cis[(!is.na(tab_gene_cis$uniq_BRCA) & (tab_gene_cis$uniq_BRCA | tab_gene_cis$uniq_OV | tab_gene_cis$uniq_CO)),]
  tab_gene_cis_anno <- tab_gene_cis_anno[!(duplicated(tab_gene_cis_anno[, c("pair", "FDR_pho_kin", "SUB_GENE.path")])),]
  tab_gene_cis2p <- rbind(tab_gene_cis2p, tab_gene_cis_anno)
  
  tab_gene_as_en <- sup_cans_tab_en_reg[sup_cans_tab_en_reg$KINASE == gene & sup_cans_tab_en_reg$SELF == "trans",]
  tab_gene_as_en_anno <- tab_gene_as_en[!is.na(tab_gene_as_en$SUB_GENE.path) | (!is.na(tab_gene_as_en$uniq_BRCA) & (tab_gene_as_en$uniq_BRCA | tab_gene_as_en$uniq_OV | tab_gene_as_en$uniq_CO)),]
  tab_gene_as_en_anno <- tab_gene_as_en_anno[!(duplicated(tab_gene_as_en_anno[, c("pair", "FDR_pho_kin", "SUB_GENE.path")])),]
  tab_gene_as_en2p <- rbind(tab_gene_as_en2p, tab_gene_as_en_anno)
  
  tab_gene_as_sub <- sup_cans_tab_en_reg[sup_cans_tab_en_reg$SUBSTRATE == gene & sup_cans_tab_en_reg$SELF == "trans",]
  tab_gene_as_sub_anno <- tab_gene_as_sub[!is.na(tab_gene_as_sub$SUB_GENE.path) | (!is.na(tab_gene_as_sub$uniq_BRCA) & (tab_gene_as_sub$uniq_BRCA | tab_gene_as_sub$uniq_OV | tab_gene_as_sub$uniq_CO)),]
  tab_gene_as_sub_anno <- tab_gene_as_sub_anno[!(duplicated(tab_gene_as_sub_anno[, c("pair", "FDR_pho_kin", "SUB_GENE.path")])),]
  tab_gene_as_sub2p <- rbind(tab_gene_as_sub2p, tab_gene_as_sub_anno)

}
sum_gene_as_en <- data.frame(table(unique(tab_gene_as_en2p[, c("SUBSTRATE", "KINASE")])[, "KINASE"]))
colnames(sum_gene_as_en)  <- c("Gene", "num_downstream")
sum_gene_as_sub <- data.frame(table(unique(tab_gene_as_sub2p[, c("SUBSTRATE", "KINASE")])[, "SUBSTRATE"]))
colnames(sum_gene_as_sub)  <- c("Gene", "num_upstream")
sum_gene <- merge(sum_gene_as_sub, sum_gene_as_en, all = T)
sum_gene <- sum_gene[order(sum_gene$num_downstream, decreasing = T),]
sum_gene <- merge(sum_gene, sub_genes2pathways, by.x = c("Gene"), by.y = c("SUB_GENE"), all.x = T)
sum_gene <- sum_gene[order(sum_gene$SUB_GENE.path),]

tab_gene_cis2p$log10_FDR <- -log10(tab_gene_cis2p$FDR_pro_kin); tab_gene_cis2p$coef <- tab_gene_cis2p$coef_pro_kin
tab_gene_as_en2p$log10_FDR <- -log10(tab_gene_as_en2p$FDR_pho_kin); tab_gene_as_en2p$coef <- tab_gene_as_en2p$coef_pho_kin
tab_gene_as_sub2p$log10_FDR <- -log10(tab_gene_as_sub2p$FDR_pho_kin); tab_gene_as_sub2p$coef <- tab_gene_as_sub2p$coef_pho_kin

tab_genes <- rbind(tab_gene_cis2p, tab_gene_as_en2p, tab_gene_as_sub2p)

# write table to cytoscape ------------------------------------------------
for (path_name in unique(sum_gene$SUB_GENE.path)) {
  path_genes <- sum_gene$Gene[sum_gene$SUB_GENE.path == path_name]
  network <- tab_genes[tab_genes$SUB_GENE %in% path_genes | tab_genes$GENE %in% path_genes,]
  network <- unique(network)
  network$site <- paste0(network$SUB_GENE, "_", network$SUB_MOD_RSD)
  network$target_arrow_shape <- 1
  network$target_arrow_shape[!is.na(network$sig_BRCA) & (network$uniq_BRCA | network$uniq_OV | network$uniq_CO | network$shared3can)] <- 2
  
  ## make another data table for showing phosphoprotein
  network2add <- network[!duplicated(network$site),]
  network2add$GENE <- as.vector(network2add$site)
  network2add$site <- as.vector(network2add$SUB_GENE)
  network2add$target_arrow_shape <- 0
  network2add$log10_FDR <- 1
  network2add$coef <- 0
  
  network4cyto <- rbind(network, network2add)
  network4cyto <- network4cyto[, c("GENE", "site", "log10_FDR", "coef", "target_arrow_shape", "Cancer", "SELF")]
  write.table(x = network4cyto, file = paste0(makeOutDir(resultD = resultD), path_name, "_4cytoscape_network.txt"),
              quote = F, row.names = F, col.names = T, sep = "\t")
  
  node_tab_kin <- data.frame(shared_name = unique(network$GENE), text = unique(network$GENE))
  node_tab_kin$type <- ifelse(node_tab_kin$shared_name %in% tab_pp$KINASE, "phosphatase", "kinase")
  node_tab_site <- data.frame(shared_name = network$site[!duplicated(network$site)], text = network$SUB_MOD_RSD[!duplicated(network$site)])
  node_tab_site$type <- "phosphosite"
  node_tab_sub <- data.frame(shared_name = unique(network$SUB_GENE), text = unique(network$SUB_GENE))
  node_tab_sub$type <- "substrate"
  node_tab_sub <- node_tab_sub[!(node_tab_sub$shared_name %in% node_tab_kin$shared_name),]
  node_tab <- rbind(node_tab_kin, node_tab_site, node_tab_sub)
  write.table(x = node_tab, file = paste0(makeOutDir(resultD = resultD), path_name, "_4cytoscape_node.txt"), 
              quote = F, row.names = F, col.names = T, sep = "\t")
}

for (path_name in unique(sum_gene$SUB_GENE.path)) {
  path_genes <- sum_gene$Gene[sum_gene$SUB_GENE.path == path_name]
  network <- tab_genes[tab_genes$SUB_GENE %in% path_genes | tab_genes$GENE %in% path_genes,]
  network <- network[!is.na(network$uniq_BRCA) & (network$uniq_BRCA | network$uniq_OV | network$uniq_CO | network$shared3can),]
  if (nrow(network) > 0) {
    network <- unique(network)
    network$site <- paste0(network$SUB_GENE, "_", network$SUB_MOD_RSD)
    network$target_arrow_shape <- 2
    
    ## make another data table for showing phosphoprotein
    network2add <- network[!duplicated(network$site),]
    network2add$GENE <- as.vector(network2add$site)
    network2add$site <- as.vector(network2add$SUB_GENE)
    network2add$target_arrow_shape <- 0
    network2add$log10_FDR <- 1
    network2add$coef <- 0
    network4cyto <- rbind(network, network2add)
    network4cyto <- network4cyto[, c("GENE", "site", "log10_FDR", "coef", "target_arrow_shape", "Cancer","SELF")]
    write.table(x = network4cyto, file = paste0(makeOutDir(resultD = resultD), path_name, "_commonsite_4cytoscape_network.txt"),
                quote = F, row.names = F, col.names = T, sep = "\t")
    
    node_tab_kin <- data.frame(shared_name = unique(network$GENE), text = unique(network$GENE))
    node_tab_kin$type <- ifelse(node_tab_kin$shared_name %in% tab_pp$KINASE, "phosphatase", "kinase")
    node_tab_site <- data.frame(shared_name = network$site[!duplicated(network$site)], text = network$SUB_MOD_RSD[!duplicated(network$site)])
    node_tab_site$type <- "phosphosite"
    node_tab_sub <- data.frame(shared_name = unique(network$SUB_GENE), text = unique(network$SUB_GENE))
    node_tab_sub$type <- "substrate"
    node_tab_sub <- node_tab_sub[!(node_tab_sub$shared_name %in% node_tab_kin$shared_name),]
    node_tab <- rbind(node_tab_kin, node_tab_site, node_tab_sub)
    write.table(x = node_tab, file = paste0(makeOutDir(resultD = resultD), path_name, "_commonsite_4cytoscape_node.txt"), 
                quote = F, row.names = F, col.names = T, sep = "\t")
  }
}


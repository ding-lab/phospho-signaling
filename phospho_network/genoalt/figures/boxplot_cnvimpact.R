# Yige Wu @ WashU 2018 July
## overview figure for mutational impact of kinases on substrates

# source ------------------------------------------------------------------
source("./cptac2p_analysis/phospho_network/phospho_network_plotting.R")
source("./cptac2p_analysis/phospho_network/phospho_network_shared.R")
source('/Users/yigewu/Box Sync/cptac2p_analysis/dependencies/tables/load_TCGA_pathways.R')

library(ggrepel)
# set variables -----------------------------------------------------------
sig_thres <- 0.05
id_vars <- c("GENE", "SUB_GENE", "SUB_MOD_RSD", "cancer")
prefix_vars <- c("p", "meddiff_bottom", "meddiff_upper", "meddiff")
freq_thres <- 1
top_num <- 15
networkin_thres <- 5
top_substrate2show <- 1

# inputs ------------------------------------------------------------------
## input enzyme-substrate table
# ptms_site_pairs_sup <- read_csv(paste0(ppnD, "compile_enzyme_substrate/tables/compile_omnipath_networkin_depod/omnipath_networkin_enzyme_substrate_site_level_union_source_cleaned_addDEPOD_extended.csv"))
ptms_site_pairs_sup <- read_csv(paste0(ppnD, "compile_enzyme_substrate/tables/compile_omnipath_networkin_depod_signor/omnipath_networkin_DEPOD_SignorNotSiteMapped.csv"))
# ptms_site_pairs_sup <- ptms_site_pairs_sup[ptms_site_pairs_sup$Source == "NetKIN" | (ptms_site_pairs_sup$Source != "NetKIN" & ptms_site_pairs_sup$networkin_score >= 5),]
ptms_site_pairs_sup <- ptms_site_pairs_sup[ptms_site_pairs_sup$Source != "NetKIN" | (ptms_site_pairs_sup$Source == "NetKIN" & ptms_site_pairs_sup$networkin_score >= networkin_thres),]
ptms_site_pairs_sup <- data.frame(ptms_site_pairs_sup)

## input drive gene list
driver_genes <- read_excel("./Ding_Lab/Projects_Current/TCGA_data/gene_lists/mmc1.xlsx", 
                           sheet = "Table S1", skip = 3)
driver_genes <- data.frame(driver_genes)
oncogenes <- driver_genes$Gene[grepl(x = driver_genes$Tumor.suppressor.or.oncogene.prediction..by.20.20.., pattern = "oncogene")]
tsgs <- driver_genes$Gene[grepl(x = driver_genes$Tumor.suppressor.or.oncogene.prediction..by.20.20.., pattern = "tsg")]
oncogenes <- c(oncogenes, "CTNND1", "PIK3R1")
tsgs <- tsgs[!(tsgs %in% c("CDH1", "CTNND1", "PIK3R1"))]
oncogenes <- unique(oncogenes)
tsgs <- tsgs(oncogenes)

# count per enzyme  ------------------------
count_enzyme_cans <- NULL
mut_cnv_sig_cans <- NULL
for (cancer in cancers_sort) {
  mut_cnv_tab <- fread(input = paste0(ppnD, "genoalt/tables/mut_cnv_impact/", cancer, "_mut_cnv_tab.txt"), data.table = F)
  mut_cnv_tab$pair_pro <- paste0(mut_cnv_tab$GENE, ":", mut_cnv_tab$SUB_GENE)
  mut_cnv_tab <- mut_cnv_tab[mut_cnv_tab$pair_pro %in% ptms_site_pairs_sup$pair_pro,]
  mut_cnv_tab$pair <- paste0(mut_cnv_tab$GENE, ":", mut_cnv_tab$SUB_GENE, ":", mut_cnv_tab$SUB_MOD_RSD)
  mut_cnv_tab <- merge(mut_cnv_tab, unique(ptms_site_pairs_sup[, c("GENE", "enzyme_type")]), all.x = T)
  ## filter mutation-impacted kinase-substrate pairs
  mut_cnv_mut_tab <- mut_cnv_tab
  mut_cnv_mut_tab <- mut_cnv_tab[mut_cnv_tab$p_mut > 0 & mut_cnv_tab$p_mut < sig_thres,]
  mut_cnv_mut_tab <- mut_cnv_mut_tab[(mut_cnv_mut_tab$meddiff_bottom_mut > 0) == (mut_cnv_mut_tab$meddiff_upper_mut > 0),]
  mut_cnv_mut_tab <- mut_cnv_mut_tab[order(mut_cnv_mut_tab$p_mut),]
  count_enzyme_mut <- data.frame(table(mut_cnv_mut_tab$GENE))
  mut_sig_can <- mut_cnv_mut_tab[, c(id_vars, paste0(paste0(prefix_vars, "", sep = "_"), "mut"))]
  colnames(mut_sig_can) <- c(id_vars, prefix_vars)
  if (nrow(count_enzyme_mut) > 0){
    count_enzyme_mut$genoalt_type <- "mutation"
    mut_sig_can$genoalt_type <- "mutation"
  }

  ## filter amplification-impacted kinase-substrate pairs
  mut_cnv_amp_tab <- mut_cnv_tab
  mut_cnv_amp_tab <- mut_cnv_tab[mut_cnv_tab$p_amp > 0 & mut_cnv_tab$p_amp < sig_thres,]
  mut_cnv_amp_tab <- mut_cnv_amp_tab[((mut_cnv_amp_tab$meddiff_bottom_amp > 0) & (mut_cnv_amp_tab$meddiff_upper_amp > 0 ) & mut_cnv_amp_tab$enzyme_type == "kinase") | ((mut_cnv_amp_tab$meddiff_bottom_amp < 0) & (mut_cnv_amp_tab$meddiff_upper_amp < 0 ) & mut_cnv_amp_tab$enzyme_type == "phosphatase"),]
  mut_cnv_amp_tab <- mut_cnv_amp_tab[order(mut_cnv_amp_tab$p_amp),]
  count_enzyme_amp <- data.frame(table(mut_cnv_amp_tab$GENE))
  count_enzyme_amp$genoalt_type <- "amplification"
  amp_sig_can <- mut_cnv_amp_tab[, c(id_vars, paste0(paste0(prefix_vars, "", sep = "_"), "amp"))]
  colnames(amp_sig_can) <- c(id_vars, prefix_vars)
  amp_sig_can$genoalt_type <- "amplification"
  
  ## filter amplification-impacted kinase-substrate pairs
  mut_cnv_del_tab <- mut_cnv_tab
  mut_cnv_del_tab <- mut_cnv_tab[mut_cnv_tab$p_del > 0 & mut_cnv_tab$p_del < sig_thres,]
  mut_cnv_del_tab <- mut_cnv_del_tab[((mut_cnv_del_tab$meddiff_bottom_del > 0) & (mut_cnv_del_tab$meddiff_upper_del > 0 ) & mut_cnv_del_tab$enzyme_type == "phosphatase") | ((mut_cnv_del_tab$meddiff_bottom_del < 0) & (mut_cnv_del_tab$meddiff_upper_del < 0 ) & mut_cnv_del_tab$enzyme_type == "kinase"),]
  mut_cnv_del_tab <- mut_cnv_del_tab[order(mut_cnv_del_tab$p_del),]
  count_enzyme_del <- data.frame(table(mut_cnv_del_tab$GENE))
  count_enzyme_del$genoalt_type <- "deletion"
  del_sig_can <- mut_cnv_del_tab[, c(id_vars, paste0(paste0(prefix_vars, "", sep = "_"), "del"))]
  colnames(del_sig_can) <- c(id_vars, prefix_vars)
  del_sig_can$genoalt_type <- "deletion"
  
  count_enzyme_can <- rbind(count_enzyme_mut, count_enzyme_amp, count_enzyme_del)
  colnames(count_enzyme_can) <- c("GENE", "Freq", "genoalt_type")
  count_enzyme_can <- merge(count_enzyme_can, unique(ptms_site_pairs_sup[, c("GENE", "enzyme_type")]), all.x = T)
  count_enzyme_can$cancer <- cancer
  count_enzyme_cans <- rbind(count_enzyme_cans, count_enzyme_can)
  
  mut_cnv_sig_can <- rbind(mut_sig_can, amp_sig_can, del_sig_can)
  mut_cnv_sig_cans <- rbind(mut_cnv_sig_cans, mut_cnv_sig_can)
}
count_enzyme_cans <- data.frame(count_enzyme_cans)
count_enzyme_cans$driver_type <- ""
count_enzyme_cans$driver_type[count_enzyme_cans$GENE %in% oncogenes] <- "oncogene"
count_enzyme_cans$driver_type[count_enzyme_cans$GENE %in% tsgs] <- "tsg"

## annotate enzyme and substrate
mut_cnv_sig_cans$enzyme_driver_type <- ""
mut_cnv_sig_cans$enzyme_driver_type[mut_cnv_sig_cans$GENE %in% oncogenes] <- "oncogene"
mut_cnv_sig_cans$enzyme_driver_type[mut_cnv_sig_cans$GENE %in% tsgs] <- "tsg"
mut_cnv_sig_cans$substrate_driver_type <- ""
mut_cnv_sig_cans$substrate_driver_type[mut_cnv_sig_cans$SUB_GENE %in% oncogenes] <- "oncogene"
mut_cnv_sig_cans$substrate_driver_type[mut_cnv_sig_cans$SUB_GENE %in% tsgs] <- "tsg"
write.table(x = mut_cnv_sig_cans, file = paste0(makeOutDir(resultD = resultD), "mut_cnv_sig_cans.txt"), row.names = F, quote = F, sep = "\t")

# amp & del  ----------------------------------------
for (genoalt_type in c("amplification", "deletion")) {#
  tab2p <- count_enzyme_cans
  tab2p <- tab2p[tab2p$genoalt_type == genoalt_type,]
  tab2p <- tab2p[tab2p$Freq > freq_thres,]
  ## divide genes into 2 types: only in one cancer type or >1 cancer types
  count_can_gene <- data.frame(table(tab2p$GENE))
  tab2p$cancer_type_specific <- ""
  tab2p$cancer_type_specific[tab2p$GENE %in% count_can_gene$Var1[count_can_gene$Freq == 1]] <- as.vector(tab2p$cancer[tab2p$GENE %in% count_can_gene$Var1[count_can_gene$Freq == 1]])
  tab2p$cancer_type_specific[tab2p$GENE %in% count_can_gene$Var1[count_can_gene$Freq > 1]] <- "multi-cancer"
  tab2p$cancer_type_specific <- factor(tab2p$cancer_type_specific, levels = c(cancers_sort, "multi-cancer"))
  tab2p$x <- paste0(tab2p$GENE, ":", tab2p$cancer)

  ## sort genes that only appear in one cancer type
  tab2p_sorted <- NULL
  for (cancer in cancers_sort) {
    tab2p_tmp <- tab2p[tab2p$cancer_type_specific == cancer,]
    tab2p_tmp <- tab2p_tmp[order(tab2p_tmp$Freq, decreasing = T),]
    tab2p_tmp <- tab2p_tmp[1:min(top_num, nrow(tab2p_tmp)),]
    tab2p_sorted <- rbind(tab2p_sorted, tab2p_tmp)
  }
  genes_cancer_type_specific <- as.vector(tab2p_sorted$GENE)
  
  ## sort genes that appear in multiple cancer type
  tab2p <- tab2p[tab2p$cancer_type_specific == "multi-cancer",]
  tab2p_sorted <- NULL
  top_genes <- NULL
  for (cancer in cancers_sort) {
    tab2p_tmp <- tab2p[tab2p$cancer_type_specific == "multi-cancer" & tab2p$cancer == cancer,]
    tab2p_tmp <- tab2p_tmp[order(tab2p_tmp$Freq, decreasing = T),]
    top_genes <- c(top_genes, as.vector(tab2p_tmp$GENE[1:min(top_num, nrow(tab2p_tmp))]))
    tab2p_sorted <- rbind(tab2p_sorted, tab2p_tmp)
  }
  top_genes_tab <- data.frame(table(top_genes))
  top_genes_tab <- top_genes_tab[order(top_genes_tab$Freq, decreasing = T),]
  tab2p_sorted <- tab2p_sorted[tab2p_sorted$GENE %in% top_genes_tab$top_genes[1:min(nrow(top_genes_tab), top_num)],]
  genes_cancer_type_shared <- as.vector(tab2p_sorted$GENE[!duplicated(tab2p_sorted$GENE)])
  genes_sort <- c(genes_cancer_type_shared, genes_cancer_type_specific)
  
  genoalt_sig_cans <- mut_cnv_sig_cans[mut_cnv_sig_cans$GENE %in% genes_sort & mut_cnv_sig_cans$genoalt_type == genoalt_type,]
  genoalt_sig_cans <- merge(genoalt_sig_cans, unique(ptms_site_pairs_sup[, c("GENE", "enzyme_type")]), all.x = T)
  
  tab4plot <- NULL
  for (gene in genes_sort) {
    tab_tmp <- genoalt_sig_cans[genoalt_sig_cans$GENE == gene,]
    tab4plot <- rbind(tab4plot, unique(tab_tmp))
  }
  tab4plot$log10_pvalue <- -log10(tab4plot$p)
  tab4plot$GENE <- factor(tab4plot$GENE, levels = genes_sort)
  tab4plot$text <- paste0(tab4plot$SUB_GENE, "_", tab4plot$SUB_MOD_RSD)
  if (genoalt_type == "amplification") {
    driver_color <- ifelse(genes_sort %in% oncogenes, set1[1], "black")
    meddiff_cap <- quantile(x = genoalt_sig_cans$meddiff, probs = 0.99)
    tab4plot <- tab4plot[tab4plot$meddiff <= meddiff_cap, ]
  }
  if (genoalt_type == "deletion") {
    driver_color <- ifelse(genes_sort %in% tsgs, set1[2], "black")
    meddiff_cap <- quantile(x = genoalt_sig_cans$meddiff, probs = 0.01)
    tab4plot <- tab4plot[tab4plot$meddiff >= meddiff_cap, ]
  }
  enzyme_face <- ifelse(genes_sort %in% ptms_site_pairs_sup$GENE[ptms_site_pairs_sup$enzyme_type == "phosphatase"], "bold", "plain")
  p <- ggplot(tab4plot)
  p <- p + geom_point(data=tab4plot, aes(x = GENE, y = meddiff, color = cancer, size = log10_pvalue), 
                      shape = 16, stroke = 0, alpha = 0.6)
  if (genoalt_type == "amplification") {
    p <- p + geom_text_repel(data = tab4plot[((tab4plot$enzyme_driver_type == "oncogene") & ((tab4plot$meddiff > 0 & tab4plot$substrate_driver_type == "oncogene" & tab4plot$enzyme_type == "kinase"))) | tab4plot$enzyme_type == "phosphatase",],
                             aes(x = GENE, y = meddiff, label = text, color = cancer), size = 3, force = 2)
  }
  if (genoalt_type == "deletion") {
    p <- p + geom_text_repel(data = tab4plot[tab4plot$enzyme_driver_type == "tsg" & ((tab4plot$meddiff < 0 & tab4plot$substrate_driver_type == "tsg" & tab4plot$enzyme_type == "kinase")) | (tab4plot$enzyme_type == "phosphatase"),],
                             aes(x = GENE, y = meddiff, label = text, color = cancer), size = 3, force = 2)
  }
  p <- p + geom_hline(yintercept = 0, linetype = 2, color = "grey")
  p <- p + scale_color_manual(values = color_cancers2)
  p <- p + theme_minimal()
  p <- p + xlab('enzyme')+ylab("(median substrate phosphorylation difference)")
  p <- p + theme(axis.title=element_text(size=10))
  p <- p + theme(axis.text.y = element_text(colour="black", size=10))
  p <- p + theme(axis.text.x = element_text(size = 8, angle = 90, hjust = 0.5, vjust = 0.5, colour = as.vector(tab4plot$driver_color), face = as.vector(tab4plot$enzyme_face)))
  p <- p + theme(axis.text.x = element_text(size = 10, angle = 90, hjust = 0.5, vjust = 0.5, colour = driver_color, face = enzyme_face))
  p
  fn = paste(makeOutDir(resultD = resultD), 'enzyme_', genoalt_type, '_affect_phosphosites_p_value', sig_thres,'_phosphosites_scatterplot.pdf',sep ="")
  ggsave(filename = fn, width = 15, height = 6)
}


# plot mutation-impacted substrate pairs ----------------------------------
for (genoalt_type in c("mutation")) {
  tab2p <- count_enzyme_cans
  tab2p <- tab2p[tab2p$genoalt_type == genoalt_type,]
  tab2p <- tab2p[tab2p$Freq > 1,]
  
  
  ## divide genes into 2 types: only in one cancer type or >1 cancer types
  count_can_gene <- data.frame(table(tab2p$GENE))
  tab2p$cancer_type_specific <- ""
  tab2p$cancer_type_specific[tab2p$GENE %in% count_can_gene$Var1[count_can_gene$Freq == 1]] <- as.vector(tab2p$cancer[tab2p$GENE %in% count_can_gene$Var1[count_can_gene$Freq == 1]])
  tab2p$cancer_type_specific[tab2p$GENE %in% count_can_gene$Var1[count_can_gene$Freq > 1]] <- "multi-cancer"
  tab2p$cancer_type_specific <- factor(tab2p$cancer_type_specific, levels = c(cancers_sort, "multi-cancer"))
  tab2p$x <- paste0(tab2p$GENE, ":", tab2p$cancer)

  genes_sort <- unique(tab2p$GENE)
  genoalt_sig_cans <- mut_cnv_sig_cans[mut_cnv_sig_cans$GENE %in% genes_sort & mut_cnv_sig_cans$genoalt_type == genoalt_type,]
  genoalt_sig_cans <- merge(genoalt_sig_cans, unique(ptms_site_pairs_sup[, c("GENE", "enzyme_type")]), all.x = T)
  
  tab4plot <- NULL
  for (gene in genes_sort) {
    tab_tmp <- genoalt_sig_cans[genoalt_sig_cans$GENE == gene,]
    tab4plot <- rbind(tab4plot, unique(tab_tmp))
  }
  ## annotate the substrates
  sub_genes2pathways <- map2TCGApathwaways(gene_list = unique(tab4plot$SUB_GENE), pathway_list = tcga_pathways_pluskegg_and_pathway)
  sub_genes2pathways <- data.frame(SUB_GENE = rep(x = names(sub_genes2pathways), sapply(X = sub_genes2pathways, FUN = function(x) length(x))), 
                                   SUB_GENE.path = unlist(sub_genes2pathways, use.names = F))
  tab4plot <- merge(tab4plot, sub_genes2pathways, all.x = T)
  tmp <- as.vector(tab4plot$SUB_GENE.path)
  tmp[is.na(tmp)] <- "NA"
  tab4plot$SUB_GENE.path <- tmp
  
  ## annotate driver genes to TCGA pathways
  genes2pathways <- map2TCGApathwaways(gene_list = unique(as.vector(tab4plot$GENE)), pathway_list = tcga_pathways)
  genes2pathways <- data.frame(GENE = rep(x = names(genes2pathways), sapply(X = genes2pathways, FUN = function(x) length(x))), 
                               GENE.path = unlist(genes2pathways, use.names = F))
  tab4plot <- merge(tab4plot, genes2pathways, all.x = T)
  
  ## annotate features for ploting
  tab4plot$log10_pvalue <- -log10(tab4plot$p)
  tab4plot$text <- paste0(tab4plot$SUB_GENE, "_", tab4plot$SUB_MOD_RSD)
  tab4plot$same_path <- (as.vector(tab4plot$GENE.path) == as.vector(tab4plot$SUB_GENE.path))
  tab4plot$pair_pro <- paste0(tab4plot$GENE, ":", tab4plot$SUB_GENE)
  tab4plot$pair <- paste0(tab4plot$pair_pro, "_", tab4plot$SUB_MOD_RSD)
  tab4plot <- tab4plot[order(tab4plot$pair, tab4plot$same_path, decreasing = T),]
  tab4plot <- tab4plot[!duplicated(tab4plot$pair),]
  tab4plot$coef <- as.vector(tab4plot$meddiff)
  
  # annotate the substrate with top effect size by pathway ------------------
  top_effect_size_path <- vector(mode = "logical", length = nrow(tab4plot))
  for (gene in unique(tab4plot$GENE)) {
    tab_gene <- tab4plot[tab4plot$GENE == gene,]
    for (cancer in unique(tab_gene$Cancer)) {
      tab_can <- tab_gene[tab_gene$Cancer == cancer,]
      for (pathway_name in unique(tab_can$SUB_GENE.path)) {
        tab_path <- tab_can[tab_can$SUB_GENE.path == pathway_name,]
        abs_coefs <- abs(tab_path$coef)
        abs_coefs <- abs_coefs[order(abs_coefs, decreasing = T)]
        top_coefs <- abs_coefs[1:min(length(abs_coefs), top_substrate2show)]
        top_effect_size_path[abs(tab4plot$coef) %in%  top_coefs] <- TRUE
      }
    }
  }
  tab4plot$top_effect_size_path <- top_effect_size_path
  
  # hide the text label for not annotated substrates -----------------------------------------------------
  tab4plot$text2p <- as.vector(tab4plot$text)
  tab4plot$text2p[tab4plot$SUB_GENE.path == "NA" & (tab4plot$GENE %in% as.vector(tab4plot$GENE[tab4plot$SUB_GENE.path != "NA"]))] <- ""
  tab4plot$text2p[tab4plot$SUB_GENE.path == "NA" & !(tab4plot$GENE %in% as.vector(tab4plot$GENE[tab4plot$SUB_GENE.path != "NA"])) & !tab4plot$top_effect_size_path] <- ""
  
  enzyme_face <- ifelse(genes_sort %in% ptms_site_pairs_sup$GENE[ptms_site_pairs_sup$enzyme_type == "phosphatase"], "bold", "plain")
  names(enzyme_face) <- genes_sort
  driver_color <- rep("black", length(genes_sort)) 
  driver_color[genes_sort %in% oncogenes] <- set1[1]
  driver_color[genes_sort %in%  tsgs] <- set1[2]
  names(driver_color) <- genes_sort
  
  tab4plot$GENE <- factor(tab4plot$GENE, levels = genes_sort)

  tab4plot_driver <- tab4plot[tab4plot$GENE %in% driver_genes$Gene,]
  tab4plot_driver$GENE <- factor(tab4plot_driver$GENE, levels = genes_sort[genes_sort %in% as.vector(driver_genes$Gene)])
  
  p <- ggplot(tab4plot_driver)
  p <- p + geom_point(data=tab4plot_driver, aes(x = GENE, y = meddiff, color = cancer, size = log10_pvalue), 
                      shape = 16, stroke = 0, alpha = 0.6)
  p <- p + geom_text_repel(data = tab4plot_driver, 
                           aes(x = GENE, y = meddiff, label = text2p, color = SUB_GENE.path, alpha = ifelse(top_effect_size_path, 0.8, 0.2)), 
                           size = 2, force = 2)
  p <- p + geom_hline(yintercept = 0, linetype = 2, color = "grey")
  p <- p + scale_color_manual(values = c(color_cancers2, tcga_path_colors))
  p <- p + theme_minimal()
  p <- p + xlab('enzyme')+ylab("(median substrate phosphorylation difference)")
  p <- p + theme(axis.title=element_text(size=10))
  p <- p + theme(axis.text.y = element_text(colour="black", size=10))
  p <- p + theme(axis.text.x = element_text(size = 18, angle = 90, hjust = 0.5, vjust = 0.5, colour = driver_color[genes_sort %in% as.vector(driver_genes$Gene)], face = enzyme_face[genes_sort %in% as.vector(driver_genes$Gene)]))
  p
  fn = paste(makeOutDir(resultD = resultD), 'driver_enzyme_', genoalt_type, '_affect_phosphosites_p_value', sig_thres,'_phosphosites_scatterplot.pdf',sep ="")
  ggsave(filename = fn, width = 9, height = 3)
}



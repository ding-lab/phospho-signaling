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
## input enzyme-substrate table
ptms_site_pairs_sup <- read_csv(paste0(ppnD, "compile_enzyme_substrate/tables/compile_omnipath_networkin_depod_signor_manual/Omnipath_NetworKin_DEPOD_SignorNotSiteMapped_manualAdded.csv"))
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
tsgs <- unique(tsgs)

## input significant mutation impact
mut_cnv_sig_cans <- fread(input = paste0(ppnD, "genoalt/figures/barplot_mut_cnv_impact/mut_cnv_sig_cans_sig_thres", sig_thres, ".txt"), data.table = F)
count_enzyme_cans <- fread(input = paste0(ppnD, "genoalt/figures/barplot_mut_cnv_impact/count_enzyme_cans_sig_thres", sig_thres, ".txt"), data.table = F)

# amp & del  ----------------------------------------
for (genoalt_type in c("amplification", "deletion")) {#
  tab2p <- mut_cnv_sig_cans
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
  p <- p + geom_errorbar(aes(x = GENE, ymin=meddiff_bottom, ymax=meddiff_upper), width=.2, position=position_dodge(0.05))
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
  for (self in c("trans")) {
    tab2p <- count_enzyme_cans
    tab2p <- tab2p[tab2p$genoalt_type == genoalt_type & tab2p$SELF == self,]
    tab2p <- tab2p[tab2p$Freq > 0,]
    tab2p <- tab2p[order(tab2p$Freq, decreasing = T),]
    
    genes_sort <- unique(tab2p$GENE)
    genoalt_sig_cans <- mut_cnv_sig_cans[mut_cnv_sig_cans$GENE %in% genes_sort & mut_cnv_sig_cans$genoalt_type == genoalt_type & mut_cnv_sig_cans$SELF == self,]
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
    tab4plot$text <- paste0(tab4plot$SUB_GENE, "\n", tab4plot$SUB_MOD_RSD)
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
      for (cancer in unique(tab_gene$cancer)) {
        tab_can <- tab_gene[tab_gene$cancer == cancer,]
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
    
    top_effect_size <- vector(mode = "logical", length = nrow(tab4plot))
    for (gene in unique(tab4plot$GENE)) {
      tab_gene <- tab4plot[tab4plot$GENE == gene,]
      for (cancer in unique(tab_gene$cancer)) {
        tab_can <- tab_gene[tab_gene$cancer == cancer,]
        abs_coefs <- abs(tab_can$coef)
        abs_coefs <- abs_coefs[order(abs_coefs, decreasing = T)]
        top_coefs <- abs_coefs[1:min(length(abs_coefs), 1+top_substrate2show)]
        top_effect_size[abs(tab4plot$coef) %in%  top_coefs] <- TRUE
      }
    }
    tab4plot$top_effect_size <- top_effect_size
    tab4plot$GENE_can <- paste0(tab4plot$GENE, ":", tab4plot$cancer)
    
    # hide the text label for not annotated substrates -----------------------------------------------------
    tab4plot$text2p <- as.vector(tab4plot$text)
    tab4plot$text2p[tab4plot$SUB_GENE.path == "NA" & !tab4plot$top_effect_size_path] <- ""
    tab4plot$text2p[tab4plot$SUB_GENE.path == "NA" & (tab4plot$GENE_can %in% as.vector(tab4plot$GENE_can[tab4plot$SUB_GENE.path != "NA"]))] <- ""
    
    # add alpha colunn according to significance ------------------------------
    tab4plot$alpha <- tab4plot$log10_pvalue/max(tab4plot$log10_pvalue, na.rm = T)
    
    enzyme_face <- ifelse(genes_sort %in% ptms_site_pairs_sup$GENE[ptms_site_pairs_sup$enzyme_type == "phosphatase"], "bold", "plain")
    names(enzyme_face) <- genes_sort
    driver_color <- rep("black", length(genes_sort)) 
    driver_color[genes_sort %in% oncogenes] <- set1[1]
    driver_color[genes_sort %in%  tsgs] <- set1[2]
    names(driver_color) <- genes_sort
    
    tab4plot$GENE <- factor(tab4plot$GENE, levels = genes_sort)
    
    tab4plot_driver <- tab4plot[(tab4plot$GENE %in% driver_genes$Gene) | tab4plot$cancer == "OV",]
    tab4plot_driver$cancer <- factor(tab4plot_driver$cancer, levels = cancers_sort)
    
    p <- ggplot(tab4plot_driver)
    p <- p + geom_point(data=tab4plot_driver, aes(x = GENE, y = meddiff,size = log10_pvalue, alpha = alpha), 
                       stroke = 0, color = "#6A3D9A", shape = 16)
    p <- p + geom_errorbar(data = tab4plot_driver[tab4plot_driver$top_effect_size,], 
                           aes(x = GENE, ymin=meddiff_bottom, ymax=meddiff_upper), 
                           width=.2, position=position_dodge(0.05), alpha = 0.1)
    
    if (self == "trans") {
      p <- p + geom_text_repel(data = tab4plot_driver[tab4plot_driver$top_effect_size_path,], 
                               aes(x = GENE, y = meddiff, label = text2p, color = SUB_GENE.path), 
                               size = 2.5, force = 2, alpha = 0.8)
    } else {
      p <- p + geom_text_repel(data = tab4plot_driver[tab4plot_driver$top_effect_size,], 
                               aes(x = GENE, y = meddiff, label = text2p, color = SUB_GENE.path), 
                               size = 2, force = 2, alpha = 0.8)
    }
    p <- p + facet_grid(.~cancer, scales = "free", space = "free")
    p <- p + theme(panel.background = element_rect(fill = brewer.pal(9, "Greys")[2], size = 0))
    p <- p + theme(strip.background.x = element_rect(color = "black"))
    p <- p + theme(panel.spacing.x = unit(0, "lines"), strip.text.x = element_text(size = 10))
    p <- p + geom_hline(yintercept = 0, linetype = 2, color = "grey")
    p <- p + scale_color_manual(values = c(color_cancers2, tcga_path_colors))
    p <- p + xlab('enzyme')+ylab("(median substrate phosphorylation difference)")
    p <- p + theme(axis.title.y = element_text(size=8))
    p <- p + theme(axis.title.x = element_blank())
    p <- p + theme(axis.text.y = element_text(colour="black", size=10))
    p <- p + theme(axis.text.x = element_text(size = 15, angle = 90, hjust = 0.5, vjust = 0.5, face = enzyme_face[genes_sort %in% as.vector(driver_genes$Gene)]))
    p <- p + guides(size = FALSE, alpha = FALSE)
    # p <- p + theme(axis.text.x = element_text(size = 12, angle = 90, hjust = 0.5, vjust = 0.5, colour = driver_color[genes_sort %in% as.vector(driver_genes$Gene)], face = enzyme_face[genes_sort %in% as.vector(driver_genes$Gene)]))
    p
    fn = paste(makeOutDir(resultD = resultD), 'driver_enzyme_', genoalt_type, '_affect_', self, '_phosphosites_p_value', sig_thres,'_phosphosites_scatterplot.pdf',sep ="")
    ggsave(filename = fn, width = ifelse(self == "cis", 4, 14), height = ifelse(self == "cis", 3, 6.5), useDingbats = F)
  }
}


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
reg_nonNA2test <- c(15,20)

top_kinase2show <- 10
top_substrate2show <- 2
pdf_sizes <- list(kinase = list(trans = list(BRCA = list(width = 7, height = 10),
                                             OV = list(width = 5, height = 8.5),
                                             UCEC = list(width = 5, height = 4),
                                             CO = list(width = 6, height = 10)),
                                cis = list(BRCA = list(width = 6, height = 4),
                                           OV = list(width = 4, height = 4),
                                           UCEC = list(width = 5, height = 4),
                                           CO = list(width = 4, height = 4))),
                  phosphatase = list(trans = list(BRCA = list(width = 5, height = 4.5),
                                                  OV = list(width = 3, height = 3),
                                                  UCEC = list(width = 4, height = 4),
                                                  CO = list(width = 3, height = 3))))
# cptac_phase2process <- "cptac3"
# cancers2process <- "UCEC"
cptac_phase2process <- "cptac2p_3can"
cancers2process <- cancers_sort
fdr_thres2process <- c(0.05)
pdf_sizes_by_facetx <- list("1" = 5, "3" = 12)

# inputs -------------------------------------------------------------------
## input driver genes
# driver_genes <- read.delim(file = "./TCGA_data/reference_files/Consensus.Genelist.full.txt")
driver_pancan <- loadGeneList(gene_type = "driver", cancer = "PANCAN", is.soft.limit = "")

## limit phosphatase implicated in each cancer
pp_cans <- list()
pp_cans[["CO"]] <- c("PTPRF", "PTPRG", "PTPRT", "PTPN3", "PTPN13", "PTPN14", "PTPRA", "PTPRS", "PTPN5", "PTPN21", "PTPN23", "PTPN1", "PTPRD", "PTPRH", "DUSP1", "PTPRJ")
pp_cans[["OV"]] <- c("PTEN", "PPP2CA", "CIP2A", "DUSP1", "PTPN13", "CDC25A", "CDC25B", "PTPN6", "PPM1D", "PP2C", "PTP4A3", "PPAP2A")
pp_cans[["BRCA"]] <- c("PTPN1", "PTPN13", "PTEN", "PTPN12", "PTPN9", "DUSP1", "DUSP3", "DUSP4", "DUSP5", "DUSP23", "CDC25A", "CDC25B", "CDC25C", "PPP1CA", "PPP1CB", "PPP1CC", "PPM1D", "EYA2")
pp_cans[["other"]] <- c("PPP3CB", "PPP3CA", "PPP3CC", "PTPN11")
# plot SMGs/driver as kinase/phosphatase--------------------------------------------------------------------
for (enzyme_type in c("kinase")) {
  enzyme_type <- "kinase"
  subdir1 <- paste0(makeOutDir(resultD = resultD), enzyme_type, "/")
  dir.create(subdir1)
  for (reg_nonNA in reg_nonNA2test) {
    # reg_nonNA <- 20
    subdir2 <- paste0(subdir1, "reg_nonNA", reg_nonNA, "/")
    dir.create(subdir2)
    
    sup_cans_tab_en <- fread(input = paste0(ppnD, "regression/tables/", "change_regression_nonNA/", 
                                            enzyme_type, "_substrate_regression_", cptac_phase2process, "_tumor",
                                            "_reg_nonNA", reg_nonNA, ".txt"), data.table = F)
    for (fdr_thres in fdr_thres2process) {
      # fdr_thres <- 0.05
      subdir3 <- paste0(subdir2, "fdr_thres", fdr_thres, "/")
      dir.create(subdir3)
      
      sup_cans_tab_en <- markSigKS(regression = sup_cans_tab_en, sig_thres = fdr_thres, enzyme_type = enzyme_type)
      sup_cans_tab_en$regulated <- (sup_cans_tab_en$fdr_sig & sup_cans_tab_en$coef_sig)
      sup_cans_tab_en_reg <- sup_cans_tab_en[sup_cans_tab_en$regulated,]
      
      for (self in c("trans")) {
        subdir4 <- paste0(subdir3, self, "/")
        dir.create(subdir4)
        
        tab4plot_cans <- NULL
        for (cancer in cancers2process) {
          subdir5 <- paste0(subdir4, cancer, "/")
          dir.create(subdir5)
          driver_can <- loadGeneList(gene_type = "driver", cancer = cancer, is.soft.limit = "")
          genes2show <- unique(c(driver_can, SMGs[[cancer]]))
          
          tab_can <- sup_cans_tab_en_reg[sup_cans_tab_en_reg$SELF == self & sup_cans_tab_en_reg$Cancer == cancer,]
          if (nrow(tab_can) > 0) {
            if (self == "cis") {
              # tab4plot <- unique(tab_can[(tab_can$GENE %in% c(driver_can, top_kinases)), c("GENE", "SUB_GENE", "SUB_MOD_RSD","Cancer", "FDR_pro_kin", "coef_pro_kin")])
              tab4plot <- unique(tab_can[, c("GENE", "SUB_GENE", "SUB_MOD_RSD","Cancer", "FDR_pro_kin", "coef_pro_kin")])
            }
            if (self == "trans") {
              ## filter cis by showing only the SMGs
              tab4plot <- unique(tab_can[(tab_can$GENE %in% genes2show | tab_can$SUB_GENE %in% genes2show), c("GENE", "SUB_GENE", "SUB_MOD_RSD","Cancer", "FDR_pho_kin", "coef_pho_kin")])
              
              if (enzyme_type == "phosphatase") {
                tab4plot <- unique(tab_can[, c("GENE", "SUB_GENE", "SUB_MOD_RSD","Cancer", "FDR_pho_kin", "coef_pho_kin")])
              } 
            }
            colnames(tab4plot) <- c("GENE", "SUB_GENE", "SUB_MOD_RSD", "Cancer", "FDR", "coef")
            
            ## annotate driver genes to TCGA pathways
            genes2pathways <- map2TCGApathwaways(gene_list = unique(as.vector(tab4plot$GENE)), pathway_list = tcga_pathways)
            genes2pathways <- data.frame(GENE = rep(x = names(genes2pathways), sapply(X = genes2pathways, FUN = function(x) length(x))), 
                                         GENE.path = unlist(genes2pathways, use.names = F))
            tab4plot <- merge(tab4plot, genes2pathways, all.x = T)
            
            ## annoate substrate genes to TCGA pathways
            if (self == "trans" & cancer == "BRCA" & enzyme_type == "kinase") {
              sub_genes2pathways <- map2TCGApathwaways(gene_list = unique(tab4plot$SUB_GENE), pathway_list = tcga_pathways)
            } else {
              sub_genes2pathways <- map2TCGApathwaways(gene_list = unique(tab4plot$SUB_GENE), pathway_list = tcga_pathways_pluskegg)
            }
            sub_genes2pathways <- data.frame(SUB_GENE = rep(x = names(sub_genes2pathways), sapply(X = sub_genes2pathways, FUN = function(x) length(x))), 
                                             SUB_GENE.path = unlist(sub_genes2pathways, use.names = F))
            tab4plot <- merge(tab4plot, sub_genes2pathways, all.x = T)
            tmp <- as.vector(tab4plot$SUB_GENE.path)
            tmp[is.na(tmp)] <- "NA"
            tmp[tmp == "NA" & tmp %in% c(driver_can, SMGs[[cancer]])] <- "driver"
            tab4plot$SUB_GENE.path <- tmp
            
            if (enzyme_type == "kinase") {
              tab4plot <- tab4plot[(tab4plot$GENE %in% tab4plot$GENE[!is.na(tab4plot$SUB_GENE.path)]) | (tab4plot$GENE %in% driver_can),]
            }
            
            # add more annotation columns to table ------------------------------------
            tab4plot$same_path <- (as.vector(tab4plot$GENE.path) == as.vector(tab4plot$SUB_GENE.path))
            tab4plot$pair_pro <- paste0(tab4plot$GENE, ":", tab4plot$SUB_GENE)
            tab4plot$pair <- paste0(tab4plot$pair_pro, "_", tab4plot$SUB_MOD_RSD)
            tab4plot <- tab4plot[order(tab4plot$pair, tab4plot$same_path, decreasing = T),]
            tab4plot <- tab4plot[!duplicated(tab4plot$pair),]
            
            tab4plot$log10FDR <- -log10(tab4plot$FDR)
            tab4plot$site <- paste0(tab4plot$SUB_GENE, "\n", tab4plot$SUB_MOD_RSD)
            tab4plot$x <- paste0(tab4plot$GENE, "-", tab4plot$Cancer)
            tab4plot$site_coef <- paste0(tab4plot$site, "\n(coef=", signif(tab4plot$coef, digits = 2 ), ")")
            tab4plot$text <- as.vector(tab4plot$site_coef)
            if (self == "trans" & cptac_phase2process == "cptac2p_3can") {
              facet_y <- as.vector(tab4plot$SUB_GENE.path)
              tab4plot$facet_y <- facet_y
              tab4plot$facet_y <- factor(tab4plot$facet_y, levels = c(pathway_names, "NA"))
            }
            
            
            ## sort x-axis genes
            genes_sort <- as.vector(unique(tab4plot$GENE))
            genes_sort <- c(genes_sort[(genes_sort %in% driver_can)], genes_sort[!(genes_sort %in% driver_can)])
            
            ## sort gene on x-axis again
            genes2pathways <- genes2pathways[genes2pathways$GENE %in% tab4plot$GENE,]
            genes_sort2 <- NULL
            for (path in unique(tab4plot$GENE.path[!is.na(tab4plot$GENE.path)])) {
              path_genes <- as.vector(unique(genes2pathways$GENE[genes2pathways$GENE.path == path & !is.na(genes2pathways$GENE.path)]))
              genes_sort2 <- c(genes_sort2, path_genes)
            }
            genes_sort2 <- c(genes_sort2, genes_sort[!(genes_sort %in% genes_sort2)])
            
            # assign color for x-axis genes
            color_x_axis <- NULL
            for (gene in genes_sort2) {
              gene_path <- unique(as.vector(tab4plot$GENE.path[tab4plot$GENE == gene]))
              gene_path <- gene_path[!is.na(gene_path)]
              if (length(gene_path) > 0){
                color_x_axis <- c(color_x_axis, tcga_path_colors[gene_path])
              } else {
                color_x_axis <- c(color_x_axis, "black")
              }
            }
            names(color_x_axis) <- genes_sort2
            
            ## cap y-axis effect size for better visualization
            if (enzyme_type == "kinase") {
              lim = quantile(x = tab4plot$coef, probs = 0.99)
              cap <- min(lim, 2)
              tab4plot$coef_capped <- tab4plot$coef
              tab4plot$coef_capped[tab4plot$coef > cap] <- cap
            } else {
              lim = quantile(x = tab4plot$coef, probs = 0.01)
              cap <- max(lim, -2)
              tab4plot$coef_capped <- tab4plot$coef
              tab4plot$coef_capped[tab4plot$coef < cap] <- cap
            }
            tab4plot$text[abs(tab4plot$coef) < abs(cap)] <- as.vector(tab4plot$site[abs(tab4plot$coef) < abs(cap)])
            
            ## annotate the substrates with top effect size
            tab4plot <- tab4plot[order(tab4plot$pair_pro, tab4plot$coef, decreasing = T),]
            top_effect_size <- vector(mode = "logical", length = nrow(tab4plot))
            for (gene in unique(tab4plot$GENE)) {
              for (cancer in unique(tab4plot$Cancer)) {
                tab_gene <- tab4plot[tab4plot$GENE == gene & tab4plot$Cancer == cancer,]
                abs_coefs <- abs(tab_gene$coef)
                abs_coefs <- abs_coefs[order(abs_coefs, decreasing = T)]
                top_coefs <- abs_coefs[1:min(length(abs_coefs), top_substrate2show)]
                top_effect_size[abs(tab4plot$coef) %in% top_coefs] <- TRUE
              }
            }
            tab4plot$top_effect_size <- top_effect_size
            
            
            # annotate the substrate with top effect size by pathway ------------------
            top_effect_size_path <- vector(mode = "logical", length = nrow(tab4plot))
            for (gene in unique(tab4plot$GENE)) {
              tab_gene <- tab4plot[tab4plot$GENE == gene,]
              for (cancer in unique(tab4plot$Cancer)) {
                tab_can <- tab_gene[tab_gene$Cancer == cancer,]
                for (pathway_name in unique(tab_gene$SUB_GENE.path)) {
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
            if (self == "trans")  {
              tab4plot$text2p <- as.vector(tab4plot$text)
              tab4plot$text2p[tab4plot$SUB_GENE.path == "NA" & (tab4plot$GENE %in% as.vector(tab4plot$GENE[tab4plot$SUB_GENE.path != "NA"]))] <- ""
              tab4plot$text2p[tab4plot$SUB_GENE.path == "NA" & !(tab4plot$GENE %in% as.vector(tab4plot$GENE[tab4plot$SUB_GENE.path != "NA"])) & !tab4plot$top_effect_size] <- ""
            } else {
              tab4plot$text2p <- as.vector(tab4plot$text)
              tab4plot$text2p[!tab4plot$top_effect_size] <- ""
            }
            
            tab4plot$GENE <- factor(tab4plot$GENE, levels = genes_sort2)
            tab4plot$Cancer <- factor(tab4plot$Cancer, levels = cancers2process)
            # p <- ggplot()
            # if (enzyme_type == "kinase") {
            #   p <- p + geom_point(data=tab4plot, aes(x = GENE, y = coef_capped, size = log10FDR, alpha = ifelse(!is.na(SUB_GENE.path), 0.4, 0.3)), 
            #                       shape = 16, stroke = 0, color = color_cancers2[cancer])
            #   p <- p + geom_text_repel(data = tab4plot, mapping = aes(x = GENE, y = coef_capped, label = text2p, color = SUB_GENE.path, alpha = ifelse(top_effect_size_path, 0.8, 0.2)), 
            #                            size = 2, force = 1)
            #   if (self == "trans" & cptac_phase2process == "cptac2p_3can" & enzyme_type == "kinase") {
            #     p <- p + facet_grid(facet_y~., scales = "free_x", space = "fixed")
            #     p <- p + theme(panel.spacing.y = unit(0, "lines"))
            #   }
            # } else {
            #   p <- p + geom_point(data=tab4plot, aes(x = GENE, y = coef_capped, size = log10FDR), 
            #                       shape = 16, stroke = 0, color = color_cancers2[cancer], alpha = 0.8)
            #   p <- p + geom_text_repel(data = tab4plot, mapping = aes(x = GENE, y = coef_capped, label = text2p, color = SUB_GENE.path, alpha = ifelse(top_effect_size_path, 0.8, 0.2)), 
            #                            size = 2, force = 1, alpha = 0.8)
            #   if (self == "trans" & cptac_phase2process == "cptac2p_3can" & enzyme_type == "kinase") {
            #     p <- p + facet_grid(facet_y~., scales = "free_x", space = "fixed")
            #     p <- p + theme(panel.spacing.y = unit(0, "lines"))
            #   }
            # }
            # p <- p + geom_hline(yintercept = 0, linetype = 2, color = "grey")
            # p <- p + scale_color_manual(values = tcga_path_colors)
            # p <- p + theme_minimal() + theme(panel.background = element_rect(fill = brewer.pal(9, "Greys")[2], size = 0))
            # p <- p + xlab('kinase')+ylab("effect size")
            # p <- p + theme(axis.text.x = element_text(angle = 90, face = "bold", hjust = 0.5, vjust = 0.5, size = 10))
            # p <- p + theme(axis.title=element_text(size=10))
            # p <- p + theme(axis.text.y = element_text(colour="black", size=10))
            # p <- p + theme(title = element_text(size = 6, face = "italic"))
            # p <- p + theme(panel.spacing.x = unit(0, "lines"), strip.text.x = element_text(size = 9, face = "bold"))
            # p <- p + theme(strip.background.x = element_rect(size = 1, color = "black", linetype = "solid"))
            # p
            # fn = paste0(subdir5, cancer, 'substrates_of_', self,'_regulated_driver_gene_', enzyme_type, '_pairs_in_TCGA_pathways_reg_nonNA',reg_nonNA,'.png')
            # ggsave(file = fn, height=pdf_sizes[[enzyme_type]][[self]][[cancer]]$height, width = pdf_sizes[[enzyme_type]][[self]][[cancer]]$width, device = png())
            # dev.off()
            tab4plot_cans <- rbind(tab4plot_cans, tab4plot)
          }
        }
        
        tab4plot <- tab4plot_cans
        tab4plot <- tab4plot[!is.na(tab4plot$SUB_GENE.path) | (is.na(tab4plot$SUB_GENE.path) & tab4plot$SUB_GENE %in% c(driver_can, SMGs[[cancer]])),] 
        tab4plot <- tab4plot[tab4plot$Cancer %in% c("BRCA", "CO"),]
        tab4plot$cancer_full_name <- cancer_full_names[as.vector(tab4plot$Cancer)]
        
        p <- ggplot()
        if (enzyme_type == "kinase") {
          p <- p + geom_point(data=tab4plot, aes(x = GENE, y = coef_capped, size = log10FDR, alpha = ifelse(!is.na(SUB_GENE.path), 0.4, 0.3)), 
                              shape = 16, stroke = 0, color = color_cancers2[cancer])
          p <- p + geom_text_repel(data = tab4plot, mapping = aes(x = GENE, y = coef_capped, label = text2p, color = SUB_GENE.path, alpha = ifelse(top_effect_size_path, 0.8, 0.2)), 
                                   size = 2, force = 1)
        } else {
          p <- p + geom_point(data=tab4plot, aes(x = GENE, y = coef_capped, size = log10FDR), 
                              shape = 16, stroke = 0, color = color_cancers2[cancer], alpha = 0.8)
          p <- p + geom_text_repel(data = tab4plot, mapping = aes(x = GENE, y = coef_capped, label = text2p, color = SUB_GENE.path, alpha = ifelse(top_effect_size_path, 0.8, 0.2)), 
                                   size = 2, force = 1)
        }
        if (self == "trans" & cptac_phase2process == "cptac2p_3can") {
          p <- p + facet_grid(facet_y~cancer_full_name, scales = "free_x", space = "free_x")
          p <- p + theme(panel.spacing.y = unit(0, "lines"))
        }
        p <- p + scale_color_manual(values = tcga_path_colors)
        p <- p + xlab(enzyme_type)+ylab("beta coefficient")
        p <- p + theme_bw()
        p <- p + theme(axis.text.x = element_text(angle = 90, face = "bold", hjust = 0.5, vjust = 0.5, size = 10))
        p <- p + theme(axis.title=element_text(size=10))
        p <- p + theme(axis.text.y = element_text(colour="black", size=10))
        p <- p + theme(title = element_text(size = 6, face = "italic"))
        p <- p + theme(panel.spacing.x = unit(0, "lines"), strip.text.x = element_text(size = 9, face = "bold"))
        p <- p + theme(panel.spacing.y = unit(0, "lines"), strip.text.y = element_text(size = 9, face = "bold"))
        p <- p + theme(strip.background.x = element_rect(size = 1, color = "black", linetype = "solid"))
        p <- p + guides(alpha = F)
        p
        fn = paste0(subdir4, 'substrates_of_', self,'_regulated_driver_gene_', enzyme_type, '_pairs_in_TCGA_pathways_reg_nonNA',reg_nonNA, '_fdr_thres', fdr_thres, '_BRCA_CO_only_faceted.png')
        ggsave(filename = fn, height=8, width = 10, device = png())
        dev.off()
        
        p <- ggplot()
        if (enzyme_type == "kinase") {
          p <- p + geom_point(data=tab4plot, aes(x = GENE, y = coef_capped, size = log10FDR, alpha = ifelse(!is.na(SUB_GENE.path), 0.4, 0.3)), 
                              shape = 16, stroke = 0, color = color_cancers2[cancer])
          p <- p + geom_text_repel(data = tab4plot, mapping = aes(x = GENE, y = coef_capped, label = text2p, color = SUB_GENE.path, alpha = ifelse(top_effect_size_path, 0.8, 0.2)), 
                                   size = 2, force = 1)
        } else {
          p <- p + geom_point(data=tab4plot, aes(x = GENE, y = coef_capped, size = log10FDR), 
                              shape = 16, stroke = 0, color = color_cancers2[cancer], alpha = 0.8)
          p <- p + geom_text_repel(data = tab4plot, mapping = aes(x = GENE, y = coef_capped, label = text2p, color = SUB_GENE.path, alpha = ifelse(top_effect_size_path, 0.8, 0.2)), 
                                   size = 2, force = 1)
        }
        if (self == "trans" & cptac_phase2process == "cptac2p_3can") {
          p <- p + facet_grid(.~cancer_full_name, scales = "free_x", space = "free_x")
          p <- p + theme(panel.spacing.x = unit(0, "lines"))
        }
        p <- p + scale_color_manual(values = tcga_path_colors)
        p <- p + xlab(enzyme_type)+ylab("beta coefficient")
        p <- p + theme_bw()
        p <- p + theme(axis.text.x = element_text(angle = 90, face = "bold", hjust = 0.5, vjust = 0.5, size = 10))
        p <- p + theme(axis.title=element_text(size=10))
        p <- p + theme(axis.text.y = element_text(colour="black", size=10))
        p <- p + theme(title = element_text(size = 6, face = "italic"))
        p <- p + theme(panel.spacing.x = unit(0, "lines"), strip.text.x = element_text(size = 6, face = "bold"))
        p <- p + theme(panel.spacing.y = unit(0, "lines"), strip.text.y = element_text(size = 9, face = "bold"))
        p <- p + theme(strip.background.x = element_rect(size = 1, color = "black", linetype = "solid"))
        p <- p + guides(alpha = F, guide_legend(title="New Legend Title"))
        p
        fn = paste0(subdir4, 'substrates_of_', self,'_regulated_driver_gene_', enzyme_type, '_pairs_in_TCGA_pathways_reg_nonNA',reg_nonNA, '_fdr_thres', fdr_thres, '_BRCA_CO_only.png')
        ggsave(filename = fn, height=4, width = 6, device = png())
        dev.off()
      }
    }
  }
}


# Phosphatase curated -----------------------------------------------------
for (enzyme_type in c("phosphatase")) {
  subdir1 <- paste0(makeOutDir(resultD = resultD), enzyme_type, "/")
  dir.create(subdir1)
  for (reg_nonNA in reg_nonNA2test) {
    subdir2 <- paste0(subdir1, "reg_nonNA", reg_nonNA, "/")
    dir.create(subdir2)
    
    sup_cans_tab_en <- fread(input = paste0(ppnD, "regression/tables/", "change_regression_nonNA/", 
                                            enzyme_type, "_substrate_regression_", cptac_phase2process, "_tumor",
                                            "_reg_nonNA", reg_nonNA, ".txt"), data.table = F)
    for (fdr_thres in fdr_thres2process) {
      subdir3 <- paste0(subdir2, "fdr_thres", fdr_thres, "/")
      dir.create(subdir3)
      
      sup_cans_tab_en <- markSigKS(regression = sup_cans_tab_en, sig_thres = fdr_thres, enzyme_type = enzyme_type)
      sup_cans_tab_en$regulated <- (sup_cans_tab_en$fdr_sig & sup_cans_tab_en$coef_sig)
      sup_cans_tab_en_reg <- sup_cans_tab_en[sup_cans_tab_en$regulated,]
      
      for (self in c("trans")) {
        subdir4 <- paste0(subdir3, self, "/")
        dir.create(subdir4)
        
        tab4plot_cans <- NULL
        for (cancer in cancers2process) {
          subdir5 <- paste0(subdir4, cancer, "/")
          dir.create(subdir5)
          driver_can <- loadGeneList(gene_type = "driver", cancer = cancer, is.soft.limit = "")
          
          tab_can <- sup_cans_tab_en_reg[sup_cans_tab_en_reg$SELF == self & sup_cans_tab_en_reg$Cancer == cancer,]
          if (nrow(tab_can) > 0) {
            if (self == "cis") {
              # tab4plot <- unique(tab_can[(tab_can$GENE %in% c(driver_can, top_kinases)), c("GENE", "SUB_GENE", "SUB_MOD_RSD","Cancer", "FDR_pro_kin", "coef_pro_kin")])
              tab4plot <- unique(tab_can[, c("GENE", "SUB_GENE", "SUB_MOD_RSD","Cancer", "FDR_pro_kin", "coef_pro_kin")])
            }
            if (self == "trans") {
              ## filter cis by showing only the SMGs
              tab4plot <- unique(tab_can[(tab_can$GENE %in% c(driver_can, driver_pancan)), c("GENE", "SUB_GENE", "SUB_MOD_RSD","Cancer", "FDR_pho_kin", "coef_pho_kin")])
              
              if (enzyme_type == "phosphatase") {
                tab4plot <- unique(tab_can[, c("GENE", "SUB_GENE", "SUB_MOD_RSD","Cancer", "FDR_pho_kin", "coef_pho_kin")])
              } 
            }
            colnames(tab4plot) <- c("GENE", "SUB_GENE", "SUB_MOD_RSD", "Cancer", "FDR", "coef")
            
            ## annotate driver genes to TCGA pathways
            genes2pathways <- map2TCGApathwaways(gene_list = unique(as.vector(tab4plot$GENE)), pathway_list = tcga_pathways)
            genes2pathways <- data.frame(GENE = rep(x = names(genes2pathways), sapply(X = genes2pathways, FUN = function(x) length(x))), 
                                         GENE.path = unlist(genes2pathways, use.names = F))
            tab4plot <- merge(tab4plot, genes2pathways, all.x = T)
            
            ## annoate substrate genes to TCGA pathways
            if (self == "trans" & cancer == "BRCA" & enzyme_type == "kinase") {
              sub_genes2pathways <- map2TCGApathwaways(gene_list = unique(tab4plot$SUB_GENE), pathway_list = tcga_pathways)
            } else {
              sub_genes2pathways <- map2TCGApathwaways(gene_list = unique(tab4plot$SUB_GENE), pathway_list = tcga_pathways_pluskegg)
            }
            sub_genes2pathways <- data.frame(SUB_GENE = rep(x = names(sub_genes2pathways), sapply(X = sub_genes2pathways, FUN = function(x) length(x))), 
                                             SUB_GENE.path = unlist(sub_genes2pathways, use.names = F))
            tab4plot <- merge(tab4plot, sub_genes2pathways, all.x = T)
            tmp <- as.vector(tab4plot$SUB_GENE.path)
            tmp[is.na(tmp)] <- "NA"
            tab4plot$SUB_GENE.path <- tmp
            
            if (enzyme_type == "kinase") {
              tab4plot <- tab4plot[(tab4plot$GENE %in% tab4plot$GENE[!is.na(tab4plot$SUB_GENE.path)]) | (tab4plot$GENE %in% driver_can),]
            }
            
            # add more annotation columns to table ------------------------------------
            tab4plot$same_path <- (as.vector(tab4plot$GENE.path) == as.vector(tab4plot$SUB_GENE.path))
            tab4plot$pair_pro <- paste0(tab4plot$GENE, ":", tab4plot$SUB_GENE)
            tab4plot$pair <- paste0(tab4plot$pair_pro, "_", tab4plot$SUB_MOD_RSD)
            tab4plot <- tab4plot[order(tab4plot$pair, tab4plot$same_path, decreasing = T),]
            tab4plot <- tab4plot[!duplicated(tab4plot$pair),]
            
            tab4plot$log10FDR <- -log10(tab4plot$FDR)
            tab4plot$site <- paste0(tab4plot$SUB_GENE, "\n", tab4plot$SUB_MOD_RSD)
            tab4plot$x <- paste0(tab4plot$GENE, "-", tab4plot$Cancer)
            tab4plot$site_coef <- paste0(tab4plot$site, "\n(coef=", signif(tab4plot$coef, digits = 2 ), ")")
            tab4plot$text <- as.vector(tab4plot$site_coef)
            if (self == "trans" & cptac_phase2process == "cptac2p_3can") {
              facet_y <- as.vector(tab4plot$SUB_GENE.path)
              tab4plot$facet_y <- facet_y
              tab4plot$facet_y <- factor(tab4plot$facet_y, levels = c(pathway_names, "NA"))
            }
            
            
            ## sort x-axis genes
            genes_sort <- as.vector(unique(tab4plot$GENE))
            genes_sort <- c(genes_sort[(genes_sort %in% driver_can)], genes_sort[!(genes_sort %in% driver_can)])
            
            ## sort gene on x-axis again
            genes2pathways <- genes2pathways[genes2pathways$GENE %in% tab4plot$GENE,]
            genes_sort2 <- NULL
            for (path in unique(tab4plot$GENE.path[!is.na(tab4plot$GENE.path)])) {
              path_genes <- as.vector(unique(genes2pathways$GENE[genes2pathways$GENE.path == path & !is.na(genes2pathways$GENE.path)]))
              genes_sort2 <- c(genes_sort2, path_genes)
            }
            genes_sort2 <- c(genes_sort2, genes_sort[!(genes_sort %in% genes_sort2)])
            
            # assign color for x-axis genes
            color_x_axis <- NULL
            for (gene in genes_sort2) {
              gene_path <- unique(as.vector(tab4plot$GENE.path[tab4plot$GENE == gene]))
              gene_path <- gene_path[!is.na(gene_path)]
              if (length(gene_path) > 0){
                color_x_axis <- c(color_x_axis, tcga_path_colors[gene_path])
              } else {
                color_x_axis <- c(color_x_axis, "black")
              }
            }
            names(color_x_axis) <- genes_sort2
            
            ## cap y-axis effect size for better visualization
            if (enzyme_type == "kinase") {
              lim = quantile(x = tab4plot$coef, probs = 0.99)
              cap <- min(lim, 2)
              tab4plot$coef_capped <- tab4plot$coef
              tab4plot$coef_capped[tab4plot$coef > cap] <- cap
            } else {
              lim = quantile(x = tab4plot$coef, probs = 0.01)
              cap <- max(lim, -2)
              tab4plot$coef_capped <- tab4plot$coef
              tab4plot$coef_capped[tab4plot$coef < cap] <- cap
            }
            tab4plot$text[abs(tab4plot$coef) < abs(cap)] <- as.vector(tab4plot$site[abs(tab4plot$coef) < abs(cap)])
            
            ## annotate the substrates with top effect size
            tab4plot <- tab4plot[order(tab4plot$pair_pro, tab4plot$coef, decreasing = T),]
            top_effect_size <- vector(mode = "logical", length = nrow(tab4plot))
            for (gene in unique(tab4plot$GENE)) {
              for (cancer in unique(tab4plot$Cancer)) {
                tab_gene <- tab4plot[tab4plot$GENE == gene & tab4plot$Cancer == cancer,]
                abs_coefs <- abs(tab_gene$coef)
                abs_coefs <- abs_coefs[order(abs_coefs, decreasing = T)]
                top_coefs <- abs_coefs[1:min(length(abs_coefs), top_substrate2show)]
                top_effect_size[abs(tab4plot$coef) %in% top_coefs] <- TRUE
              }
            }
            tab4plot$top_effect_size <- top_effect_size
            
            
            # annotate the substrate with top effect size by pathway ------------------
            top_effect_size_path <- vector(mode = "logical", length = nrow(tab4plot))
            for (gene in unique(tab4plot$GENE)) {
              tab_gene <- tab4plot[tab4plot$GENE == gene,]
              for (cancer in unique(tab4plot$Cancer)) {
                tab_can <- tab_gene[tab_gene$Cancer == cancer,]
                for (pathway_name in unique(tab_gene$SUB_GENE.path)) {
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
            if (self == "trans")  {
              tab4plot$text2p <- as.vector(tab4plot$text)
              tab4plot$text2p[tab4plot$SUB_GENE.path == "NA" & (tab4plot$GENE %in% as.vector(tab4plot$GENE[tab4plot$SUB_GENE.path != "NA"]))] <- ""
              tab4plot$text2p[tab4plot$SUB_GENE.path == "NA" & !(tab4plot$GENE %in% as.vector(tab4plot$GENE[tab4plot$SUB_GENE.path != "NA"])) & !tab4plot$top_effect_size] <- ""
            } else {
              tab4plot$text2p <- as.vector(tab4plot$text)
              tab4plot$text2p[!tab4plot$top_effect_size] <- ""
            }
            
            tab4plot$GENE <- factor(tab4plot$GENE, levels = genes_sort2)
            tab4plot$Cancer <- factor(tab4plot$Cancer, levels = cancers2process)
            p <- ggplot()
            if (enzyme_type == "kinase") {
              p <- p + geom_point(data=tab4plot, aes(x = GENE, y = coef_capped, size = log10FDR, alpha = ifelse(!is.na(SUB_GENE.path), 0.4, 0.3)), 
                                  shape = 16, stroke = 0, color = color_cancers2[cancer])
              p <- p + geom_text_repel(data = tab4plot, mapping = aes(x = GENE, y = coef_capped, label = text2p, color = SUB_GENE.path, alpha = ifelse(top_effect_size_path, 0.8, 0.2)), 
                                       size = 2, force = 1)
              if (self == "trans" & cptac_phase2process == "cptac2p_3can" & enzyme_type == "kinase") {
                p <- p + facet_grid(facet_y~., scales = "free_x", space = "fixed")
                p <- p + theme(panel.spacing.y = unit(0, "lines"))
              }
            } else {
              p <- p + geom_point(data=tab4plot, aes(x = GENE, y = coef_capped, size = log10FDR), 
                                  shape = 16, stroke = 0, color = color_cancers2[cancer], alpha = 0.8)
              p <- p + geom_text_repel(data = tab4plot, mapping = aes(x = GENE, y = coef_capped, label = text2p, color = SUB_GENE.path, alpha = ifelse(top_effect_size_path, 0.8, 0.2)), 
                                       size = 2, force = 1, alpha = 0.8)
              if (self == "trans" & cptac_phase2process == "cptac2p_3can" & enzyme_type == "kinase") {
                p <- p + facet_grid(facet_y~., scales = "free_x", space = "fixed")
                p <- p + theme(panel.spacing.y = unit(0, "lines"))
              }
            }
            p <- p + geom_hline(yintercept = 0, linetype = 2, color = "grey")
            p <- p + scale_color_manual(values = tcga_path_colors)
            p <- p + theme_minimal() + theme(panel.background = element_rect(fill = brewer.pal(9, "Greys")[2], size = 0))
            p <- p + xlab('kinase')+ylab("effect size")
            p <- p + theme(axis.text.x = element_text(angle = 90, face = "bold", hjust = 0.5, vjust = 0.5, size = 10, color = color_x_axis))
            p <- p + theme(axis.title=element_text(size=10))
            p <- p + theme(axis.text.y = element_text(colour="black", size=10))
            p <- p + theme(title = element_text(size = 6, face = "italic"))
            p <- p + theme(panel.spacing.x = unit(0, "lines"), strip.text.x = element_text(size = 9, face = "bold"))
            p <- p + theme(strip.background.x = element_rect(size = 1, color = "black", linetype = "solid"))
            # fn = paste0(subdir5, cancer, 'substrates_of_', self,'_regulated_driver_gene_', enzyme_type, '_pairs_in_TCGA_pathways_reg_nonNA',reg_nonNA,'.pdf')
            # pdf(file = fn, height=pdf_sizes[[enzyme_type]][[self]][[cancer]]$height, width = pdf_sizes[[enzyme_type]][[self]][[cancer]]$width, useDingbats = F)
            # print(p)
            # dev.off()
            
            tab4plot_cans <- rbind(tab4plot_cans, tab4plot)
          }
        }
        
        tab4plot <- tab4plot_cans
        tab4plot <- tab4plot[tab4plot$GENE %in% unlist(pp_cans),]
        tab4plot <- tab4plot[!(duplicated(tab4plot$FDR) & duplicated(tab4plot$coef)),]
        tab4plot <- tab4plot[tab4plot$Cancer %in% c("BRCA", "CO"),]
        p <- ggplot()
        if (enzyme_type == "kinase") {
          p <- p + geom_point(data=tab4plot, aes(x = GENE, y = coef_capped, size = log10FDR, alpha = ifelse(!is.na(SUB_GENE.path), 0.4, 0.3)), 
                              shape = 16, stroke = 0, color = color_cancers2[cancer])
          p <- p + geom_text_repel(data = tab4plot, mapping = aes(x = GENE, y = coef_capped, label = text2p, color = SUB_GENE.path, alpha = ifelse(top_effect_size_path, 0.8, 0.2)), 
                                   size = 2, force = 1)
        } else {
          p <- p + geom_point(data=tab4plot, aes(x = GENE, y = coef_capped, size = log10FDR), 
                              shape = 16, stroke = 0, color = color_cancers2[cancer], alpha = 0.8)
          p <- p + geom_text_repel(data = tab4plot, mapping = aes(x = GENE, y = coef_capped, label = text2p, color = SUB_GENE.path, alpha = ifelse(top_effect_size_path, 0.8, 0.2)), 
                                   size = 2, force = 1)
        }
        if (self == "trans" & cptac_phase2process == "cptac2p_3can") {
          p <- p + facet_grid(facet_y~Cancer, scales = "free_x", space = "free_x")
          p <- p + theme(panel.spacing.y = unit(0, "lines"))
        }
        p <- p + scale_color_manual(values = tcga_path_colors)
        p <- p + theme_bw()
        p <- p + xlab(enzyme_type)+ylab(paste0("standardized coefficient of ", enzyme_type, " phosphoprotein abundance"))
        p <- p + theme(axis.text.x = element_text(angle = 90, face = "bold", hjust = 0.5, vjust = 0.5, size = 10))
        p <- p + theme(axis.title=element_text(size=12))
        p <- p + theme(axis.text.y = element_text(colour="black", size=10))
        p <- p + theme(title = element_text(size = 6, face = "italic"))
        p <- p + theme(panel.spacing.x = unit(0, "lines"), strip.text.x = element_text(size = 9, face = "bold"))
        p <- p + theme(panel.spacing.y = unit(0, "lines"))
        p
        fn = paste0(subdir4, 'substrates_of_', self,'_regulated_driver_gene_', enzyme_type, '_pairs_in_TCGA_pathways_reg_nonNA',reg_nonNA, '_fdr_thres', fdr_thres, '_curated.pdf')
        pdf(file = fn, height=6, width = 8, useDingbats = F)
        print(p)
        dev.off()
      }
    }
  }
}





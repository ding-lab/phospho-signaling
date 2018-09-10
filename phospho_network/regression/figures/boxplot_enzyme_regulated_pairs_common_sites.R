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
reg_nonNA2test <- seq(from = 5, to = 30, by = 5)
reg_nonNA2test <- c(20, 25)

top_kinase2show <- 10
top_substrate2show <- 2
top_site2show <- 1

pdf_sizes <- list(kinase = list(trans = list(shared3can = list(width = 8, height = 3),
                                             uniq_BRCA = list(width = 9, height = 3),
                                             uniq_OV = list(width = 9, height = 3),
                                             uniq_CO = list(width = 9, height = 3)),
                                cis = list(shared3can = list(width = 9, height = 3),
                                             uniq_BRCA = list(width = 9, height = 3),
                                             uniq_OV = list(width = 5, height = 3),
                                             uniq_CO = list(width = 5, height = 3))))
text_sizes <- list(kinase = list(trans = list(shared3can = 2,
                                             uniq_BRCA = 2,
                                             uniq_OV = 2,
                                             uniq_CO = 2),
                                cis = list(shared3can = 2,
                                           uniq_BRCA = 2,
                                           uniq_OV = 2,
                                           uniq_CO = 2)))

# cptac_phase2process <- "cptac3"
# cancers2process <- "UCEC"
cptac_phase2process <- "cptac2p_3can"
cancers2process <- cancers_sort
fdr_thres2process <- c(0.05, 0.1, 0.2, 0.3)
fdr_thres2process <- c(0.05, 0.1)

pdf_sizes_by_facetx <- list("1" = 5, "3" = 12)

# inputs -------------------------------------------------------------------
## input driver genes
# driver_genes <- read.delim(file = "./TCGA_data/reference_files/Consensus.Genelist.full.txt")

driver_pancan <- loadGeneList(gene_type = "driver", cancer = "PANCAN", is.soft.limit = "")
driver_can <- NULL
for (cancer in cancers2process) {
  driver_can <- unique(c(driver_can, as.vector(loadGeneList(gene_type = "driver", cancer = cancer, is.soft.limit = ""))))
}
## annotate kinase substrate regulation
for (enzyme_type in c("kinase")) {
  enzyme_type <- "kinase"
  subdir1 <- paste0(makeOutDir(resultD = resultD), enzyme_type, "/")
  dir.create(subdir1)
  for (reg_nonNA in reg_nonNA2test) {
    subdir2 <- paste0(subdir1, "reg_nonNA", reg_nonNA, "/")
    dir.create(subdir2)
    
    sup_cans_tab_en <- fread(input = paste0(ppnD, "regression/tables/", "change_regression_nonNA/", 
                                            enzyme_type, "_substrate_regression_", cptac_phase2process, "_tumor",
                                            "_reg_nonNA", reg_nonNA, "_commonsite_nonNA", reg_nonNA, ".txt"), data.table = F)
    for (fdr_thres in fdr_thres2process) {
      subdir3 <- paste0(subdir2, "fdr_thres", fdr_thres, "/")
      dir.create(subdir3)
      
      sup_cans_tab_en <- markSigSiteCan(regression = sup_cans_tab_en, sig_thres = fdr_thres, enzyme_type = enzyme_type)
      sup_cans_tab_en$regulated <- (sup_cans_tab_en$fdr_sig & sup_cans_tab_en$coef_sig)
      sup_cans_tab_en_reg <- sup_cans_tab_en[sup_cans_tab_en$regulated,]
      
      cross_can_list <- list()
      # for (pair_cat in c("uniq_BRCA", "uniq_OV", "uniq_CO")) {
        for (pair_cat in c("shared3can", "uniq_BRCA", "uniq_OV", "uniq_CO")) {
          
        subdir4 <- paste0(subdir3, pair_cat, "/")
        dir.create(subdir4)
        
        cross_can_list[[pair_cat]] <- sup_cans_tab_en[!is.na(sup_cans_tab_en[, pair_cat]) & sup_cans_tab_en[, pair_cat],]
        for (self in c("trans")) {
          subdir5 <- paste0(subdir4, self, "/")
          dir.create(subdir5)
          
          tab4plot_cans <- cross_can_list[[pair_cat]]
          tab4plot_cans <- tab4plot_cans[tab4plot_cans$SELF == self,]

          if (nrow(tab4plot_cans) > 0) {
            if (self == "cis") {
              # tab4plot <- unique(tab_can[(tab_can$GENE %in% c(driver_can, top_kinases)), c("GENE", "SUB_GENE", "SUB_MOD_RSD","Cancer", "FDR_pro_kin", "coef_pro_kin")])
              tab4plot <- unique(tab4plot_cans[, c("GENE", "SUB_GENE", "SUB_MOD_RSD","Cancer", "FDR_pro_kin", "coef_pro_kin")])
            }
            if (self == "trans") {
              if (pair_cat != "shared3can") {
                ## filter cis by showing only the SMGs
                tab4plot <- unique(tab4plot_cans[(tab4plot_cans$GENE %in% c(driver_can, driver_pancan)), c("GENE", "SUB_GENE", "SUB_MOD_RSD","Cancer", "FDR_pho_kin", "coef_pho_kin")])
              } else {
                tab4plot <- unique(tab4plot_cans[, c("GENE", "SUB_GENE", "SUB_MOD_RSD","Cancer", "FDR_pho_kin", "coef_pho_kin")])
              }
            }
            colnames(tab4plot) <- c("GENE", "SUB_GENE", "SUB_MOD_RSD", "Cancer", "FDR", "coef")
            
            # annotate driver genes to TCGA pathways ----------------------------------
            genes2pathways <- map2TCGApathwaways(gene_list = unique(as.vector(tab4plot$GENE)), pathway_list = tcga_pathways)
            genes2pathways <- data.frame(GENE = rep(x = names(genes2pathways), sapply(X = genes2pathways, FUN = function(x) length(x))), 
                                         GENE.path = unlist(genes2pathways, use.names = F))
            tab4plot <- merge(tab4plot, genes2pathways, all.x = T)
            
            # annotate substrate genes to TCGA pathways ----------------------------------
            sub_genes2pathways <- map2TCGApathwaways(gene_list = unique(tab4plot$SUB_GENE), pathway_list = tcga_pathways_pluskegg)
            sub_genes2pathways <- data.frame(SUB_GENE = rep(x = names(sub_genes2pathways), sapply(X = sub_genes2pathways, FUN = function(x) length(x))), 
                                             SUB_GENE.path = unlist(sub_genes2pathways, use.names = F))
            tab4plot <- merge(tab4plot, sub_genes2pathways, all.x = T)
          
            if (enzyme_type == "kinase" & (!(self == "trans" & pair_cat == "shared3can" ))) {
              tab4plot <- tab4plot[(tab4plot$GENE %in% tab4plot$GENE[!is.na(tab4plot$SUB_GENE.path)]) | (tab4plot$GENE %in% driver_can),]
            }
            
            # add more annotation columns to table ------------------------------------
            tab4plot$same_path <- (as.vector(tab4plot$GENE.path) == as.vector(tab4plot$SUB_GENE.path))
            tab4plot$pair_pro <- paste0(tab4plot$GENE, ":", tab4plot$SUB_GENE)
            tab4plot$pair <- paste0(tab4plot$pair_pro, "_", tab4plot$SUB_MOD_RSD)
            tab4plot$pair_can <- paste0(tab4plot$pair, "-", tab4plot$Cancer)
            tab4plot <- tab4plot[order(tab4plot$pair_can, tab4plot$same_path, decreasing = T),]
            tab4plot <- tab4plot[!(duplicated(tab4plot$pair_can)),]
            
            tab4plot$log10FDR <- -log10(tab4plot$FDR)
            tab4plot$site <- paste0(tab4plot$SUB_GENE, "\n", tab4plot$SUB_MOD_RSD)
            tab4plot$x <- paste0(tab4plot$GENE, "-", tab4plot$Cancer)
            tab4plot$site_coef <- paste0(tab4plot$site, "\n(coef=", signif(tab4plot$coef, digits = 2 ), ")")
            tab4plot$text <- as.vector(tab4plot$site_coef)
            
            # cap y-axis and edit text according to y-axis cap ---------------------------------------
            tab4plot$coef_capped <- tab4plot$coef
            if (pair_cat != "shared3can") {
              lim = quantile(x = tab4plot$coef, probs = 0.99)
              cap_upper <- min(lim, 2)
              tab4plot$coef_capped[tab4plot$coef > cap_upper] <- cap_upper
              lim = quantile(x = tab4plot$coef, probs = 0.01)
              cap_bottom <- max(lim, -2)
              tab4plot$coef_capped[tab4plot$coef < cap_bottom] <- cap_bottom
              tab4plot$text[tab4plot$coef < cap_upper & tab4plot$coef > cap_bottom] <- as.vector(tab4plot$site[tab4plot$coef < cap_upper & tab4plot$coef > cap_bottom])
            } else {
              if (enzyme_type == "kinase") {
                lim = quantile(x = tab4plot$coef, probs = 0.99)
                cap_upper <- min(lim, 2)
                tab4plot$coef_capped[tab4plot$coef > cap_upper] <- cap_upper
                tab4plot$text[tab4plot$coef < cap_upper] <- as.vector(tab4plot$site[tab4plot$coef < cap_upper])
              } else {
                lim = quantile(x = tab4plot$coef, probs = 0.01)
                cap_bottom <- max(lim, -2)
                tab4plot$coef_capped[tab4plot$coef < cap_bottom] <- cap_bottom
                tab4plot$text[tab4plot$coef > cap_bottom] <- as.vector(tab4plot$site[tab4plot$coef > cap_bottom])
              }
            }
           
            if (self == "trans" & cptac_phase2process == "cptac2p_3can") {
              facet_y <- as.vector(tab4plot$SUB_GENE.path)
              tab4plot$facet_y <- facet_y
              tab4plot$facet_y <- factor(tab4plot$facet_y, levels = c(pathway_names, "NA"))
            }
            tmp <- as.vector(tab4plot$SUB_GENE.path)
            tmp[is.na(tmp)] <- "NA"
            tab4plot$SUB_GENE.path <- tmp
            
            ## sort x-axis genes
            genes_sort <- as.vector(unique(tab4plot$GENE))
            genes_sort <- c(genes_sort[(genes_sort %in% driver_can)], genes_sort[!(genes_sort %in% driver_can)])
            
            ## sort gene on x-axis again
            genes2pathways <- genes2pathways[genes2pathways$GENE %in% tab4plot$GENE,]

            genes_sort2 <- NULL
            for (path in c(pathway_names, "NA")) {
              path_genes <- as.vector(unique(genes2pathways$GENE[genes2pathways$GENE.path == path & !is.na(genes2pathways$GENE.path)]))
              if (length(path_genes) >0 ) {
                genes_sort2 <- c(genes_sort2, path_genes)
              }
            }
            genes_sort2 <- c(genes_sort2, genes_sort[!(genes_sort %in% genes_sort2)])
            
            # assign color for x-axis genes
            color_enzyme <- NULL
            for (gene in genes_sort2) {
              gene_path <- unique(as.vector(tab4plot$GENE.path[tab4plot$GENE == gene]))
              gene_path <- gene_path[!is.na(gene_path)]
              if (length(gene_path) > 0){
                color_enzyme <- c(color_enzyme, tcga_path_colors[gene_path])
              } else {
                color_enzyme <- c(color_enzyme, "black")
              }
            }
            names(color_enzyme) <- genes_sort2
            

            ## annotate the substrates with top effect size
            tab4plot <- tab4plot[order(tab4plot$pair_pro, tab4plot$coef, decreasing = T),]
            top_effect_size <- vector(mode = "logical", length = nrow(tab4plot))
            for (gene in unique(tab4plot$GENE)) {
              tab_gene <- tab4plot[tab4plot$GENE == gene,]
              for (site in unique(tab_gene$site)) {
                tab_site <- tab_gene[tab_gene$site == site,]
                order_type <- (enzyme_type == "kinase")
                coefs <- as.vector(tab_site$coef)
                coefs <- coefs[order(coefs, decreasing = order_type)]
                top_coefs <- coefs[1:min(length(coefs), top_site2show)]
                top_effect_size[tab4plot$coef %in% top_coefs] <- TRUE
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
              ## take out for the substrates whose kinase already is annotated with some substrates
              tab4plot$text2p[tab4plot$SUB_GENE.path == "NA" & (tab4plot$GENE %in% as.vector(tab4plot$GENE[tab4plot$SUB_GENE.path != "NA"]))] <- ""
              tab4plot$text2p[tab4plot$SUB_GENE.path == "NA" & !(tab4plot$GENE %in% as.vector(tab4plot$GENE[tab4plot$SUB_GENE.path != "NA"])) & !tab4plot$top_effect_size] <- ""
              if (pair_cat != "shared3can") {
                tab4plot <- tab4plot[tab4plot$text2p != "",]
              }
            } else {
              tab4plot$text2p <- as.vector(tab4plot$text)
              tab4plot$text2p[!tab4plot$top_effect_size] <- ""
            }
            tab4plot$GENE <- factor(tab4plot$GENE, levels = genes_sort2)
            tab4plot$Cancer <- factor(tab4plot$Cancer, levels = cancers2process)
            tab4plot$site <- factor(tab4plot$site, levels = unique(as.vector(tab4plot$site)[order(tab4plot$SUB_GENE.path)]))
            # ggplot -----------------------------------------------------------------
            
            p <- ggplot()
            if (enzyme_type == "kinase") {
              p <- p + geom_point(data=tab4plot, aes(x = site, y = coef_capped, size = log10FDR, alpha = ifelse(!is.na(SUB_GENE.path), 0.8, 0.6), color = Cancer), 
                                  shape = 16, stroke = 0)
              p <- p + geom_text_repel(data = tab4plot[tab4plot$top_effect_size,], mapping = aes(x = site, y = coef_capped, label = text2p, color = SUB_GENE.path),
                                       size = text_sizes[[enzyme_type]][[self]][[pair_cat]], force = 1, alpha = 0.8, fontface = "bold")
            } else {
              p <- p + geom_point(data=tab4plot, aes(x = SUB_GENE, y = coef_capped, size = log10FDR, color = Cancer), 
                                  shape = 16, stroke = 0, alpha = 0.8)
              p <- p + geom_text_repel(data = tab4plot, mapping = aes(x = SUB_GENE, y = coef_capped, label = text2p, color = SUB_GENE.path, alpha = ifelse(top_effect_size_path, 0.8, 0.2)), 
                                       size = text_sizes[[enzyme_type]][[self]][[pair_cat]], force = 1)
            }
            if (cptac_phase2process == "cptac2p_3can") {
              p <- p + facet_grid(.~GENE, scales = "free_x", space = "free_x")
              p <- p + theme(panel.background = element_rect(fill = brewer.pal(9, "Greys")[2], size = 0))
              p <- p + theme(strip.background.x = element_rect(color = "black"))
              p <- p + theme(panel.spacing.x = unit(0, "lines"), strip.text.x = element_text(size = 4))
            }
            p <- p + geom_hline(yintercept = 0, linetype = 2, color = "grey")
            p <- p + scale_color_manual(values = c(tcga_path_colors, color_cancers2))
            p <- p + xlab('kinase')+ylab("effect size")
            p <- p + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank())
            p <- p + theme(axis.text.y = element_text(colour="black", size=10))
            fn = paste0(subdir5, 'substrates_of_', self,'_regulated_driver_gene_', pair_cat, "_", enzyme_type, '_pairs_in_TCGA_pathways_reg_nonNA',reg_nonNA, '_fdr_thres', fdr_thres, '.pdf')
            pdf(file = fn, 
                height= pdf_sizes[[enzyme_type]][[self]][[pair_cat]]$height, 
                width = pdf_sizes[[enzyme_type]][[self]][[pair_cat]]$width, useDingbats = F)
            print(p)
            dev.off()
            
          }
        }
      }
      
    }
  }
}









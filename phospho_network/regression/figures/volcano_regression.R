# Yige Wu @ WashU 2018 Jan
# draw volcano plots for regression result (3 cancer types together)

# source ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_plotting.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
library(ggrepel)

# inputs ------------------------------------------------------------------
# set variables -----------------------------------------------------------
reg_nonNA2test <- c(20)
top_num <- 5
# cptac_phase2process <- "cptac3"
# cancers2process <- "UCEC"
cptac_phase2process <- "cptac2p_3can"
cancers2process <- cancers_sort
pdf_sizes_by_facetx <- list("cptac3" = 5, "cptac2p_3can" = 12)
cap_by_cptac_phase <- list("cptac2p_3can" = list(trans = 15, cis = 30),
                      "cptac3" = list(trans = 5, cis = 30))
fdr_thres <- c(0.05, 0.2); names(fdr_thres) <- c("kinase", "phosphatase")

# kinases ## unique in each cancer -----------------------------------------------------------------
for (reg_nonNA in reg_nonNA2test) {
  for (enzyme_type in c("kinase")) {
    # regression <- fread(input = paste0(ppnD, "regression/tables/", "change_regression_nonNA/", 
    #                                         enzyme_type, "_substrate_regression_cptac2p_3can_tumor",
    #                                         "_reg_nonNA", reg_nonNA, ".txt"), data.table = F)
    regression <- fread(input = paste0(ppnD, "regression/tables/", "change_regression_nonNA/", 
                                       enzyme_type, "_substrate_regression_", cptac_phase2process, "_tumor",
                                       "_reg_nonNA", reg_nonNA, ".txt"), data.table = F)
    regression <- markSigKS(regression = regression, sig_thres = fdr_thres[enzyme_type], enzyme_type = enzyme_type)
    regression$regulated <- (regression$fdr_sig & regression$coef_sig)
    
    subdir1 <- paste0(makeOutDir(resultD = resultD), enzyme_type, "/")
    dir.create(subdir1)
    for (self in c("trans")) {
      tab2p <- regression
      tab2p <- tab2p[tab2p$SELF == self,]
      tab2p <- unique(tab2p[, c("coef_pho_kin", "FDR_pho_kin", "pair", "Cancer","regulated", "GENE", "SUB_GENE")])
      tab2p <- tab2p[sample(x = 1:nrow(tab2p), size = nrow(tab2p), replace = F),]
      tab2p$coef_capped <- remove_outliers(x = as.vector(tab2p$coef_pho_kin), out_thres = 2, na.rm = T)
      tab2p$log_FDR <- -log10(tab2p$FDR_pho_kin)
      tab2p$y <- as.vector(tab2p$log_FDR)
      
      cap <- cap_by_cptac_phase[[cptac_phase2process]][[self]]
      cap <- min(max(tab2p$log_FDR), cap)
      tab2p$y[tab2p$y > cap] <- cap
      tab2p$y[tab2p$y < (-cap)] <- -cap
      tab2p <- tab2p[order(tab2p$log_FDR, decreasing = T),]

      tab2p_regulated <- tab2p[tab2p$regulated,]
      tab2p_unregulated <- tab2p[!(tab2p$regulated),]
      tab2p_regulated$Cancer <- factor(tab2p_regulated$Cancer, levels = cancers2process)
      tab2p_unregulated$Cancer <- factor(tab2p_unregulated$Cancer, levels = cancers2process)
      
      sub_genes2pathways <- map2TCGApathwaways(gene_list = unique(tab2p_regulated$SUB_GENE), pathway_list = tcga_pathways_pluskegg_and_pathway)
      sub_genes2pathways <- data.frame(SUB_GENE = rep(x = names(sub_genes2pathways), sapply(X = sub_genes2pathways, FUN = function(x) length(x))), 
                                       SUB_GENE.path = unlist(sub_genes2pathways, use.names = F))
      
      
      tab2p_text <- NULL
      for (cancer in cancers_sort) {
        genes2show <- unique(c(SMGs[[cancer]], driver_genes$Gene[driver_genes$Cancer == ifelse(cancer == "CO", "COADREAD", cancer)]))
        tab2p_text2add <- tab2p[tab2p$pair %in% head(tab2p$pair[tab2p$regulated & tab2p$Cancer == cancer & (tab2p$GENE %in% genes2show & tab2p$SUB_GENE %in% sub_genes2pathways$SUB_GENE) ], n = 6) & tab2p$Cancer == cancer,]
        if (nrow(tab2p_text2add) == 0) {
          tab2p_text2add <- tab2p[tab2p$pair %in% head(tab2p$pair[tab2p$regulated & tab2p$Cancer == cancer & (tab2p$GENE %in% genes2show | tab2p$SUB_GENE %in% genes2show) ], n = 6) & tab2p$Cancer == cancer,]
          
        }
        tab2p_text <- rbind(tab2p_text, tab2p_text2add)
      }
      
      for (text_cutoff in c(-log10(fdr_thres[enzyme_type]))) {
        subdir2 <- paste0(subdir1, "textcutoff", text_cutoff, "/")
        dir.create(subdir2)
        p = ggplot()
        p = p + geom_point(data = tab2p_regulated, mapping = aes(x=coef_capped, y= y, color = Cancer), 
                           stroke = 0, alpha = 0.3, shape = 16, size = 2)
        p = p + geom_point(data = tab2p_unregulated, mapping = aes(x=coef_capped, y= y), 
                           stroke = 0, alpha = 0.1, shape = 16, size = 2, color = "grey")
        p <- p + geom_text_repel(data = tab2p_text, mapping = aes(x=coef_capped, y = y, label = pair), force = 2, color = "black")
        p = p + geom_hline(yintercept = text_cutoff, color = 'grey', linetype = 2)
        p = p + geom_vline(xintercept = 0, color = 'grey', linetype = 2)
        p <- p + facet_grid(.~Cancer, scales = "fixed", space = "fixed")
        p <- p + theme_minimal() + theme(panel.background = element_rect(fill = brewer.pal(9, "Greys")[2], size = 0))
        p = p + scale_color_manual(values = color_cancers2)
        p = p + scale_fill_manual(values = color_cancers2)
        p = p + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
        p = p + theme(strip.text.x = element_text(size = 20, face = "bold")) # panel.spacing.x=unit(0, "lines"), 
        p = p + labs(x = paste0("standardized coefficient for ", enzyme_type, " phosphorylation level"), y="-log10(FDR)")
        p
        # fn = paste(subdir2, "volcano_trans_uniq_", cancer, "_", enzyme_type, '_substrate_cap', cap, '.pdf',sep ="")
        fn = paste(subdir2, "volcano_", cptac_phase2process, "_", self, "_", enzyme_type, '_substrate_cap', cap, '_reg_nonNA', reg_nonNA, '.jpeg',sep ="")
        ggsave(file=fn, height=4, width=pdf_sizes_by_facetx[[cptac_phase2process]], device = "jpeg")
      }
    }
    
    for (self in c("cis")) {
      tab2p <- regression
      tab2p <- tab2p[tab2p$SELF == self,]
      tab2p <- unique(tab2p[, c("coef_pro_kin", "FDR_pro_kin", "pair", "Cancer","regulated", "GENE", "SUB_GENE")])
      tab2p <- tab2p[sample(x = 1:nrow(tab2p), size = nrow(tab2p), replace = F),]
      tab2p$coef_capped <- remove_outliers(x = as.vector(tab2p$coef_pro_kin), out_thres = 2, na.rm = T)
      tab2p$log_FDR <- -log10(tab2p$FDR_pro_kin)
      tab2p$y <- as.vector(tab2p$log_FDR)
      cap <- cap_by_cptac_phase[[cptac_phase2process]][[self]]
      cap <- min(max(tab2p$log_FDR), cap)
      tab2p$y[tab2p$y > cap] <- cap
      tab2p$y[tab2p$y < (-cap)] <- -cap
      tab2p <- tab2p[order(tab2p$log_FDR, decreasing = T),]
      tab2p$Cancer <- factor(tab2p$Cancer, levels = cancers2process)
      
      sub_genes2pathways <- map2TCGApathwaways(gene_list = unique(tab2p$SUB_GENE), pathway_list = tcga_pathways_pluskegg_and_pathway)
      sub_genes2pathways <- data.frame(SUB_GENE = rep(x = names(sub_genes2pathways), sapply(X = sub_genes2pathways, FUN = function(x) length(x))), 
                                       SUB_GENE.path = unlist(sub_genes2pathways, use.names = F))
      
      tab2p_text <- NULL
      for (cancer in cancers_sort) {
        genes2show <- unique(c(SMGs[[cancer]], driver_genes$Gene[driver_genes$Cancer == ifelse(cancer == "CO", "COADREAD", cancer)]))
        tab2p_text2add <- tab2p[tab2p$pair %in% head(tab2p$pair[tab2p$regulated & tab2p$Cancer == cancer & (tab2p$GENE %in% genes2show) ], n = 6) & tab2p$Cancer == cancer,]
        if (nrow(tab2p_text2add) == 0) {
          tab2p_text2add <- tab2p[tab2p$pair %in% head(tab2p$pair[tab2p$regulated & tab2p$Cancer == cancer & (tab2p$GENE %in% sub_genes2pathways$SUB_GENE) ], n = 6) & tab2p$Cancer == cancer,]
        }
        tab2p_text <- rbind(tab2p_text, tab2p_text2add)
      }
      
      
      
      for (text_cutoff in c(-log10(fdr_pk))) {
        subdir2 <- paste0(subdir1, "textcutoff", text_cutoff, "/")
        dir.create(subdir2)
        p = ggplot()
        p = p + geom_point(data = tab2p[tab2p$regulated,], mapping = aes(x=coef_capped, y= y, color = Cancer), 
                           stroke = 0, alpha = 0.3, shape = 16, size = 2)
        p = p + geom_point(data = tab2p[!(tab2p$regulated),], mapping = aes(x=coef_capped, y= y), 
                           stroke = 0, alpha = 0.1, shape = 16, size = 2, color = "grey")
        p <- p + geom_text_repel(data = tab2p_text, mapping = aes(x=coef_capped, y = y, label = pair), force = 2, color = "black")
        p = p + geom_hline(yintercept = text_cutoff, color = 'grey', linetype = 2)
        p = p + geom_vline(xintercept = 0, color = 'grey', linetype = 2)
        p <- p + facet_grid(.~Cancer, scales = "fixed", space = "fixed")
        p <- p + theme_minimal() + theme(panel.background = element_rect(fill = brewer.pal(9, "Greys")[2], size = 0))
        p = p + scale_color_manual(values = color_cancers2)
        p = p + scale_fill_manual(values = color_cancers2)
        p = p + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
        p = p + theme(strip.text.x = element_text(size = 20, face = "bold")) #panel.spacing.x=unit(0, "lines"), 
        p = p + labs(x = paste0("standardized coefficient for ", enzyme_type, " protein level"), y="-log10(FDR)")
        p
        # fn = paste(subdir2, "volcano_trans_uniq_", cancer, "_", enzyme_type, '_substrate_cap', cap, '.pdf',sep ="")
        fn = paste(subdir2, "volcano_", cptac_phase2process, "_", self, "_", enzyme_type, '_substrate_cap', cap,'_reg_nonNA', reg_nonNA, '.jpeg',sep ="")
        ggsave(file=fn, height=4, width=pdf_sizes_by_facetx[[cptac_phase2process]], device = "jpeg")
      }
    }
  }
}




# kinase all cancers -----------------------------------------------------------------
# enzyme_type <- "kinase"
# regression <- fread(input = paste0("./cptac2p/analysis_results/phospho_network/regression/tables/regression_omnipath&newworkin_protein_level/kinase_substrate_regression_cptac2p_3can_tumor.txt"), data.table = F)
# regression$GENE <- as.vector(regression$KINASE)
# regression$SUB_GENE <- as.vector(regression$SUBSTRATE)
# regression <- markSigSiteCan(regression = regression, sig_thres = fdr_pk, enzyme_type = enzyme_type)
# regression$regulated <- (regression$fdr_sig & regression$coef_sig)
# 
# subdir1 <- paste0(makeOutDir(resultD = resultD), enzyme_type, "/")
# dir.create(subdir1) 
# 
# for (text_cutoff in c(-log10(fdr_pk))) {
#   ## get KSEA scores integrated with differential phosphorylation results
#   tab2p <- regression[regression$SELF == "trans",]
#   tab2p <- tab2p[sample(x = 1:nrow(tab2p), size = nrow(tab2p), replace = F),]
#   tab2p$coef_capped <- remove_outliers(x = as.vector(tab2p$coef_pho_kin), out_thres = 2, na.rm = T)
#   tab2p$log_FDR <- -log10(tab2p$FDR_pho_kin)
#   tab2p$y <- as.vector(tab2p$log_FDR)
#   cap <- 15
#   tab2p$y[tab2p$y > cap] <- cap
#   tab2p$y[tab2p$y < (-cap)] <- -cap
#   tab2p <- tab2p[order(tab2p$log_FDR, decreasing = T),]
#   
#   p = ggplot()
#   p = p + geom_point(data = tab2p[tab2p$regulated,], mapping = aes(x=coef_capped, y= y, color = Cancer), stroke = 0, shape = 16, alpha = 0.5)
#   p = p + geom_point(data = tab2p[!(tab2p$regulated),], mapping = aes(x=coef_capped, y= y), color = "grey", stroke = 0, shape = 16, alpha = 0.1)
#   p = p + geom_hline(yintercept = text_cutoff, color = 'grey', linetype = 2)
#   p = p + geom_vline(xintercept = 0, color = 'grey', linetype = 2)
#   p = p + theme_bw() + theme_nogrid()
#   p = p + scale_color_manual(values = color_cancers2)
#   p = p + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
#   p = p + theme(panel.spacing.y=unit(0, "lines"), strip.text.y = element_text(size = 12))
#   p = p + labs(x = paste0("Coefficient for ", enzyme_type, " phosphorylation level"), y="-log10(FDR)")
#   p = p + theme(legend.position = "none")
#   p
#   subdir2 <- paste0(subdir1, "textcutoff", text_cutoff, "/")
#   dir.create(subdir2)
#   fn = paste(subdir2, "volcano_trans_", "cancers" , "_", enzyme_type, '_substrate_cap', cap, '.pdf',sep ="")
#   ggsave(file=fn, height=4, width=5, useDingbats=FALSE)
#   
#   p = ggplot()
#   p = p + geom_point(data = tab2p[tab2p$regulated,], mapping = aes(x=coef_capped, y= y, color = Cancer, alpha = ifelse(Cancer == "BRCA", 0.4, 0.5)), stroke = 0, shape = 16)
#   p = p + geom_point(data = tab2p[!(tab2p$regulated),], mapping = aes(x=coef_capped, y= y), color = "grey", stroke = 0, shape = 16, alpha = 0.1)
#   p = p + geom_hline(yintercept = text_cutoff, color = 'grey', linetype = 2)
#   p = p + geom_vline(xintercept = 0, color = 'grey', linetype = 2)
#   p = p + theme_bw() + theme_nogrid()
#   p = p + scale_color_manual(values = color_cancers2)
#   p = p + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
#   p = p + theme(panel.spacing.y=unit(0, "lines"), strip.text.y = element_text(size = 12))
#   p = p + labs(x = paste0("Coefficient for ", enzyme_type, " phosphorylation level"), y="-log10(FDR)")
#   p = p + theme(legend.position = "none")
#   p
#   subdir2 <- paste0(subdir1, "textcutoff", text_cutoff, "/")
#   dir.create(subdir2)
#   fn = paste(subdir2, "volcano_trans_", "cancers" , "_", enzyme_type, '_substrate_cap', cap, '_BRCA_lightened.pdf',sep ="")
#   ggsave(file=fn, height=4, width=5, useDingbats=FALSE)
#   
#   p = ggplot()
#   p = p + geom_point(data = tab2p[tab2p$regulated,], mapping = aes(x=coef_capped, y= y, color = Cancer), stroke = 0, alpha = 0.5, shape = 16)
#   p = p + geom_point(data = tab2p[!(tab2p$regulated),], mapping = aes(x=coef_capped, y= y), color = "grey", stroke = 0, shape = 16, alpha = 0.1)
#   p = p + geom_hline(yintercept = text_cutoff, color = 'grey', linetype = 2)
#   p = p + geom_vline(xintercept = 0, color = 'grey', linetype = 2)
#   p = p + theme_minimal() + theme_nogrid()
#   p = p + facet_grid(.~Cancer, scales = "fixed", space = "fixed")
#   p = p + geom_text_repel(data = tab2p[1:top_num,], mapping = aes(x=coef_capped, y=y,
#                                                                   label= as.character(pair),
#                                                                   color = Cancer), size= 2, alpha = 0.7, force = 1)
#   p = p + scale_color_manual(values = color_cancers2)
#   p = p + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
#   p = p + theme(panel.spacing.y=unit(0, "lines"), strip.text.y = element_text(size = 12))
#   p = p + labs(x = paste0("Coefficient for ", enzyme_type, " phosphorylation level"), y="-log10(FDR)")
#   p = p + theme(legend.position = "none")
#   p
#   subdir2 <- paste0(subdir1, "textcutoff", text_cutoff, "/")
#   dir.create(subdir2)
#   fn = paste(subdir2, "volcano_trans_", "cancers" , "_", enzyme_type, '_substrate_cap', cap, '_top_labeled.pdf',sep ="")
#   ggsave(file=fn, height=3, width=8, useDingbats=FALSE)
# }
# 
# for (text_cutoff in c(-log10(fdr_pk))) {
#   ## get KSEA scores integrated with differential phosphorylation results
#   tab2p <- regression[regression$SELF == "cis",]
#   tab2p <- tab2p[sample(x = 1:nrow(tab2p), size = nrow(tab2p), replace = F),]
#   tab2p$coef_capped <- remove_outliers(x = as.vector(tab2p$coef_pro_kin), out_thres = 2, na.rm = T)
#   tab2p$log_FDR <- -log10(tab2p$FDR_pro_kin)
#   tab2p$y <- as.vector(tab2p$log_FDR)
#   cap <- 50
#   tab2p$y[tab2p$y > cap] <- cap
#   tab2p$y[tab2p$y < (-cap)] <- -cap
#   tab2p <- tab2p[order(tab2p$log_FDR, decreasing = T),]
#   tab2p$Cancer <- factor(tab2p$Cancer, levels = cancers2process)
#   
#   p = ggplot()
#   p = p + geom_point(data = tab2p[tab2p$regulated,], mapping = aes(x=coef_capped, y= y, color = Cancer), stroke = 0, shape = 16, alpha = 0.5)
#   p = p + geom_point(data = tab2p[!(tab2p$regulated),], mapping = aes(x=coef_capped, y= y), color = "grey", stroke = 0, shape = 16, alpha = 0.1)
#   p = p + geom_hline(yintercept = text_cutoff, color = 'grey', linetype = 2)
#   p = p + geom_vline(xintercept = 0, color = 'grey', linetype = 2)
#   p = p + theme_bw() + theme_nogrid()
#   p = p + scale_color_manual(values = color_cancers2)
#   p = p + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
#   p = p + theme(panel.spacing.y=unit(0, "lines"), strip.text.y = element_text(size = 12))
#   p = p + labs(x = paste0("Coefficient for ", enzyme_type, " phosphorylation level"), y="-log10(FDR)")
#   p = p + theme(legend.position = "none")
#   p
#   subdir2 <- paste0(subdir1, "textcutoff", text_cutoff, "/")
#   dir.create(subdir2)
#   fn = paste(subdir2, "volcano_cis_", "cancers" , "_", enzyme_type, '_substrate_cap', cap, '.pdf',sep ="")
#   ggsave(file=fn, height=4, width=5, useDingbats=FALSE)
#   
#   p = ggplot()
#   p = p + geom_point(data = tab2p[tab2p$regulated,], mapping = aes(x=coef_capped, y= y, color = Cancer, alpha = ifelse(Cancer == "BRCA", 0.4, 0.5)), stroke = 0, shape = 16)
#   p = p + geom_point(data = tab2p[!(tab2p$regulated),], mapping = aes(x=coef_capped, y= y), color = "grey", stroke = 0, shape = 16, alpha = 0.1)
#   p = p + geom_hline(yintercept = text_cutoff, color = 'grey', linetype = 2)
#   p = p + geom_vline(xintercept = 0, color = 'grey', linetype = 2)
#   p = p + theme_bw() + theme_nogrid()
#   p = p + scale_color_manual(values = color_cancers2)
#   p = p + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
#   p = p + theme(panel.spacing.y=unit(0, "lines"), strip.text.y = element_text(size = 12))
#   p = p + labs(x = paste0("Coefficient for ", enzyme_type, " phosphorylation level"), y="-log10(FDR)")
#   p = p + theme(legend.position = "none")
#   p
#   subdir2 <- paste0(subdir1, "textcutoff", text_cutoff, "/")
#   dir.create(subdir2)
#   fn = paste(subdir2, "volcano_cis_", "cancers" , "_", enzyme_type, '_substrate_cap', cap, '_BRCA_lightened.pdf',sep ="")
#   ggsave(file=fn, height=4, width=5, useDingbats=FALSE)
#   
#   p = ggplot()
#   p = p + geom_point(data = tab2p[tab2p$regulated,], mapping = aes(x=coef_capped, y= y, color = Cancer), stroke = 0, alpha = 0.5, shape = 16)
#   p = p + geom_point(data = tab2p[!(tab2p$regulated),], mapping = aes(x=coef_capped, y= y), color = "grey", stroke = 0, shape = 16, alpha = 0.1)
#   p = p + geom_hline(yintercept = text_cutoff, color = 'grey', linetype = 2)
#   p = p + geom_vline(xintercept = 0, color = 'grey', linetype = 2)
#   p = p + theme_minimal() + theme_nogrid()
#   p = p + facet_grid(.~Cancer, scales = "fixed", space = "fixed")
#   p = p + geom_text_repel(data = tab2p[1:top_num,], mapping = aes(x=coef_capped, y=y,
#                                                                   label= as.character(pair),
#                                                                   color = Cancer), size= 2, alpha = 0.7, force = 1)
#   p = p + scale_color_manual(values = color_cancers2)
#   p = p + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
#   p = p + theme(panel.spacing.y=unit(0, "lines"), strip.text.y = element_text(size = 12))
#   p = p + labs(x = paste0("Coefficient for ", enzyme_type, " phosphorylation level"), y="-log10(FDR)")
#   p = p + theme(legend.position = "none")
#   p
#   subdir2 <- paste0(subdir1, "textcutoff", text_cutoff, "/")
#   dir.create(subdir2)
#   fn = paste(subdir2, "volcano_cis_", "cancers" , "_", enzyme_type, '_substrate_cap', cap, '_top_labeled.pdf',sep ="")
#   ggsave(file=fn, height=3, width=8, useDingbats=FALSE)
# }
# 
# # phosphatase per cancer -----------------------------------------------------------------
# for (enzyme_type in c("phosphatase")) {
#   regression <- fread(input = paste0("./cptac2p/analysis_results/phospho_network/regression/tables/regression_omnipath&newworkin&depod_protein_level/phosphatase_substrate_regression_cptac2p_3can_tumor.txt"), data.table = F)
#   subdir1 <- paste0(makeOutDir(resultD = resultD), enzyme_type, "/")
#   dir.create(subdir1)  
#   for (cancer in cancers2process) {
#     for (self in c("trans")) {
#       df0 <- regression[regression$SELF == self & !is.na(regression$FDR_pho_kin) & regression$Cancer == cancer,]
#       # df1 <- df0[-log10(df0$FDR_pho_kin) < 30,]
#       df1 <- df0
#       df1$coef_pho_kin_filtered = remove_outliers(df1$coef_pho_kin, out_thres = 2)
#       
#       df <- unique(df1[,c("coef_pho_kin_filtered", "FDR_pho_kin", "P_pho_kin", "pair")])
#       df <- df[!is.na(df$coef_pho_kin_filtered) & !is.na(df$FDR_pho_kin),]
#       df$pass_FDR <- ifelse(df$coef_pho_kin_filtered < 0 & df$FDR_pho_kin < fdr_pp, TRUE, FALSE)
#       df$pass_pvalue <- ifelse(df$coef_pho_kin_filtered < 0 & df$P_pho_kin < fdr_pp, TRUE, FALSE)
#       
#       for (text_cutoff in c(-log10(fdr_pp))) {
#         ## get KSEA scores integrated with differential phosphorylation results
#         p = ggplot(df)
#         p = p + geom_point(aes(x=coef_pho_kin_filtered, y=-log10(FDR_pho_kin), 
#                                color = pass_FDR), stroke = 0, alpha = 0.5)
#         p = p + theme_bw() + theme_nogrid()
#         p = p + geom_hline(yintercept = text_cutoff, color = 'grey', linetype = 2)
#         p = p + geom_vline(xintercept = 0, color = 'grey', linetype = 2)
#         p = p + geom_text_repel(aes(x=coef_pho_kin_filtered, y=-log10(FDR_pho_kin), 
#                                     label= ifelse(-log10(FDR_pho_kin) > text_cutoff & coef_pho_kin_filtered > 0, as.character(pair), NA)), 
#                                 size=2, color = "grey", alpha = 0.5)
#         p = p + geom_text_repel(aes(x=coef_pho_kin_filtered, y=-log10(FDR_pho_kin), 
#                                     label= ifelse(-log10(FDR_pho_kin) > text_cutoff & coef_pho_kin_filtered < 0, as.character(pair), NA)),
#                                 size=3, color = "black", alpha = 1, force = 2)
#         p = p + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
#         p = p + labs(x = paste0("Coefficient for ", enzyme_type, " phosphorylation level"), y="-log10(FDR)")
#         p = p + scale_color_manual(values = c("TRUE" = "#377EB8", "FALSE" = "grey", "NA" = "grey")) 
#         p = p + theme(legend.position="none")
#         p
#         subdir2 <- paste0(subdir1, "textcutoff", text_cutoff, "/")
#         dir.create(subdir2)
#         fn = paste(subdir2, "volcano_trans_", cancer , "_", enzyme_type, '_substrate_textcutoff', text_cutoff, '.pdf',sep ="")
#         ggsave(file=fn, height=4, width=4, useDingbats=FALSE)
#       }
#     }
#   }
# }
# 
# # phosphatase all cancer -----------------------------------------------------------------
# regression <- fread(input = paste0("./cptac2p/analysis_results/phospho_network/regression/tables/regression_omnipath&newworkin&depod_protein_level/phosphatase_substrate_regression_cptac2p_3can_tumor.txt"), data.table = F)
# regression$GENE <- as.vector(regression$KINASE)
# regression$SUB_GENE <- as.vector(regression$SUBSTRATE)
# regression <- markSigSiteCan(regression = regression, sig_thres = fdr_pp, enzyme_type = "phosphatase")
# regression$regulated <- (regression$fdr_sig & regression$coef_sig)
# 
# for (text_cutoff in c(-log10(fdr_pp))) {
#   ## get KSEA scores integrated with differential phosphorylation results
#   tab2p <- regression[regression$SELF == "trans",]
#   tab2p <- tab2p[sample(x = 1:nrow(tab2p), size = nrow(tab2p), replace = F),]
#   tab2p$coef_capped <- remove_outliers(x = as.vector(tab2p$coef_pho_kin), out_thres = 2, na.rm = T)
#   p = ggplot()
#   p = p + geom_point(data = tab2p[tab2p$regulated,], mapping = aes(x=coef_capped, y=-log10(FDR_pho_kin), 
#                          color = Cancer), stroke = 0, alpha = 0.8, shape = 16)
#   p = p + geom_point(data = tab2p[!(tab2p$regulated),], mapping = aes(x=coef_capped, y=-log10(FDR_pho_kin)), color = "grey", stroke = 0, alpha = 0.2, shape = 16)
#   p = p + geom_hline(yintercept = text_cutoff, color = 'grey', linetype = 2)
#   p = p + geom_vline(xintercept = 0, color = 'grey', linetype = 2)
#   p = p + theme_bw() + theme_nogrid()
#   # p = p + geom_text_repel(data = tab2p[tab2p$regulated,], mapping = aes(x=coef_capped, y=-log10(FDR_pho_kin),
#   #                             label= as.character(pair),
#   #                             color = Cancer), size= 1, alpha = 0.7, force = 1)
#   p = p + scale_color_manual(values = color_cancers2)
#   p = p + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
#   p = p + theme(panel.spacing.y=unit(0, "lines"), strip.text.y = element_text(size = 12))
#   p = p + labs(x = paste0("Coefficient for ", enzyme_type, " phosphorylation level"), y="-log10(FDR)")
#   p
#   subdir2 <- paste0(subdir1, "textcutoff", text_cutoff, "/")
#   dir.create(subdir2)
#   fn = paste(subdir2, "volcano_trans_", "cancers" , "_", enzyme_type, '_substrate_textcutoff', text_cutoff, '.pdf',sep ="")
#   ggsave(file=fn, height=4, width=5, useDingbats=FALSE)
# 
# }
# 
# regression <- fread(input = paste0("./cptac2p/analysis_results/phospho_network/regression/tables/regression_omnipath&newworkin&depod_protein_level/phosphatase_substrate_regression_cptac2p_3can_tumor.txt"), data.table = F)
# regression$GENE <- as.vector(regression$KINASE)
# regression$SUB_GENE <- as.vector(regression$SUBSTRATE)
# regression$site <- paste0(regression$SUBSTRATE, ":", regression$SUB_MOD_RSD)
# regression <- regression[!is.na(regression$FDR_pho_kin),]
# regression <- regression[order(regression$FDR_pho_kin),]
# regression <- markSigSiteCan(regression = regression, sig_thres = fdr_pp, enzyme_type = "phosphatase")
# regression$regulated <- (regression$fdr_sig & regression$coef_sig)
# table(regression$regulated)
# length(unique(regression$site[regression$regulated]))
# table(regression$regulated, regression$Cancer)

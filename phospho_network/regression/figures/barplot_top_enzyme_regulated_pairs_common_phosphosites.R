# Yige Wu @ WashU 2018 Jul
# barplot for the top kinases with high and low kinase-substrate pairs

# source ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_plotting.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/dependencies/tables/load_TCGA_pathways.R')
# variables ---------------------------------------------------------------
top_enzyme2show <- 10
nudge_ys <- c(0.3, 0.1); names(nudge_ys) <- c("trans", "cis")
ylims <- list(cis = c(1, 150),
             trans = c(1,10000))
# cptac_phase2process <- "cptac3"
# cancers2process <- "UCEC"
cptac_phase2process <- "cptac2p_3can"
cancers2process <- cancers_sort
pdf_sizes <- list(kinase = list(trans = list(width = 5.5, height = 3), cis = list(width = 6.5, height = 3)))
pdf_sizes_by_facetx <- list("1" = 5, "3" = 12)
reg_nonNA2test <- seq(from = 5, to = 25, by = 5)
reg_nonNA2test <- c(20, 25)
fdr_thres2process <- c(0.05, 0.1)
pair_cat2process <- c("shared3can", "uniq_BRCA", "uniq_OV", "uniq_CO")

# inputs -------------------------------------------------------------------
driver_pancan <- loadGeneList(gene_type = "driver", cancer = "PANCAN", is.soft.limit = "")

# 3 cancers separate-faceted, sorted by the number of substrates in TCGA pathways------------------------------------------------
for (enzyme_type in c("kinase")) {
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
      
      # cancer separate pdf -----------------------------------------------------
      cross_can_list <- list()
      pairs_per_en <- NULL
      for (pair_cat in pair_cat2process) {
      # for (pair_cat in c("uniq_BRCA")) {
        subdir4 <- paste0(subdir3, pair_cat, "/")
        dir.create(subdir4)
        
        for (self in c("trans")) {
          subdir5 <- paste0(subdir4, self, "/")
          dir.create(subdir5)

          tab_self <- sup_cans_tab_en
          tab_self <- tab_self[tab_self$SELF == self,]
          tab_self_regulated <- tab_self[tab_self$regulated & !is.na(tab_self[, pair_cat]) & tab_self[, pair_cat],]
          
          pairs_per_en_tmp <- data.frame(table(unique(tab_self_regulated[, c("pair", "GENE")])[, c("GENE")]))
          pairs_per_en_tmp$pair_cat <- pair_cat
          pairs_per_en <- rbind(pairs_per_en, pairs_per_en_tmp)

          # annotate driver genes to TCGA pathways ----------------------------------
          # genes2pathways <- map2TCGApathwaways(gene_list = unique(as.vector(tab_self_regulated$GENE)), pathway_list = tcga_pathways)
          # genes2pathways <- data.frame(GENE = rep(x = names(genes2pathways), sapply(X = genes2pathways, FUN = function(x) length(x))), 
          #                              GENE.path = unlist(genes2pathways, use.names = F))
          # tab_self_regulated <- merge(tab_self_regulated, genes2pathways, all.x = T)
          # 
          # # annotate substrate genes to TCGA pathways ----------------------------------
          # sub_genes2pathways <- map2TCGApathwaways(gene_list = unique(tab_self_regulated$SUB_GENE), pathway_list = tcga_pathways_pluskegg)
          # sub_genes2pathways <- data.frame(SUB_GENE = rep(x = names(sub_genes2pathways), sapply(X = sub_genes2pathways, FUN = function(x) length(x))), 
          #                                  SUB_GENE.path = unlist(sub_genes2pathways, use.names = F))
          # tab_self_regulated <- merge(tab_self_regulated, sub_genes2pathways, all.x = T)
          # 
          # ## fill in the NA values with string
          # tmp <- as.vector(tab_self_regulated$SUB_GENE.path)
          # tmp[is.na(tmp)] <- "other_regulated"
          # tab_self_regulated$SUB_GENE.path <- tmp
          # 
          # ## REMOVE replicate substrate annotation
          # tab_self_regulated$same_path <- (as.vector(tab_self_regulated$GENE.path) == as.vector(tab_self_regulated$SUB_GENE.path))
          # tab_self_regulated$pair_pro <- paste0(tab_self_regulated$GENE, ":", tab_self_regulated$SUB_GENE)
          # tab_self_regulated$pair <- paste0(tab_self_regulated$pair_pro, "_", tab_self_regulated$SUB_MOD_RSD)
          # tab_self_regulated$pair_can <- paste0(tab_self_regulated$pair, "-", tab_self_regulated$Cancer)
          # tab_self_regulated <- tab_self_regulated[order(tab_self_regulated$pair_can, tab_self_regulated$same_path, decreasing = T),]
          # tab_self_regulated <- tab_self_regulated[!(duplicated(tab_self_regulated$pair_can)),]
          # 
          # ## reformat the unregulated
          # tab_self <- tab_self_regulated
          # pairs_per_en_per_path <- data.frame(table(unique(tab_self[, c("pair", "GENE", "regulated", "Cancer", "SUB_GENE.path")])[, c("regulated", "GENE","Cancer", "SUB_GENE.path")]))
          # tab4sort_enzyme <- data.frame(table(unique(tab_self[tab_self$SUB_GENE.path != "unregulated" & tab_self$SUB_GENE.path != "other_regulated", c("pair", "GENE", "regulated", "Cancer")])[, c("regulated", "GENE","Cancer")]))
          # 
          # ## count of regulated pairs per enzymef
          # ## get the top enzymes per cancer
          # for (cancer in strsplit(pair_cat, split = "uniq_")[[1]][2]) {
          #   tmp <- tab4sort_enzyme[tab4sort_enzyme$Freq > 0 & tab4sort_enzyme$regulated == "TRUE" & tab4sort_enzyme$Cancer == cancer,]
          #   if (nrow(tmp) > 0) {
          #     tmp <- tmp[order(tmp$Freq, decreasing = T),]
          #     genes_sorted <- as.vector(tmp$GENE[1:min(nrow(tmp), top_enzyme2show)])
          #     
          #     tab2p <- pairs_per_en_per_path[pairs_per_en_per_path$GENE %in% as.vector(tmp$GENE[1:min(nrow(tmp), top_enzyme2show)]) & pairs_per_en_per_path$Cancer == cancer,]
          #     tab2p <- tab2p[tab2p$Freq>0,]
          #     tab2p$color_cat <- as.vector(tab2p$SUB_GENE.path)
          #     tab2p$Cancer <- factor(tab2p$Cancer, levels = cancers2process)
          #     tab2p <- rbind(tab2p[tab2p$regulated == "TRUE" &  tab2p$SUB_GENE.path != "other_regulated",], tab2p[tab2p$regulated == "TRUE" &  tab2p$SUB_GENE.path == "other_regulated",], tab2p[tab2p$regulated == "FALSE",])
          #     tab2p$GENE <- factor(tab2p$GENE, levels = genes_sorted)
          #     ## edit colors
          #     color_manual <- c(tcga_path_colors, color_cancers2[cancer], "grey")
          #     names(color_manual) <- c(names(tcga_path_colors), "other_regulated", "unregulated")
          #     
          #     ## sort by pathways and calculate the y axis
          #     tab2p_sorted <- NULL
          #     for (gene in genes_sorted) {
          #       y <- 0
          #       for (color_cat in names(color_manual)) {
          #         tab2p_tmp <- tab2p[tab2p$GENE == gene & tab2p$color_cat == color_cat,]
          #         if (nrow(tab2p_tmp) > 0){
          #           tab2p_tmp$y <- (0.5*(tab2p_tmp$Freq) + y)
          #           tab2p_tmp$y_stack <- (1*(tab2p_tmp$Freq) + y)
          #           y <- 1*(tab2p_tmp$y_stack)
          #         }
          #         tab2p_sorted <- rbind(tab2p_sorted, tab2p_tmp)
          #       }
          #     }
          #     
          #     tab2p_sorted_regulated <- tab2p_sorted[tab2p_sorted$color_cat != "unregulated",]
          #     tab2p_sorted_regulated_annotated <- tab2p_sorted_regulated[tab2p_sorted_regulated$color_cat != "other_regulated",]
          #     
          #     p <- ggplot()
          #     p <- p + geom_bar(data=tab2p_sorted_regulated, aes(y = Freq, x = GENE, fill = color_cat, group = Cancer),
          #                       stat="identity", position='stack', color = NA)
          #     p <- p + scale_fill_manual(values = color_manual)
          #     p <- p + geom_text(data = tab2p_sorted_regulated, mapping = aes(x = GENE, y = y, label = Freq), nudge_y = nudge_ys[self], size = 2, color = "black")
          #     p
          #     p <- p + theme_minimal() + theme(panel.background = element_rect(fill = brewer.pal(9, "Greys")[2], size = 0))
          #     p <- p + xlab(enzyme_type)+ylab("number of substrate phosphosites")
          #     p <- p + theme(axis.title=element_text(size=10))
          #     p <- p + theme(axis.text.y = element_text(colour="black", size=10), axis.title.y = element_text(size = 6),
          #                    axis.text.x = element_text(size = 10, angle = 90, hjust = 0.5, vjust = 0.5, face = "bold"),  axis.title.x = element_blank(),
          #                    axis.ticks.x = element_blank(),
          #                    panel.spacing.x = unit(0, "lines"), strip.text.x = element_text(size = 9, face = "bold"))
          #     fn = paste0(subdir5, 'substrates_of_', self,'_regulated_', pair_cat, '_top', top_enzyme2show, "_", enzyme_type, '_pairs_in_TCGA_pathways_reg_nonNA',reg_nonNA, '_fdr_thres', fdr_thres, '_wregulated.pdf')
          #     ggsave(file=fn, height=3, width=4)
          #     
          #     p <- ggplot()
          #     p <- p + geom_bar(data=tab2p_sorted_regulated_annotated, aes(y = Freq, x = GENE, fill = color_cat, group = Cancer),
          #                       stat="identity", position='stack', color = NA)
          #     p <- p + scale_fill_manual(values = color_manual)
          #     p <- p + geom_text(data = tab2p_sorted_regulated_annotated,
          #                        mapping = aes(x = GENE, y = y, label = Freq), nudge_y = nudge_ys[self], size = 2, color = "black")
          #     p
          #     p <- p + theme_minimal() + theme(panel.background = element_rect(fill = brewer.pal(9, "Greys")[2], size = 0))
          #     p <- p + xlab(enzyme_type)+ylab("number of substrate phosphosites")
          #     p <- p + theme(axis.title=element_text(size=10))
          #     p <- p + theme(axis.text.y = element_text(colour="black", size=10), axis.title.y = element_text(size = 6),
          #                    axis.text.x = element_text(size = 10, angle = 90, hjust = 0.5, vjust = 0.5, face = "bold"),  axis.title.x = element_blank(),
          #                    axis.ticks.x = element_blank(),
          #                    panel.spacing.x = unit(0, "lines"), strip.text.x = element_text(size = 9, face = "bold"))
          #     fn = paste0(subdir5, 'substrates_of_', self,'_regulated_', pair_cat, '_top', top_enzyme2show, "_", enzyme_type, '_pairs_in_TCGA_pathways_reg_nonNA',reg_nonNA, '_fdr_thres', fdr_thres, '.pdf')
          #     ggsave(file=fn, height=3, width=4)
            }
          }
        # }
      # }

      # 3 cancers together-stack ------------------------------------------------
      for (self in c("trans")) {
        tab2p <- pairs_per_en[pairs_per_en$Freq > 0,]
        colnames(tab2p) <- c("GENE", "Freq", "pair_cat")
        ## only take kinases in top of either of each category
        enzymes2show <- NULL
        for (pair_cat in pair_cat2process) {
          tab_pair_cat <- tab2p[tab2p$pair_cat == pair_cat,]
          tab_pair_cat <- tab_pair_cat[order(tab_pair_cat$Freq, decreasing = T),]
          enzymes2show <- unique(c(enzymes2show, as.vector(tab_pair_cat$GENE[1:min(nrow(tab_pair_cat), top_enzyme2show)])))
        }
        
        ## sort by pair category
        tab2p_sorted <- NULL
        count <- 1
        tab2add <- tab2p[tab2p$GENE %in% enzymes2show,]
        while (nrow(tab2add) > 0 & count <= length(pair_cat2process)) {
          pair_cat <- pair_cat2process[count]
          tab_can <- tab2add[tab2add$pair_cat == pair_cat,]
          if (nrow(tab_can)>0) {
            tab_can <- tab_can[order(tab_can$Freq, decreasing = T),]
            for (gene in tab_can$GENE) {
              for (cancer2 in pair_cat2process) {
                tab_can2 <- tab2add[tab2add$GENE == gene & tab2add$pair_cat == cancer2,]
                if (nrow(tab_can2) > 0) {
                  tab2p_sorted <- rbind(tab2p_sorted, tab_can2)
                }
              }
            }
          }
          tab2add <- tab2add[!(tab2add$GENE %in% tab2p_sorted$GENE),]
          count <- count + 1
        }
        tab2p_sorted$pair_cat <- factor(tab2p_sorted$pair_cat, levels = rev(pair_cat2process))
        tab2p_sorted$GENE <- factor(tab2p_sorted$GENE, levels = unique(as.vector(tab2p_sorted$GENE)))
        
        ## give color to each pair category
        color_pair <- c("#6A3D9A", color_cancers2[cancers2process])
        names(color_pair) <- pair_cat2process
        
        p <- ggplot()
        p <- p + geom_bar(data=tab2p_sorted, aes(y = Freq, x = GENE, fill = pair_cat),
                          stat="identity", position='stack', color = "#000000")
        # p <- p + geom_text(data=tab2p[tab2p$regulated == "TRUE" & tab2p$Freq > 0,],
        #                    aes(y = Freq, x = Cancer, color = color_cat, label = Freq), nudge_y = 0.2)
        p <- p + scale_fill_manual(values = color_pair)
        p <- p + scale_color_manual(values = color_pair)
        p <- p + theme_bw()
        p <- p + theme_nogrid()
        p <- p + xlab(enzyme_type)+ylab("number of substrate phosphosites")
        p <- p + theme(axis.title=element_text(size=10))
        p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5, size = 8, face = "bold"))
        p <- p + theme(axis.text.y = element_text(colour="black", size=10), axis.ticks.x = element_blank(),
                       panel.spacing.x = unit(0, "lines"), strip.text.x = element_text(size = 9, face = "bold"))
        p
        fn = paste0(subdir3, '3can_', self,'_', enzyme_type, '_top_uniq_regulated_cross_cancers.pdf')
        ggsave(file=fn, height=pdf_sizes[[enzyme_type]][[self]]$height, width=pdf_sizes[[enzyme_type]][[self]]$width)
      }
    }
  }
}

# 3 cancers separate-faceted ------------------------------------------------
# for (reg_nonNA in reg_nonNA2test) {
#   for (enzyme_type in c("kinase")) {
#     regression <- fread(input = paste0(ppnD, "regression/tables/", "change_regression_nonNA/", 
#                                        enzyme_type, "_substrate_regression_", cptac_phase2process, "_tumor",
#                                        "_reg_nonNA", reg_nonNA, ".txt"), data.table = F)
#     for (self in c("cis", "trans")) {
#       tab_self <- regression
#       tab_self <- tab_self[tab_self$SELF == self,]
#       sub_gene.path <- map2TCGApathwaways(gene_list = unique(tab_self$SUB_GENE), tcga_pathways = tcga_pathways)
#       tab_self$SUB_GENE.path <- tab_self
#       pairs_per_en <- data.frame(table(unique(tab_self[, c("pair", "GENE", "regulated", "Cancer")])[, c("regulated", "GENE","Cancer")]))
# 
#       ## count of regulated pairs per enzymef
#       ## get the top enzymes per cancer
#       for (cancer in cancers2process) {
#         tmp <- pairs_per_en[pairs_per_en$Freq > 0 & pairs_per_en$regulated == "TRUE" & pairs_per_en$Cancer == cancer,]
#         tmp <- tmp[order(tmp$Freq, decreasing = T),]
#         tab2p <- pairs_per_en[pairs_per_en$GENE %in% as.vector(tmp$GENE[1:min(nrow(tmp), top_enzyme2show)]) & pairs_per_en$Cancer == cancer,]
#         
#         color_cat <- vector(mode = "character", length = nrow(tab2p))
#         color_cat[tab2p$regulated == "TRUE"] <- as.vector(tab2p$Cancer)[tab2p$regulated == "TRUE"]
#         color_cat[tab2p$regulated == "FALSE"]  <- "other"
#         tab2p$color_cat <- color_cat
#         tab2p$Cancer <- factor(tab2p$Cancer, levels = cancers2process)
#         tab2p$color_cat <- factor(tab2p$color_cat, levels = c("other", cancers2process))
#         tab2p <- rbind(tab2p[tab2p$regulated == "TRUE",], tab2p[tab2p$regulated == "FALSE",])
#         tab2p$GENE <- factor(tab2p$GENE, levels = as.vector(tmp$GENE[1:min(nrow(tmp), top_enzyme2show)]))
#         
#         p <- ggplot()
#         p <- p + geom_bar(data=tab2p, aes(y = Freq, x = GENE, fill = color_cat, group = Cancer), 
#                           stat="identity", position='stack', color = NA)
#         p <- p + geom_text(data=tab2p[tab2p$regulated == "TRUE" & tab2p$Freq > 0,], 
#                            aes(y = Freq, x = GENE, label = Freq), nudge_y = nudge_ys[self], size = 2, color = "black")
#         p <- p + scale_fill_manual(values = color_cancers2)
#         p <- p + scale_color_manual(values = color_cancers2)
#         p <- p + theme_minimal() + theme(panel.background = element_rect(fill = brewer.pal(9, "Greys")[2], size = 0))
#         if (self == "cis") {
#           p <- p + scale_y_log10(limits = ylims[[self]])
#         } else {
#           p <- p + scale_y_log10()
#         }
#         p <- p + xlab(enzyme_type)+ylab("number of substrate phosphosites")
#         p <- p + theme(axis.title=element_text(size=10))
#         p <- p + theme(axis.text.y = element_text(colour="black", size=10), axis.title.y = element_text(size = 6),
#                        axis.text.x = element_text(size = 10, angle = 90, hjust = 0.5, vjust = 0.5, face = "bold"),  axis.title.x = element_blank(),
#                        axis.ticks.x = element_blank(),
#                        panel.spacing.x = unit(0, "lines"), strip.text.x = element_text(size = 9, face = "bold"))
#         p
#         fn = paste0(makeOutDir(resultD = resultD), cancer, '_', self,'_', enzyme_type, '_regulated_pairs_top', top_enzyme2show, '_across_cancers.pdf')
#         ggsave(file=fn, height=3, width=4)
#       }
#     }
#     
#   }
# }




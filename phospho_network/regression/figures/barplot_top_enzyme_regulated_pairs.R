# Yige Wu @ WashU 2018 Jul
# barplot for the top kinases with high and low kinase-substrate pairs

# source ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_plotting.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/dependencies/tables/load_TCGA_pathways.R')
# variables ---------------------------------------------------------------
top_num <- 10
nudge_ys <- c(0.3, 0.1); names(nudge_ys) <- c("trans", "cis")
ylims <- list(cis = c(1, 150),
             trans = c(1,10000))
cptac_phase2process <- "cptac3"
cancers2process <- "UCEC"
# cptac_phase2process <- "cptac2p_3can"
# cancers2process <- cancers_sort
pdf_sizes_by_facetx <- list("1" = 5, "3" = 12)
reg_nonNA2test <- c(25)
# inputs -------------------------------------------------------------------

# 3 cancers separate-faceted, sorted by the number of substrates in TCGA pathways------------------------------------------------
for (reg_nonNA in reg_nonNA2test) {
  for (enzyme_type in c("kinase")) {
    regression <- fread(input = paste0(ppnD, "regression/tables/", "change_regression_nonNA/", 
                                       enzyme_type, "_substrate_regression_", cptac_phase2process, "_tumor",
                                       "_reg_nonNA", reg_nonNA, ".txt"), data.table = F)
    for (self in c("cis")) {
      tab_self <- regression
      tab_self <- tab_self[tab_self$SELF == self,]
      tab_self_regulated <- tab_self[tab_self$regulated,]
      sub_gene.path <- unlist(map2TCGApathwaways(gene_list = unique(tab_self_regulated$SUB_GENE), pathway_list = tcga_pathways_pluskegg))
      tab_self_regulated$SUB_GENE.path <- sub_gene.path[as.vector(tab_self_regulated$SUB_GENE)]
      ## fill in the NA values with string
      tmp <- as.vector(tab_self_regulated$SUB_GENE.path)
      tmp[is.na(tmp)] <- "other_regulated"
      tab_self_regulated$SUB_GENE.path <- tmp
      ## reformat the unregulated
      tab_self_unregulated <- tab_self[!tab_self$regulated,]; tab_self_unregulated$SUB_GENE.path <- "unregulated"
      tab_self <- rbind(tab_self_regulated, tab_self_unregulated)
      sum_tab3 <- data.frame(table(unique(tab_self[, c("pair", "GENE", "regulated", "Cancer", "SUB_GENE.path")])[, c("regulated", "GENE","Cancer", "SUB_GENE.path")]))
      tab4sort_enzyme <- data.frame(table(unique(tab_self[tab_self$SUB_GENE.path != "unregulated" & tab_self$SUB_GENE.path != "other_regulated", c("pair", "GENE", "regulated", "Cancer")])[, c("regulated", "GENE","Cancer")]))

      ## count of regulated pairs per enzymef
      ## get the top enzymes per cancer
      for (cancer in cancers2process) {
        tmp <- tab4sort_enzyme[tab4sort_enzyme$Freq > 0 & tab4sort_enzyme$regulated == "TRUE" & tab4sort_enzyme$Cancer == cancer,]
        tmp <- tmp[order(tmp$Freq, decreasing = T),]
        genes_sorted <- as.vector(tmp$GENE[1:min(nrow(tmp), top_num)])
        
        tab2p <- sum_tab3[sum_tab3$GENE %in% as.vector(tmp$GENE[1:min(nrow(tmp), top_num)]) & sum_tab3$Cancer == cancer,]
        tab2p <- tab2p[tab2p$Freq>0,]
        tab2p$color_cat <- as.vector(tab2p$SUB_GENE.path)
        tab2p$Cancer <- factor(tab2p$Cancer, levels = cancers2process)
        tab2p <- rbind(tab2p[tab2p$regulated == "TRUE" &  tab2p$SUB_GENE.path != "other_regulated",], tab2p[tab2p$regulated == "TRUE" &  tab2p$SUB_GENE.path == "other_regulated",], tab2p[tab2p$regulated == "FALSE",])
        tab2p$GENE <- factor(tab2p$GENE, levels = genes_sorted)
        ## edit colors
        color_manual <- c(tcga_path_colors, color_cancers2[cancer], "grey")
        names(color_manual) <- c(names(tcga_path_colors), "other_regulated", "unregulated")
        
        ## sort by pathways and calculate the y axis
        tab2p_sorted <- NULL
        for (gene in genes_sorted) {
          y <- 0
          for (color_cat in names(color_manual)) {
            tab2p_tmp <- tab2p[tab2p$GENE == gene & tab2p$color_cat == color_cat,]
            if (nrow(tab2p_tmp) > 0){
              tab2p_tmp$y <- (0.5*(tab2p_tmp$Freq) + y)
              tab2p_tmp$y_stack <- (1*(tab2p_tmp$Freq) + y)
              y <- 1*(tab2p_tmp$y_stack)
            }
            tab2p_sorted <- rbind(tab2p_sorted, tab2p_tmp)
          }
        }
           
        tab2p_sorted_regulated <- tab2p_sorted[tab2p_sorted$color_cat != "unregulated",]
        tab2p_sorted_regulated_annotated <- tab2p_sorted_regulated[tab2p_sorted_regulated$color_cat != "other_regulated",]
        p <- ggplot()
        p <- p + geom_bar(data=tab2p_sorted_regulated, aes(y = Freq, x = GENE, fill = color_cat, group = Cancer), 
                          stat="identity", position='stack', color = NA)
        p <- p + scale_fill_manual(values = color_manual)
        p <- p + geom_text(data = tab2p_sorted_regulated, mapping = aes(x = GENE, y = y, label = Freq), nudge_y = nudge_ys[self], size = 2, color = "black")
        p
        p <- p + theme_minimal() + theme(panel.background = element_rect(fill = brewer.pal(9, "Greys")[2], size = 0))
        p <- p + xlab(enzyme_type)+ylab("number of substrate phosphosites")
        p <- p + theme(axis.title=element_text(size=10))
        p <- p + theme(axis.text.y = element_text(colour="black", size=10), axis.title.y = element_text(size = 6),
                       axis.text.x = element_text(size = 10, angle = 90, hjust = 0.5, vjust = 0.5, face = "bold"),  axis.title.x = element_blank(),
                       axis.ticks.x = element_blank(),
                       panel.spacing.x = unit(0, "lines"), strip.text.x = element_text(size = 9, face = "bold"))
        fn = paste0(makeOutDir(resultD = resultD), cancer, '_', self,'_', enzyme_type, '_regulated_pairs_top', top_num, '_tcga_pathways_across_cancers_wregulated.pdf')
        ggsave(file=fn, height=3, width=4)
        
        p <- ggplot()
        p <- p + geom_bar(data=tab2p_sorted_regulated_annotated, aes(y = Freq, x = GENE, fill = color_cat, group = Cancer), 
                          stat="identity", position='stack', color = NA)
        p <- p + scale_fill_manual(values = color_manual)
        p <- p + geom_text(data = tab2p_sorted_regulated_annotated,
                           mapping = aes(x = GENE, y = y, label = Freq), nudge_y = nudge_ys[self], size = 2, color = "black")
        p
        p <- p + theme_minimal() + theme(panel.background = element_rect(fill = brewer.pal(9, "Greys")[2], size = 0))
        p <- p + xlab(enzyme_type)+ylab("number of substrate phosphosites")
        p <- p + theme(axis.title=element_text(size=10))
        p <- p + theme(axis.text.y = element_text(colour="black", size=10), axis.title.y = element_text(size = 6),
                       axis.text.x = element_text(size = 10, angle = 90, hjust = 0.5, vjust = 0.5, face = "bold"),  axis.title.x = element_blank(),
                       axis.ticks.x = element_blank(),
                       panel.spacing.x = unit(0, "lines"), strip.text.x = element_text(size = 9, face = "bold"))
        fn = paste0(makeOutDir(resultD = resultD), cancer, '_', self,'_', enzyme_type, '_regulated_pairs_top', top_num, '_tcga_pathways_across_cancers_reg_nonNA', reg_nonNA, '.pdf')
        ggsave(file=fn, height=3, width=4)
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
#       sum_tab3 <- data.frame(table(unique(tab_self[, c("pair", "GENE", "regulated", "Cancer")])[, c("regulated", "GENE","Cancer")]))
# 
#       ## count of regulated pairs per enzymef
#       ## get the top enzymes per cancer
#       for (cancer in cancers2process) {
#         tmp <- sum_tab3[sum_tab3$Freq > 0 & sum_tab3$regulated == "TRUE" & sum_tab3$Cancer == cancer,]
#         tmp <- tmp[order(tmp$Freq, decreasing = T),]
#         tab2p <- sum_tab3[sum_tab3$GENE %in% as.vector(tmp$GENE[1:min(nrow(tmp), top_num)]) & sum_tab3$Cancer == cancer,]
#         
#         color_cat <- vector(mode = "character", length = nrow(tab2p))
#         color_cat[tab2p$regulated == "TRUE"] <- as.vector(tab2p$Cancer)[tab2p$regulated == "TRUE"]
#         color_cat[tab2p$regulated == "FALSE"]  <- "other"
#         tab2p$color_cat <- color_cat
#         tab2p$Cancer <- factor(tab2p$Cancer, levels = cancers2process)
#         tab2p$color_cat <- factor(tab2p$color_cat, levels = c("other", cancers2process))
#         tab2p <- rbind(tab2p[tab2p$regulated == "TRUE",], tab2p[tab2p$regulated == "FALSE",])
#         tab2p$GENE <- factor(tab2p$GENE, levels = as.vector(tmp$GENE[1:min(nrow(tmp), top_num)]))
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
#         fn = paste0(makeOutDir(resultD = resultD), cancer, '_', self,'_', enzyme_type, '_regulated_pairs_top', top_num, '_across_cancers.pdf')
#         ggsave(file=fn, height=3, width=4)
#       }
#     }
#     
#   }
# }



# 3 cancers together-dodge ------------------------------------------------
# for (self in c("trans")) {
#   sup_cans_tab_en_self <- sup_cans_tab_en[sup_cans_tab_en$SELF == self,]
#   tab_self <- sup_cans_tab_en_self
#   sum_tab1 <- data.frame(table(unique(tab_self[, c("pair", "GENE", "Cancer")])[, c("GENE")]))
#   colnames(sum_tab1) <- c("GENE", "num_pairs_cans")
#   sum_tab2 <- data.frame(table(unique(tab_self[tab_self$regulated, c("pair", "GENE", "Cancer")])[, c("GENE")]))
#   sum_tab2 <- merge(sum_tab2, sum_tab1, all.x = T)
#   sum_tab2$ratio <- sum_tab2$Freq/sum_tab2$num_pairs_cans
#   sum_tab2 <- sum_tab2[order(sum_tab2$ratio, decreasing = T),]
#   sum_tab3 <- data.frame(table(unique(tab_self[, c("pair", "GENE", "regulated", "Cancer")])[, c("regulated", "GENE","Cancer")]))
#   sum_tab4 <- data.frame(table(unique(tab_self[ c("pair", "GENE", "Cancer")])[, c("GENE","Cancer")]))
#   colnames(sum_tab4) <- c("GENE", "Cancer", "num_pairs_can")
#   sum_tab3 <- merge(sum_tab3, sum_tab4, all.x = T)
#   sum_tab3$ratio_can <- sum_tab3$Freq/sum_tab3$num_pairs_can
#   
#   tab2p <- sum_tab3[sum_tab3$Freq > 0 & sum_tab3$regulated == "TRUE",]
#   tab2p_count2 <- data.frame(table(tab2p$GENE))
#   # tab2p <- tab2p[tab2p$GENE %in% tab2p_count2$Var1[tab2p_count2$Freq == length(cancers2process)],]
#   tab2p_count <- group_by(tab2p, GENE)
#   tab2p_count <- data.frame(summarise(tab2p_count, ave_ratio = sum(ratio_can)/3))
#   tab2p <- sum_tab3[sum_tab3$GENE %in% tab2p$GENE ,]
#   tab2p_count <- tab2p_count[order(tab2p_count$ave_ratio, decreasing = T),]
#   tab2p <- tab2p[tab2p$GENE %in% tab2p_count$GENE[1:20],]
#   
#   color_cat <- vector(mode = "character", length = nrow(tab2p))
#   color_cat[tab2p$regulated == "TRUE"] <- as.vector(tab2p$Cancer)[tab2p$regulated == "TRUE"]
#   color_cat[tab2p$regulated == "FALSE"]  <- "other"
#   tab2p$color_cat <- color_cat
#   tab2p$Cancer <- factor(tab2p$Cancer, levels = cancers2process)
#   tab2p$color_cat <- factor(tab2p$color_cat, levels = c("other", cancers2process))
#   tab2p <- rbind(tab2p[tab2p$regulated == "TRUE",], tab2p[tab2p$regulated == "FALSE",])
#   tab2p$GENE <- factor(tab2p$GENE, levels = as.vector(tab2p_count$GENE)[order(tab2p_count$ave_ratio, decreasing = T)])
#   p <- ggplot()
#   p <- p + geom_bar(data=tab2p, aes(y = Freq, x = Cancer, fill = color_cat, group = Cancer), 
#                     stat="identity", position='stack', color = "#000000")
#   p <- p + facet_grid(.~GENE, drop=T, space = "free",scales = "free")#, space = "free", scales = "free")
#   p <- p + geom_text(data=tab2p[tab2p$regulated == "TRUE" & tab2p$Freq > 0,], 
#                      aes(y = Freq, x = Cancer, color = color_cat, label = Freq), nudge_y = 0.2)
#   p <- p + scale_fill_manual(values = color_cat_man)
#   p <- p + scale_color_manual(values = color_cat_man)
#   p <- p + theme_bw()
#   p <- p + theme_nogrid()
#   p <- p + scale_y_log10()
#   p <- p + xlab(enzyme_type)+ylab("number of substrate phosphosites")
#   p <- p + theme(axis.title=element_text(size=10))
#   p <- p + theme(axis.text.y = element_text(colour="black", size=10), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
#                   panel.spacing.x = unit(0, "lines"), strip.text.x = element_text(size = 9, face = "bold"))
#   p
#   fn = paste0(makeOutDir(resultD = resultD),'3can_', self,'_', enzyme_type, '_top_regulated_average_ratio_across_cancers.ptab2p')
#   ggsave(file=fn, height=3, width=8)
# }
# 

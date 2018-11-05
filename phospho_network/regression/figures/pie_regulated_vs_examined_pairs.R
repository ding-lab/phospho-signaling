# Yige Wu @ WashU 2018 Aug
# pie chart showing proportion of regulated pairs by each substrate by each enzyme, also showing regulated ratio
# source ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_plotting.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
if (!require("rPython")) install.packages("rPython")
if (!require("ggsunburst")) install.packages("http://genome.crg.es/~didac/ggsunburst/ggsunburst_0.0.9.tar.gz", repos=NULL, type="source")
library(ggsunburst)

# set variables -----------------------------------------------------------
reg_nonNA2test <- c(20, 25)
top_kinase2show <- 10
top_substrate2show <- 2
fdr_thres <- c(0.2, 0.05); names(fdr_thres) <- c("phosphatase", "kinase")
# inputs ------------------------------------------------------------------

# pie chart of top regulating enzymes within all regulated pairs-------------------------------------
for (reg_nonNA in reg_nonNA2test) {
  for (enzyme_type in c("kinase")) {
    sup_cans_tab_en <- fread(input = paste0(ppnD, "regression/tables/", "change_regression_nonNA/", 
                                            enzyme_type, "_substrate_regression_cptac2p_3can_tumor",
                                            "_reg_nonNA", reg_nonNA, ".txt"), data.table = F)
    
    sup_cans_tab_en_reg <- sup_cans_tab_en[sup_cans_tab_en$regulated,]
    for (cancer in unique(sup_cans_tab_en_reg$Cancer)) {
      # cancer <- "BRCA"
      driver_can <- loadGeneList(gene_type = "driver", cancer = cancer, is.soft.limit = "")
      driver_pancan <- loadGeneList(gene_type = "driver", cancer = "PANCAN", is.soft.limit = "")
      for (self in c("trans")) {
        # self <- "trans"
        tab_can <- sup_cans_tab_en_reg[sup_cans_tab_en_reg$Cancer == cancer & sup_cans_tab_en_reg$SELF == self,]
        if (nrow(tab_can) > 0) {
          tab_can$pair_pro <- paste0(tab_can$GENE, ":", tab_can$SUB_GENE)
          tab4orderKinase <- data.frame(table(tab_can[, c("GENE")]))
          tab4orderKinase <- tab4orderKinase[order(tab4orderKinase$Freq, decreasing = T),]
          tab4orderKinase$kinase_order <- 1:nrow(tab4orderKinase)
          top_kinases <- tab4orderKinase$Var1[1:min(nrow(tab4orderKinase), top_kinase2show)]
          
          tab4sunburst <- data.frame(table(tab_can[, c("pair_pro")]))
          tab4sunburst <- tab4sunburst[tab4sunburst$Freq > 0,]
          tab4sunburst$kinase <- str_split_fixed(string = tab4sunburst$Var1, pattern = ":", 2)[,1]
          tab4sunburst$substrate <- str_split_fixed(string = tab4sunburst$Var1, pattern = ":", 2)[,2]
          colnames(tab4sunburst)<- c("node", "size", "kinase", "substrate")
          tab4sunburst$parent <- tab4sunburst$kinase
          
          ## keep top substrates, condense others
          tab4sunburst <- tab4sunburst[order(tab4sunburst$kinase, tab4sunburst$size, decreasing = T),]
          tab4sunburst$substrate_order <- unlist(sapply(tab4sunburst$kinase[!duplicated(tab4sunburst$kinase)], FUN = function(kinase, tab) 1:length(which(tab$kinase == kinase)), tab = tab4sunburst))
          # substrate_show_condition <- (tab4sunburst$substrate_order <= top_substrate2show & tab4sunburst$size >= 5)
          if (self == "trans") {
            if (enzyme_type == "kinase") {
              substrate_show_condition <- (tab4sunburst$substrate_order <= top_substrate2show & tab4sunburst$size >= 5) | ((tab4sunburst$substrate %in% driver_can) & tab4sunburst$size >= 5)
              
            }
            if (enzyme_type == "phosphatase") {
              substrate_show_condition <- (tab4sunburst$size >= 1)
              
            }
          }
          if (self == "cis") {
            substrate_show_condition <- rep(x = c(TRUE), nrow(tab4sunburst))
          }
          
          tab4sunburst_sub_top <- tab4sunburst[substrate_show_condition,]
          tab4sunburst_sub_other <- tab4sunburst[!(substrate_show_condition),]
          tab4sunburst_short <- tab4sunburst
          
          ## add kinase order
          tab4sunburst_short <- merge(tab4sunburst_short, tab4orderKinase[, c("Var1", "kinase_order")], 
                                      by.x = c("kinase"), by.y = c("Var1"), all.x =T)
          tab4sunburst_short <- tab4sunburst_short[order(tab4sunburst_short$kinase_order, tab4sunburst_short$substrate_order),]
          
          ## write to file
          file <- paste0(makeOutDir(resultD = resultD), 'tab4sunburst.csv')
          write.table(tab4sunburst_short, file = file, row.names = F, sep = ",")
          sb <- sunburst_data(input = file, type = 'node_parent', sep = ",", node_attributes = c("kinase","size", "substrate"))
          
          sb$rects[!sb$rects$leaf,]$kinase <- sb$rects[!sb$rects$leaf,]$name
          ## make not top kinase not showing kinase name
          tmp <- as.vector(sb$node_labels$label)
          tmp[!(tmp %in% top_kinases)] <- ""
          sb$node_labels$label <- tmp
          ## all substrate white and not top kinae grey
          tmp <- as.vector(sb$rects$kinase)
          tmp[!(sb$rects$kinase %in% tab4orderKinase$Var1[1:min(nrow(tab4orderKinase), top_kinase2show)])] <- "other"
          tmp[!is.na(sb$rects$substrate)] <- "substrate"
          sb$rects$fill <- tmp
          ## make not showing substrate name
          sb$leaf_labels$label <- ""
          
          color4enzyme <- c(rainbow(n = length(top_kinases)), "grey", "white")
          names(color4enzyme) <- c(as.vector(tab4orderKinase$Var1[1:min(nrow(tab4orderKinase), top_kinase2show)]), "other", "substrate")
          
          p <- sunburst(sb, rects.fill.aes = "fill", leaf_labels = F,
                        node_labels = T, node_labels.min = 1, leaf_labels.size = 2, node_labels.size = 2.5)
          p <- p + scale_fill_manual(values = color4enzyme)
          p <- p + theme(legend.position = "none")
          p
          ggsave(filename = paste0(makeOutDir(resultD = resultD), 
                                   cancer, "_", enzyme_type, "substrate_pairs_", self, "_regulated_reg_nonNA", reg_nonNA, 
                                   "_AllEnzyme_showtop", top_kinase2show,  ".pdf"), 
                 width = 6, height = 6)
        }
        
      }
    }
  }
}

# pie chart of each top regulating enzyme -------------------------------------
for (reg_nonNA in reg_nonNA2test) {
  for (enzyme_type in c("phosphatase")) {
    sup_cans_tab_en <- fread(input = paste0(ppnD, "regression/tables/", "change_regression_nonNA/", 
                                            enzyme_type, "_substrate_regression_cptac2p_3can_tumor",
                                            "_reg_nonNA", reg_nonNA, ".txt"), data.table = F)
    sup_cans_tab_en <- markSigKS(regression = sup_cans_tab_en, sig_thres = fdr_thres[enzyme_type], enzyme_type = enzyme_type)
    sup_cans_tab_en$regulated <- (sup_cans_tab_en$fdr_sig & sup_cans_tab_en$coef_sig)
    sup_cans_tab_en_reg <- sup_cans_tab_en[sup_cans_tab_en$regulated,]
    for (cancer in unique(sup_cans_tab_en_reg$Cancer)) {
      # cancer <- "BRCA"
      driver_can <- loadGeneList(gene_type = "driver", cancer = cancer, is.soft.limit = "")
      driver_pancan <- loadGeneList(gene_type = "driver", cancer = "PANCAN", is.soft.limit = "")
      for (self in c("trans")) {
        # self <- "trans"
        tab_can <- sup_cans_tab_en_reg[sup_cans_tab_en_reg$Cancer == cancer & sup_cans_tab_en_reg$SELF == self,]
        if (nrow(tab_can) > 0) {
          tab_can$pair_pro <- paste0(tab_can$GENE, ":", tab_can$SUB_GENE)
          tab4orderKinase <- data.frame(table(tab_can[, c("GENE")]))
          tab4orderKinase <- tab4orderKinase[order(tab4orderKinase$Freq, decreasing = T),]
          tab4orderKinase$kinase_order <- 1:nrow(tab4orderKinase)
          
          tab4sunburst <- data.frame(table(tab_can[, c("pair_pro")]))
          tab4sunburst <- tab4sunburst[tab4sunburst$Freq > 0,]
          tab4sunburst$kinase <- str_split_fixed(string = tab4sunburst$Var1, pattern = ":", 2)[,1]
          tab4sunburst$substrate <- str_split_fixed(string = tab4sunburst$Var1, pattern = ":", 2)[,2]
          colnames(tab4sunburst)<- c("node", "size", "kinase", "substrate")
          tab4sunburst$parent <- tab4sunburst$kinase
          
          ## keep top kinases
          tab4sunburst <- tab4sunburst[tab4sunburst$kinase %in% tab4orderKinase$Var1[1:min(nrow(tab4orderKinase), top_kinase2show)],]
          
          ## keep top substrates, condense others
          tab4sunburst <- tab4sunburst[order(tab4sunburst$kinase, tab4sunburst$size, decreasing = T),]
          tab4sunburst$substrate_order <- unlist(sapply(tab4sunburst$kinase[!duplicated(tab4sunburst$kinase)], FUN = function(kinase, tab) 1:length(which(tab$kinase == kinase)), tab = tab4sunburst))
          # substrate_show_condition <- (tab4sunburst$substrate_order <= top_substrate2show & tab4sunburst$size >= 5)
          if (self == "trans") {
            if (enzyme_type == "kinase") {
              substrate_show_condition <- (tab4sunburst$substrate_order <= top_substrate2show & tab4sunburst$size >= 5) | ((tab4sunburst$substrate %in% driver_can) & tab4sunburst$size >= 5)
              
            }
            if (enzyme_type == "phosphatase") {
              substrate_show_condition <- (tab4sunburst$size >= 1)
              
            }
          }
          if (self == "cis") {
            substrate_show_condition <- rep(x = c(TRUE), nrow(tab4sunburst))
          }
          
          tab4sunburst_sub_top <- tab4sunburst[substrate_show_condition,]
          tab4sunburst_sub_other <- tab4sunburst[!(substrate_show_condition),]
          tab4sunburst_short <- tab4sunburst
          
          ## add kinase order
          tab4sunburst_short <- merge(tab4sunburst_short, tab4orderKinase[, c("Var1", "kinase_order")], 
                                      by.x = c("kinase"), by.y = c("Var1"), all.x =T)
          tab4sunburst_short <- tab4sunburst_short[order(tab4sunburst_short$kinase_order, tab4sunburst_short$substrate_order),]
          
          ## write to file
          file <- paste0(makeOutDir(resultD = resultD), 'tab4sunburst.csv')
          write.table(tab4sunburst_short, file = file, row.names = F, sep = ",")
          sb <- sunburst_data(input = file, type = 'node_parent', sep = ",", node_attributes = c("kinase","size", "substrate"))
          sb$rects[!sb$rects$leaf,]$kinase <- sb$rects[!sb$rects$leaf,]$name
          tmp <- as.vector(sb$rects$kinase)
          tmp[!is.na(sb$rects$substrate) & sb$rects$substrate != "other"] <- "substrate"
          tmp[!is.na(sb$rects$substrate) & sb$rects$substrate == "other"] <- "other"
          sb$rects$fill <- tmp
          ## make not top substrate not showing name
          tmp <- as.vector(sb$leaf_labels$substrate)
          tmp[!(tmp %in% tab4sunburst_sub_top$substrate)] <- ""
          sb$leaf_labels$label <- tmp
          
          color4enzyme <- c(rainbow(n = length(unique(tab4sunburst_short$kinase))), "grey", "white")
          names(color4enzyme) <- c(as.vector(unique(tab4sunburst_short$kinase)), "other", "substrate")
          
          p <- sunburst(sb, rects.fill.aes = "fill",
                        node_labels = T, node_labels.min = 1, leaf_labels.size = 2, node_labels.size = 2.5)
          if (self == "cis") {
            p <- p + geom_text(data = sb$leaf_labels,
                               aes(x=x, y=0.1, label=paste(size), angle=angle, hjust=hjust), size = 5)
          } else {
            p <- p + geom_text(data = sb$leaf_labels,
                               aes(x=x, y=0.1, label=paste(size), angle=angle, hjust=hjust), size = 2)
            
          }
          p <- p + scale_fill_manual(values = color4enzyme)
          p <- p + theme(legend.position = "none")
          p
          ggsave(filename = paste0(makeOutDir(resultD = resultD), 
                                   cancer, "_", enzyme_type, "substrate_pairs_", self, "_regulated_reg_nonNA", reg_nonNA, 
                                   "_topenzyme", top_kinase2show, "_topsubstrate", top_substrate2show, ".pdf"), 
                 width = 6, height = 6)
        }

      }
    }
  }
}

# pie chart of pairs involving driver gene substrate -------------------------------------
top_kinase2show <- 15
for (reg_nonNA in reg_nonNA2test) {
  for (enzyme_type in c("phosphatase")) {
    sup_cans_tab_en <- fread(input = paste0(ppnD, "regression/tables/", "change_regression_nonNA/", 
                                            enzyme_type, "_substrate_regression_cptac2p_3can_tumor",
                                            "_reg_nonNA", reg_nonNA, ".txt"), data.table = F)
    sup_cans_tab_en <- markSigKS(regression = sup_cans_tab_en, sig_thres = fdr_thres[enzyme_type], enzyme_type = enzyme_type)
    sup_cans_tab_en$regulated <- (sup_cans_tab_en$fdr_sig & sup_cans_tab_en$coef_sig)
    sup_cans_tab_en_reg <- sup_cans_tab_en[sup_cans_tab_en$regulated,]
    for (cancer in unique(sup_cans_tab_en_reg$Cancer)) {
      # cancer <- "BRCA"
      driver_can <- loadGeneList(gene_type = "driver", cancer = cancer, is.soft.limit = "")
      driver_pancan <- loadGeneList(gene_type = "driver", cancer = "PANCAN", is.soft.limit = "")
      for (self in c("trans")) {
        # self <- "trans"
        tab_can <- sup_cans_tab_en_reg[sup_cans_tab_en_reg$Cancer == cancer & sup_cans_tab_en_reg$SELF == self,]
        tab_can <- tab_can[tab_can$SUB_GENE %in% driver_can,]
        if (nrow(tab_can) > 0) {
          tab_can$pair_pro <- paste0(tab_can$GENE, ":", tab_can$SUB_GENE)
          tab4orderKinase <- data.frame(table(tab_can[, c("GENE")]))
          tab4orderKinase <- tab4orderKinase[order(tab4orderKinase$Freq, decreasing = T),]
          tab4orderKinase$kinase_order <- 1:nrow(tab4orderKinase)
          
          tab4sunburst <- data.frame(table(tab_can[, c("pair_pro")]))
          tab4sunburst <- tab4sunburst[tab4sunburst$Freq > 0,]
          tab4sunburst$kinase <- str_split_fixed(string = tab4sunburst$Var1, pattern = ":", 2)[,1]
          tab4sunburst$substrate <- str_split_fixed(string = tab4sunburst$Var1, pattern = ":", 2)[,2]
          colnames(tab4sunburst)<- c("node", "size", "kinase", "substrate")
          tab4sunburst$parent <- tab4sunburst$kinase
          
          ## keep top kinases
          tab4sunburst <- tab4sunburst[tab4sunburst$kinase %in% tab4orderKinase$Var1[1:min(nrow(tab4orderKinase), top_kinase2show)],]
          
          ## keep top substrates, condense others
          tab4sunburst <- tab4sunburst[order(tab4sunburst$kinase, tab4sunburst$size, decreasing = T),]
          tab4sunburst$substrate_order <- unlist(sapply(tab4sunburst$kinase[!duplicated(tab4sunburst$kinase)], FUN = function(kinase, tab) 1:length(which(tab$kinase == kinase)), tab = tab4sunburst))
          # substrate_show_condition <- (tab4sunburst$substrate_order <= top_substrate2show & tab4sunburst$size >= 5)
          if (self == "trans") {
            if (enzyme_type == "kinase") {
              substrate_show_condition <- (tab4sunburst$substrate_order <= top_substrate2show & tab4sunburst$size >= 5) | ((tab4sunburst$substrate %in% driver_can) & tab4sunburst$size >= 5)
              
            }
            if (enzyme_type == "phosphatase") {
              substrate_show_condition <- (tab4sunburst$size >= 1)
              
            }
          }
          if (self == "cis") {
            substrate_show_condition <- rep(x = c(TRUE), nrow(tab4sunburst))
          }
          
          tab4sunburst_sub_top <- tab4sunburst[substrate_show_condition,]
          tab4sunburst_sub_other <- tab4sunburst[!(substrate_show_condition),]
          tab4sunburst_short <- tab4sunburst
          
          ## add kinase order
          tab4sunburst_short <- merge(tab4sunburst_short, tab4orderKinase[, c("Var1", "kinase_order")], 
                                      by.x = c("kinase"), by.y = c("Var1"), all.x =T)
          tab4sunburst_short <- tab4sunburst_short[order(tab4sunburst_short$kinase_order, tab4sunburst_short$substrate_order),]
          
          ## write to file
          file <- paste0(makeOutDir(resultD = resultD), 'tab4sunburst.csv')
          write.table(tab4sunburst_short, file = file, row.names = F, sep = ",")
          sb <- sunburst_data(input = file, type = 'node_parent', sep = ",", node_attributes = c("kinase","size", "substrate"))
          sb$rects[!sb$rects$leaf,]$kinase <- sb$rects[!sb$rects$leaf,]$name
          tmp <- as.vector(sb$rects$kinase)
          tmp[!is.na(sb$rects$substrate) & sb$rects$substrate != "other"] <- "substrate"
          tmp[!is.na(sb$rects$substrate) & sb$rects$substrate == "other"] <- "other"
          sb$rects$fill <- tmp

          color4enzyme <- c(rainbow(n = length(unique(tab4sunburst_short$kinase))), "grey", "white")
          names(color4enzyme) <- c(as.vector(unique(tab4sunburst_short$kinase)), "other", "substrate")
          
          if (self == "cis") {
            p <- sunburst(sb, rects.fill.aes = "fill",
                          node_labels = T, node_labels.min = 1, leaf_labels.size = 2, node_labels.size = 8)
            p <- p + geom_text(data = sb$leaf_labels,
                               aes(x=x, y=0.1, label=paste(size), angle=angle, hjust=hjust), size = 10)
          } else {
            p <- sunburst(sb, rects.fill.aes = "fill",
                          node_labels = T, node_labels.min = 1, leaf_labels.size = 2, node_labels.size = 2.5)
            p <- p + geom_text(data = sb$leaf_labels,
                               aes(x=x, y=0.1, label=paste(size), angle=angle, hjust=hjust), size = 2)
          }

          p <- p + scale_fill_manual(values = color4enzyme)
          p <- p + theme(legend.position = "none")
          p
          ggsave(filename = paste0(makeOutDir(resultD = resultD), 
                                   cancer, "_", enzyme_type, "_driver_substrate_pairs_", self, "_regulated_reg_nonNA", reg_nonNA, 
                                   "_topenzyme", top_kinase2show, "_topsubstrate", top_substrate2show, ".pdf"), 
                 width = 6, height = 6)
        }
        
      }
    }
  }
}

# pie chart of pairs involving driver gene enzyme -------------------------------------
for (reg_nonNA in reg_nonNA2test) {
  for (enzyme_type in c("kinase", "phosphatase")) {
    sup_cans_tab_en <- fread(input = paste0(ppnD, "regression/tables/", "change_regression_nonNA/", 
                                            enzyme_type, "_substrate_regression_cptac2p_3can_tumor",
                                            "_reg_nonNA", reg_nonNA, ".txt"), data.table = F)
    sup_cans_tab_en <- markSigKS(regression = sup_cans_tab_en, sig_thres = fdr_thres[enzyme_type], enzyme_type = enzyme_type)
    sup_cans_tab_en$regulated <- (sup_cans_tab_en$fdr_sig & sup_cans_tab_en$coef_sig)
    sup_cans_tab_en_reg <- sup_cans_tab_en[sup_cans_tab_en$regulated,]
    for (cancer in unique(sup_cans_tab_en_reg$Cancer)) {
      # cancer <- "BRCA"
      driver_can <- loadGeneList(gene_type = "driver", cancer = cancer, is.soft.limit = "")
      driver_pancan <- loadGeneList(gene_type = "driver", cancer = "PANCAN", is.soft.limit = "")
      for (self in c("trans")) {
        # self <- "trans"
        tab_can <- sup_cans_tab_en_reg[sup_cans_tab_en_reg$Cancer == cancer & sup_cans_tab_en_reg$SELF == self,]
        tab_can <- tab_can[tab_can$GENE %in% driver_can,]
        # tab_can <- tab_can[tab_can$GENE %in% c(driver_can, driver_pancan),]
        
        if (nrow(tab_can) > 0) {
          tab_can$pair_pro <- paste0(tab_can$GENE, ":", tab_can$SUB_GENE)
          tab4orderKinase <- data.frame(table(tab_can[, c("GENE")]))
          tab4orderKinase <- tab4orderKinase[order(tab4orderKinase$Freq, decreasing = T),]
          tab4orderKinase$kinase_order <- 1:nrow(tab4orderKinase)
          
          tab4sunburst <- data.frame(table(tab_can[, c("pair_pro")]))
          tab4sunburst <- tab4sunburst[tab4sunburst$Freq > 0,]
          tab4sunburst$kinase <- str_split_fixed(string = tab4sunburst$Var1, pattern = ":", 2)[,1]
          tab4sunburst$substrate <- str_split_fixed(string = tab4sunburst$Var1, pattern = ":", 2)[,2]
          colnames(tab4sunburst)<- c("node", "size", "kinase", "substrate")
          tab4sunburst$parent <- tab4sunburst$kinase
          
          ## keep top kinases
          tab4sunburst <- tab4sunburst[tab4sunburst$kinase %in% tab4orderKinase$Var1[1:min(nrow(tab4orderKinase), top_kinase2show)],]
          
          ## keep top substrates, condense others
          tab4sunburst <- tab4sunburst[order(tab4sunburst$kinase, tab4sunburst$size, decreasing = T),]
          tab4sunburst$substrate_order <- unlist(sapply(tab4sunburst$kinase[!duplicated(tab4sunburst$kinase)], FUN = function(kinase, tab) 1:length(which(tab$kinase == kinase)), tab = tab4sunburst))
          # substrate_show_condition <- (tab4sunburst$substrate_order <= top_substrate2show & tab4sunburst$size >= 5)
          if (self == "trans") {
            if (enzyme_type == "kinase") {
              substrate_show_condition <- (tab4sunburst$substrate_order <= top_substrate2show & tab4sunburst$size >= 5) | ((tab4sunburst$substrate %in% driver_can) & tab4sunburst$size >= 5)
              
            }
            if (enzyme_type == "phosphatase") {
              substrate_show_condition <- (tab4sunburst$size >= 1)
              
            }
          }
          if (self == "cis") {
            substrate_show_condition <- rep(x = c(TRUE), nrow(tab4sunburst))
          }
          
          tab4sunburst_sub_top <- tab4sunburst[substrate_show_condition,]
          tab4sunburst_sub_other <- tab4sunburst[!(substrate_show_condition),]
          tab4sunburst_short <- tab4sunburst
          
          ## add kinase order
          tab4sunburst_short <- merge(tab4sunburst_short, tab4orderKinase[, c("Var1", "kinase_order")], 
                                      by.x = c("kinase"), by.y = c("Var1"), all.x =T)
          tab4sunburst_short <- tab4sunburst_short[order(tab4sunburst_short$kinase_order, tab4sunburst_short$substrate_order),]
          
          ## write to file
          file <- paste0(makeOutDir(resultD = resultD), 'tab4sunburst.csv')
          write.table(tab4sunburst_short, file = file, row.names = F, sep = ",")
          sb <- sunburst_data(input = file, type = 'node_parent', sep = ",", node_attributes = c("kinase","size", "substrate"))
          sb$rects[!sb$rects$leaf,]$kinase <- sb$rects[!sb$rects$leaf,]$name
          tmp <- as.vector(sb$rects$kinase)
          tmp[!is.na(sb$rects$substrate) & sb$rects$substrate != "other"] <- "substrate"
          tmp[!is.na(sb$rects$substrate) & sb$rects$substrate == "other"] <- "other"
          sb$rects$fill <- tmp
          
          color4enzyme <- c(rainbow(n = length(unique(tab4sunburst_short$kinase))), "grey", "white")
          names(color4enzyme) <- c(as.vector(unique(tab4sunburst_short$kinase)), "other", "substrate")
          
          p <- sunburst(sb, rects.fill.aes = "fill",
                        node_labels = T, node_labels.min = 1, leaf_labels.size = 2, node_labels.size = 2.5)
          p <- p + geom_text(data = sb$leaf_labels,
                             aes(x=x, y=0.1, label=paste(size), angle=angle, hjust=hjust), size = 2)
          p <- p + scale_fill_manual(values = color4enzyme)
          p <- p + theme(legend.position = "none")
          p
          ggsave(filename = paste0(makeOutDir(resultD = resultD), 
                                   cancer, "_driver_", enzyme_type, "_substrate_pairs_", self, "_regulated_reg_nonNA", reg_nonNA, 
                                   "_topenzyme", top_kinase2show, "_topsubstrate", top_substrate2show, ".pdf"), 
                 width = 5, height = 6)
        }
        
      }
    }
  }
}


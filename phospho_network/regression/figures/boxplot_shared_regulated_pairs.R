# Yige Wu @ WashU 2018 Jul
# barplot for the top kinases with high and low kinase-substrate pairs

# source ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/pan3can_aes.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_plotting.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
library(dplyr)
library(ggrepel)

# variables ---------------------------------------------------------------
cancers_sort <- c("BRCA", "OV", "CO")
reg_sig <- 0.05
diff_sig <- 0.2
diff_log2fc <- 1
diff_log2fc2p <- 2
color_cat_man <- c(colors['BRCA'], colors['OV'], colors['COAD'], "#bdbdbd"); names(color_cat_man) <- c("BRCA", "OV", "CO", "other")

# inputs -------------------------------------------------------------------
## input driver genes
driver_genes <- read.delim(file = "./TCGA_data/reference_files/Consensus.Genelist.full.txt")

## input regression table
sup_cans_tab <- NULL
for (cancer in cancers_sort) {
  sup_tab <- fread(paste0(ppnD, "kinase_activity/tables/fisher_es_pairs/", 
                          cancer, "_KSEA_enzyme_diffexp_substrate.txt"), data.table = F)
  sup_tab$Cancer <- cancer
  sup_cans_tab <- rbind(sup_cans_tab, sup_tab)
}

# summarize kinase-substrate pairs ----------------------------------------
enzyme_type <- "kinase"
sup_cans_tab_pk <- sup_cans_tab[sup_cans_tab$enzyme_type == enzyme_type,]
sup_cans_tab_pk <- data.frame(sup_cans_tab_pk)
sup_cans_tab_pk <- markSigSiteCan(sup_cans_tab_pk, sig_thres = reg_sig, enzyme_type = enzyme_type)
## annotate kinase substrate regulation
sup_cans_tab_pk$regulated <- (sup_cans_tab_pk$coef_sig & sup_cans_tab_pk$fdr_sig)



# plot trans coef all--------------------------------------------------------------
for (self in c("trans")) {
  sup_cans_tab_pk_self <- sup_cans_tab_pk[sup_cans_tab_pk$SELF == self,]
  tab_tmp <- sup_cans_tab_pk_self
  sum_tab1 <- data.frame(table(unique(tab_tmp[, c("pair", "GENE", "Cancer")])[, c("GENE")]))
  colnames(sum_tab1) <- c("GENE", "num_pairs_cans")
  sum_tab2 <- data.frame(table(unique(tab_tmp[tab_tmp$regulated, c("pair", "GENE", "Cancer")])[, c("GENE")]))
  sum_tab2 <- merge(sum_tab2, sum_tab1, all.x = T)
  sum_tab2$ratio <- sum_tab2$Freq/sum_tab2$num_pairs_cans
  sum_tab2 <- sum_tab2[order(sum_tab2$ratio, decreasing = T),]
  sum_tab3 <- data.frame(table(unique(tab_tmp[, c("pair", "GENE", "regulated", "Cancer")])[, c("regulated", "GENE","Cancer")]))
  sum_tab4 <- data.frame(table(unique(tab_tmp[ c("pair", "GENE", "Cancer")])[, c("GENE","Cancer")]))
  colnames(sum_tab4) <- c("GENE", "Cancer", "num_pairs_can")
  sum_tab3 <- merge(sum_tab3, sum_tab4, all.x = T)
  sum_tab3$ratio_can <- sum_tab3$Freq/sum_tab3$num_pairs_can
  
  df <- sum_tab3[sum_tab3$Freq > 1 & sum_tab3$regulated == "TRUE",]
  df_count2 <- data.frame(table(df$GENE))
  df <- df[df$GENE %in% df_count2$Var1[df_count2$Freq == length(cancers_sort)],]
  df_count <- group_by(df, GENE)
  df_count <- data.frame(summarise(df_count, ave_ratio = sum(ratio_can)/3))
  df <- sum_tab3[sum_tab3$GENE %in% df$GENE ,]
  df_count <- df_count[order(df_count$ave_ratio, decreasing = T),]
  df <- df[df$GENE %in% df_count$GENE[1:20],]
  
  tab4plot <- unique(sup_cans_tab_pk[(sup_cans_tab_pk$GENE %in% df$GENE) & !is.na(sup_cans_tab_pk$GENE) & (!is.na(sup_cans_tab_pk$regulated) & sup_cans_tab_pk$fdr_sig) & sup_cans_tab_pk$SELF == self, 
                                     c("GENE", "SUB_GENE", "SUB_MOD_RSD", "Cancer", "FDR_pho_kin", "coef_pho_kin")])
  colnames(tab4plot) <- c("GENE", "SUB_GENE", "SUB_MOD_RSD", "Cancer", "FDR", "coef")
  tab4plot$log10FDR <- -log10(tab4plot$FDR)
  tab4plot$site <- paste0(tab4plot$SUB_GENE, ":", tab4plot$SUB_MOD_RSD)
  tab4plot$x <- paste0(tab4plot$GENE, "-", tab4plot$Cancer)
  lim = median(tab4plot$coef) + 1.5*IQR(tab4plot$coef)
  cap <- min(lim, 2)
  tab4plot$coef_capped <- tab4plot$coef
  tab4plot$coef_capped[tab4plot$coef > cap] <- cap
  tab4plot$coef_capped[tab4plot$coef < -cap] <- (-cap)
  tab4plot$Cancer <- factor(tab4plot$Cancer, levels = cancers_sort)
  tab4plot$GENE <- factor(tab4plot$GENE, levels = as.vector(df_count$GENE)[order(df_count$ave_ratio, decreasing = T)])
  
  minmax_can <- vector(mode = "logical", length = nrow(tab4plot))
  for (gene in unique(tab4plot$GENE)) {
    for (cancer in unique(tab4plot$Cancer)) {
      tab_gene <- tab4plot[tab4plot$GENE == gene & tab4plot$Cancer == cancer,]
      min_fc <- min(tab_gene$coef)
      max_fc <- max(tab_gene$coef)
      minmax_can[tab4plot$coef == min_fc] <- TRUE
      minmax_can[tab4plot$coef == max_fc] <- TRUE
    }
  }
  tab4plot$minmax_can <- minmax_can
  
  minmax_cans <- vector(mode = "logical", length = nrow(tab4plot))
  for (gene in unique(tab4plot$GENE)) {
    tab_gene <- tab4plot[tab4plot$GENE == gene,]
    min_fc <- min(tab_gene$coef)
    max_fc <- max(tab_gene$coef)
    minmax_cans[tab4plot$coef == min_fc] <- TRUE
    minmax_cans[tab4plot$coef == max_fc] <- TRUE
  }
  tab4plot$minmax_cans <- minmax_cans
  
  tab4plot <- tab4plot[tab4plot$minmax_cans,]
  
  p <- ggplot()
  p <- p + geom_point(data=tab4plot, aes(x = GENE, y = coef, fill = Cancer, size = log10FDR), 
                      shape = 16, stroke = 0, alpha = 0.6)
  p <- p + geom_hline(yintercept = 0, linetype = 2, color = "grey")
  p <- p + geom_text_repel(data=tab4plot, 
                           aes(x = GENE, y = coef, label = site, color = Cancer), size = 2.5, force = 1)
  p <- p + ggtitle(label = paste0("phosphosites with the largest/lowest effect size correlated with each kinase across cancers (FDR<", reg_sig, ")"))
  p <- p + scale_fill_manual(values = color_cat_man)
  p <- p + scale_color_manual(values = color_cat_man)
  p <- p + theme_bw()
  p <- p + theme_nogrid()
  p <- p + xlab('kinase')+ylab("effect size")
  p <- p + theme(axis.text.x = element_text(angle = 90, face = "bold", hjust = 0.5, vjust = 0.5, size = 15))
  p <- p + theme(axis.title=element_text(size=10))
  p <- p + theme(axis.text.y = element_text(colour="black", size=15))
  p <- p + theme(title = element_text(size = 6, face = "italic"))
  fn = paste0(makeOutDir(resultD = resultD),'substrates_of_', self,'_kinases_top_regulated_average_ratio_across_cancers.pdf')
  pdf(file = fn, height=5, width = 10)
  print(p)
  dev.off()
}

# plot cis ----------------------------------------------------------------
for (self in c("cis")) {
  sup_cans_tab_pk_self <- sup_cans_tab_pk[sup_cans_tab_pk$SELF == self,]
  tab_tmp <- sup_cans_tab_pk_self
  sum_tab1 <- data.frame(table(unique(tab_tmp[, c("pair", "GENE", "Cancer")])[, c("GENE")]))
  colnames(sum_tab1) <- c("GENE", "num_pairs_cans")
  sum_tab2 <- data.frame(table(unique(tab_tmp[tab_tmp$regulated, c("pair", "GENE", "Cancer")])[, c("GENE")]))
  sum_tab2 <- merge(sum_tab2, sum_tab1, all.x = T)
  sum_tab2$ratio <- sum_tab2$Freq/sum_tab2$num_pairs_cans
  sum_tab2 <- sum_tab2[order(sum_tab2$ratio, decreasing = T),]
  sum_tab3 <- data.frame(table(unique(tab_tmp[, c("pair", "GENE", "regulated", "Cancer")])[, c("regulated", "GENE","Cancer")]))
  sum_tab4 <- data.frame(table(unique(tab_tmp[ c("pair", "GENE", "Cancer")])[, c("GENE","Cancer")]))
  colnames(sum_tab4) <- c("GENE", "Cancer", "num_pairs_can")
  sum_tab3 <- merge(sum_tab3, sum_tab4, all.x = T)
  sum_tab3$ratio_can <- sum_tab3$Freq/sum_tab3$num_pairs_can
  
  df <- sum_tab3[sum_tab3$Freq > 1 & sum_tab3$regulated == "TRUE",]
  df_count2 <- data.frame(table(df$GENE))
  df <- df[df$GENE %in% df_count2$Var1[df_count2$Freq == length(cancers_sort)],]
  df_count <- group_by(df, GENE)
  df_count <- data.frame(summarise(df_count, ave_ratio = sum(ratio_can)/3))
  df <- sum_tab3[sum_tab3$GENE %in% df$GENE ,]
  df_count <- df_count[order(df_count$ave_ratio, decreasing = T),]
  df <- df[df$GENE %in% df_count$GENE[1:20],]
  
  tab4plot <- unique(sup_cans_tab_pk[(sup_cans_tab_pk$GENE %in% df$GENE) & !is.na(sup_cans_tab_pk$GENE) & (!is.na(sup_cans_tab_pk$regulated) & sup_cans_tab_pk$fdr_sig) & sup_cans_tab_pk$SELF == self, 
                                     c("GENE", "SUB_GENE", "SUB_MOD_RSD","Cancer", "FDR_pro_kin", "coef_pro_kin")])
  colnames(tab4plot) <- c("GENE", "SUB_GENE", "SUB_MOD_RSD", "Cancer", "FDR", "coef")
  tab4plot$log10FDR <- -log10(tab4plot$FDR)
  tab4plot$site <- paste0(tab4plot$SUB_GENE, ":", tab4plot$SUB_MOD_RSD)
  tab4plot$x <- paste0(tab4plot$GENE, "-", tab4plot$Cancer)
  lim = median(tab4plot$coef) + 1.5*IQR(tab4plot$coef)
  cap <- min(lim, 2)
  tab4plot$coef_capped <- tab4plot$coef
  tab4plot$coef_capped[tab4plot$coef > cap] <- cap
  tab4plot$coef_capped[tab4plot$coef < -cap] <- (-cap)
  tab4plot$Cancer <- factor(tab4plot$Cancer, levels = cancers_sort)
  tab4plot$GENE <- factor(tab4plot$GENE, levels = as.vector(df_count$GENE)[order(df_count$ave_ratio, decreasing = T)])
  minmax_can <- vector(mode = "logical", length = nrow(tab4plot))
  for (gene in unique(tab4plot$GENE)) {
    for (cancer in unique(tab4plot$Cancer)) {
      tab_gene <- tab4plot[tab4plot$GENE == gene & tab4plot$Cancer == cancer,]
      min_fc <- min(tab_gene$coef)
      max_fc <- max(tab_gene$coef)
      minmax_can[tab4plot$coef == min_fc] <- TRUE
      minmax_can[tab4plot$coef == max_fc] <- TRUE
    }
  }
  tab4plot$minmax_can <- minmax_can
  
  minmax_cans <- vector(mode = "logical", length = nrow(tab4plot))
  for (gene in unique(tab4plot$GENE)) {
    tab_gene <- tab4plot[tab4plot$GENE == gene,]
    min_fc <- min(tab_gene$coef)
    max_fc <- max(tab_gene$coef)
    minmax_cans[tab4plot$coef == min_fc] <- TRUE
    minmax_cans[tab4plot$coef == max_fc] <- TRUE
  }
  tab4plot$minmax_cans <- minmax_cans
  
  tab4plot <- tab4plot[tab4plot$minmax_cans,]  
  p <- ggplot()
  p <- p + geom_point(data=tab4plot, aes(x = GENE, y = coef, fill = Cancer, size = log10FDR), 
                      shape = 16, stroke = 0, alpha = 0.6)
  p <- p + geom_hline(yintercept = 0, linetype = 2, color = "grey")
  p <- p + geom_text_repel(data=tab4plot, 
                           aes(x = GENE, y = coef, label = site, color = Cancer), size = 2.5, force = 1)
  p <- p + ggtitle(label = paste0("phosphosites with the largest/lowest effect size correlated with each kinase across cancers (FDR<", reg_sig, ")"))
  p <- p + scale_fill_manual(values = color_cat_man)
  p <- p + scale_color_manual(values = color_cat_man)
  p <- p + theme_bw()
  p <- p + theme_nogrid()
  p <- p + xlab('kinase')+ylab("effect size")
  p <- p + theme(axis.text.x = element_text(angle = 90, face = "bold", hjust = 0.5, vjust = 0.5, size = 15))
  p <- p + theme(axis.title=element_text(size=10))
  p <- p + theme(axis.text.y = element_text(colour="black", size=10))
  p <- p + theme(title = element_text(size = 6, face = "italic"))
  fn = paste0(makeOutDir(resultD = resultD),'substrates_of_', self,'_kinases_top_regulated_average_ratio_across_cancers.pdf')
  pdf(file = fn, height=4, width = 6)
  print(p)
  dev.off()
}

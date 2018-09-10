# Yige Wu @ WashU 2018 Feb
# overlap differentially phosphorylated sites and trans-regulated phosphosites, find common kinases


# input -------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
resultDnow <- makeOutDir()

for (protein in c("kinase")) {
  tn = paste0(resultD, "signal_hub/tables/overlap_diffexp_regression/", protein,"_substrate_regression_overlap_cptac2p_3can_diffPhosphositeAnnotated.txt")
  table_3can <- fread(tn, data.table = F)
  table_3can <- markSigSiteCan(regression = table_3can, sig_thres = sig, protein_type = protein)
  table_3can_sig <- table_3can[table_3can$fdr_sig & table_3can$coef_sig,]
  table_3can_sig$diffexp_type[is.na(table_3can_sig$diffexp_type)] <- "undiff"
  sum_table <- data.frame(table(table_3can_sig[,c("KINASE", "Cancer", "diffexp_type")]))
  sum_table$diffexp_type <- factor(sum_table$diffexp_type, levels = c("undiff", "up", "down"))
  # sum_table$Cancer <- factor(sum_table$Cancer, levels = c("BRCA", "OV", "CO"))
  diff_sites = sum_table[sum_table$diffexp_type %in% c("up", "down"),] %>%
    group_by(Cancer, KINASE) %>%
    summarise(diff_sites = sum(Freq))
  diff_sites <- data.frame(diff_sites)
  sum_table <- merge(sum_table, diff_sites, all.x =T); sum_table <- sum_table[sum_table$diff_sites > 0,]
  sum_table$kinase <- reorder(sum_table$KINASE, -sum_table$diff_sites)
  sum_table$Freq_p <- NA
  sum_table$Freq_p[sum_table$Freq > 0] <- sum_table$Freq[sum_table$Freq > 0]
  sum_table$Freq_p[sum_table$diffexp_type == "undiff"] <- NA
  
  ## make barplot
  p <- ggplot(data=sum_table, aes(y = Freq, x = kinase, fill = diffexp_type, label = Freq_p))
  p <- p + geom_bar(stat="identity", position='stack')
  p <- p + geom_text(position = position_stack(vjust = 0.5), size = 2)
  p <- p + scale_fill_manual(values = c("up" = "#E41A1C", "down" = "#377EB8", "undiff" = "grey"))
  p <- p + facet_grid(.~Cancer, drop=T, space = "free",scales = "free")#, space = "free", scales = "free")
  p <- p + theme_bw() + theme_nogrid() + scale_y_log10()
  p <- p + ylab("number of trans regulated phosphosites") +
    ggtitle("Trans-regulated phosphosites (FDR<.05, beta>0) overlap with differentially phophorylated sites (FDR<.05, >2 fold change)")
  p <- p + theme(axis.title.x = element_blank())
  p <- p + theme(axis.text.x = element_text(colour="black", size= 5.5, angle= 30, hjust=0.7), 
                 axis.text.y = element_text(colour="black", size=10))
  p
  ggsave(filename = paste0(resultDnow, 'cptac2p_3can_trans_', protein,'_substrate_overlap_diffPhosphosites.pdf'),
         width = 8, height = 4)
}

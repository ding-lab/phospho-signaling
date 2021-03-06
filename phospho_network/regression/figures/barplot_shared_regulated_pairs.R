# Yige Wu @ WashU 2018 Jul
# barplot for the top kinases with high and low kinase-substrate pairs

# source ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/pan3can_aes.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_plotting.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
library(dplyr)
# variables ---------------------------------------------------------------
cancers_sort <- c("BRCA", "OV", "CO")
reg_sig <- 0.05
diff_sig <- 0.2
diff_log2fc <- 1
color_cat_man <- c(colors['BRCA'], colors['OV'], colors['COAD'], "#bdbdbd"); names(color_cat_man) <- c("BRCA", "OV", "CO", "other")

# inputs -------------------------------------------------------------------
enzyme_type <- "kinase"
sup_cans_tab_en <- fread(input = paste0("./cptac2p/analysis_results/phospho_network/regression/tables/regression_omnipath&newworkin&depod_protein_level/", enzyme_type, "_substrate_regression_cptac2p_3can_tumor.txt"), data.table = F)
sup_cans_tab_en$pair <- paste0(sup_cans_tab_en$KINASE, ":", sup_cans_tab_en$SUBSTRATE, ":", sup_cans_tab_en$SUB_MOD_RSD)
colnames(sup_cans_tab_en)[1:2] <- c("GENE", "SUB_GENE")

sup_cans_tab_en <- markSigSiteCan(sup_cans_tab_en, sig_thres = reg_sig, enzyme_type = enzyme_type)
## annotate kinase substrate regulation
sup_cans_tab_en$regulated <- (sup_cans_tab_en$coef_sig & sup_cans_tab_en$fdr_sig)


# summarize kinase-substrate pairs ----------------------------------------
enzyme_type <- "kinase"
sup_cans_tab_en <- sup_cans_tab[sup_cans_tab$enzyme_type == enzyme_type,]
sup_cans_tab_en <- data.frame(sup_cans_tab_en)
sup_cans_tab_en <- markSigSiteCan(sup_cans_tab_en, sig_thres = reg_sig, enzyme_type = enzyme_type)
## annotate kinase substrate regulation
sup_cans_tab_en$regulated <- (sup_cans_tab_en$coef_sig & sup_cans_tab_en$fdr_sig)

## get table
for (self in c("cis", "trans")) {
  sup_cans_tab_en_self <- sup_cans_tab_en[sup_cans_tab_en$SELF == self,]
  tab_tmp <- sup_cans_tab_en_self
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
  
  color_cat <- vector(mode = "character", length = nrow(df))
  color_cat[df$regulated == "TRUE"] <- as.vector(df$Cancer)[df$regulated == "TRUE"]
  color_cat[df$regulated == "FALSE"]  <- "other"
  df$color_cat <- color_cat
  df$Cancer <- factor(df$Cancer, levels = cancers_sort)
  df$color_cat <- factor(df$color_cat, levels = c("other", cancers_sort))
  df <- rbind(df[df$regulated == "TRUE",], df[df$regulated == "FALSE",])
  df$GENE <- factor(df$GENE, levels = as.vector(df_count$GENE)[order(df_count$ave_ratio, decreasing = T)])
  p <- ggplot()
  p <- p + geom_bar(data=df, aes(y = Freq, x = Cancer, fill = color_cat, group = Cancer), stat="identity", position='stack', color = "#000000")
  p <- p + facet_grid(.~GENE, drop=T, space = "free",scales = "free")#, space = "free", scales = "free")
  p <- p + scale_fill_manual(values = color_cat_man)
  p <- p + theme_bw()
  p <- p + theme_nogrid()
  p <- p + scale_y_log10()
  p <- p + xlab('kinase')+ylab("number of substrate phosphosites")
  p <- p + theme(axis.title=element_text(size=10))
  p <- p + theme( axis.text.y = element_text(colour="black", size=10), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                  panel.spacing.x = unit(0, "lines"), strip.text.x = element_text(size = 5))
  p
  fn = paste0(makeOutDir(resultD = resultD),'3can_', self,'_kinases_top_regulated_average_ratio_across_cancers.pdf')
  ggsave(file=fn, height=3, width=10)
}


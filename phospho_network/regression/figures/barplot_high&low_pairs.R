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

## annotate up and down regulation
sup_cans_tab_pk$diffexp_type <- "other"
diffexp_sig <- ((!is.na(sup_cans_tab_pk$KSEA_pvalue) & sup_cans_tab_pk$KSEA_pvalue < diff_sig) | (!is.na(sup_cans_tab_pk$diffexp_log2FC & abs(sup_cans_tab_pk$diffexp_log2FC) > diff_log2fc)))
sup_cans_tab_pk$diffexp_type[diffexp_sig & sup_cans_tab_pk$enzyme_direction == "up" & sup_cans_tab_pk$substrate_direction == "up"] <- "up"
sup_cans_tab_pk$diffexp_type[diffexp_sig & sup_cans_tab_pk$enzyme_direction == "down" & sup_cans_tab_pk$substrate_direction == "down"] <- "down"

## get table
sum_tab1 <- data.frame(table(unique(sup_cans_tab_pk[, c("pair", "GENE", "diffexp_type", "Cancer")])[, c("GENE")]))
colnames(sum_tab1) <- c("GENE", "num_pairs_cans")
sum_tab2 <- data.frame(table(unique(sup_cans_tab_pk[, c("pair", "GENE", "diffexp_type", "Cancer")])[, c("GENE","diffexp_type")]))
sum_tab2 <- sum_tab2[sum_tab2$diffexp_type != "other",]
sum_tab2 <- merge(sum_tab2, sum_tab1, all.x = T)
sum_tab2$ratio <- sum_tab2$Freq/sum_tab2$num_pairs_cans
sum_tab2 <- sum_tab2[order(sum_tab2$ratio, decreasing = T),]
sum_tab3 <- data.frame(table(unique(sup_cans_tab_pk[, c("pair", "GENE", "diffexp_type", "Cancer")])[, c("GENE","diffexp_type", "Cancer")]))
sum_tab4 <- data.frame(table(unique(sup_cans_tab_pk[, c("pair", "GENE", "diffexp_type", "Cancer")])[, c("GENE","Cancer")]))
colnames(sum_tab4) <- c("GENE", "Cancer", "num_pairs_can")
sum_tab3 <- merge(sum_tab3, sum_tab4, all.x = T)
sum_tab3$ratio_can <- sum_tab3$Freq/sum_tab3$num_pairs_can

for (diffexp_type in c("up", "down")) {
  df <- sum_tab3[sum_tab3$diffexp_type == diffexp_type & sum_tab3$Freq > 1,]
  df_count2 <- data.frame(table(df$GENE))
  df <- df[df$GENE %in% df_count2$Var1[df_count2$Freq == length(cancers_sort)],]
  df_count <- group_by(df, GENE)
  df_count <- data.frame(summarise(df_count, ave_ratio = sum(ratio_can)/3))
  df <- sum_tab3[sum_tab3$GENE %in% df$GENE & (sum_tab3$diffexp_type == diffexp_type | sum_tab3$diffexp_type == "other"),]
  df_count <- df_count[order(df_count$ave_ratio, decreasing = T),]
  df <- df[df$GENE %in% df_count$GENE[1:50],]
  color_cat <- vector(mode = "character", length = nrow(df))
  color_cat[df$diffexp_type != "other"] <- as.vector(df$Cancer)[df$diffexp_type != "other"]
  color_cat[df$diffexp_type == "other"]  <- "other"
  df$color_cat <- color_cat
  df$Cancer <- factor(df$Cancer, levels = cancers_sort)
  df$color_cat <- factor(df$color_cat, levels = c("other", cancers_sort))
  df <- rbind(df[df$diffexp_type == diffexp_type,], df[df$diffexp_type == "other",])
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
  fn = paste0(makeOutDir(resultD = resultD),'3can_', diffexp_type,'_kinases_top_average_ratio_across_cancers.pdf')
  ggsave(file=fn, height=3, width=10)
}

# Yige Wu @ WashU 2018 Jul
# 

# source ------------------------------------------------------------------
source("./cptac2p_analysis/phospho_network/phospho_network_plotting.R")
source('./cptac2p_analysis/phospho_network/phospho_network_shared.R')
library(ggrepel)


# set variables -----------------------------------------------------------
color_cat_man <- c(colors['BRCA'], colors['OV'], colors['COAD'], "#bdbdbd"); names(color_cat_man) <- c("BRCA", "OV", "CO", "other")



# inputs ------------------------------------------------------------------
graph_stat_cans <- fread(input = paste0(ppnD, "signal_hub/tables/calculate_BC_CC/igraph_stats.txt"), data.table = F)

sup_cans_tab <- NULL
for (cancer in cancers_sort) {
  sup_tab <- fread(paste0(ppnD, "kinase_activity/tables/integrate_enzyme_KSEA_diffexp_sub_FC_plus_cistrans_regression_plus_genoalt/", 
                          cancer, "_KSEA_enzyme_diffexp_substrate.txt"), data.table = F)
  sup_tab$Cancer <- cancer
  sup_cans_tab <- rbind(sup_cans_tab, sup_tab)
}

# classify ----------------------------------------------------------------
sup_cans_tab_pk <- sup_cans_tab[sup_cans_tab$enzyme_type == "kinase",]
sup_cans_tab_pk <- markSigSiteCan(sup_cans_tab_pk, sig_thres = reg_sig, enzyme_type = "kinase")
sup_cans_tab_pp <- sup_cans_tab[sup_cans_tab$enzyme_type == "phosphotase",]
sup_cans_tab_pp <- markSigSiteCan(sup_cans_tab_pp, sig_thres = reg_sig, enzyme_type = "phosphotase")
sup_cans_tab <- rbind(sup_cans_tab_pk, sup_cans_tab_pp)
## annotate kinase substrate regulation
sup_cans_tab$regulated <- (sup_cans_tab$coef_sig & sup_cans_tab$fdr_sig)
## annotate up and down regulation
sup_cans_tab$diffexp_type <- "other"
diffexp_sig <- ((!is.na(sup_cans_tab$KSEA_pvalue) & sup_cans_tab$KSEA_pvalue < diff_sig) | (!is.na(sup_cans_tab$diffexp_log2FC & abs(sup_cans_tab$diffexp_log2FC) > diff_log2fc)))
sup_cans_tab$diffexp_type[diffexp_sig & sup_cans_tab$enzyme_direction == "up" & sup_cans_tab$substrate_direction == "up"] <- "up"
sup_cans_tab$diffexp_type[diffexp_sig & sup_cans_tab$enzyme_direction == "down" & sup_cans_tab$substrate_direction == "down"] <- "down"
sup_cans_tab_pk <- sup_cans_tab[sup_cans_tab$enzyme_type == "kinase",]
sup_cans_tab_pp <- sup_cans_tab[sup_cans_tab$enzyme_type == "phosphotase",]


# betweenness~degree_total ------------------------------------------------
tab2p <- graph_stat_cans
tab2p <- tab2p[tab2p$name %in% sup_cans_tab_pk$GENE,]
tab2p$Cancer <- factor(tab2p$Cancer, levels = cancers_sort)
# p <- ggplot()
# p <- p + geom_point(data = tab2p, mapping = aes(x = betweenness, y = degree_total))
# p <- p + geom_text_repel(data = tab2p[tab2p$betweenness > 0 & tab2p$degree_total > quantile(x = tab2p$degree_total, probs = 0.90),], mapping = aes(x = betweenness, y = degree_total, label = name))
# p <- p + facet_grid(.~Cancer, scales = "free", space = "fixed")
# p <- p + theme_bw()
# p <- p + theme(panel.spacing.y = unit(0, "lines"),
#                panel.spacing.x = unit(0, "lines"))
# p
# fn = paste0(makeOutDir(resultD = resultD),'kinase_betweenness~degree_total_faceted.pdf')
# ggsave(file=fn, height=3, width=5)

p <- ggplot()
p <- p + geom_point(data = tab2p, mapping = aes(x = betweenness, y = degree_total, color = Cancer))
p <- p + geom_text_repel(data = tab2p[tab2p$betweenness > quantile(x = tab2p$betweenness, probs = 0.95) | tab2p$degree_total > quantile(x = tab2p$degree_total, probs = 0.95),], 
                         mapping = aes(x = betweenness, y = degree_total, label = name, color = Cancer))
p <- p + theme_bw()
p <- p + theme(panel.spacing.y = unit(0, "lines"),
               panel.spacing.x = unit(0, "lines"))
p <- p + scale_color_manual(values = color_cat_man)
p
fn = paste0(makeOutDir(resultD = resultD),'kinase_betweenness~degree_total.pdf')
ggsave(file=fn, height=5, width=6)

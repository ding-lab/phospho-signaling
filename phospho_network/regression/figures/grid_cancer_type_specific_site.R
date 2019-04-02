# Yige Wu @ WashU 2019 Feb
## show the different sites are regulated in the same kinase-substrate protein pair in different cancers


# source-------------------------
wd <- getwd()
if (wd != "/Users/yigewu/Box Sync") {
  setwd("/Users/yigewu/Box Sync")
}
source('./cptac2p_analysis/phospho_network/phospho_network_shared.R')
source('./cptac2p_analysis/phospho_network/phospho_network_plotting.R')


# set variables -----------------------------------------------------------
reg_nonNA <- 20
size <- 83
cancers2process <- c("BRCA", "OV", "CO", "UCEC", "CCRCC")


# gather regression results ----------------------------------------------------------------
fdr_thres <- 0.1
file_path_tmp <- paste0("./cptac2p/analysis_results/phospho_network/regression/tables/generate_regression_regulated_uniq_marked/",  "regression_size", size, "_FDR", fdr_thres, "_detected_in_", paste0(cancers2process, collapse = "_"),".txt")
sup_tab <- fread(input = file_path_tmp, data.table = F)
sup_tab$pair_pro <- paste0(sup_tab$GENE, ":", sup_tab$SUB_GENE)

# bussiness ---------------------------------------------------------------

for (SELF in c("trans")) {
  tab2p_sup <- NULL
  
  for (col2evaluate in c("GENE", "SUB_GENE")) {
    ## take out kinase-substrate protein pairs that purely based on motif prediction
    # tab_tmp <- sup_tab[sup_tab$SELF == SELF,]
    tab_tmp <- tab_tmp %>%
      filter(pair_pro %in% omnipath_tab$pair_pro[!(omnipath_tab$Source %in% c("MIMP", "NetKIN", "PhosphoNetworks"))])
    
    tab_tmp$is.smg <- get_SMG_by_cancer(gene_vector = tab_tmp[, col2evaluate], cancer_vector = tab_tmp$Cancer)
    
    ## get 3 substrate phosphosites from each protein pairs but also make sure to include the correlated ones
    tab_pairs <- tab_tmp %>%
      arrange(-regulated, FDR_pho_kin) %>%
      select(pair, regulated, FDR_pho_kin, Cancer, pair_pro, SUB_MOD_RSD, SUB_GENE, is.smg)
    tab_pairs <- tab_pairs[!duplicated(tab_pairs$pair),]
    tab_pairs$row_id <- 1:nrow(tab_pairs)
    tab_pairs$is_top <- get_top_by_id(value_vector = -(tab_pairs$row_id), id_vector = tab_pairs$pair_pro, num_top = 3)
    
    pair_cancer_summary <- tab_tmp %>%
      filter(pair_pro %in% tab_pairs$pair_pro[tab_pairs$is_top]) %>%
      # filter(SUB_MOD_RSD %in% tab_pairs$SUB_MOD_RSD[tab_pairs$is_top]) %>%
      filter(SUB_MOD_RSD %in% tab_pairs$SUB_MOD_RSD[tab_pairs$is_top & tab_pairs$is.smg]) %>%
      # select(pair_pro, Cancer, SUB_MOD_RSD, pair, GENE, SUB_GENE, regulated) %>%
      filter(!(SUB_GENE %in% c("FLNA", "LMNA", "SRRM2", "YAP1"))) %>%
      filter(!(SUB_MOD_RSD %in% "S37" & SUB_GENE %in% "RB1")) %>%
      unique()
    
    pair_pro_cancer_summary <- pair_cancer_summary %>%
      select(pair_pro, pair) %>%
      unique() %>%
      select(pair_pro) %>%
      table() %>%
      data.frame()
    colnames(pair_pro_cancer_summary) <- c("pair_pro", "Freq")
    
    pair_cancer_summary <- merge(pair_cancer_summary, pair_pro_cancer_summary, by = c("pair_pro"), all.x = T)
    pair_cancer_summary %>% head()
    
    pair_cancer_summary <- pair_cancer_summary %>%
      filter(pair_pro %in% pair_cancer_summary$pair_pro[pair_cancer_summary$regulated])
    
    tab2p <- pair_cancer_summary
    ## filter for cancer SMGs in substrate
    tab2p$is.smg <- get_SMG_by_cancer(gene_vector = tab2p[, col2evaluate], cancer_vector = tab2p$Cancer)
    
    
    tab2p <- tab2p %>%
      filter(pair %in% tab2p$pair[tab2p$is.smg])
    tab2p$bar_len <- 1/tab2p$Freq
    tab2p$id <- paste0(tab2p$pair_pro, ":", tab2p$Cancer)
    tab2p <- tab2p[order(tab2p$id),]
    tab2p <- tab2p[order(tab2p$id, tab2p$SUB_MOD_RSD),]
    
    y_print <- vector(mode = "numeric", length = nrow(tab2p))
    for (i in 1:nrow(tab2p)) {
      if (duplicated(tab2p$id)[i]) {
        y_print[i] <- y_print[i-1] + tab2p$bar_len[i]
      } else {
        y_print[i] <- tab2p$bar_len[i]/2
      }
    }
    tab2p$y_print <- y_print
    tab2p$Cancer <- order_cancer_rev(tab2p$Cancer)
    tab2p$SUB_MOD_RSD <- factor(tab2p$SUB_MOD_RSD, levels = unique(tab2p$SUB_MOD_RSD))
    
    p <- ggplot()
    p <- p + geom_bar(data = tab2p, mapping = aes(x = pair_pro, y = bar_len, group = SUB_MOD_RSD), stat = "identity", position = "stack", color = "white",  fill = "grey90")
    # p <- p + geom_bar(data = tab2p, mapping = aes(x = pair_pro, y = bar_len, fill = regulated, group = SUB_MOD_RSD), stat = "identity", position = "stack", color = "white")
    p <- p + geom_text(data = tab2p, mapping = aes(x = pair_pro, y = y_print, label = SUB_MOD_RSD, alpha = regulated, color = regulated), size = 2.5)
    # p <- p + scale_fill_manual(values = c("TRUE" = "#FB9A99", "FALSE" = "grey90"))
    p <- p + scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.5))
    p <- p + scale_color_manual(values = c("TRUE" = set1[1], "FALSE" = "grey50"))
    p <- p + facet_grid(SUB_GENE~Cancer, drop=T,shrink = T, space = "free",scales = "free")#, space = "free", scales = "free")
    p <- p + coord_flip()
    p <- p + theme_nogrid()
    p <- p + theme(axis.text.x = element_blank(), axis.ticks = element_blank(),
                   axis.title.y = element_blank(), axis.title.x = element_blank(),
                   strip.text.y = element_text(size = 10, angle = 0),
                   strip.background = element_rect(fill = "white", color = "white"),
                   panel.spacing.y = unit(0, "lines"),
                   panel.spacing.x = unit(0, "lines"))
    p <- p + guides(fill = F)
    p
    fn <- paste0(makeOutDir(resultD = resultD), "cancer_specific_", SELF, "_regulated_pairs_", col2evaluate, "_is_SMG.pdf")
    ggsave(filename = fn, width = 8, height = 8)
    
    tab2p_sup <- rbind(tab2p_sup, tab2p)
  }
  
  tab2p <- tab2p_sup
  tab2p$facet_y <- tab2p$SUB_GENE
  row_tmp <- (tab2p$GENE %in% unlist(SMGs))
  tab2p$facet_y[row_tmp] <- tab2p$GENE[row_tmp]
  
  p <- ggplot()
  p <- p + geom_bar(data = tab2p, mapping = aes(x = pair_pro, y = bar_len, group = SUB_MOD_RSD), stat = "identity", position = "stack", color = "white",  fill = "grey90")
  # p <- p + geom_bar(data = tab2p, mapping = aes(x = pair_pro, y = bar_len, fill = regulated, group = SUB_MOD_RSD), stat = "identity", position = "stack", color = "white")
  p <- p + geom_text(data = tab2p, mapping = aes(x = pair_pro, y = y_print, label = SUB_MOD_RSD, alpha = regulated, color = regulated), size = 2.5)
  # p <- p + scale_fill_manual(values = c("TRUE" = "#FB9A99", "FALSE" = "grey90"))
  p <- p + scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.5))
  p <- p + scale_color_manual(values = c("TRUE" = set1[1], "FALSE" = "grey50"))
  p <- p + facet_grid(facet_y~Cancer, drop=T,shrink = T, space = "free",scales = "free")#, space = "free", scales = "free")
  p <- p + coord_flip()
  p <- p + theme_nogrid()
  p <- p + theme(axis.text.x = element_blank(), axis.ticks = element_blank(),
                 axis.title.y = element_blank(), axis.title.x = element_blank(),
                 strip.text.y = element_text(size = 10, angle = 0),
                 strip.background = element_rect(fill = "white", color = "white"),
                 panel.spacing.y = unit(0, "lines"),
                 panel.spacing.x = unit(0, "lines"))
  p <- p + guides(fill = F)
  p
  fn <- paste0(makeOutDir(resultD = resultD), "cancer_specific_", SELF, "_regulated_pairs_SMG_related.pdf")
  ggsave(filename = fn, width = 8, height = 8)
}

stop()
# mut_cnv_cans %>%
#   filter(GENE == "CTNNB1", SUB_GENE %in% c("PRKD1", "PRKCD", "PRKACB", "PRKACA", "PAK4", "GSK3B"), cancer == "CO", p < 0.1) %>%
#   arrange(p)
# ## PRKD1
# 
# mut_cnv_cans %>%
#   filter(SUB_GENE %in% c("GSK3B"), cancer == "CO") %>%
#   arrange(p)
# 
# mut_cnv_cans %>%
#   filter(GENE == "CTNNB1", SUB_GENE %in% c("PRKD1", "PRKCD", "PRKACB", "PRKACA", "PAK4", "GSK3B"), cancer == "UCEC", p < 0.05) %>%
#   arrange(p)
# ## confusing down-regualtion
# 
# mut_cnv_cans %>%
#   filter(GENE == "RB1", SUB_GENE %in% c("PRKAA1", "PPP1CB", "CDK18", "CDK1"), cancer == "OV") %>%
#   arrange(p)
# 
# mut_cnv_cans %>%
#   filter(GENE == "RB1", SUB_GENE %in% c("PRKAA1", "PPP1CB", "CDK18", "CDK1"), cancer == "BRCA") %>%
#   arrange(p)
# 
# mut_cnv_cans %>%
#   filter(GENE == "TP53", SUB_GENE %in% c("IKBKB"), cancer == "OV") %>%
#   arrange(p)
# 
# mut_cnv_cans %>%
#   filter(GENE == "PTEN", SUB_GENE %in% c("CREB1"), cancer == "CCRCC") %>%
#   arrange(p)
# 
# mut_cnv_cans %>%
#   filter(GENE == "MAP2K4", SUB_GENE %in% c("PAK1"), cancer == "BRCA") %>%
#   arrange(p)

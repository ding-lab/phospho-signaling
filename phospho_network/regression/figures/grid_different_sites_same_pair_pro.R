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

for (SELF in c("cis", "trans")) {
  tab_tmp <- sup_tab[sup_tab$SELF == SELF,]
  
  pair_pro_cancer_summary <- data.frame(table(unique(tab_tmp[tab_tmp$regulated,c("pair", "pair_pro", "Cancer", "regulated", "SELF")])[,c("pair_pro", "Cancer")]))
  pair_pro_cancer_summary <- pair_pro_cancer_summary[pair_pro_cancer_summary$Freq > 0,]
  pair_pro_cancer_summary %>% 
    arrange(-Freq) %>%
    head()
  
  pair_cancer_summary <- unique(tab_tmp[tab_tmp$regulated & tab_tmp$pair_pro %in% pair_pro_cancer_summary$pair_pro, c("pair_pro", "Cancer", "SUB_MOD_RSD", "pair", "GENE", "SUB_GENE")])
  pair_cancer_summary <- merge(pair_cancer_summary, pair_pro_cancer_summary, by = c("pair_pro", "Cancer"), all.x = T)
  pair_cancer_summary %>% head()
  
  tab2p <- pair_cancer_summary
  ## filter for cancer SMGs in substrate
  tab2p$is.SUB_GENE.smg <- get_SMG_by_cancer(gene_vector = tab2p[, "SUB_GENE"], cancer_vector = tab2p$Cancer)
  tab2p$is.GENE.smg <- get_SMG_by_cancer(gene_vector = tab2p[, "GENE"], cancer_vector = tab2p$Cancer)
  
  tab2p <- tab2p[tab2p$is.GENE.smg | tab2p$is.SUB_GENE.smg,]
  tab2p$bar_len <- 1/tab2p$Freq
  tab2p$id <- paste0(tab2p$pair_pro, ":", tab2p$Cancer)
  tab2p <- tab2p[order(tab2p$id),]
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
    
  p <- ggplot()
  p <- p + geom_bar(data = tab2p, mapping = aes(x = pair_pro, y = bar_len, fill = Cancer), stat = "identity", position = "stack", color = "black")
  p <- p + geom_text(data = tab2p, mapping = aes(x = pair_pro, y = y_print, label = SUB_MOD_RSD), size = 2.5)
  p <- p + scale_fill_manual(values = color_cancers2)
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
  fn <- paste0(makeOutDir(resultD = resultD), "cancer_specific_", SELF, "_regulated_pairs_in_SMGs.pdf")
  ggsave(filename = fn, width = 4.5, height = 5)

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

# Yige Wu @ WashU 2018 Jul
# barplot for the top kinases with high and low kinase-substrate pairs

# source ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/pan3can_aes.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_plotting.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
library(dplyr)
# variables ---------------------------------------------------------------
cancers_sort <- c("BRCA", "OV", "CO")
fdr_pk <- 0.05
color_cat_man <- c(colors['BRCA'], colors['OV'], colors['COAD'], "#bdbdbd"); names(color_cat_man) <- c("BRCA", "OV", "CO", "other")
top_num <- 10
# inputs -------------------------------------------------------------------
enzyme_type <- "kinase"
regression <- fread(input = paste0("./cptac2p/analysis_results/phospho_network/regression/tables/regression_omnipath&newworkin&depod_protein_level/kinase_substrate_regression_cptac2p_3can_tumor.txt"), data.table = F)
regression$GENE <- as.vector(regression$KINASE)
regression$SUB_GENE <- as.vector(regression$SUBSTRATE)
regression <- markSigSiteCan(regression = regression, sig_thres = fdr_pk, enzyme_type = enzyme_type)
regression$regulated <- (regression$fdr_sig & regression$coef_sig)
pair_cans_count  <- data.frame(table(regression$pair))
pair_cans_count3 <- pair_cans_count[pair_cans_count$Freq == length(cancers_sort),]

# summarize kinase-substrate pairs ----------------------------------------
## get table
for (self in c("cis", "trans")) {
  tab_tmp <- regression
  tab_tmp <- tab_tmp[tab_tmp$SELF == self & tab_tmp$pair %in% pair_cans_count3$Var1,]
  tab_tmp <- tab_tmp[(tab_tmp$uniq_BRCA | tab_tmp$uniq_OV | tab_tmp$uniq_CO),]
  tab_tmp <- unique(tab_tmp[, c("GENE", "pair", "uniq_BRCA", "uniq_OV", "uniq_CO")])
  tab_tmp.m <- melt(tab_tmp, id.vars = c("GENE", "pair"))
  tab_tmp.m <- tab_tmp.m[tab_tmp.m$value == "TRUE",]
  sum_tab1 <- data.frame(table(tab_tmp.m[, c("GENE", "variable")]))
  sum_tab1$Cancer <- str_split_fixed(string = sum_tab1$variable, pattern = "_", 2)[,2]
  sum_tab2 <- sum_tab1[sum_tab1$Freq > 1,]
  sum_tab1 <- sum_tab1[sum_tab1$GENE %in% sum_tab2$GENE,]
  tab2p <- sum_tab1
  tab2p_sort <- NULL
  genes2p_sort <- NULL
  for (cancer in cancers_sort) {
    tab2p_tmp <- tab2p[tab2p$Cancer == cancer,]
    tab2p_tmp <- tab2p_tmp[order(tab2p_tmp$Freq, decreasing = T),]
    tab2p_tmp <- tab2p_tmp[!(tab2p_tmp$GENE %in% genes2p_sort),]
    genes2p <- as.vector(tab2p_tmp$GENE[1:top_num])
    tab2p_tmp <- sum_tab1[sum_tab1$GENE %in% genes2p,]
    tab2p_tmp$facet <- cancer
    tab2p_sort <- rbind(tab2p_sort, tab2p_tmp)
    genes2p_sort <- c(genes2p_sort, genes2p)
  }
  tab2p_sort$facet <- factor(tab2p_sort$facet, levels = cancers_sort)
  tab2p_sort$GENE <- factor(tab2p_sort$GENE, levels = unique(genes2p_sort))
  
  p <- ggplot()
  p <- p + geom_bar(data=tab2p_sort, aes(y = Freq, x = GENE, fill = Cancer, group = Cancer), 
                    stat="identity", position='dodge', color = "#000000")
  p <- p + scale_fill_manual(values = color_cat_man)
  p <- p + theme_bw()
  p <- p + theme_nogrid()
  p <- p + scale_y_sqrt()
  # p <- p + scale_y_continuous(trans='log2')
  # p <- p + scale_y_log10()
  p <- p + xlab('kinase')+ylab("number of cancer-type specific substrate phosphosites")
  p <- p + theme(axis.title.y=element_text(size=5), axis.ticks.x = element_blank())
  p <- p + theme(axis.text.x = element_text(size = 10, face = "bold", angle = 90, hjust = 0.5, vjust = 0.5))
  p <- p + theme(axis.text.y = element_text(colour="black", size=10), axis.ticks.x = element_blank(),
                  panel.spacing.x = unit(0, "lines"), strip.text.x = element_text(size = 5))
  p
  fn = paste0(makeOutDir(resultD = resultD),'3can_', self,'_kinases_cancer_type_specific_top', top_num, '.pdf')
  ggsave(file=fn, height=3, width=6)
}

table(unique(sup_cans_tab_en[sup_cans_tab_en$regulated, c("Cancer", "pair")])[, c("Cancer")])

length(unique(regression$pair [regression$regulated]))

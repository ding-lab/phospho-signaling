# Yige Wu @ WashU 2018 Jan
# draw volcano plots for regression result (3 cancer types together)

# source ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_plotting.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
library(ggrepel)


# inputs ------------------------------------------------------------------
## input enzyme substrate table
# ptms_site_pairs_sup <- read_csv(paste0(ppnD, "compile_enzyme_substrate/tables/compile_omnipath_networkin/omnipath_networkin_enzyme_substrate_site_level_union.csv"))
ptms_site_pairs_sup <- read_csv(paste0(ppnD, "compile_enzyme_substrate/tables/compile_omnipath_networkin_depod/omnipath_networkin_enzyme_substrate_site_level_union_source_cleaned_addDEPOD_extended.csv"))

# set variables -----------------------------------------------------------
color_map <- c("firebrick1", "black")
names(color_map) <- c("TRUE", "FALSE")
color_direction <- c(set1[1], set1[2], "black")
names(color_direction) <- c("up", "down", "NA")
fdr_pp <- 0.1
fdr_pk <- 0.05
color_cat_man <- c(colors['BRCA'], colors['OV'], colors['COAD'], "#bdbdbd"); names(color_cat_man) <- c("BRCA", "OV", "CO", "other")

# kinases -----------------------------------------------------------------
enzyme_type <- "kinase"
regression <- fread(input = paste0("./cptac2p/analysis_results/phospho_network/regression/tables/regression_omnipath&newworkin&depod_protein_level/", enzyme_type, "_substrate_regression_cptac2p_3can_tumor.txt"), data.table = F)
regression$GENE <- as.vector(regression$KINASE)
regression$SUB_GENE <- as.vector(regression$SUBSTRATE)
regression <- markSigSiteCan(regression = regression, sig_thres = fdr_pk, enzyme_type = enzyme_type)
regression$regulated <- (regression$fdr_sig & regression$coef_sig)
subdir1 <- paste0(makeOutDir(resultD = resultD), enzyme_type, "/")
dir.create(subdir1)

for (self in c("trans")) {
  tab2p <- regression
  tab2p <- tab2p[tab2p$SELF == self,]
  tab2p_back <- tab2p[!(tab2p$regulated),]; tab2p_back <- tab2p_back[sample(x = 1:nrow(tab2p_back), size = 10000),]
  p <- ggplot()
  p <- p + geom_point(data = tab2p_back, aes(x=sd_pho_kin, y = sd_pho_sub), color = "grey", alpha = 0.2)
  p <- p + geom_point(data = tab2p[(tab2p$regulated),], aes(x=sd_pho_kin, y = sd_pho_sub), color = "red", alpha = 0.3)
  p <- p + facet_grid(Cancer~., scales = "fixed", space = "free")
  p
  fn = paste0(subdir1,'substrates_of_', self,'_', enzyme_type, '_top_regulated_average_ratio_across_cancers.pdf')
  pdf(file = fn, height=12, width = 4)
  print(p)
  dev.off()
  
  p <- ggplot()
  p <- p + geom_boxplot(data = tab2p, aes(x = regulated, y = sd_pho_sub))
  p <- p + facet_grid(Cancer~., scales = "fixed", space = "free")
  p <- p + ylim(c(0,3))
  p
  fn = paste0(subdir1,'sd_pho_sub', self,'_', enzyme_type, '_top_regulated_average_ratio_across_cancers.pdf')
  pdf(file = fn, height=12, width = 4)
  print(p)
  dev.off()
  
  p <- ggplot()
  p <- p + geom_boxplot(data = tab2p, aes(x = regulated, y = sd_pho_kin))
  p <- p + facet_grid(Cancer~., scales = "fixed", space = "free")
  p <- p + ylim(c(0,3))
  p
  fn = paste0(subdir1,'v', self,'_', enzyme_type, '_top_regulated_average_ratio_across_cancers.pdf')
  pdf(file = fn, height=12, width = 4)
  print(p)
  dev.off()

}




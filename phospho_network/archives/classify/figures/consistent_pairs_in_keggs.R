# Yige Wu @ WashU 2018 Jan
# showing consistent pairs (kinase from KSEA, substrate from differential expression) in pathways


# source ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_plotting.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/pan3can_aes.R')


# inputs ------------------------------------------------------------------
load("~/Box Sync/pan3can_shared_data/analysis_results/2015-08-01_Gene_Set.RData")
keywords <- c("cell cycle", "Hippo", "Notch", "PI3K", "TGF", "Ras", "p53", "Wnt")


# get kegg pathways -------------------------------------------------------
all_keggs <- names(KEGG)
keggs <- NULL
for (key in keywords) {
  keggs <- c(keggs, all_keggs[grepl(key, all_keggs, ignore.case = T)])
}
keggs <- unique(keggs)

for (cancer in c("BRCA")) {
# for (cancer in c("BRCA", "OV", "CO")) {
  ptms_table <- fread(input = paste0("cptac2p/cptac_shared/analysis_results/phospho_network/classify/tables/integrate_KSEA_diffexp/", cancer, "_KSEA_enzyme_diffexp_substrate.txt"), data.table = F)
  us_all_genes <- unique(c(as.vector(ptms_table$GENE), as.vector(ptms_table$SUB_GENE)))
  for (kegg_tmp in keggs) {
    gene_tmp <- us_all_genes[us_all_genes %in% KEGG[[kegg_tmp]]]
    if (length(gene_tmp) >0) {
      print(paste0(cancer, ":", kegg_tmp, ":", length(gene_tmp)))
      ptms_table[, kegg_tmp] <- ifelse(ptms_table$GENE %in% gene_tmp | ptms_table$SUB_GENE %in% gene_tmp, TRUE, FALSE)
    }
  }
  ptms_table2m <- ptms_table[rowSums(ptms_table == TRUE) > 0,]
  ptms_table_m <- melt(ptms_table2m, id.vars = c("GENE", "SUB_GENE", "SUB_MOD_RSD", "KINASE", "KIN_ACC_ID", "KIN_ORGANISM", "SUBSTRATE", "SUB_GENE_ID", "SUB_ACC_ID", "SUB_ORGANISM", "SITE_GRP_ID", "SITE_...7_AA", "Source", "enzyme_type", "enzyme_direction", "substrate_direction",
                                               "networkin_score", "kinase_zscore", "kinase_pvalue", "substrate_FC"))
  ptms_table_m <- ptms_table_m[ptms_table_m$value == "TRUE",]
  ptms_table_m$same_direction <- (as.vector(ptms_table_m$enzyme_direction) == as.vector(ptms_table_m$substrate_direction))
  ptms_table_m$consistent <- ifelse((ptms_table_m$enzyme_type == "kinase" & ptms_table_m$same_direction) | (ptms_table_m$enzyme_type == "phosphotase" & !ptms_table_m$same_direction), TRUE, FALSE)
  ptms_table_m <- unique(ptms_table_m)
  
  ##  annotate driver genes
  drivers <- loadGeneList(id = "driver", cancer = cancer)
  ptms_table_m$enzyme.is.driver <- (ptms_table_m$GENE %in% drivers) 
  ptms_table_m$substrate.is.driver <- (ptms_table_m$SUB_GENE %in% drivers) 
  tsg <- loadGeneList(id = "tsg", cancer = cancer, is.soft.limit = "soft")
  oncogene <- loadGeneList(id = "oncogene", cancer = cancer, is.soft.limit = "soft")
  ptms_table_m$enzyme.driver_type[ptms_table_m$GENE %in% tsg] <- "tsg"
  ptms_table_m$enzyme.driver_type[ptms_table_m$GENE %in% oncogene] <- "oncogene"
  ptms_table_m$substrate.driver_type[ptms_table_m$SUB_GENE %in% tsg] <- "tsg"
  ptms_table_m$substrate.driver_type[ptms_table_m$SUB_GENE %in% oncogene] <- "oncogene"
  ptms_table_m_driver <- ptms_table_m[(ptms_table_m$enzyme.is.driver | ptms_table_m$substrate.is.driver) & ptms_table_m$consistent,]
  write.table(x = ptms_table_m_driver,
              file = paste0(makeOutDir(), cancer, "_ptms_table_m_driver.txt"),
              row.names = F, col.names = T, quote = F, sep = "\t")
  
  p <- ggplot(ptms_table_m)
  p <- p + geom_tile(aes(x = GENE, y = SUB_GENE, fill = enzyme_direction, color = substrate_direction, alpha = 0.5*((consistent == TRUE) + 1)))
  p <- p + facet_grid(variable~., drop = T, space = "free",scales = "free")
  p <- p + theme_bw()
  p = p + theme(strip.text.y = element_text(size = 7, angle = 0), strip.text.x = element_text(size = 5))
  p = p + theme(axis.text.x = element_text(size = 8, angle = -30, hjust = 0, vjust = 1))
  p <- p + scale_fill_manual(values = c("up" = set1[1], "down" = set1[2]))
  p <- p + scale_color_manual(values = c("up" = set1[1], "down" = set1[2]))
  p = p + theme(panel.spacing.x=unit(0, "lines"), panel.spacing.y=unit(0, "lines"))
  p
  ggsave(filename = paste0(makeOutDir(), "ks_pairs_within_keggs_in_", cancer, ".pdf"), 
         height = 35, width = 20)
}


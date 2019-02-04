# Yige Wu @ WashU 2018 July
## overview figure for mutational impact of kinases on substrates

# source ------------------------------------------------------------------
source("Box Sync/cptac2p_analysis/phospho_network/phospho_network_plotting.R")
source("./cptac2p_analysis/phospho_network/phospho_network_shared.R")
source("./cptac2p_analysis/dependencies/tables/load_TCGA_pathways.R")


# set variables -----------------------------------------------------------
show_p_thres <- 0.2
sig_p_thres <- 0.05
# driver_genes2show <- unique(c(as.vector(driver_genes$alteredGene[driver_genes$Cancer == "UCEC"]), SMGs[["UCEC"]]))

# inputs ------------------------------------------------------------------
mut_cnv_cans <- fread(input = paste0(ppnD, "genoalt/tables/test_mut_impact_cis_pro_separate_variant/mut_cnv_cans.txt"), data.table = F)
fdr <- vector(mode = "numeric", length = nrow(mut_cnv_cans))
for (cancer in unique(mut_cnv_cans$cancer)) {
  for (self in c("cis", "trans")) {
    rows_tmp <- which(mut_cnv_cans$cancer == cancer)
    fdr[rows_tmp] <- p.adjust(p = mut_cnv_cans$p[rows_tmp], method = "fdr")
  }
}
mut_cnv_cans$fdr <- fdr

## annotate the substrates
affectedGenes2pathways <- map2TCGApathwaways(gene_list = unique(mut_cnv_cans$affectedGene), pathway_list = tcga_pathways_pluskegg_and_pathway)
affectedGenes2pathways <- data.frame(affectedGene = rep(x = names(affectedGenes2pathways), sapply(X = affectedGenes2pathways, FUN = function(x) length(x))), 
                                 affectedGene.path = unlist(affectedGenes2pathways, use.names = F))
affectedGenes2pathways <- affectedGenes2pathways[!(affectedGenes2pathways$affectedGene == "APC" & affectedGenes2pathways$affectedGene.path != "WNT") & !(affectedGenes2pathways$affectedGene == "GSK3B" & affectedGenes2pathways$affectedGene.path != "PI3K") & !(affectedGenes2pathways$affectedGene == "IRS1" & affectedGenes2pathways$affectedGene.path != "RTK RAS"),]
mut_cnv_cans$affectedGene.path <- NULL
mut_cnv_cans <- merge(mut_cnv_cans, affectedGenes2pathways, all.x = T)

for (genoalt_type in c("mut")) {
  ## get the list of kinase/phosphatases to show
  enzymes2show <- unique(mut_cnv_cans$alteredGene[mut_cnv_cans$p < 0.3 & mut_cnv_cans$p > 0 & mut_cnv_cans$genoalt_type == genoalt_type & (mut_cnv_cans$alteredGene %in% driver_genes$Gene[driver_genes$Cancer %in% c("COADREAD", "BRCA", "OV", "CO", "KIRC", "UCEC")])])
  # enzymes2show <- enzymes2show[!(enzymes2show %in% c("ACVR2A", "ARHGAP35", "CDK12", "EPHA2", "FLNA", "HUWE1", "INPPL1", "NCOR1", "RASA1", "SED2", "STAG2", "TSC1",
  #                                                    "ZFHX3", "ZFP36L1", "NFE2L2", "PDGFRA", "PPP2R1A"))]
  
  ## get the list of substrates to show
  subtrates2show <- enzymes2show
  df <- mut_cnv_cans[mut_cnv_cans$genoalt_type == genoalt_type,]
  df <- df[(df$alteredGene %in% enzymes2show) & (df$affectedGene %in% subtrates2show),]
  df$sig <- (df$p < sig_p_thres)
  df$fdr_sig <- (df$fdr < 0.1)
  
  # df <- df[!(df$alteredGene == "CTNNB1" & df$affectedGene %in% c("APC", "AXIN1")),]

  ## remove the duplicated pathway annotation
  tmp <- as.vector(df$affectedGene.path)
  tmp[is.na(tmp)] <- "other"
  tmp[tmp == "Mismatch repair"] <- "MMR"
  tmp[df$affectedGene == "STAT3"] <- "JAK-STAT"
  tmp[df$affectedGene == "AXIN1"] <- "WNT"
  
  df$affectedGene.path <- tmp
  df <- df[order(df$affectedGene.path, df$affectedGene),]
  df$pair_can  <- paste0(df$pair, ":", df$cancer)
  df <- df[!duplicated(df$pair_can),]
  df <- df[order(df$p, decreasing = T),]
  
  df$log10_Pvalue <- -log10(df$p)
  cap <- min(max(abs(df$meddiff), na.rm = T), 1)
  df$meddiff_capped <- df$meddiff
  df$meddiff_capped[df$meddiff > cap] <- cap
  df$meddiff_capped[df$meddiff < (-cap)] <- (-cap)
  df$driver_gene_type.en[df$driver_gene_type.en == ""] <- "other"
  df$driver_gene_type.en[df$alteredGene %in% c("INSR", "PPM1D")] <- "oncogene"
  
  df$driver_gene_type.en <- factor(df$driver_gene_type.en, levels = c("oncogene", "tsg", "other"))
  df$affectedGene.path <- factor(df$affectedGene.path, levels = unique(c(names(tcga_pathways_pluskegg_and_pathway), "WNT","MMR", "JAK-STAT", "other")))
  
  # names(sig_colors) <- c(paste0('P-value<', sig_p_thres), paste0('P-value>', sig_p_thres))
  
  # df <- df[!(df$affectedGene %in% c("NCOR1", "ZFYVE16", "ROCK1", "RAF1", "MAPK1")) & (df$affectedGene.path != "other"),]
  tmp <- as.vector(unique(df$alteredGene))
  df$alteredGene <- factor(df$alteredGene, levels = c(c("CTNNB1", "APC", "PTEN", "TP53", "MSH6"), tmp[!(tmp %in% c("CTNNB1", "APC", "PTEN", "TP53", "MSH6"))]))
  df$cancer <- factor(df$cancer, levels = rev(c("BRCA", "OV", "CO", "UCEC", "CCRCC")))
  p <- ggplot()
  # p <- p + geom_point(data = df, mapping = aes(x = affectedGene, y = cancer, fill = meddiff_capped, size = log10_Pvalue,
  #                                              colour = ifelse(sig == TRUE, paste0('P-value<', sig_p_thres), paste0('P-value>', sig_p_thres))), shape = 21, alpha = 0.6)
  p <- p + geom_point(data = df, mapping = aes(x = affectedGene, y = cancer, fill = meddiff_capped, size = log10_Pvalue,
                                               colour = ifelse(fdr_sig == TRUE, paste0('FDR<', 0.1), paste0('FDR>', 0.1))), shape = 21, alpha = 0.6)
  sig_colors <- c("black", "white"); names(sig_colors) <- c(paste0('FDR<', 0.1), paste0('FDR>', 0.1)); p <- p + scale_color_manual(values = sig_colors)
  p = p + scale_fill_gradientn(name= "protein change", na.value=NA, colours=RdBu1024, limit=c(-cap,cap))
  # p <- p + facet_grid(. ~ affectedGene.path, drop=T, space = "free",scales = "free")#, space = "free", scales = "free")
  p <- p + xlab("Associated Protein") + ylab("Cancer Type")
  p <- p + guides(colour = guide_legend(title = "P-value"))
  p <- p + theme_bw()
  p <- p + theme(panel.spacing.x = unit(0, "lines"), panel.spacing.y = unit(0, "lines"))
  p <- p + theme(axis.text.y = element_text(size = 12, face = "bold"))
  p <- p + theme(axis.text.x = element_text(size = 10, face = "bold", angle = 90, hjust = 0.5, vjust = 0.5))
  p <- p + theme(panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "grey80"))
  p
  fn = paste(makeOutDir(resultD = resultD), genoalt_type, "_sig_cancer_sig_thres", sig_p_thres, '_driver_enzyme_only_cleaned.png',sep = "")
  ggsave(filename = fn, width = 8, height = 2, device = png())
  dev.off()
  fn = paste(makeOutDir(resultD = resultD), genoalt_type, "_sig_cancer_sig_thres", sig_p_thres, '_driver_enzyme_only_cleaned.pdf',sep = "")
  ggsave(filename = fn, width = 8, height = 2)
  dev.off()
  
  sig_fdr_thres <- 0.1
  df <- df[!(df$affectedGene.path %in% c("NOTCH", "TP53", "Spliceosome")),]
  p <- ggplot()
  p <- p + geom_point(data = df, mapping = aes(x = affectedGene, y = cancer, fill = meddiff_capped, size = log10_Pvalue,
                                               colour = ifelse(fdr_sig == TRUE, paste0('FDR<', 0.1), paste0('FDR>', 0.1))), 
                      shape = 21, alpha = 0.8)
  sig_colors <- c("black", "white"); names(sig_colors) <- c(paste0('FDR<', sig_fdr_thres), paste0('FDR>', sig_fdr_thres)); p <- p + scale_color_manual(values = sig_colors)
  p <- p + facet_grid(. ~ affectedGene.path, drop=T, space = "free",scales = "free")#, space = "free", scales = "free")
  p = p + scale_fill_gradientn(name= "protein change", na.value=NA, colours=RdBu1024, limit=c(-cap,cap))
  p <- p + xlab("Affected Protein") + ylab("Cancer Type")
  p <- p + guides(colour = guide_legend(title = "FDR"), size = FALSE)
  p <- p + theme_bw()
  p <- p + theme(panel.spacing.x = unit(0, "lines"), panel.spacing.y = unit(0, "lines"))
  p <- p + theme(axis.text.y = element_text(size = 12, face = "bold"))
  p <- p + theme(axis.text.x = element_text(size = 10, face = "bold", angle = 90, hjust = 0.5, vjust = 0.5))
  p <- p + theme(panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "grey80"))
  p <- p + theme(strip.text.x = element_text(angle = 90))
  p
  fn = paste(makeOutDir(resultD = resultD), genoalt_type, "_sig_cancer_FDR", sig_fdr_thres, '.png',sep = "")
  ggsave(filename = fn, width = 9, height = 2.5, device = png())
  dev.off()
  fn = paste(makeOutDir(resultD = resultD), genoalt_type, "_sig_cancer_FDR", sig_fdr_thres, '.pdf',sep = "")
  ggsave(filename = fn, width = 9, height = 2.5)
  dev.off()
}


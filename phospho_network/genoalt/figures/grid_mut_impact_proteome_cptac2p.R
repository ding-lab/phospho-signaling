# Yige Wu @ WashU 2018 July
## overview figure for mutational impact of kinases on substrates

# source ------------------------------------------------------------------
setwd(dir = "Box Sync/")
source("./cptac2p_analysis/phospho_network/phospho_network_plotting.R")
source("./cptac2p_analysis/phospho_network/phospho_network_shared.R")
source("./cptac2p_analysis/dependencies/tables/load_TCGA_pathways.R")


# set variables -----------------------------------------------------------
show_p_thres <- 0.2
sig_p_thres <- 0.05

# inputs ------------------------------------------------------------------
# mut_cnv_cans <- fread(input = paste0(ppnD, "genoalt/tables/test_mut_cnv_impact_p_adjust/mut_cnv_cans_p_adjusted.csv"), data.table = F)
mut_cnv_cans <- fread(input = paste0(ppnD, "genoalt/tables/test_mut_cna_impact_cptac2p/mut_cnv_sig_cans.txt"), data.table = F)

## annotate the substrates
sub_genes2pathways <- map2TCGApathwaways(gene_list = unique(mut_cnv_cans$SUB_GENE), pathway_list = tcga_pathways_pluskegg_and_pathway)
sub_genes2pathways <- data.frame(SUB_GENE = rep(x = names(sub_genes2pathways), sapply(X = sub_genes2pathways, FUN = function(x) length(x))), 
                                 SUB_GENE.path = unlist(sub_genes2pathways, use.names = F))
sub_genes2pathways <- sub_genes2pathways[!(sub_genes2pathways$SUB_GENE == "APC" & sub_genes2pathways$SUB_GENE.path != "WNT") & !(sub_genes2pathways$SUB_GENE == "GSK3B" & sub_genes2pathways$SUB_GENE.path != "PI3K") & !(sub_genes2pathways$SUB_GENE == "IRS1" & sub_genes2pathways$SUB_GENE.path != "RTK RAS"),]
sub_genes2pathways <- sub_genes2pathways[!(sub_genes2pathways$SUB_GENE %in% c("MAP3K5", "MAP3K8", "BRAF", "ARAF", "RAF1", "MAP2K1", "MAP2K2", "MAP2K4", "MAP3K1", "MAPK1", "MAPK3", "MAPK9", "MAPK14", "MAPK13") & sub_genes2pathways$SUB_GENE.path != "MAPK"),]
mut_cnv_cans$SUB_GENE.path <- NULL
mut_cnv_cans <- merge(mut_cnv_cans, sub_genes2pathways, all.x = T)
mut_cnv_cans_sig <- mut_cnv_cans[mut_cnv_cans$p < sig_p_thres,]
for (genoalt_type in c("mutation")) {
  ## get the list of kinase/phosphatases to show
  genes2show <- unique(c(unlist(SMGs), as.vector(driver_genes$Gene[driver_genes$Cancer %in% cancers_sort])))
  enzymes2show_brca <- unique(mut_cnv_cans$GENE[mut_cnv_cans$p < sig_p_thres & mut_cnv_cans$p > 0 & mut_cnv_cans$cancer == "BRCA" & mut_cnv_cans$genoalt_type == genoalt_type & (mut_cnv_cans$GENE %in% genes2show | mut_cnv_cans$SUB_GENE %in% genes2show)])
  enzymes2show_ov <- unique(mut_cnv_cans$GENE[mut_cnv_cans$p < sig_p_thres & mut_cnv_cans$p > 0 & mut_cnv_cans$cancer == "OV" & mut_cnv_cans$genoalt_type == genoalt_type])
  enzymes2show_ov <- enzymes2show_ov[!(enzymes2show_ov %in% c("KMT2A", "HUWE1"))]
  enzymes2show_co <- unique(mut_cnv_cans$GENE[mut_cnv_cans$p < sig_p_thres & mut_cnv_cans$p > 0 & mut_cnv_cans$cancer == "CO" & mut_cnv_cans$genoalt_type == genoalt_type & (mut_cnv_cans$GENE %in% genes2show | mut_cnv_cans$SUB_GENE %in% genes2show)])
  enzymes2show_co <- intersect(enzymes2show_co, unique(c(SMGs[["CO"]], as.vector(driver_genes$Gene[driver_genes$Cancer == "CO"]))))
  # enzymes2show_brca <- unique(mut_cnv_cans$GENE[mut_cnv_cans$p < sig_p_thres & mut_cnv_cans$p > 0 & mut_cnv_cans$cancer == "BRCA" & mut_cnv_cans$genoalt_type == genoalt_type])
  # enzymes2show_brca <- intersect(enzymes2show_brca, unique(c(SMGs[["BRCA"]], as.vector(driver_genes$Gene[driver_genes$Cancer == "BRCA"]))))
  # enzymes2show_ov <- unique(mut_cnv_cans$GENE[mut_cnv_cans$p < sig_p_thres & mut_cnv_cans$p > 0 & mut_cnv_cans$cancer == "OV" & mut_cnv_cans$genoalt_type == genoalt_type])
  # enzymes2show_co <- unique(mut_cnv_cans$GENE[mut_cnv_cans$p < sig_p_thres & mut_cnv_cans$p > 0 & mut_cnv_cans$cancer == "CO" & mut_cnv_cans$genoalt_type == genoalt_type])
  # enzymes2show_co <- intersect(enzymes2show_co, unique(c(SMGs[["CO"]], as.vector(driver_genes$Gene[driver_genes$Cancer == "CO"]))))
  # enzymes2show_co <- c(enzymes2show_co, "CDK12")
  enzymes2show <- unique(c(enzymes2show_brca, enzymes2show_ov, enzymes2show_co))
  enzymes2show <- unique(c(enzymes2show, "RAF1", "AHNAK", "AHNAK2", "CREBBP", "EIF4G1", "PLEC", "PRKDC", "ROCK1", "TSC1", "TSC2", "SOS1", "ERBB3","AR", "INSR", "MTOR", "DNMT1", "ARAF", "ATR", "PIK3C3", "RFC1", "MAP2K1", "ATM", "CDK12", "EP300", "ERCC5", "TGFBR2", "ERBB2", "MAP3K4", "MAP3K5", "HDAC4", "MARK2", "MARK4", "RPS6KA2", "RPS6KA4", "MAP3K11", "FGFR4", "JAK1", "ACVR2A", "RPTOR", "PPM1D"))
  # enzymes2show <- enzymes2show[!(enzymes2show %in% c(""))]
  
  ## get the list of substrates to show
  subtrates2show_brca <- unique(mut_cnv_cans$SUB_GENE[mut_cnv_cans$p < sig_p_thres & mut_cnv_cans$p > 0 & mut_cnv_cans$cancer == "BRCA" & mut_cnv_cans$genoalt_type == genoalt_type & mut_cnv_cans$GENE %in% enzymes2show & (!is.na(mut_cnv_cans$SUB_GENE.path) | mut_cnv_cans %in% genes2show)])
  subtrates2show_brca <- c(subtrates2show_brca, "GSK3B", "FOXO3", "MAP2K4")
  subtrates2show_ov <- unique(mut_cnv_cans$SUB_GENE[mut_cnv_cans$p < sig_p_thres & mut_cnv_cans$p > 0 & mut_cnv_cans$cancer == "OV" & mut_cnv_cans$genoalt_type == genoalt_type])
  subtrates2show_co <- unique(mut_cnv_cans$SUB_GENE[mut_cnv_cans$p < sig_p_thres & mut_cnv_cans$p > 0 & mut_cnv_cans$cancer == "CO" & mut_cnv_cans$genoalt_type == genoalt_type & (!is.na(mut_cnv_cans$SUB_GENE.path) | mut_cnv_cans %in% genes2show) & mut_cnv_cans$GENE %in% enzymes2show])
  subtrates2show_co <- subtrates2show_co[!(subtrates2show_co %in% c("ARHGAP35", "HSP90AA1", "PLA2G4A", "RALBP1", "YAP1"))]
  subtrates2show_co <- c(subtrates2show_co, "TP53", "PTK2")
  subtrates2show <- unique(c(subtrates2show_brca, subtrates2show_ov, subtrates2show_co))
  subtrates2show <- subtrates2show[!(subtrates2show %in% c("ARHGAP35", "HNRNPU" , "HNRNPA1", "PML", "SRC", "MAPK14", "MAPK13", "MAPK9", "KMT2A", "TOP2B"))]
  
  df <- mut_cnv_cans[mut_cnv_cans$genoalt_type == genoalt_type,]
  df <- df[(df$GENE %in% enzymes2show) & (df$SUB_GENE %in% subtrates2show),]
  df$sig <- (df$p < sig_p_thres)

  ## remove the duplicated pathway annotation
  tmp <- as.vector(df$SUB_GENE.path)
  tmp[is.na(tmp)] <- "other"
  tmp[tmp == "TP53"] <- "Cell Cycle"
  tmp[tmp == "TGF-Beta"] <- "other"
  tmp[tmp == "Spliceosome"] <- "other"
  tmp[tmp == "Mismatch repair"] <- "MMR"
  
  df$SUB_GENE.path <- tmp
  df <- df[order(df$SUB_GENE.path, df$SUB_GENE),]
  df$pair_can  <- paste0(df$pair, ":", df$cancer)
  df <- df[!duplicated(df$pair_can),]
  df <- df[order(df$p, decreasing = T),]
  
  df$log10_Pvalue <- -log10(df$p)
  cap <- min(max(abs(df$meddiff), na.rm = T),2)
  df$meddiff_capped <- df$meddiff
  df$meddiff_capped[df$meddiff > cap] <- cap
  df$meddiff_capped[df$meddiff < (-cap)] <- (-cap)
  df$driver_gene_type.en[df$driver_gene_type.en == ""] <- "other"
  df$driver_gene_type.en[df$GENE %in% c("INSR", "PPM1D")] <- "oncogene"
  
  df$driver_gene_type.en <- factor(df$driver_gene_type.en, levels = c("oncogene", "tsg", "other"))
  df$SUB_GENE.path <- factor(df$SUB_GENE.path, levels = c(names(tcga_pathways_pluskegg_and_pathway), "MMR", "other"))
  
  sig_colors <- c("black", "white")
  names(sig_colors) <- c(paste0('P-value<', sig_p_thres), paste0('P-value>', sig_p_thres))

  
  # p <- ggplot()
  # p <- p + geom_point(data = df, mapping = aes(x = GENE, y = SUB_GENE, fill = meddiff_capped, size = log10_Pvalue, 
  #                                              colour = ifelse(sig == TRUE, paste0('P-value<', sig_p_thres), paste0('P-value>', sig_p_thres))), shape = 21, alpha = 0.6)
  # p <- p + scale_color_manual(values = sig_colors)
  # p = p + scale_fill_gradientn(name= "phosphorylation change", na.value=NA, colours=RdBu1024, limit=c(-cap,cap))
  # p <- p + facet_grid(SUB_GENE.path ~ cancer + driver_gene_type.en, drop=T, space = "free",scales = "free")#, space = "free", scales = "free")
  # p <- p + theme_bw()
  # p <- p + theme(panel.spacing.x = unit(0, "lines"), panel.spacing.y = unit(0, "lines"))
  # p <- p + theme(axis.text.y = element_text(size = 10, face = "bold"))
  # p <- p + theme(axis.text.x = element_text(size = 10, face = "bold", angle = 90, hjust = 0.5, vjust = 0.5))
  # p
  # fn = paste(makeOutDir(resultD = resultD), genoalt_type, "_sig_cancer_sig_thres", sig_p_thres, '.pdf',sep = "")
  # ggsave(filename = fn, width = 14, height = 10, useDingbats = F)
  

  p <- ggplot()
  p <- p + geom_point(data = df, mapping = aes(x = SUB_GENE, y = GENE, fill = meddiff_capped, size = log10_Pvalue, 
                                               colour = ifelse(sig == TRUE, paste0('P-value<', sig_p_thres), paste0('P-value>', sig_p_thres))), shape = 21, alpha = 0.6)
  p <- p + scale_color_manual(values = sig_colors)
  p = p + scale_fill_gradientn(name= "phosphorylation change", na.value=NA, colours=RdBu1024, limit=c(-cap,cap))
  p <- p + facet_grid(cancer + driver_gene_type.en ~ SUB_GENE.path, drop=T, space = "free",scales = "free")#, space = "free", scales = "free")
  p <- p + xlab("Affected Phosphoprotein") + ylab("Mutated Gene")
  p <- p + guides(colour = guide_legend(title = "P-value"))
  p <- p + theme_bw()
  p <- p + theme(panel.spacing.x = unit(0, "lines"), panel.spacing.y = unit(0, "lines"))
  p <- p + theme(axis.text.y = element_text(size = 10, face = "bold"))
  p <- p + theme(axis.text.x = element_text(size = 10, face = "bold", angle = 90, hjust = 0.5, vjust = 0.5))
  p
  # fn = paste(makeOutDir(resultD = resultD), genoalt_type, "_sig_cancer_sig_thres", sig_p_thres, '.png',sep = "")
  # ggsave(filename = fn, width = 14, height = 10, device = png())
  
  
  p <- ggplot()
  p <- p + geom_point(data = df[(df$driver_gene_type.en %in% c("oncogene", "tsg") & df$cancer %in% c("BRCA", "CO")) | df$cancer == "OV",], mapping = aes(x = SUB_GENE, y = GENE, fill = meddiff_capped, size = log10_Pvalue, 
                                               colour = ifelse(sig == TRUE, paste0('P-value<', sig_p_thres), paste0('P-value>', sig_p_thres))), shape = 21, alpha = 0.6)
  p <- p + scale_color_manual(values = sig_colors)
  p = p + scale_fill_gradientn(name= "phosphorylation change", na.value=NA, colours=RdBu1024, limit=c(-cap,cap))
  p <- p + facet_grid(cancer + driver_gene_type.en ~ SUB_GENE.path, drop=T, space = "free",scales = "free")#, space = "free", scales = "free")
  p <- p + xlab("Associated Phosphoprotein") + ylab("Mutated Gene")
  p <- p + guides(colour = guide_legend(title = "P-value"))
  p <- p + theme_bw()
  p <- p + theme(panel.spacing.x = unit(0, "lines"), panel.spacing.y = unit(0, "lines"))
  p <- p + theme(axis.text.y = element_text(size = 12, face = "bold"))
  p <- p + theme(axis.text.x = element_text(size = 10, face = "bold", angle = 90, hjust = 0.5, vjust = 0.5))
  p <- p + theme(panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "grey80"))
  p
  fn = paste(makeOutDir(resultD = resultD), genoalt_type, "_sig_cancer_sig_thres", sig_p_thres, '_driver_enzyme_only.png',sep = "")
  ggsave(filename = fn, width = 12, height = 8, device = png())
  
  fn = paste(makeOutDir(resultD = resultD), genoalt_type, "_sig_cancer_sig_thres", sig_p_thres, '_driver_enzyme_only.pdf',sep = "")
  ggsave(filename = fn, width = 12, height = 8)
}
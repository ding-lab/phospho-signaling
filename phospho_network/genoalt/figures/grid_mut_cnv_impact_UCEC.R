# Yige Wu @ WashU 2018 July
## overview figure for mutational impact of kinases on substrates

# source ------------------------------------------------------------------
source("Box Sync/cptac2p_analysis/phospho_network/phospho_network_plotting.R")
source("./cptac2p_analysis/phospho_network/phospho_network_shared.R")
source("./cptac2p_analysis/dependencies/tables/load_TCGA_pathways.R")


# set variables -----------------------------------------------------------
show_p_thres <- 0.2
sig_p_thres <- 0.05
driver_genes2show <- unique(c(as.vector(driver_genes$Gene[driver_genes$Cancer == "UCEC"]), SMGs[["UCEC"]]))

# inputs ------------------------------------------------------------------
mut_cnv_cans <- fread(input = paste0(ppnD, "genoalt/tables/test_mut_cna_impact_UCEC/mut_cnv_sig_cans.txt"), data.table = F)
# mut_cnv_cans2 <- fread(input = paste0(ppnD, "genoalt/tables/test_mut_cna_impact_UCEC/UCEC_mut_cnv_tab.txt"), data.table = F)

## annotate the substrates
sub_genes2pathways <- map2TCGApathwaways(gene_list = unique(mut_cnv_cans$SUB_GENE), pathway_list = tcga_pathways_pluskegg_and_pathway)
sub_genes2pathways <- data.frame(SUB_GENE = rep(x = names(sub_genes2pathways), sapply(X = sub_genes2pathways, FUN = function(x) length(x))), 
                                 SUB_GENE.path = unlist(sub_genes2pathways, use.names = F))
sub_genes2pathways <- sub_genes2pathways[!(sub_genes2pathways$SUB_GENE == "APC" & sub_genes2pathways$SUB_GENE.path != "WNT") & !(sub_genes2pathways$SUB_GENE == "GSK3B" & sub_genes2pathways$SUB_GENE.path != "PI3K") & !(sub_genes2pathways$SUB_GENE == "IRS1" & sub_genes2pathways$SUB_GENE.path != "RTK RAS") & !(sub_genes2pathways$SUB_GENE == "TP53" & sub_genes2pathways$SUB_GENE.path != "TP53"),]
mut_cnv_cans$SUB_GENE.path <- NULL
mut_cnv_cans <- merge(mut_cnv_cans, sub_genes2pathways, all.x = T)

for (genoalt_type in c("mut")) {
  ## get the list of kinase/phosphatases to show
  enzymes2show <- unique(mut_cnv_cans$GENE[mut_cnv_cans$p < sig_p_thres & mut_cnv_cans$p > 0 & mut_cnv_cans$cancer == "UCEC" & mut_cnv_cans$genoalt_type == genoalt_type & (mut_cnv_cans$GENE %in% driver_genes2show | mut_cnv_cans$driver_gene_type.en != "")])
  enzymes2show <- enzymes2show[!(enzymes2show %in% c("ACVR2A", "ARHGAP35", "CDK12", "EPHA2", "FLNA", "HUWE1", "INPPL1", "NCOR1", "RASA1", "SED2", "STAG2", "TSC1",
                                                     "ZFHX3", "ZFP36L1", "NFE2L2", "PDGFRA", "PPP2R1A"))]
  
  ## get the list of substrates to show
  subtrates2show <- unique(mut_cnv_cans$SUB_GENE[mut_cnv_cans$p < sig_p_thres & mut_cnv_cans$p > 0 & mut_cnv_cans$cancer == "UCEC" & mut_cnv_cans$genoalt_type == genoalt_type & mut_cnv_cans$GENE %in% enzymes2show & (!is.na(mut_cnv_cans$SUB_GENE.path) | mut_cnv_cans$SUB_GENE %in% driver_genes2show)])
  subtrates2show <- c(subtrates2show, "STAT3")
  df <- mut_cnv_cans[mut_cnv_cans$genoalt_type == genoalt_type,]
  df <- df[(df$GENE %in% enzymes2show) & (df$SUB_GENE %in% subtrates2show),]
  df$sig <- (df$p < sig_p_thres)
  # df <- df[!(df$GENE == "CTNNB1" & df$SUB_GENE %in% c("APC", "AXIN1")),]

  ## remove the duplicated pathway annotation
  tmp <- as.vector(df$SUB_GENE.path)
  tmp[is.na(tmp)] <- "other"
  tmp[tmp == "Mismatch repair"] <- "MMR"
  tmp[df$SUB_GENE == "STAT3"] <- "JAK-STAT"
  tmp[df$SUB_GENE == "AXIN1"] <- "WNT"
  
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
  # df$SUB_GENE.path <- factor(df$SUB_GENE.path, levels = c(names(tcga_pathways_pluskegg_and_pathway), "other"))
  
  sig_colors <- c("black", "white")
  names(sig_colors) <- c(paste0('P-value<', sig_p_thres), paste0('P-value>', sig_p_thres))

  df <- df[!(df$SUB_GENE %in% c("NCOR1", "ZFYVE16", "ROCK1", "RAF1", "MAPK1")) & (df$SUB_GENE.path != "other"),]
  tmp <- as.vector(unique(df$GENE))
  df$GENE <- factor(df$GENE, levels = c(c("CTNNB1", "APC", "PTEN", "TP53", "MSH6"), tmp[!(tmp %in% c("CTNNB1", "APC", "PTEN", "TP53", "MSH6"))]))
  p <- ggplot()
  p <- p + geom_point(data = df, mapping = aes(x = SUB_GENE, y = GENE, fill = meddiff_capped, size = log10_Pvalue,
                                               colour = ifelse(sig == TRUE, paste0('P-value<', sig_p_thres), paste0('P-value>', sig_p_thres))), shape = 21, alpha = 0.6)
  p <- p + scale_color_manual(values = sig_colors)
  p = p + scale_fill_gradientn(name= "phosphorylation change", na.value=NA, colours=RdBu1024, limit=c(-cap,cap))
  p <- p + facet_grid(driver_gene_type.en ~ SUB_GENE.path, drop=T, space = "free",scales = "free")#, space = "free", scales = "free")
  p <- p + xlab("Associated Phosphoprotein") + ylab("Mutated Gene")
  p <- p + guides(colour = guide_legend(title = "P-value"))
  p <- p + theme_bw()
  p <- p + theme(panel.spacing.x = unit(0, "lines"), panel.spacing.y = unit(0, "lines"))
  p <- p + theme(axis.text.y = element_text(size = 12, face = "bold"))
  p <- p + theme(axis.text.x = element_text(size = 10, face = "bold", angle = 90, hjust = 0.5, vjust = 0.5))
  p <- p + theme(panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "grey80"))
  p <- p + theme(strip.text.x = element_text(angle = 90))
  p
  # fn = paste(makeOutDir(resultD = resultD), genoalt_type, "_sig_cancer_sig_thres", sig_p_thres, '_driver_enzyme_only_cleaned.png',sep = "")
  # ggsave(filename = fn, width = 12, height = 8, device = png())
  
  fn = paste(makeOutDir(resultD = resultD), genoalt_type, "_sig_cancer_sig_thres", sig_p_thres, '_driver_enzyme_only_cleaned.pdf',sep = "")
  ggsave(filename = fn, width = 12, height = 8)
 }
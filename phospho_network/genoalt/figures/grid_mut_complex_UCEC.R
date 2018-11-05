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
driver_genes2show <- c(driver_genes2show, "TP53", "JAK1", "STAT1", "STAT3", "APC", "CTNNB1")

# inputs ------------------------------------------------------------------
## input mutational association table
mut_cnv_cans <- fread(input = paste0(ppnD, "genoalt/tables/test_mut_complex_UCEC/mut_cnv_sig_cans.txt"), data.table = F)
# mut_cnv_cans <- fread(input = paste0(ppnD, "genoalt/tables/test_mut_complex_UCEC/UCEC_mut_cnv_tab.txt"), data.table = F)
mut_cnv_cans$cancer <- "UCEC"
## annotate the substrates
sub_genes2pathways <- map2TCGApathwaways(gene_list = unique(mut_cnv_cans$geneB), pathway_list = tcga_pathways_pluskegg_and_pathway)
sub_genes2pathways <- data.frame(geneB = rep(x = names(sub_genes2pathways), sapply(X = sub_genes2pathways, FUN = function(x) length(x))), 
                                 geneB.path = unlist(sub_genes2pathways, use.names = F))
sub_genes2pathways <- sub_genes2pathways[!(sub_genes2pathways$geneB == "APC" & sub_genes2pathways$geneB.path != "WNT") & !(sub_genes2pathways$geneB == "GSK3B" & sub_genes2pathways$geneB.path != "PI3K") & !(sub_genes2pathways$geneB == "IRS1" & sub_genes2pathways$geneB.path != "RTK RAS"),]
mut_cnv_cans <- merge(mut_cnv_cans, sub_genes2pathways, all.x = T)

for (genoalt_type in c("mutation")) {
  ## get the list of kinase/phosphatases to show
  # enzymes2show <- unique(mut_cnv_cans$geneA[mut_cnv_cans$p < sig_p_thres & mut_cnv_cans$p > 0 & mut_cnv_cans$genoalt_type == genoalt_type])
  enzymes2show <- unique(mut_cnv_cans$geneA[mut_cnv_cans$p < sig_p_thres & mut_cnv_cans$p > 0 & mut_cnv_cans$cancer == "UCEC" & mut_cnv_cans$genoalt_type == genoalt_type & (mut_cnv_cans$geneA %in% driver_genes2show | mut_cnv_cans$driver_gene_type.en != "")])
  # enzymes2show <- c(enzymes2show)

  ## get the list of substrates to show
  # subtrates2show <- unique(mut_cnv_cans$geneB[mut_cnv_cans$p < sig_p_thres & mut_cnv_cans$p > 0 & mut_cnv_cans$genoalt_type == genoalt_type & mut_cnv_cans$geneA %in% enzymes2show])
  subtrates2show <- unique(mut_cnv_cans$geneB[mut_cnv_cans$p < sig_p_thres & mut_cnv_cans$p > 0 & mut_cnv_cans$cancer == "UCEC" & mut_cnv_cans$genoalt_type == genoalt_type & mut_cnv_cans$geneA %in% enzymes2show & (!is.na(mut_cnv_cans$geneB.path) | mut_cnv_cans$geneB %in% driver_genes2show)])
  subtrates2show <- c(subtrates2show, "TP53BP1", "STAT3")
  
  df <- mut_cnv_cans[mut_cnv_cans$genoalt_type == genoalt_type,]
  df <- df[(df$geneA %in% enzymes2show) & (df$geneB %in% subtrates2show),]
  df$sig <- (df$p < sig_p_thres)

  ## remove the duplicated pathway annotation
  tmp <- as.vector(df$geneB.path)
  tmp[is.na(tmp)] <- "other"
  tmp[tmp == "Mismatch repair"] <- "MMR"
  tmp[df$geneB == "AXIN1"] <- "WNT"
  tmp[df$geneB == "STAT3"] <- "JAK-STAT"
  
  df$geneB.path <- tmp
  df <- df[order(df$geneB.path, df$geneB),]
  df$pair_can  <- paste0(df$pair, ":", df$cancer)
  df <- df[!duplicated(df$pair_can),]
  df <- df[order(df$p, decreasing = T),]
  
  df$log10_Pvalue <- -log10(df$p)
  cap <- min(max(abs(df$meddiff), na.rm = T),2)
  df$meddiff_capped <- df$meddiff
  df$meddiff_capped[df$meddiff > cap] <- cap
  df$meddiff_capped[df$meddiff < (-cap)] <- (-cap)
  df$driver_gene_type.en[df$driver_gene_type.en == ""] <- "other"
  df$driver_gene_type.en[df$geneA %in% c("INSR", "PPM1D")] <- "oncogene"
  
  df$driver_gene_type.en <- factor(df$driver_gene_type.en, levels = c("oncogene", "tsg", "other"))
  # df$geneB.path <- factor(df$geneB.path, levels = c(names(tcga_pathways_pluskegg_and_pathway), "other"))
  
  sig_colors <- c("black", "white")
  names(sig_colors) <- c(paste0('P-value<', sig_p_thres), paste0('P-value>', sig_p_thres))

  # df <- df[!(df$geneB %in% c("NCOR1", "ZFYVE16", "ROCK1", "RAF1", "MAPK1")) & (df$geneB.path != "other"),]
  p <- ggplot()
  p <- p + geom_point(data = df, mapping = aes(x = geneB, y = geneA, fill = meddiff_capped, size = log10_Pvalue, 
                                               colour = ifelse(sig == TRUE, paste0('P-value<', sig_p_thres), paste0('P-value>', sig_p_thres))), shape = 21, alpha = 0.6)
  p <- p + scale_color_manual(values = sig_colors)
  p = p + scale_fill_gradientn(name= "phosphorylation change", na.value=NA, colours=RdBu1024, limit=c(-cap,cap))
  p <- p + facet_grid(driver_gene_type.en ~ geneB.path, drop=T, space = "free",scales = "free")#, space = "free", scales = "free")
  p <- p + xlab("Associated Phosphoprotein") + ylab("Mutated Gene")
  p <- p + guides(colour = guide_legend(title = "P-value"))
  p <- p + theme_bw()
  p <- p + theme(panel.spacing.x = unit(0, "lines"), panel.spacing.y = unit(0, "lines"))
  p <- p + theme(axis.text.y = element_text(size = 12, face = "bold"))
  p <- p + theme(axis.text.x = element_text(size = 10, face = "bold", angle = 90, hjust = 0.5, vjust = 0.5))
  p <- p + theme(panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "grey80"))
  p
  fn = paste(makeOutDir(resultD = resultD), genoalt_type, "_sig_cancer_sig_thres", sig_p_thres, '.png',sep = "")
  ggsave(filename = fn, width = 15, height = 10, device = png())
 }
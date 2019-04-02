# Yige Wu @ WashU 2018 July
## overview figure for mutational impact of kinases on substrates

# source ------------------------------------------------------------------
setwd(dir = "Box Sync/")
source("./cptac2p_analysis/phospho_network/phospho_network_plotting.R")
source("./cptac2p_analysis/phospho_network/phospho_network_shared.R")
source("./cptac2p_analysis/dependencies/tables/load_TCGA_pathways.R")


# set variables -----------------------------------------------------------
sig_p_thres <- 0.05
num_genoalt_thres <- 4
gene_altered <- "TP53"
affected_exp_type <- "PHO"
affected_exp_type <- "PRO_PHO"
affected_exp_type <- "RNA_PRO"
cancers2process <- c("BRCA", "OV", "CO", "UCEC", "CCRCC", "LIHC")
interaction_type <- "TF"
do.smg <- F
sig_colors <- c("black", "white")
cap <- 2

mut_impact_tab <- load_mut_impact_proteome()


# draw kinases ------------------------------------------------------------
for (variant_class in c("not_silent", "missense", "truncation")) {
  tab2p <- mut_impact_tab
  tab2p <- tab2p %>%
    filter(GENE == gene_altered) %>%
    filter(SUB_GENE.is_kinase == T) %>%
    filter(fdr < 0.05) %>%
    filter(variant_class == variant_class) %>%
    mutate(log10_fdr = -log10(fdr))
  tab2p$meddiff_capped <- add_cap(x = tab2p$meddiff, cap = cap)
  
  fn = paste(makeOutDir(resultD = resultD), gene_altered, "_", variant_class, "_impact_", "kinase" , "_num_genoalt_thres_", num_genoalt_thres, '.pdf',sep = "")
  
  p <- ggplot()
  p <- p + geom_point(data = tab2p, mapping = aes(x = SUB_GENE, y = cancer, fill = meddiff_capped, size = log10_fdr,
                                                  colour = ifelse(sig == TRUE, paste0('FDR<', sig_p_thres), paste0('FDR>', sig_p_thres))), alpha = 0.6, shape = 21)
  p <- p + facet_grid(affected_exp_type ~ ., drop=T, space = "free",scales = "free")#, space = "free", scales = "free")
  p <- p + scale_color_manual(values = sig_colors)
  p = p + scale_fill_gradientn(name= "RNA/protein\nchange", na.value=NA, colours=RdBu1024, limit=c(-cap,cap))
  p <- p + xlab("Affected gene") + ylab("Mutated Gene")
  p <- p + guides(colour = guide_legend(title = "FDR"))
  p <- p + theme_bw()
  p <- p + theme(panel.spacing.x = unit(0, "lines"), panel.spacing.y = unit(0, "lines"))
  p <- p + theme(axis.text.y = element_text(size = 10, face = "bold"))
  p <- p + theme(axis.text.x = element_text(size = 10, face = "bold", angle = 90, hjust = 0.5, vjust = 0.5))
  p <- p + theme(panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "grey80"))
  p <- p + theme(strip.text.x = element_text(angle = 90))
  p
  ggsave(filename = fn, width = 8, height = 5)
}

# draw phosphatase ------------------------------------------------------------
for (variant_class in c("not_silent", "missense", "truncation")) {
  tab2p <- mut_impact_tab
  tab2p <- tab2p %>%
    filter(GENE == gene_altered) %>%
    filter(SUB_GENE.is_phosphatase == T) %>%
    filter(fdr < 0.05) %>%
    filter(variant_class == variant_class) %>%
    mutate(log10_fdr = -log10(fdr))
  tab2p$meddiff_capped <- add_cap(x = tab2p$meddiff, cap = cap)
  
  fn = paste(makeOutDir(resultD = resultD), gene_altered, "_", variant_class, "_impact_", "phosphatase" , "_num_genoalt_thres_", num_genoalt_thres, '.pdf',sep = "")
  
  p <- ggplot()
  p <- p + geom_point(data = tab2p, mapping = aes(x = SUB_GENE, y = cancer, fill = meddiff_capped, size = log10_fdr,
                                                  colour = ifelse(sig == TRUE, paste0('FDR<', sig_p_thres), paste0('FDR>', sig_p_thres))), alpha = 0.6, shape = 21)
  p <- p + facet_grid(affected_exp_type ~ ., drop=T, space = "free",scales = "free")#, space = "free", scales = "free")
  p <- p + scale_color_manual(values = sig_colors)
  p = p + scale_fill_gradientn(name= "RNA/protein\nchange", na.value=NA, colours=RdBu1024, limit=c(-cap,cap))
  p <- p + xlab("Affected gene") + ylab("Mutated Gene")
  p <- p + guides(colour = guide_legend(title = "FDR"))
  p <- p + theme_bw()
  p <- p + theme(panel.spacing.x = unit(0, "lines"), panel.spacing.y = unit(0, "lines"))
  p <- p + theme(axis.text.y = element_text(size = 10, face = "bold"))
  p <- p + theme(axis.text.x = element_text(size = 10, face = "bold", angle = 90, hjust = 0.5, vjust = 0.5))
  p <- p + theme(panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "grey80"))
  p <- p + theme(strip.text.x = element_text(angle = 90))
  p
  ggsave(filename = fn, width = 8, height = 5)
}


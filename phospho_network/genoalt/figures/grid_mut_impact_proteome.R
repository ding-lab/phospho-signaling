# Yige Wu @ WashU 2018 July
## overview figure for mutational impact of kinases on substrates

# source ------------------------------------------------------------------
baseD = "/Users/yigewu/Box\ Sync/"
wd <- getwd()
if (wd != baseD) {
  setwd(baseD)
}
source("./cptac2p_analysis/phospho_network/phospho_network_plotting.R")
source("./cptac2p_analysis/phospho_network/phospho_network_shared.R")
source("./cptac2p_analysis/dependencies/tables/load_TCGA_pathways.R")


# set variables -----------------------------------------------------------
show_p_thres <- 0.2
sig_p_thres <- 0.05
num_genoalt_thres <- 4
cancers2process <- c("BRCA", "OV", "CO", "UCEC", "CCRCC", "LIHC")
affected_exp_type <- "PRO"
do.smg <- TRUE
pathway_order <- c("TP53", "Cell Cycle","WNT", "PI3K","RTK RAS", "MAPK", "SWI/SNF", "NOTCH", "MMR", "HIPPO", "TGF-Beta", "other")

# inputs mutation impact proteome test results and reformat ------------------------------------------------------------------
mut_cnv_cans <- fread(input = paste0(ppnD, "genoalt/tables/test_mut_impact_proteome/mut_impact_proteome_RNA_cptac2p_cptac3_tab.txt"), data.table = F, sep = "\t")
mut_cnv_cans <- mut_cnv_cans[mut_cnv_cans$num >= num_genoalt_thres,]
mut_cnv_cans$pair_pro <- paste0(mut_cnv_cans$GENE, ":",  mut_cnv_cans$SUB_GENE)

## annotate the substrates to oncogenic pathways
sub_genes2pathways <- map2TCGApathwaways(gene_list = unique(mut_cnv_cans$SUB_GENE), pathway_list = tcga_pathways_pluskegg_and_pathway)
sub_genes2pathways <- data.frame(SUB_GENE = rep(x = names(sub_genes2pathways), sapply(X = sub_genes2pathways, FUN = function(x) length(x))), 
                                 SUB_GENE.path = unlist(sub_genes2pathways, use.names = F))
### correct substrate pathway annotation
tmp <- as.vector(sub_genes2pathways$SUB_GENE.path)
tmp[is.na(tmp)] <- "other"
tmp[sub_genes2pathways$SUB_GENE.path == "Mismatch repair"] <- "MMR"
tmp[sub_genes2pathways$SUB_GENE %in% c("STAT3", "JAK1")] <- "JAK-STAT"
tmp[sub_genes2pathways$SUB_GENE %in% c("ARID1A", "SMARCA4", "PBRM1", "ARID1B")] <- "SWI/SNF"
tmp[sub_genes2pathways$SUB_GENE %in% c("KRAS", "IRS1", "BRAF", "ARAF", "RAF1")] <- "RTK RAS"
tmp[sub_genes2pathways$SUB_GENE %in% c("AXIN1", "APC", "TCF7L2")] <- "WNT"
tmp[sub_genes2pathways$SUB_GENE %in% c("BCL2", "BCL2L1", "TP53")] <- "TP53"
tmp[sub_genes2pathways$SUB_GENE %in% c("MAP3K8", "MAP3K5", "MAP3K4",
                                       "MAP2K1", "MAP2K2", "MAP2K4", "MAP3K1", "MAPK1", "MAPK3", "MAPK9", "MAPK14", "MAPK13")] <- "MAPK"
tmp[sub_genes2pathways$SUB_GENE %in% c("GSK3B")] <- "PI3K"

sub_genes2pathways$SUB_GENE.path <- tmp
sub_genes2pathways <- unique(sub_genes2pathways)

mut_cnv_cans$SUB_GENE.path <- NULL
mut_cnv_cans <- merge(mut_cnv_cans, sub_genes2pathways, all.x = T)
tmp <- as.vector(mut_cnv_cans$SUB_GENE.path)
tmp[is.na(tmp)] <- "other"
tmp[mut_cnv_cans$SUB_GENE %in% c("STAT3", "JAK1")] <- "JAK-STAT"
tmp[mut_cnv_cans$SUB_GENE %in% c("ARID1A", "SMARCA4", "PBRM1", "ARID1B")] <- "SWI/SNF"
tmp[mut_cnv_cans$SUB_GENE %in% c("KRAS", "IRS1", "BRAF", "ARAF", "RAF1")] <- "RTK RAS"
tmp[mut_cnv_cans$SUB_GENE %in% c("AXIN1", "APC", "TCF7L2")] <- "WNT"
tmp[mut_cnv_cans$SUB_GENE %in% c("BCL2", "BCL2L1", "TP53")] <- "TP53"
tmp[mut_cnv_cans$SUB_GENE %in% c("MAP3K8", "MAP3K5", "MAP3K4",
                                       "MAP2K1", "MAP2K2", "MAP2K4", "MAP3K1", "MAPK1", "MAPK3", "MAPK9", "MAPK14", "MAPK13")] <- "MAPK"
tmp[mut_cnv_cans$SUB_GENE %in% c("GSK3B")] <- "PI3K"
mut_cnv_cans$SUB_GENE.path <- tmp

## annotate the altered gene to be oncogene or tsgs
mut_cnv_cans$driver_gene_type.en <- ""
mut_cnv_cans$driver_gene_type.en[mut_cnv_cans$GENE %in% oncogenes] <- "oncogene"
mut_cnv_cans$driver_gene_type.en[mut_cnv_cans$GENE %in% tsgs] <- "tsg"

## annotate with expression type
mut_cnv_cans$expression_type <- mut_cnv_cans$SUB_MOD_RSD
mut_cnv_cans$expression_type[!(mut_cnv_cans$SUB_MOD_RSD %in% c("RNA", "PRO"))] <- "PHO"

for (variant_class in c("not_silent")) {
  gene_altered2show <- unique(c(unlist(SMGs[cancers2process]), driver_genes$Gene[driver_genes$Cancer %in% c(cancers2process, "KIRC")]))
  gene_altered_had2show <- c("ARID1B:UCEC", "MAP3K1:UCEC", "SMARCA4:UCEC")
  gene_altered_not_show <- c("ACVR2A", "ARHGAP35", "CDK12", "EPHA2", "FLNA", "HUWE1", "INPPL1", "NCOR1", "RASA1", "SED2", "STAG2", "TSC1",
                         "ZFHX3", "ZFP36L1", "NFE2L2", "PDGFRA", "PPP2R1A", "ZMYM2", "PDS5B", "GNAS", "BCOR", "ATF7IP", "SMC1A", "DHX9", "IL6ST", "XPO1", "JAK1")
  gene_affected2show <- c("STAT3", "RB1", "ARID1A", "CDK1", "CDK2", "CDKN1A", "CCNA2", "TP53BP1", "MAPKAPK5", "AURKA", "LEF1", "SMARCC2", "SMARCA4", "SMARCE1")
  gene_affected_not_show <- c("NCOR1", "ZFYVE16", "ROCK1", "NCOR2", "SAP30", "ZMYM2", "WWTR1", "YAP1", "FAS", 
                              "SKP1", "SMAD3", "YWHAE", "PPP3CB", "PPP2R5C", "PPP2R5A", "SOS1", "THBS1", "TSC2", 
                              "SERPINE1", "SIN3A", "RPL22", "RBL22", "RRM2B", "HDAC1", "CREBBP", "HNRNPA1")
  
  for (genoalt_type in c("mut")) {
    tab4cytoscape <- NULL
    for (affected_exp_type in c("PRO", "PHO")) {
      df <- mut_cnv_cans
      df <- df[df$expression_type == affected_exp_type,]
      df <- df[df$cancer %in% cancers2process,]
      df <- df[df$variant_class == variant_class,]
      df$GENE_cancer <- paste0(df$GENE, ":", df$cancer)
      df_origin <- df
      
      if (do.smg) {
        df_smg <- NULL
        for (cancer in cancers2process) {
          df_new <- df[df$cancer == cancer & df$GENE %in% SMGs[[cancer]],]
          df_smg <- rbind(df_smg, df_new)
        }
        df <- df_smg
        
        enzymes2show <- unique(df$GENE[df$p < sig_p_thres &  (df$GENE %in% gene_altered2show)])
        enzymes2show <- c(enzymes2show, gene_altered_had2show)
        
        subtrates2show <- unique(df$SUB_GENE[df$p < sig_p_thres & (df$SUB_GENE.path != "other" | df$SUB_GENE %in% gene_altered2show)])
        subtrates2show <- c(subtrates2show, gene_altered_had2show)
        subtrates2show <- subtrates2show[!(subtrates2show %in% gene_affected_not_show)]
        
        ## do the data frame for plotting
        df <- df[(df$GENE %in% enzymes2show) & (df$SUB_GENE %in% subtrates2show),]
        df <- rbind(df, df_origin[df_origin$SUB_GENE %in% subtrates2show & df_origin$GENE_cancer %in% gene_altered_had2show,])
        
      } else {
        ## for altered genes to show either SMGs or driver genes for these cancer types
        enzymes2show <- unique(df$GENE[df$p < sig_p_thres &  (df$GENE %in% gene_altered2show) & df$SELF == "trans"])
        enzymes2show <- enzymes2show[!(enzymes2show %in% gene_altered_not_show)]
        ## get the list of substrates to show
        subtrates2show <- unique(df$SUB_GENE[df$p < sig_p_thres & df$GENE %in% enzymes2show & (df$SUB_GENE.path != "other" | df$SUB_GENE %in% gene_altered2show)])
        subtrates2show <- c(subtrates2show, gene_affected2show)
        subtrates2show <- subtrates2show[!(subtrates2show %in% gene_affected_not_show)]
        
        ## do the data frame for plotting
        df <- df[(df$GENE %in% enzymes2show) & (df$SUB_GENE %in% subtrates2show),]
        enzymes2show <- unique(df$GENE[df$p < sig_p_thres &  (df$GENE %in% gene_altered2show) & df$SELF == "trans"])
        df <- df[(df$GENE %in% enzymes2show) & (df$SUB_GENE %in% subtrates2show),]
        subtrates2show <- unique(df$SUB_GENE[df$p < sig_p_thres & df$GENE %in% enzymes2show & (df$SUB_GENE.path != "other"  | df$SUB_GENE %in% gene_altered2show)])
        df <- df[(df$GENE %in% enzymes2show) & (df$SUB_GENE %in% subtrates2show),]
      }
      
      
      df$sig <- (df$p < sig_p_thres)
      df$pair_can  <- paste0(df$pair, ":", df$cancer)
      df <- df[order(df$SUB_GENE.path, df$SUB_GENE),]
      df <- df[order(df$p, decreasing = T),]
      df <- df[!duplicated(df$pair_can),]
      df$log10_Pvalue <- -log10(df$p)
      
      cap <- min(max(abs(df$meddiff), na.rm = T),2)
      df$meddiff_capped <- df$meddiff
      df$meddiff_capped[df$meddiff > cap] <- cap
      df$meddiff_capped[df$meddiff < (-cap)] <- (-cap)
      
      df$driver_gene_type.en[df$driver_gene_type.en == ""] <- "other"
      df$driver_gene_type.en <- factor(df$driver_gene_type.en, levels = c("oncogene", "tsg", "other"))
      df$SUB_GENE.path <- factor(df$SUB_GENE.path, levels = pathway_order)
      
      sig_colors <- c("black", "white")
      names(sig_colors) <- c(paste0('P-value<', sig_p_thres), paste0('P-value>', sig_p_thres))
      
      tmp <- as.vector(unique(df$GENE))
      # df$GENE <- factor(df$GENE, levels = c(c("CTNNB1", "APC", "PTEN", "TP53", "MSH6"), tmp[!(tmp %in% c("CTNNB1", "APC", "PTEN", "TP53", "MSH6"))]))
      
      df$shape <- ifelse(df$SELF == "cis", "a", "b")
      df$shape <- as.factor(df$shape)
      
      if (affected_exp_type == "PRO") {
        p <- ggplot()
        p <- p + geom_point(data = df, mapping = aes(x = SUB_GENE, y = GENE, fill = meddiff_capped, size = log10_Pvalue, shape = SELF,
                                                     colour = ifelse(sig == TRUE, paste0('P-value<', sig_p_thres), paste0('P-value>', sig_p_thres))), alpha = 0.6)
        p <- p  + scale_shape_manual(values = c("cis" = 21, "trans" = 22))
        p <- p + scale_color_manual(values = sig_colors)
        p = p + scale_fill_gradientn(name= "protein\nchange", na.value=NA, colours=RdBu1024, limit=c(-cap,cap))
        p <- p + facet_grid(cancer + driver_gene_type.en ~ SUB_GENE.path, drop=T, space = "free",scales = "free")#, space = "free", scales = "free")
        p <- p + xlab("Affected Protein") + ylab("Mutated Gene")
        p <- p + guides(colour = guide_legend(title = "P-value"))
        p <- p + theme_bw()
        p <- p + theme(panel.spacing.x = unit(0, "lines"), panel.spacing.y = unit(0, "lines"))
        p <- p + theme(axis.text.y = element_text(size = 10, face = "bold"))
        p <- p + theme(axis.text.x = element_text(size = 10, face = "bold", angle = 90, hjust = 0.5, vjust = 0.5))
        p <- p + theme(panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "grey80"))
        p <- p + theme(strip.text.x = element_text(angle = 90))
        p
        fn = paste(makeOutDir(resultD = resultD), variant_class, "_", genoalt_type, "_impact_", affected_exp_type, "_num_genoalt_thres_", num_genoalt_thres, "_do.smg_", do.smg, '.pdf',sep = "")
        if (do.smg) {
          ggsave(filename = fn, width = 15, height = 8)
        } else {
          ggsave(filename = fn, width = 20, height = 15)
        }
      }
      if (affected_exp_type == "PHO") {
        p <- ggplot()
        p <- p + geom_point(data = df, mapping = aes(x = SUB_GENE, y = GENE, fill = meddiff_capped, size = log10_Pvalue, shape = SELF,
                                                     colour = ifelse(sig == TRUE, paste0('P-value<', sig_p_thres), paste0('P-value>', sig_p_thres))), alpha = 0.6)
        p <- p  + scale_shape_manual(values = c("cis" = 21, "trans" = 22))
        p <- p + scale_color_manual(values = sig_colors)
        p = p + scale_fill_gradientn(name= "phosphoprotein\nchange", na.value=NA, colours=RdBu1024, limit=c(-cap,cap))
        p <- p + facet_grid(cancer + driver_gene_type.en ~ SUB_GENE.path, drop=T, space = "free",scales = "free")#, space = "free", scales = "free")
        p <- p + xlab("Affected Phosphorotein") + ylab("Mutated Gene")
        p <- p + guides(colour = guide_legend(title = "P-value"))
        p <- p + theme_bw()
        p <- p + theme(panel.spacing.x = unit(0, "lines"), panel.spacing.y = unit(0, "lines"))
        p <- p + theme(axis.text.y = element_text(size = 10, face = "bold"))
        p <- p + theme(axis.text.x = element_text(size = 10, face = "bold", angle = 90, hjust = 0.5, vjust = 0.5))
        p <- p + theme(panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "grey80"))
        p <- p + theme(strip.text.x = element_text(angle = 90))
        p
        fn = paste(makeOutDir(resultD = resultD), variant_class, "_", genoalt_type, "_impact_", affected_exp_type, "_num_genoalt_thres_", num_genoalt_thres, "_do.smg_", do.smg, '.pdf',sep = "")
        if (do.smg) {
          ggsave(filename = fn, width = 15, height = 8)
        } else {
          ggsave(filename = fn, width = 20, height = 15)
        }
      }
      tab4cytoscape <- rbind(tab4cytoscape, df)
    }
    source("./cptac2p_analysis/phospho_network/genoalt/tables/generate_table4cytoscape.R")
    
  }
}


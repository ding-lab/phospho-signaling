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
genes_altered <- "TP53"
affected_exp_type <- "PHO"
affected_exp_type <- "PRO_PHO"
affected_exp_type <- "RNA_PRO"
cancers2process <- c("BRCA", "OV", "CO", "UCEC", "CCRCC", "LIHC")
interaction_type <- "TF"
do.smg <- F

# inputs mutation impact proteome test results and reformat ------------------------------------------------------------------
mut_cnv_cans <- fread(input = paste0(ppnD, "genoalt/tables/test_mut_impact_proteome/mut_impact_proteome_RNA_cptac2p_cptac3_tab.txt"), data.table = F, sep = "\t")
## https://wustl.app.box.com/folder/64828665245
mut_cnv_cans <- mut_cnv_cans[mut_cnv_cans$num >= num_genoalt_thres,]
mut_cnv_cans$pair_pro <- paste0(mut_cnv_cans$GENE, ":",  mut_cnv_cans$SUB_GENE)

# limit to specific gene alteractions -------------------------------------
mut_cnv_cans <- mut_cnv_cans[mut_cnv_cans$GENE %in% genes_altered,]

# limit to interactions just to corum -------------------------------------
# ## input corum
# corum_tab <- fread(input = "./cptac2p/analysis_results/phospho_network/compile_enzyme_substrate/tables/compile_corum_complex/corum_complex_pair_tab.txt", data.table = F)
# ## https://wustl.app.box.com/folder/64826080258
# corum_tab$pair_pro <- paste0(corum_tab$geneA, ":",  corum_tab$geneB)
# 
# ## filter
# mut_cnv_cans <- mut_cnv_cans[mut_cnv_cans$pair_pro %in% corum_tab$pair_pro,]


# limit to interactions just to TF ----------------------------------------
if (interaction_type == "TF") {
  TF_tab <- fread(input = "./Ding_Lab/Projects_Current/TP53_shared_data/resources/PPI/TF_interactions.txt", data.table = F)
  TF_tab %>% head()
  TF_tab <- TF_tab[TF_tab$source_genesymbol == genes_altered,]
  TF_tab %>% head()
  TF_tab %>% nrow()
  # TF_tab <- TF_tab[grepl(x = TF_tab$tfregulons_level, pattern = "A"),]
  ## filter
  TF_tab$pair_pro <- paste0(TF_tab$source_genesymbol, ":", TF_tab$target_genesymbol)
  mut_cnv_cans <- mut_cnv_cans[mut_cnv_cans$pair_pro %in% TF_tab$pair_pro,]
  mut_cnv_cans <- mut_cnv_cans[!(mut_cnv_cans$pair_pro %in% c("TP53:CHEK1", "TP53:CHEK2", "TP53:CDK1")),]
  mut_cnv_cans %>% nrow()
}

# Annotate affected protein to pathways -----------------------------------
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


# appoint an order for pathway --------------------------------------------
pathway_order <- c("TP53", "Cell Cycle", "PI3K", "WNT", "NOTCH", "MMR", "HIPPO", "TGF-Beta", "other")


for (variant_class in c("not_silent", "missense", "truncation")) {
  gene_altered2show <- genes_altered
  gene_altered_not_show <- c("")
  gene_affected2show <- c(unique(mut_cnv_cans$SUB_GENE[mut_cnv_cans$p < sig_p_thres]))
  write.table(x = gene_affected2show, file = "./genes2order.txt", row.names = F, col.names = F, quote = F)
  gene_affected_not_show <- c("")
  
  for (genoalt_type in c("mut")) {
    df <- mut_cnv_cans
    df$expression_type <- df$SUB_MOD_RSD
    df$expression_type[!(df$SUB_MOD_RSD %in% c("RNA", "PRO"))] <- "PHO"
    
    if (affected_exp_type == "PRO") {
      df <- df[df$variant_class == variant_class & df$SUB_MOD_RSD == "PRO" & df$cancer %in% cancers2process & df$genoalt_type == genoalt_type & df$p > 0,]
    }
    if (affected_exp_type == "PHO") {
      df <- df[df$variant_class == variant_class & df$SUB_MOD_RSD != "PRO" & df$cancer %in% cancers2process & df$genoalt_type == genoalt_type & df$p > 0,]
    }
    if (affected_exp_type == "PRO_PHO") {
      df <- df[df$variant_class == variant_class & df$SUB_MOD_RSD != "RNA" & df$cancer %in% cancers2process & df$genoalt_type == genoalt_type & df$p > 0,]
      df$expression_type[!(df$SUB_MOD_RSD %in% c("RNA", "PRO"))] <- "PHO"
    }
    if (affected_exp_type == "RNA_PRO") {
      df <- df[df$variant_class == variant_class & df$expression_type != "PHO" & df$cancer %in% cancers2process & df$genoalt_type == genoalt_type & df$p > 0,]
    }
    
    if (do.smg) {
      df_smg <- NULL
      for (cancer in cancers2process) {
        df_new <- df[df$cancer == cancer & df$GENE %in% SMGs[[cancer]],]
        df_smg <- rbind(df_smg, df_new)
      }
      df <- df_smg
      enzymes2show <- unique(df$GENE[df$p < sig_p_thres &  (df$GENE %in% gene_altered2show) & df$SELF == "trans"])
      subtrates2show <- unique(df$SUB_GENE[df$p < sig_p_thres & (df$SUB_GENE.path != "other" | df$SUB_GENE %in% gene_altered2show)])
      subtrates2show <- c(subtrates2show, gene_affected2show)
      subtrates2show <- subtrates2show[!(subtrates2show %in% gene_affected_not_show)]
      
      ## do the data frame for plotting
      df <- df[(df$GENE %in% enzymes2show) & (df$SUB_GENE %in% subtrates2show),]
      
    } else {
      ## for altered genes to show either SMGs or driver genes for these cancer types
      enzymes2show <- unique(df$GENE[df$p < sig_p_thres &  (df$GENE %in% gene_altered2show) & df$SELF == "trans"])
      enzymes2show <- enzymes2show[!(enzymes2show %in% gene_altered_not_show)]
      ## get the list of substrates to show
      subtrates2show <- unique(df$SUB_GENE[df$p < sig_p_thres & df$GENE %in% enzymes2show])
      subtrates2show <- c(subtrates2show, gene_affected2show)
      subtrates2show <- subtrates2show[!(subtrates2show %in% gene_affected_not_show)]
    }
    ## do the data frame for plotting
    df <- df[(df$GENE %in% enzymes2show) & (df$SUB_GENE %in% subtrates2show),]
    
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
    
    sig_colors <- c("black", "white")
    names(sig_colors) <- c(paste0('P-value<', sig_p_thres), paste0('P-value>', sig_p_thres))
    
    tmp <- as.vector(unique(df$GENE))
    # df$GENE <- factor(df$GENE, levels = c(c("CTNNB1", "APC", "PTEN", "TP53", "MSH6"), tmp[!(tmp %in% c("CTNNB1", "APC", "PTEN", "TP53", "MSH6"))]))
    
    df$shape <- ifelse(df$SELF == "cis", "a", "b")
    df$shape <- as.factor(df$shape)
    
    ## order the cancer types
    df$cancer <- factor(as.vector(df$cancer), levels = rev(cancers2process))
    
    ## order the pathway
    df$SUB_GENE.path <- factor(df$SUB_GENE.path, levels = pathway_order)
    
    
    if (affected_exp_type == "PRO") {
      p <- ggplot()
      p <- p + geom_point(data = df, mapping = aes(x = SUB_GENE, y = cancer, fill = meddiff_capped, size = log10_Pvalue, shape = SELF,
                                                   colour = ifelse(sig == TRUE, paste0('P-value<', sig_p_thres), paste0('P-value>', sig_p_thres))), alpha = 0.6)
      p <- p + facet_grid(. ~ SUB_GENE.path, drop=T, space = "free",scales = "free")#, space = "free", scales = "free")
      
      p <- p  + scale_shape_manual(values = c("cis" = 21, "trans" = 22))
      p <- p + scale_color_manual(values = sig_colors)
      p = p + scale_fill_gradientn(name= "protein\nchange", na.value=NA, colours=RdBu1024, limit=c(-cap,cap))
      p <- p + xlab("Affected Protein") + ylab("Mutated Gene")
      p <- p + guides(colour = guide_legend(title = "P-value"))
      p <- p + theme_bw()
      p <- p + theme(panel.spacing.x = unit(0, "lines"), panel.spacing.y = unit(0, "lines"))
      p <- p + theme(axis.text.y = element_text(size = 10, face = "bold"))
      p <- p + theme(axis.text.x = element_text(size = 10, face = "bold", angle = 90, hjust = 0.5, vjust = 0.5))
      p <- p + theme(panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "grey80"))
      p <- p + theme(strip.text.x = element_text(angle = 90))
      p
      fn = paste(makeOutDir(resultD = resultD), variant_class, "_", genoalt_type, "_impact_", interaction_type, "_", affected_exp_type, "_num_genoalt_thres_", num_genoalt_thres, "_do.smg_", do.smg, '.pdf',sep = "")
      if (do.smg) {
        ggsave(filename = fn, width = 15, height = 8)
      } else {
        ggsave(filename = fn, width = 8, height = 3)
      }
    }
    if (affected_exp_type == "PHO") {
      p <- ggplot()
      p <- p + geom_point(data = df, mapping = aes(x = SUB_GENE, y = cancer, fill = meddiff_capped, size = log10_Pvalue, shape = SELF,
                                                   colour = ifelse(sig == TRUE, paste0('P-value<', sig_p_thres), paste0('P-value>', sig_p_thres))), alpha = 0.6)
      p <- p + facet_grid(. ~ SUB_GENE.path, drop=T, space = "free",scales = "free")#, space = "free", scales = "free")
      p <- p  + scale_shape_manual(values = c("cis" = 21, "trans" = 22))
      p <- p + scale_color_manual(values = sig_colors)
      p = p + scale_fill_gradientn(name= "phosphoprotein\nchange", na.value=NA, colours=RdBu1024, limit=c(-cap,cap))
      p <- p + xlab("Affected Phosphorotein") + ylab("Mutated Gene")
      p <- p + guides(colour = guide_legend(title = "P-value"))
      p <- p + theme_bw()
      p <- p + theme(panel.spacing.x = unit(0, "lines"), panel.spacing.y = unit(0, "lines"))
      p <- p + theme(axis.text.y = element_text(size = 10, face = "bold"))
      p <- p + theme(axis.text.x = element_text(size = 10, face = "bold", angle = 90, hjust = 0.5, vjust = 0.5))
      p <- p + theme(panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "grey80"))
      p <- p + theme(strip.text.x = element_text(angle = 90))
      p
      fn = paste(makeOutDir(resultD = resultD), variant_class, "_", genoalt_type, "_impact_", interaction_type, "_", affected_exp_type, "_num_genoalt_thres_", num_genoalt_thres, "_do.smg_", do.smg, '.pdf',sep = "")
      if (do.smg) {
        ggsave(filename = fn, width = 15, height = 8)
      } else {
        ggsave(filename = fn, width = 8, height = 3)
      }
    }
    if (affected_exp_type == "PRO_PHO") {
      df$expression_type <- factor(as.vector(df$expression_type), levels = c("PRO", "PHO"))
      p <- ggplot()
      p <- p + geom_point(data = df, mapping = aes(x = SUB_GENE, y = cancer, fill = meddiff_capped, size = log10_Pvalue, shape = SELF,
                                                   colour = ifelse(sig == TRUE, paste0('P-value<', sig_p_thres), paste0('P-value>', sig_p_thres))), alpha = 0.6)
      p <- p + facet_grid(expression_type ~ SUB_GENE.path, drop=T, space = "free",scales = "free")#, space = "free", scales = "free")
      p <- p  + scale_shape_manual(values = c("cis" = 21, "trans" = 22))
      p <- p + scale_color_manual(values = sig_colors)
      p = p + scale_fill_gradientn(name= "protein/phosphoprotein\nchange", na.value=NA, colours=RdBu1024, limit=c(-cap,cap))
      p <- p + xlab("Affected Phosphorotein") + ylab("Mutated Gene")
      p <- p + guides(colour = guide_legend(title = "P-value"))
      p <- p + theme_bw()
      p <- p + theme(panel.spacing.x = unit(0, "lines"), panel.spacing.y = unit(0, "lines"))
      p <- p + theme(axis.text.y = element_text(size = 10, face = "bold"))
      p <- p + theme(axis.text.x = element_text(size = 10, face = "bold", angle = 90, hjust = 0.5, vjust = 0.5))
      p <- p + theme(panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "grey80"))
      p <- p + theme(strip.text.x = element_text(angle = 90))
      p
      fn = paste(makeOutDir(resultD = resultD), variant_class, "_", genoalt_type, "_impact_", interaction_type, "_", affected_exp_type, "_num_genoalt_thres_", num_genoalt_thres, "_do.smg_", do.smg, '.pdf',sep = "")
      if (do.smg) {
        ggsave(filename = fn, width = 15, height = 8)
      } else {
        ggsave(filename = fn, width = 10, height = 4)
      }
    }
    
    if (affected_exp_type == "RNA_PRO") {
      df$expression_type <- factor(as.vector(df$expression_type), levels = c("RNA", "PRO"))
      p <- ggplot()
      p <- p + geom_point(data = df, mapping = aes(x = SUB_GENE, y = cancer, fill = meddiff_capped, size = log10_Pvalue, shape = SELF,
                                                   colour = ifelse(sig == TRUE, paste0('P-value<', sig_p_thres), paste0('P-value>', sig_p_thres))), alpha = 0.6)
      p <- p + facet_grid(expression_type ~ SUB_GENE.path, drop=T, space = "free",scales = "free")#, space = "free", scales = "free")
      p <- p  + scale_shape_manual(values = c("cis" = 21, "trans" = 22))
      p <- p + scale_color_manual(values = sig_colors)
      p = p + scale_fill_gradientn(name= "RNA/protein\nchange", na.value=NA, colours=RdBu1024, limit=c(-cap,cap))
      p <- p + xlab("Affected gene") + ylab("Mutated Gene")
      p <- p + guides(colour = guide_legend(title = "P-value"))
      p <- p + theme_bw()
      p <- p + theme(panel.spacing.x = unit(0, "lines"), panel.spacing.y = unit(0, "lines"))
      p <- p + theme(axis.text.y = element_text(size = 10, face = "bold"))
      p <- p + theme(axis.text.x = element_text(size = 10, face = "bold", angle = 90, hjust = 0.5, vjust = 0.5))
      p <- p + theme(panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "grey80"))
      p <- p + theme(strip.text.x = element_text(angle = 90))
      p
      fn = paste(makeOutDir(resultD = resultD), variant_class, "_", genoalt_type, "_impact_", interaction_type, "_", affected_exp_type, "_num_genoalt_thres_", num_genoalt_thres, "_do.smg_", do.smg, '.pdf',sep = "")
      if (do.smg) {
        ggsave(filename = fn, width = 15, height = 8)
      } else {
        ggsave(filename = fn, width = 10, height = 5)
      }
    }
    
  }
}


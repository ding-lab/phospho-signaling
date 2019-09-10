# Yige Wu @ WashU 2019 Aug
## overview figure for mutational impact of kinases on substrates

# source ------------------------------------------------------------------
code_top_dir <- "~/Box/Ding_Lab/Projects_Current/PanCan_Phospho-signaling/phospho-signaling_analysis/"
path2phospho_network_shared <- paste0(code_top_dir, "phospho_network/phospho_network_plotting.R")
source(path2phospho_network_shared)
path2pathway_shared <- paste0(code_top_dir, "dependencies/tables/load_TCGA_pathways.R")
source(path2pathway_shared)

# input version number ----------------------------------------------------
version_num <- 1
datastamp <- "20190813"

# set variables -----------------------------------------------------------
show_p_thres <- 0.2
sig_fdr_thres <- 0.05
cancer <- "GBM"
driver_genes2show <- unique(c(as.vector(driver_genes$Gene[driver_genes$Cancer == "GBM"]), SMGs[["GBM"]]))
num_genoalt_thres <- 5
pair_tab_annotated <- fread(input = "./Ding_Lab/Projects_Current/PanCan_Phospho-signaling/analysis_results/dependencies/tables/compile_protein_pair_table/protein_pair_table_v2.txt", data.table = F)

# for (variant_class in c("not_silent", "missense", "truncation")) {
for (variant_class in c("not_silent")) {
    
  # plot protein section ----------------------------------------------------
  affected_exp_type <- "PRO"
  do.smg <- TRUE
  genes_had2show <- c("")
  if (affected_exp_type == "PRO") {
    for (genoalt_type in c("mut")) {
      mut_cnv_cans <- fread(input = "./Ding_Lab/Projects_Current/PanCan_Phospho-signaling/analysis_results/phospho_network/genoalt/tables/test_SMG_mut_impact_proteome/GBM_not_silent_mut_impact_PRO_tab.txt", data.table = F)
      mut_cnv_cans <- mut_cnv_cans %>%
        filter(num >= num_genoalt_thres)
      mut_cnv_cans$variant_class <- variant_class
      mut_cnv_cans$affected_exp_type <- affected_exp_type
      mut_cnv_cans <- merge(mut_cnv_cans, pair_tab_annotated, by = c( "GENE", "SUB_GENE"), all.x = T)
      mut_cnv_cans$fdr <- FDR_by_id_columns(p_vector = mut_cnv_cans$p, id_columns = c("SELF", "cancer", "variant_class", "affected_exp_type", colnames(pair_tab_annotated)[!(colnames(pair_tab_annotated) %in% c("GENE", "SUB_GENE", "pair_pro"))]), df = mut_cnv_cans)
      mut_cnv_cans <- mut_cnv_cans %>%
        mutate(sig = (fdr < sig_fdr_thres))
      ## annotate the substrates
      sub_genes2pathways <- map2TCGApathwaways(gene_list = unique(mut_cnv_cans$SUB_GENE), pathway_list = tcga_pathways_pluskegg_and_pathway)
      sub_genes2pathways <- data.frame(SUB_GENE = rep(x = names(sub_genes2pathways), sapply(X = sub_genes2pathways, FUN = function(x) length(x))), 
                                       SUB_GENE.path = unlist(sub_genes2pathways, use.names = F))
      sub_genes2pathways <- sub_genes2pathways[!(sub_genes2pathways$SUB_GENE == "APC" & sub_genes2pathways$SUB_GENE.path != "WNT") & !(sub_genes2pathways$SUB_GENE == "GSK3B" & sub_genes2pathways$SUB_GENE.path != "PI3K") & !(sub_genes2pathways$SUB_GENE == "IRS1" & sub_genes2pathways$SUB_GENE.path != "RTK RAS") & !(sub_genes2pathways$SUB_GENE == "TP53" & sub_genes2pathways$SUB_GENE.path != "TP53"),]
      mut_cnv_cans$SUB_GENE.path <- NULL
      mut_cnv_cans <- merge(mut_cnv_cans, sub_genes2pathways, all.x = T)
      
      ## annotate the altered gene to be oncogene or tsgs
      mut_cnv_cans$driver_gene_type.en <- ""
      mut_cnv_cans$driver_gene_type.en[mut_cnv_cans$GENE %in% oncogenes] <- "oncogene"
      mut_cnv_cans$driver_gene_type.en[mut_cnv_cans$GENE %in% tsgs] <- "tsg"
      mut_cnv_cans$sig <- (mut_cnv_cans$fdr < sig_fdr_thres)
      ## get the list of kinase/phosphatases to show
      if (do.smg) {
        enzymes2show <- intersect(SMGs[[cancer]], unique(mut_cnv_cans$GENE[mut_cnv_cans$fdr < sig_fdr_thres]))
      } else {
        stop("meow1")
        # enzymes2show <- unique(mut_cnv_cans$GENE[mut_cnv_cans$p < sig_fdr_thres & mut_cnv_cans$p > 0 & mut_cnv_cans$cancer == cancer & mut_cnv_cans$genoalt_type == genoalt_type & (mut_cnv_cans$GENE %in% driver_genes2show | mut_cnv_cans$driver_gene_type.en != "")])
        # enzymes2show <- enzymes2show[!(enzymes2show %in% c("ACVR2A", "ARHGAP35", "CDK12", "EPHA2", "FLNA", "HUWE1", "INPPL1", "NCOR1", "RASA1", "SED2", "STAG2", "TSC1",
        #                                                    "ZFHX3", "ZFP36L1", "NFE2L2", "PDGFRA", "PPP2R1A"))]
      }
      ## get the list of substrates to show
      subtrates2show <- unique(mut_cnv_cans$SUB_GENE[mut_cnv_cans$sig & mut_cnv_cans$p > 0 & mut_cnv_cans$cancer == "GBM" & mut_cnv_cans$genoalt_type == genoalt_type & mut_cnv_cans$GENE %in% enzymes2show & (!is.na(mut_cnv_cans$SUB_GENE.path) | mut_cnv_cans$SUB_GENE %in% driver_genes2show)])
      # subtrates2show <- c(subtrates2show, smg_cancer, "STAT3", "RB1", "ARID1A", "TP53BP1")
      
      ## do the data frame for plotting
      df <- mut_cnv_cans[mut_cnv_cans$genoalt_type == genoalt_type,]
      df <- df[(df$GENE %in% enzymes2show) & (df$SUB_GENE %in% subtrates2show),]
      
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
      
      df$log10_FDR <- -log10(df$fdr)
      cap <- min(max(abs(df$meddiff), na.rm = T),2)
      df$meddiff_capped <- df$meddiff
      df$meddiff_capped[df$meddiff > cap] <- cap
      df$meddiff_capped[df$meddiff < (-cap)] <- (-cap)
      df$driver_gene_type.en[df$driver_gene_type.en == ""] <- "other"
      df$driver_gene_type.en[df$GENE %in% c("INSR", "PPM1D")] <- "oncogene"
      
      df$driver_gene_type.en <- factor(df$driver_gene_type.en, levels = c("oncogene", "tsg", "other"))
      # df$SUB_GENE.path <- factor(df$SUB_GENE.path, levels = c(names(tcga_pathways_pluskegg_and_pathway), "other"))
      
      sig_colors <- c("black", "white")
      names(sig_colors) <- c(paste0('FDR<', sig_fdr_thres), paste0('FDR>', sig_fdr_thres))
      
      # df <- df[!(df$SUB_GENE %in% c("NCOR1", "ZFYVE16", "ROCK1", "RAF1", "MAPK1")) & ((df$SUB_GENE.path != "other") | (df$SUB_GENE %in% c(smg_cancer, genes_had2show))),]
      # tmp <- as.vector(unique(df$GENE))
      # df$GENE <- factor(df$GENE, levels = c(c("CTNNB1", "APC", "PTEN", "TP53", "MSH6"), tmp[!(tmp %in% c("CTNNB1", "APC", "PTEN", "TP53", "MSH6"))]))
      df$shape <- ifelse(df$SELF == "cis", "a", "b")
      df$shape <- as.factor(df$shape)
      
      p <- ggplot()
      p <- p + geom_point(data = df, mapping = aes(x = SUB_GENE, y = GENE, fill = meddiff_capped, size = log10_FDR, shape = SELF,
                                                   colour = ifelse(sig == TRUE, paste0('FDR<', sig_fdr_thres), paste0('FDR>', sig_fdr_thres))), alpha = 0.6)
      p <- p  + scale_shape_manual(values = c("cis" = 21, "trans" = 22))
      p <- p + scale_color_manual(values = sig_colors)
      p = p + scale_fill_gradientn(name= "protein\nchange", na.value=NA, colours=RdBu1024, limit=c(-cap,cap))
      p <- p + facet_grid(driver_gene_type.en ~ SUB_GENE.path, drop=T, space = "free",scales = "free")#, space = "free", scales = "free")
      p <- p + xlab("Affected Protein") + ylab("Mutated Gene")
      p <- p + guides(colour = guide_legend(title = "FDR"))
      p <- p + theme_bw()
      p <- p + theme(panel.spacing.x = unit(0, "lines"), panel.spacing.y = unit(0, "lines"))
      p <- p + theme(axis.text.y = element_text(size = 10, face = "bold"))
      p <- p + theme(axis.text.x = element_text(size = 10, face = "bold", angle = 90, hjust = 0.5, vjust = 0.5))
      p <- p + theme(panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "grey80"))
      p <- p + theme(strip.text.x = element_text(angle = 90))
      p
      fn = paste(makeOutDir(), variant_class, "_", genoalt_type, "_impact_", affected_exp_type, "_", num_genoalt_thres, "mutations_do.smg_", do.smg, ".", format(Sys.Date(), "%Y%m%d") , ".v", version_num, '.pdf',sep = "")
      if (do.smg) {
        ggsave(filename = fn, width = 10, height = 4)
      } else {
        ggsave(filename = fn, width = 13, height = 8)
      }
    }
  }
  
  
  # plot phospho ------------------------------------------------------------
  affected_exp_type <- "PHO"
  do.smg <- TRUE
  genes_had2show <- c("")
  if (affected_exp_type == "PHO") {
    for (genoalt_type in c("mut")) {
      mut_cnv_cans <- fread(input = "./Ding_Lab/Projects_Current/PanCan_Phospho-signaling/analysis_results/phospho_network/genoalt/tables/test_SMG_mut_impact_proteome/GBM_not_silent_mut_impact_PHO_tab.txt", data.table = F)
      
      mut_cnv_cans <- mut_cnv_cans %>%
        filter(num >= num_genoalt_thres)
      
      mut_cnv_cans <- mut_cnv_cans[mut_cnv_cans$affected_exp_type == "PHO",]      
      ## annotate the substrates
      sub_genes2pathways <- map2TCGApathwaways(gene_list = unique(mut_cnv_cans$SUB_GENE), pathway_list = tcga_pathways_pluskegg_and_pathway)
      sub_genes2pathways <- data.frame(SUB_GENE = rep(x = names(sub_genes2pathways), sapply(X = sub_genes2pathways, FUN = function(x) length(x))), 
                                       SUB_GENE.path = unlist(sub_genes2pathways, use.names = F))
      sub_genes2pathways <- sub_genes2pathways[!(sub_genes2pathways$SUB_GENE == "APC" & sub_genes2pathways$SUB_GENE.path != "WNT") & !(sub_genes2pathways$SUB_GENE == "GSK3B" & sub_genes2pathways$SUB_GENE.path != "PI3K") & !(sub_genes2pathways$SUB_GENE == "IRS1" & sub_genes2pathways$SUB_GENE.path != "RTK RAS") & !(sub_genes2pathways$SUB_GENE == "TP53" & sub_genes2pathways$SUB_GENE.path != "TP53"),]
      mut_cnv_cans$SUB_GENE.path <- NULL
      mut_cnv_cans <- merge(mut_cnv_cans, sub_genes2pathways, all.x = T)
      
      ## annotate the altered gene to be oncogene or tsgs
      mut_cnv_cans$driver_gene_type.en <- ""
      mut_cnv_cans$driver_gene_type.en[mut_cnv_cans$GENE %in% oncogenes] <- "oncogene"
      mut_cnv_cans$driver_gene_type.en[mut_cnv_cans$GENE %in% tsgs] <- "tsg"
      mut_cnv_cans$sig <- (mut_cnv_cans$fdr < sig_fdr_thres)
      
      ## get the list of kinase/phosphatases to show
      if (do.smg) {
        enzymes2show <- intersect(SMGs[[cancer]], mut_cnv_cans$GENE[mut_cnv_cans$sig])
      } else {
        enzymes2show <- unique(mut_cnv_cans$GENE[mut_cnv_cans$sig & mut_cnv_cans$p > 0 & mut_cnv_cans$cancer == cancer & mut_cnv_cans$genoalt_type == genoalt_type & (mut_cnv_cans$GENE %in% driver_genes2show | mut_cnv_cans$driver_gene_type.en != "")])
        enzymes2show <- enzymes2show[!(enzymes2show %in% c("ACVR2A", "ARHGAP35", "CDK12", "EPHA2", "FLNA", "HUWE1", "INPPL1", "NCOR1", "RASA1", "SED2", "STAG2", "TSC1",
                                                           "ZFHX3", "ZFP36L1", "NFE2L2", "PDGFRA", "PPP2R1A"))]
      }
      ## get the list of substrates to show
      subtrates2show <- unique(mut_cnv_cans$SUB_GENE[mut_cnv_cans$sig & mut_cnv_cans$p > 0 & mut_cnv_cans$cancer == "GBM" & mut_cnv_cans$genoalt_type == genoalt_type & mut_cnv_cans$GENE %in% enzymes2show & (!is.na(mut_cnv_cans$SUB_GENE.path) | mut_cnv_cans$SUB_GENE %in% driver_genes2show)])
      subtrates2show <- c(subtrates2show, smg_cancer, genes_had2show)
      
      ## do the data frame for plotting
      df <- mut_cnv_cans[mut_cnv_cans$genoalt_type == genoalt_type,]
      df <- df[(df$GENE %in% enzymes2show) & (df$SUB_GENE %in% subtrates2show),]
      enzymes2show <- intersect(enzymes2show, df$GENE[df$sig])
      df <- df[(df$GENE %in% enzymes2show) & (df$SUB_GENE %in% subtrates2show),]
      df$sig <- (df$fdr < sig_fdr_thres)
      # df <- df[!(df$GENE == "CTNNB1" & df$SUB_GENE %in% c("APC", "AXIN1")),]
      
      ## remove the duplicated pathway annotation
      tmp <- as.vector(df$SUB_GENE.path)
      tmp[is.na(tmp)] <- "other"
      tmp[tmp == "Mismatch repair"] <- "MMR"
      tmp[df$SUB_GENE == "STAT3"] <- "JAK-STAT"
      tmp[df$SUB_GENE == "AXIN1"] <- "WNT"
      
      df$SUB_GENE.path <- tmp
      df <- df[order(df$SUB_GENE.path, df$SUB_GENE),]
      df <- df[order(df$p, decreasing = T),]
      df <- df[order(df$p, decreasing = F),]
      df$pair_can  <- paste0(df$GENE, ":", df$SUB_GENE, ":", df$cancer)
      df <- df[!duplicated(df$pair_can),]
      
      df$log10_FDR <- -log10(df$fdr)
      cap <- min(max(abs(df$meddiff), na.rm = T),2)
      df$meddiff_capped <- df$meddiff
      df$meddiff_capped[df$meddiff > cap] <- cap
      df$meddiff_capped[df$meddiff < (-cap)] <- (-cap)
      df$driver_gene_type.en[df$driver_gene_type.en == ""] <- "other"
      df$driver_gene_type.en[df$GENE %in% c("INSR", "PPM1D")] <- "oncogene"
      
      df$driver_gene_type.en <- factor(df$driver_gene_type.en, levels = c("oncogene", "tsg", "other"))
      # df$SUB_GENE.path <- factor(df$SUB_GENE.path, levels = c(names(tcga_pathways_pluskegg_and_pathway), "other"))
      
      sig_colors <- c("black", "white")
      names(sig_colors) <- c(paste0('FDR<', sig_fdr_thres), paste0('FDR>', sig_fdr_thres))
      
      # df <- df[!(df$SUB_GENE %in% c("NCOR1", "ZFYVE16", "ROCK1", "RAF1", "MAPK1")) & ((df$SUB_GENE.path != "other") | (df$SUB_GENE %in% c(smg_cancer, genes_had2show))),]
      # tmp <- as.vector(unique(df$GENE))
      # df$GENE <- factor(df$GENE, levels = c(c("CTNNB1", "APC", "PTEN", "TP53", "MSH6"), tmp[!(tmp %in% c("CTNNB1", "APC", "PTEN", "TP53", "MSH6"))]))
      
      p <- ggplot()
      p <- p + geom_point(data = df, mapping = aes(x = SUB_GENE, y = GENE, fill = meddiff_capped, size = log10_FDR, shape = SELF,
                                                   colour = ifelse(sig == TRUE, paste0('FDR<', sig_fdr_thres), paste0('FDR>', sig_fdr_thres))), alpha = 0.6)
      p <- p  + scale_shape_manual(values = c("cis" = 21, "trans" = 22))
      p <- p + scale_color_manual(values = sig_colors)
      p = p + scale_fill_gradientn(name= "phosphoprotein\nchange", na.value=NA, colours=RdBu1024, limit=c(-cap,cap))
      p <- p + facet_grid(driver_gene_type.en ~ SUB_GENE.path, drop=T, space = "free",scales = "free")#, space = "free", scales = "free")
      p <- p + xlab("Affected Phosphorotein") + ylab("Mutated Gene")
      p <- p + guides(colour = guide_legend(title = "FDR"))
      p <- p + theme_bw()
      p <- p + theme(panel.spacing.x = unit(0, "lines"), panel.spacing.y = unit(0, "lines"))
      p <- p + theme(axis.text.y = element_text(size = 10, face = "bold"))
      p <- p + theme(axis.text.x = element_text(size = 10, face = "bold", angle = 90, hjust = 0.5, vjust = 0.5))
      p <- p + theme(panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "grey80"))
      p <- p + theme(strip.text.x = element_text(angle = 90))
      p
      fn = paste(makeOutDir(), variant_class, "_", genoalt_type, "_impact_", affected_exp_type, "_", num_genoalt_thres, "mutations_do.smg_", do.smg, ".", format(Sys.Date(), "%Y%m%d") , ".v", version_num, '.pdf',sep = "")
      
      if (do.smg) {
        ggsave(filename = fn, width = 10, height = 4)
      } else {
        ggsave(filename = fn, width = 13, height = 8)
      }
    }
  }
  
}

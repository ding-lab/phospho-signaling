# Yige Wu @ WashU 2018 July
## overview figure for mutational impact of kinases on substrates

# source ------------------------------------------------------------------
setwd("~/Box Sync/")
source("./cptac2p_analysis/phospho_network/phospho_network_shared.R")
source("./cptac2p_analysis/phospho_network/phospho_network_plotting.R")
source("./cptac2p_analysis/dependencies/tables/load_TCGA_pathways.R")


# set variables -----------------------------------------------------------
show_p_thres <- 0.2
sig_p_thres <- 0.05
cancer <- "UCEC"
driver_genes2show <- unique(c(as.vector(driver_genes$Gene[driver_genes$Cancer == "UCEC"]), SMGs[["UCEC"]]))
affected_exp_type <- "PRO"
num_genoalt_thres <- 5

## SMGs to include
smg_cancer <- c("PTEN", 
                "PIK3CA",
                "ARID1A",
                "PIK3R1",
                "CTNNB1",
                "TP53",
                "KRAS",
                "CTCF",
                "FBXW7",
                "RPL22",
                "MUC5B",
                "FGFR2",
                "SI",
                "ANKRD30B",
                "ANKRD30A",
                "ZFHX3",
                "NCAM2",
                "MSH4",
                "ARID5B")
smg_cancer <- c("FLNA",
                "MAP3K4",
                "HUWE1",
                "NSD1",
                "JAK1",
                "RPL22",
                "SCAF4",
                "FBXW7",
                "KMT2D",
                "INPPL1",
                "ZFHX3",
                "KMT2B",
                "TP53",
                "CTCF",
                "CTNNB1",
                "KRAS",
                "PIK3R1",
                "ARID1A",
                "PIK3CA",
                "PTEN")


# for (variant_class in c("not_silent", "missense", "truncation")) {
for (variant_class in c("not_silent")) {
    
  # plot protein section ----------------------------------------------------
  affected_exp_type <- "PRO"
  do.smg <- TRUE
  genes_had2show <- c("STAT3", "RB1", "ARID1A", "CDK1", "CDK2", "CDKN1A", "CCNA2", "TP53BP1")
  if (affected_exp_type == "PRO") {
    for (genoalt_type in c("mut")) {
      # mut_cnv_cans <- fread(input = paste0(ppnD, "genoalt/tables/test_mut_impact_proteome/", cancer,"_mut_impact_", affected_exp_type , "_tab.txt"), data.table = F, sep = "\t")
      # mut_cnv_cans <- fread(input = paste0(ppnD, "genoalt/tables/test_mut_impact_proteome/", cancer, "_", variant_class, "_mut_impact_", affected_exp_type , "_tab.txt"), data.table = F, sep = "\t")
      mut_cnv_cans <- fread(input = "./cptac2p/analysis_results/phospho_network/genoalt/markdowns/mut_impact_proteome_UCEC.txt", data.table = F)
      mut_cnv_cans <- mut_cnv_cans %>%
        filter(num >= num_genoalt_thres)
      
      mut_cnv_cans <- mut_cnv_cans[mut_cnv_cans$affected_exp_type == "PRO",]
      
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
      mut_cnv_cans$sig <- (mut_cnv_cans$fdr < sig_p_thres)
      ## get the list of kinase/phosphatases to show
      if (do.smg) {
        enzymes2show <- intersect(smg_cancer, mut_cnv_cans$GENE[mut_cnv_cans$sig])
      } else {
        enzymes2show <- unique(mut_cnv_cans$GENE[mut_cnv_cans$p < sig_p_thres & mut_cnv_cans$p > 0 & mut_cnv_cans$cancer == cancer & mut_cnv_cans$genoalt_type == genoalt_type & (mut_cnv_cans$GENE %in% driver_genes2show | mut_cnv_cans$driver_gene_type.en != "")])
        enzymes2show <- enzymes2show[!(enzymes2show %in% c("ACVR2A", "ARHGAP35", "CDK12", "EPHA2", "FLNA", "HUWE1", "INPPL1", "NCOR1", "RASA1", "SED2", "STAG2", "TSC1",
                                                           "ZFHX3", "ZFP36L1", "NFE2L2", "PDGFRA", "PPP2R1A"))]
      }
      ## get the list of substrates to show
      subtrates2show <- unique(mut_cnv_cans$SUB_GENE[mut_cnv_cans$sig & mut_cnv_cans$p > 0 & mut_cnv_cans$cancer == "UCEC" & mut_cnv_cans$genoalt_type == genoalt_type & mut_cnv_cans$GENE %in% enzymes2show & (!is.na(mut_cnv_cans$SUB_GENE.path) | mut_cnv_cans$SUB_GENE %in% driver_genes2show)])
      subtrates2show <- c(subtrates2show, smg_cancer, "STAT3", "RB1", "ARID1A", "TP53BP1")
      
      ## do the data frame for plotting
      df <- mut_cnv_cans[mut_cnv_cans$genoalt_type == genoalt_type,]
      df <- df[(df$GENE %in% enzymes2show) & (df$SUB_GENE %in% subtrates2show),]
      # df$sig <- (df$p < sig_p_thres)
      df$sig <- (df$fdr < sig_p_thres)
      
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
      
      df <- df[!(df$SUB_GENE %in% c("NCOR1", "ZFYVE16", "ROCK1", "RAF1", "MAPK1")) & ((df$SUB_GENE.path != "other") | (df$SUB_GENE %in% c(smg_cancer, genes_had2show))),]
      tmp <- as.vector(unique(df$GENE))
      df$GENE <- factor(df$GENE, levels = c(c("CTNNB1", "APC", "PTEN", "TP53", "MSH6"), tmp[!(tmp %in% c("CTNNB1", "APC", "PTEN", "TP53", "MSH6"))]))
      df$shape <- ifelse(df$SELF == "cis", "a", "b")
      df$shape <- as.factor(df$shape)
      
      p <- ggplot()
      p <- p + geom_point(data = df, mapping = aes(x = SUB_GENE, y = GENE, fill = meddiff_capped, size = log10_Pvalue, shape = SELF,
                                                   colour = ifelse(sig == TRUE, paste0('P-value<', sig_p_thres), paste0('P-value>', sig_p_thres))), alpha = 0.6)
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
      fn = paste(makeOutDir(resultD = resultD), variant_class, "_", genoalt_type, "_impact_", affected_exp_type, "_num_genoalt_thres_", num_genoalt_thres, "_do.smg_", do.smg, '.pdf',sep = "")
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
  genes_had2show <- c("STAT3", "RB1", "ARID1A", "CDK1", "CDK2", "CDKN1A", "CCNA2", "TP53BP1")
  if (affected_exp_type == "PHO") {
    for (genoalt_type in c("mut")) {
      mut_cnv_cans <- fread(input = "./cptac2p/analysis_results/phospho_network/genoalt/markdowns/mut_impact_proteome_UCEC.txt", data.table = F)
      
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
      mut_cnv_cans$sig <- (mut_cnv_cans$fdr < sig_p_thres)
      
      ## get the list of kinase/phosphatases to show
      if (do.smg) {
        enzymes2show <- intersect(smg_cancer, mut_cnv_cans$GENE[mut_cnv_cans$sig])
      } else {
        enzymes2show <- unique(mut_cnv_cans$GENE[mut_cnv_cans$sig & mut_cnv_cans$p > 0 & mut_cnv_cans$cancer == cancer & mut_cnv_cans$genoalt_type == genoalt_type & (mut_cnv_cans$GENE %in% driver_genes2show | mut_cnv_cans$driver_gene_type.en != "")])
        enzymes2show <- enzymes2show[!(enzymes2show %in% c("ACVR2A", "ARHGAP35", "CDK12", "EPHA2", "FLNA", "HUWE1", "INPPL1", "NCOR1", "RASA1", "SED2", "STAG2", "TSC1",
                                                           "ZFHX3", "ZFP36L1", "NFE2L2", "PDGFRA", "PPP2R1A"))]
      }
      ## get the list of substrates to show
      subtrates2show <- unique(mut_cnv_cans$SUB_GENE[mut_cnv_cans$sig & mut_cnv_cans$p > 0 & mut_cnv_cans$cancer == "UCEC" & mut_cnv_cans$genoalt_type == genoalt_type & mut_cnv_cans$GENE %in% enzymes2show & (!is.na(mut_cnv_cans$SUB_GENE.path) | mut_cnv_cans$SUB_GENE %in% driver_genes2show)])
      subtrates2show <- c(subtrates2show, smg_cancer, genes_had2show)
      
      ## do the data frame for plotting
      df <- mut_cnv_cans[mut_cnv_cans$genoalt_type == genoalt_type,]
      df <- df[(df$GENE %in% enzymes2show) & (df$SUB_GENE %in% subtrates2show),]
      enzymes2show <- intersect(enzymes2show, df$GENE[df$sig])
      df <- df[(df$GENE %in% enzymes2show) & (df$SUB_GENE %in% subtrates2show),]
      df$sig <- (df$fdr < sig_p_thres)
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
      
      df <- df[!(df$SUB_GENE %in% c("NCOR1", "ZFYVE16", "ROCK1", "RAF1", "MAPK1")) & ((df$SUB_GENE.path != "other") | (df$SUB_GENE %in% c(smg_cancer, genes_had2show))),]
      tmp <- as.vector(unique(df$GENE))
      df$GENE <- factor(df$GENE, levels = c(c("CTNNB1", "APC", "PTEN", "TP53", "MSH6"), tmp[!(tmp %in% c("CTNNB1", "APC", "PTEN", "TP53", "MSH6"))]))
      
      p <- ggplot()
      p <- p + geom_point(data = df, mapping = aes(x = SUB_GENE, y = GENE, fill = meddiff_capped, size = log10_Pvalue, shape = SELF,
                                                   colour = ifelse(sig == TRUE, paste0('P-value<', sig_p_thres), paste0('P-value>', sig_p_thres))), alpha = 0.6)
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
      fn = paste(makeOutDir(resultD = resultD), variant_class, "_", genoalt_type, "_impact_", affected_exp_type, "_num_genoalt_thres_", num_genoalt_thres, "_do.smg_", do.smg, '.pdf',sep = "")
      if (do.smg) {
        ggsave(filename = fn, width = 10, height = 4)
      } else {
        ggsave(filename = fn, width = 13, height = 8)
      }
    }
  }
  
}

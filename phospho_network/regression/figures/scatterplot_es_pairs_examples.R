# Yige Wu @ WashU 2018 Apr
## check regulated pairs with mutational impact from kinases

# source ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_plotting.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')

library(ggrepel)
library(ggplot2)


# set variables -----------------------------------------------------------
cancers2process <- c("BRCA", "OV", "CO", "UCEC", "CCRCC")

pairs2plot <- c("PRKD1:CTNNB1:S191", "PAK4:CTNNB1:S191", "GSK3B:CTNNB1:S552")
# bussiness ------------------------------------------------------------------
for (pair in pairs2plot) {
  geneA <- str_split(string = pair, pattern = ":")[[1]][1]
  geneB <- str_split(string = pair, pattern = ":")[[1]][2]
  phosphosite <- str_split(string = pair, pattern = ":")[[1]][3]
  fn = paste0(makeOutDir(resultD = resultD), geneA, "_", geneB, "_", phosphosite, ".pdf")
  
  # for (cancer in c("BRCA")) {
  if (!file.exists(fn)) {
    sup_tab <- NULL
    sup_tab %>% head()
    for (cancer in cancers2process) {
      
      # input data first because different for each cancer type --------------------------------------------------------------
      ## input protein data
      if (cancer %in% c("BRCA", "OV", "CO")) {
        pho_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "PHO", sample_type = "tumor", pipeline_type = "CDAP", norm_type = "scaled")
        phog_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "collapsed_PHO", sample_type = "tumor", pipeline_type = "CDAP", norm_type = "scaled")
        
      } else if (cancer == "UCEC") {
        pho_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "PHO", sample_type = "tumor", pipeline_type = "PGDAC", norm_type = "median_polishing")
        phog_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "collapsed_PHO", sample_type = "tumor", pipeline_type = "PGDAC", norm_type = "median_polishing")
        
      } else if (cancer == "CCRCC") {
        pho_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "PHO", sample_type = "tumor", pipeline_type = "PGDAC", norm_type = "MD_MAD")
        phog_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "collapsed_PHO", sample_type = "tumor", pipeline_type = "PGDAC", norm_type = "MD_MAD")
        
      } else if (cancer == "LIHC") {
        pho_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "PHO", sample_type = "tumor", pipeline_type = "PGDAC", norm_type = "MD")
        phog_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "collapsed_PHO", sample_type = "tumor", pipeline_type = "PGDAC", norm_type = "MD")
        
      }
      pro_tab <- pro_tab[pro_tab$Gene %in% c(geneB),]
      
      phog_tab <- phog_tab[phog_tab$Gene %in% c(geneA),]
      pho_tab <- pho_tab[pho_tab$Gene == geneB & pho_tab$Phosphosite == phosphosite,]
      
      # make the annotation columns for each sample -----------------------------
      partIDs <- colnames(pho_tab)[!(colnames(pho_tab) %in% c("Gene", "Phosphosite", "Peptide_ID"))]
      col_anno <- data.frame(partID = partIDs)
      
      # make the matrix of values showing in heatmap ----------------------------
      sup_tab_can <- col_anno
      
      if (nrow(pho_tab) > 0) {
        pho_tab <- pho_tab[!(colnames(pho_tab) %in% c("Peptide_ID"))]
        pho_tab.m <- melt(pho_tab, id.vars = c("Gene", "Phosphosite"))
        pho_tab.m %>% head()
        colnames(pho_tab.m) <- c("Gene", "Phosphosite", "partID", "pho_sub")
        sup_tab_can <- merge(sup_tab_can, pho_tab.m[,c("partID", "pho_sub")], by = c("partID"), all.x = T)
      }
      
      if (!is.null(phog_tab)) {
        if (nrow(phog_tab) > 0) {
          phog_tab.m <- melt(phog_tab, id.vars = "Gene")
          phog_tab.m %>% head()
          colnames(phog_tab.m) <- c("Gene", "partID", "pho_kin")
          sup_tab_can <- merge(sup_tab_can, phog_tab.m[,c("partID", "pho_kin")], by = c("partID"), all.x = T)
        }
      }
      
      if (!is.null(pro_tab)) {
        if (nrow(pro_tab) > 0) {
          pro_tab.m <- melt(pro_tab, id.vars = "Gene")
          pro_tab.m %>% head()
          colnames(pro_tab.m) <- c("Gene", "partID", "pro_sub")
          sup_tab_can <- merge(sup_tab_can, pro_tab.m[,c("partID", "pro_sub")], by = c("partID"), all.x = T)
        }
      }

      if (!is.null(sup_tab_can)) {
        sup_tab_can <- unique(sup_tab_can)

        sup_tab_can$cancer <- cancer
        sup_tab <- rbind(sup_tab, sup_tab_can)
      }
      
    }
  }
  
  cancer
  sup_tab %>%
    tail()
  
  tab2p <- sup_tab
  tab2p$cancer <- order_cancer(tab2p$cancer)
  
  p = ggplot(tab2p, aes(x=pho_kin, y=pho_sub))
  p = p + geom_point(stroke = 0, alpha = 0.6, color = 'black', fill = 'black', shape = 16, size = 1.5)
  p = p + geom_vline(xintercept = 0, color = 'grey', alpha = 0.5, linetype = 2)
  p = p + geom_hline(yintercept = 0, color = 'grey', alpha = 0.5, linetype = 2)
  p = p + geom_smooth(mapping = aes(x = pho_kin, y = pho_sub, color= ifelse(cancer %in% c("UCEC", "CO"), "red" ,"grey50")), method = "glm", se=FALSE, 
                       formula = y ~ x, linetype = 2, size = 1)
  p = p + scale_color_manual(values = c("red" = set1[1], "grey50" = "grey50"))
  p <- p + facet_grid(cancer~., scales = "free_y", space = "fixed")
  p = p + labs(y = paste0(geneA, " phosphorylation abundance(log2 ratio)"), 
               x = paste0(geneB, " ", phosphosite, "\nphosphorylation abundance\n(log2 ratio"))
  p = p + theme_nogrid()
  p = p + theme(axis.title.x = element_text(size=6), axis.title.y = element_text(size=12, face = "bold"),
                legend.position = 'none',
                axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), 
                axis.text.y = element_text(colour="black", size=6))
  p = p + theme(title = element_text(size = 8))
  p = p + theme(panel.spacing.y = unit(0, "line"), panel.background = element_rect(color = "white"), 
                strip.background.y = element_rect(colour = "white", fill = "white"))
  p
  ggsave(filename = fn, width = 2.5, height = 6)
}



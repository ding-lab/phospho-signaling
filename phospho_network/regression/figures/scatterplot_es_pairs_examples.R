# Yige Wu @ WashU 2018 Apr
## check regulated pairs with mutational impact from kinases

# source ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_plotting.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')

library(ggrepel)
library(ggplot2)


# set variables -----------------------------------------------------------
cancers2process <- c("BRCA", "OV", "CO", "UCEC", "CCRCC")

# pairs2plot <- c("PRKD1:CTNNB1:S191", "PAK4:CTNNB1:S191", "GSK3B:CTNNB1:S552")
# pairs2plot <- c("AKT2:CTNNB1:S675", "GSK3B:CTNNB1:S552", "PRKCD:CTNNB1:S552", "PRKD1:CTNNB1:S191", "PRKD1:CTNNB1:S552", "PRKACA:CTNNB1:S191", "PAK4:CTNNB1:S191")

# pairs2plot <- c("MET:PTK2:S843", "PTK2:PTK2:S843", "PTEN:PTK2:S843")
pairs2plot <- c("MAPK1:RAF1:S29", "PPP1CC:RAF1:S29")

# input regression --------------------------------------------------------------
regression <- fread(input = paste0(ppnD, "regression/tables/annotate_regression_with_mut_impact/", 
                                   "regression_cptac2p_cptac3_tumor_reg_nonNA20_mut_impact_cancer_specificity_annotated.txt"),
                    data.table = F)
reg_sig <- c(0.05, 0.05); names(reg_sig) <- c("kinase", "phosphatase")
regression %>% nrow()
regression <- regression %>%
  filter(pair_pro %in% omnipath_tab$pair_pro[!(omnipath_tab$Source %in% c("NetKIN", "PhosphoNetworks", "MIMP"))] | pair_pro %in% psp_tab$pair_pro)
regression %>% nrow()
regression <- adjust_regression_by_nonNA(regression = regression, reg_nonNA = 20, reg_sig = reg_sig)
table(regression$regulated)

# bussiness ------------------------------------------------------------------
for (pair in pairs2plot) {
  geneA <- str_split(string = pair, pattern = ":")[[1]][1]
  geneB <- str_split(string = pair, pattern = ":")[[1]][2]
  phosphosite <- str_split(string = pair, pattern = ":")[[1]][3]
  
  fn = paste0(makeOutDir(resultD = resultD), geneA, "_", geneB, "_", phosphosite, ".pdf")
  if (!file.exists(fn)) {
    sup_tab <- NULL
    sup_tab %>% head()
    for (cancer in cancers2process) {
      
      # for (cancer in c("BRCA")) {

      # input data first because different for each cancer type --------------------------------------------------------------
      ## input protein data
      if (cancer %in% c("BRCA", "OV", "CO")) {
        pro_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "PRO", sample_type = "tumor", pipeline_type = "CDAP", norm_type = "scaled")
        pho_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "PHO", sample_type = "tumor", pipeline_type = "CDAP", norm_type = "scaled")
        phog_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "collapsed_PHO", sample_type = "tumor", pipeline_type = "CDAP", norm_type = "scaled")
        
      } else if (cancer == "UCEC") {
        pro_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "PRO", sample_type = "tumor", pipeline_type = "PGDAC", norm_type = "median_polishing")
        pho_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "PHO", sample_type = "tumor", pipeline_type = "PGDAC", norm_type = "median_polishing")
        phog_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "collapsed_PHO", sample_type = "tumor", pipeline_type = "PGDAC", norm_type = "median_polishing")
        
      } else if (cancer == "CCRCC") {
        pro_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "PRO", sample_type = "tumor", pipeline_type = "PGDAC", norm_type = "MD_MAD")
        pho_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "PHO", sample_type = "tumor", pipeline_type = "PGDAC", norm_type = "MD_MAD")
        phog_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "collapsed_PHO", sample_type = "tumor", pipeline_type = "PGDAC", norm_type = "MD_MAD")
      } else if (cancer == "LIHC") {
        pro_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "PRO", sample_type = "tumor", pipeline_type = "PGDAC", norm_type = "MD")
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
      } else {
        sup_tab_can$pho_sub <- NA
      }
      
      if (!is.null(phog_tab)) {
        if (nrow(phog_tab) > 0) {
          phog_tab.m <- melt(phog_tab, id.vars = "Gene")
          phog_tab.m %>% head()
          colnames(phog_tab.m) <- c("Gene", "partID", "pho_kin")
          sup_tab_can <- merge(sup_tab_can, phog_tab.m[,c("partID", "pho_kin")], by = c("partID"), all.x = T)
        } else {
          sup_tab_can$pho_kin <- NA
        }
      } else {
        sup_tab_can$pho_kin <- NA
      }
      
      if (!is.null(pro_tab)) {
        if (nrow(pro_tab) > 0) {
          pro_tab.m <- melt(pro_tab, id.vars = "Gene")
          pro_tab.m %>% head()
          colnames(pro_tab.m) <- c("Gene", "partID", "pro_sub")
          sup_tab_can <- merge(sup_tab_can, pro_tab.m[,c("partID", "pro_sub")], by = c("partID"), all.x = T)
        } else {
          sup_tab_can$pro_sub <- NA
        }
      } else {
        sup_tab_can$pro_sub <- NA
      }
      
      if (!is.null(sup_tab_can)) {
        sup_tab_can <- unique(sup_tab_can)
        sup_tab_can$cancer <- cancer
        sup_tab_can$regulated <- (paste0(pair, "_", cancer) %in% regression$pair_cancer[regression$regulated])
        sup_tab <- rbind(sup_tab, sup_tab_can)
      }
    }
    cancer
    sup_tab %>%
      tail()
    
    tab2p <- sup_tab
    tab2p$cancer <- order_cancer(tab2p$cancer)
    tab2p <- tab2p %>%
      filter(pho_kin < 2.5 & pho_kin > -2.5)
    
    p = ggplot(tab2p, aes(x=pho_kin, y=pho_sub))
    p = p + geom_point(stroke = 0, alpha = 0.6, color = 'black', fill = 'black', shape = 16, size = 1.5)
    p = p + geom_vline(xintercept = 0, color = 'grey', alpha = 0.5, linetype = 2)
    p = p + geom_hline(yintercept = 0, color = 'grey', alpha = 0.5, linetype = 2)
    p = p + geom_smooth(mapping = aes(x = pho_kin, y = pho_sub, color= ifelse(regulated == T, "red" ,"grey50")), method = "glm", se=FALSE, 
                        formula = y ~ x, linetype = 2, size = 1)
    p = p + scale_color_manual(values = c("red" = set1[1], "grey50" = "grey50"))
    p <- p + facet_grid(cancer~., scales = "free", space = "fixed")
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
    
    for (cancer_tmp in unique(sup_tab$cancer)) {
      tab2p <- sup_tab
      tab2p <- tab2p %>%
        filter(cancer == cancer_tmp) %>%
        filter(!is.na(pho_kin) & !is.na(pho_sub)) %>%
        filter(pho_kin < 2.5 & pho_kin > -2.5)
      if (nrow(tab2p) == 0) {
        next()
      }
      p = ggplot(tab2p, aes(x=pho_kin, y=pho_sub))
      p = p + geom_point(stroke = 0, alpha = 0.6, color = 'black', fill = 'black', shape = 16, size = 1.5)
      p = p + geom_vline(xintercept = 0, color = 'grey', alpha = 0.5, linetype = 2)
      p = p + geom_hline(yintercept = 0, color = 'grey', alpha = 0.5, linetype = 2)
      p = p + geom_smooth(mapping = aes(x = pho_kin, y = pho_sub, color= ifelse(regulated == T, "red" ,"grey50")), method = "glm", se=FALSE, 
                          formula = y ~ x, linetype = 2, size = 1)
      p = p + scale_color_manual(values = c("red" = set1[1], "grey50" = "grey50"))
      p <- p + facet_grid(.~cancer, scales = "free", space = "fixed")
      p = p + labs(x = paste0(geneA, " phosphoprotein level\n(log2 ratio)"), 
                   y = paste0(geneB, " ", phosphosite, " phosphorylation\n(log2 ratio)"))
      p = p + theme_nogrid()
      p = p + theme(axis.title.x = element_text(size=10, face = "bold"), 
                    axis.title.y = element_text(size=8, face = "bold"),
                    legend.position = 'none',
                    axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), 
                    axis.text.y = element_text(colour="black", size=6))
      p = p + theme(title = element_text(size = 8))
      p = p + theme(panel.spacing.y = unit(0, "line"), panel.background = element_rect(color = "white"), 
                    strip.background.x = element_rect(colour = "white", fill = color_5cancers[cancer_tmp]),
                    strip.text.x = element_text(color = "white"))
      p
      fn_cancer <- paste0(makeOutDir(resultD = resultD), cancer_tmp, "_", pair, ".pdf")
      ggsave(filename = fn_cancer, width = 2.5, height = 2.5)
    }
  }
}



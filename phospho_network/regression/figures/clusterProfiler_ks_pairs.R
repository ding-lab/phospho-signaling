# Yige Wu @ WashU 2018 Jan
# take the cancer-type specific regulated kinase-substrate pairs and do KEGG over-representation test
## TODO: download the image at the KEGG website

# input -------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
resultDnow <- makeOutDir()
source(paste0(baseD, "cptac2p_analysis/pan3can_aes.R")) # aes for general purposes
library(org.Hs.eg.db)
library(clusterProfiler)

browseKEGG <- function(x, pathID) {
  url <- paste0("http://www.kegg.jp/kegg-bin/show_pathway?", pathID, '/', x@result[pathID, "geneID"])
  browseURL(url)
  invisible(url)
}

# gene set enrichment for shared/unique ks-pairs up to substrate phosphosites --------
enrichKEGG_results <- list()
for (protein in c("kinase")) {
  tn = paste0(resultD,"regression/tables/",protein,"_substrate_regression_cptac2p_3can.txt")
  table_3can <- fread(tn, data.table = F)
  table_3can <- markSigSiteCan(table_3can, sig_thres = sig, protein_type = protein)
  for (self in c('cis', 'trans')) {
    enrichKEGG_results[[self]] <- list()
    table_self <- table_3can[table_3can$SELF == self & !is.na(table_3can$shared3can),]
    
    ## get background genes
    k_s_table <- load_ks_table(protein)
    if (self == 'cis') {
      k_s_table_self <- k_s_table[as.vector(k_s_table$GENE) == k_s_table$SUB_GENE,]
    } else {
      k_s_table_self <- k_s_table[as.vector(k_s_table$GENE) != k_s_table$SUB_GENE,]
    }
    background <- unique(c(as.character(k_s_table_self$GENE), as.character(k_s_table_self$SUB_GENE)))
    background_ids <- bitr(background, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
    
    for (geneset in c('shared3can', 'uniq_BRCA', 'uniq_OV', 'uniq_CO')) {
      ## get genes to be tested
      table_test <- table_self[table_self[, geneset],]
      genes_test <- unique(c(as.character(table_test$KINASE), as.character(table_test$SUBSTRATE)))
      if (length(genes_test) > 10) {
        genes_test_ids <- bitr(genes_test, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
        kk <- enrichKEGG(gene         = genes_test_ids$ENTREZID,
                         organism     = 'hsa',
                         universe = background_ids$ENTREZID,
                         pvalueCutoff = 0.05)
        enrichKEGG_results[[self]][[geneset]] <- kk
        if (length(kk@result$ID) > 0) {
          fn = paste(resultDnow, "cptac2p_",geneset,"_", protein,'_',self, "_SiteWise_KEGG_over-representaton_FDR_", sig, '.pdf',sep ="")
          pdf(fn, height = 4, width = 8)
          print(barplot(kk, colorBy = 'p.adjust', showCategory=12, title = paste0(self, '-regulated ks pairs(', geneset, ')')))
          dev.off()
        }
      }
    }
  }
}


# find the unique and shared enriched pathways for shared/unique k --------
uniq_BRCA_keggs <- enrichKEGG_results[['trans']][['uniq_BRCA']]@result$ID
uniq_OV_keggs <- enrichKEGG_results[['trans']][['uniq_OV']]@result$ID
uniq_CO_keggs <- enrichKEGG_results[['trans']][['uniq_CO']]@result$ID
setdiff(uniq_BRCA_keggs, union(uniq_CO_keggs, uniq_OV_keggs))
intersect(uniq_BRCA_keggs, intersect(uniq_CO_keggs, uniq_OV_keggs))

browseKEGG(enrichKEGG_results[['trans']][['uniq_BRCA']], 'hsa04722')
browseKEGG(enrichKEGG_results[['trans']][['uniq_OV']], 'hsa04722')
browseKEGG(enrichKEGG_results[['trans']][['uniq_CO']], 'hsa04722')

# gene set enrichment for shared/unique ks-pairs up to substrate protein --------
for (protein in c("kinase")) {
  tn = paste0(resultD,"regression/tables/",protein,"_substrate_regression_cptac2p_3can.txt")
  table_3can <- fread(tn, data.table = F)
  table_3can <- markSigCan(table_3can, sig_thres = sig)
  for (self in c('cis', 'trans')) {
    table_self <- table_3can[table_3can$SELF == self,]
    
    ## get background genes
    k_s_table <- load_ks_table(protein)
    if (self == 'cis') {
      k_s_table_self <- k_s_table[as.vector(k_s_table$GENE) == k_s_table$SUB_GENE,]
    } else {
      k_s_table_self <- k_s_table[as.vector(k_s_table$GENE) != k_s_table$SUB_GENE,]
    }
    background <- unique(c(as.character(k_s_table_self$GENE), as.character(k_s_table_self$SUB_GENE)))
    background_ids <- bitr(background, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
    
    for (geneset in c('shared3can', 'uniq_BRCA', 'uniq_OV', 'uniq_CO')) {
      ## get genes to be tested
      table_test <- table_self[table_self[, geneset],]
      genes_test <- unique(c(as.character(table_test$KINASE), as.character(table_test$SUBSTRATE)))
      if (length(genes_test) > 10) {
        genes_test_ids <- bitr(genes_test, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
        kk <- enrichKEGG(gene         = genes_test_ids$ENTREZID,
                         organism     = 'hsa',
                         universe = background_ids$ENTREZID,
                         pvalueCutoff = 0.05)
        if (length(kk@result$ID) > 0) {
          fn = paste(resultDnow, "cptac2p_",geneset,"_", protein,'_',self, "_KSpairs_KEGG_over-representaton_FDR_", sig, '.pdf',sep ="")
          pdf(fn, height = 4, width = 8)
          print(barplot(kk, colorBy = 'p.adjust', showCategory=12, title = paste0(self, '-regulated ks pairs(', geneset, ')')))
          dev.off()
        }
      }
    }
  }
}



# cnetplot(kk, categorySize = "p.adjust")
browseKEGG(kk, 'hsa04510')

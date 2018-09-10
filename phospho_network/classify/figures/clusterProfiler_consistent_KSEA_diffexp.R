# Yige Wu @ WashU 2018 Jan
# take the KSEA and differential consistent directed enzyme-SUB_GENE pairs and do KEGG over-representation test


# input -------------------------------------------------------------------
library(org.Hs.eg.db)
library(clusterProfiler)
library(GSA)

source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
resultDnow <- makeOutDir()
source(paste0(baseD, "cptac2p_analysis/pan3can_aes.R")) # aes for general purposes

browseKEGGresult <- function(x, pathID) {
  url <- paste0("http://www.kegg.jp/kegg-bin/show_pathway?", pathID, '/', x@result[pathID, "geneID"])
  browseURL(url)
  invisible(url)
}

browseKEGG <- function(x, pathID) {
  url <- paste0("http://www.kegg.jp/kegg-bin/show_pathway?", pathID, '/', paste(x@geneSets[[pathID]], collapse = "/"))
  cat(url)
  browseURL(url)
  invisible(url)
}

ptms_site_pairs_sup <- read_csv("~/Box Sync/cptac2p/cptac_shared/analysis_results/phospho_network/compile_enzyme_substrate/compile_omnipath/omnipath_networkin_enzyme_substrate_site_level_union.csv")
ptms_site_pairs_sup <- ptms_site_pairs_sup[ptms_site_pairs_sup$networkin_score >= 5,]

# find enriched pathways for differential phosphorylated sites --------
enrichKEGG_results <- list()
enrichMKEGG_results <- list()
for (sig in c(1)) {
  for (cancer in c("BRCA", "OV", "CO")) {
    enrichKEGG_results[[cancer]] <- list()
    enrichMKEGG_results[[cancer]] <- list()
    
    ptms_table <- fread(input = paste0("./cptac2p/cptac_shared/analysis_results/phospho_network/classify/tables/integrate_KSEA_diffexp/", cancer, "_KSEA_enzyme_diffexp_substrate_consistent.txt"), data.table = F)
    ## get background gene list
    Pho <- fread(paste0(inputD, cancer, "/", prefix[cancer], "_PHO_formatted_normalized_replicate_averaged_imputed.txt"),
                 data.table = F)
    background <- intersect(unique(Pho$Gene), unique(c(as.vector(ptms_site_pairs_sup$GENE), as.vector(ptms_site_pairs_sup$SUB_GENE))))
    background_ids <- bitr(background, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
    
    for (diffexp_type in c("up", "down")) {
      ## get gene list to test
      genes_test <- unique(c(as.vector(ptms_table$SUB_GENE[ptms_table$substrate_direction == diffexp_type]),
                             as.vector(ptms_table$GENE[ptms_table$substrate_direction == diffexp_type])))
      if (length(genes_test) > 0) {
        genes_test_ids <- bitr(genes_test, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
        
        ## KEGG over-representation test
        kk <- enrichKEGG(gene         = genes_test_ids$ENTREZID,
                         organism     = 'hsa',
                         universe = background_ids$ENTREZID,
                         pvalueCutoff = 0.1)
        enrichKEGG_results[[cancer]][[diffexp_type]] <- kk
        if (length(kk@result$ID) > 0) {
          fn = paste(resultDnow, "cptac2p_",diffexp_type,"_", cancer, "_KEGG_over-representaton_FDR_", sig, '.pdf',sep ="")
          pdf(fn, height = 4, width = 8)
          print(barplot(kk, colorBy = 'p.adjust', showCategory=12, title = paste0(diffexp_type, '-regulated phosphosites in ', cancer, ' dataset')))
          dev.off()
        }
        
        ## KEGG Module over-representation test
        mkk <- enrichMKEGG(gene = genes_test_ids$ENTREZID,
                           universe = background_ids$ENTREZID,
                           organism = 'hsa',
                           pvalueCutoff = 0.1)
        enrichMKEGG_results[[cancer]][[diffexp_type]] <- mkk
        if (length(mkk@result$ID) > 0) {
          fn = paste(resultDnow, "cptac2p_",diffexp_type,"_", cancer, "_moduleKEGG_over-representaton_FDR_", sig, '.pdf',sep ="")
          pdf(fn, height = 4, width = 8)
          print(barplot(mkk, colorBy = 'p.adjust', showCategory=12, title = paste0(diffexp_type, '-regulated phosphosites in ', cancer, ' dataset')))
          dev.off()
        }
      }
    }
  }
}

enrichKEGG_results[['BRCA']][['up']]@geneSets$hsa03040

browseKEGG(x = enrichKEGG_results[['BRCA']][['up']], pathID = 'hsa03040')





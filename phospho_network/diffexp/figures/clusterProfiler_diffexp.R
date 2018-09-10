# Yige Wu @ WashU 2018 Jan
# take the differentially expressed protein and phosphosits and do KEGG over-representation test
## TODO: download the image at the KEGG website

# input -------------------------------------------------------------------
library(org.Hs.eg.db)
library(clusterProfiler)
library(GSA)

source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
resultDnow <- makeOutDir()
source(paste0(baseD, "cptac2p_analysis/pan3can_aes.R")) # aes for general purposes


# c4.cm <- read.gmt("./cptac2p/cptac_shared/analysis_results/utils/MSigDB/c4.cm.v6.1.entrez.gmt")
# c4.cgn <- read.gmt("./cptac2p/cptac_shared/analysis_results/utils/MSigDB/c4.cgn.v6.1.entrez.gmt")
# c4.cgn <- GSA.read.gmt("./cptac2p/cptac_shared/analysis_results/utils/MSigDB/c4.cgn.v6.1.entrez.gmt")
# 
# egmt <- enricher(gene = genes_test_ids$ENTREZID, TERM2GENE=c4.cm, universe = background_ids$ENTREZID, TERM2NAME = )
# head(egmt@result)
# barplot(egmt)

browseKEGG <- function(x, pathID) {
  url <- paste0("http://www.kegg.jp/kegg-bin/show_pathway?", pathID, '/', x@result[pathID, "geneID"])
  browseURL(url)
  invisible(url)
}
Pho.diffexp3can <- fread(paste0(resultD, "diffexp/tables/differential_expression/", "Pho.diffexp3can.txt"),
                         data.table = F)
Pho.diffexp3can <- cbind(Pho.diffexp3can, formatPhosphosite(phosphosite_vector = Pho.diffexp3can$Phosphosite, Pho.diffexp3can$Gene))
Pho.diffexp3can <- Pho.diffexp3can[, !(colnames(Pho.diffexp3can) %in% c("Gene", "Phosphosite"))]

# find enriched pathways for differential phosphorylated sites --------
enrichKEGG_results <- list()
enrichMKEGG_results <- list()
for (cancer in c("BRCA", "OV", "CO")) {
  enrichKEGG_results[[cancer]] <- list()
  enrichMKEGG_results[[cancer]] <- list()
  
  Pho.diffexp_can <- Pho.diffexp3can[Pho.diffexp3can$Cancer == cancer,]
  ## get background gene list
  Pho <- fread(paste0(inputD, cancer, "/", prefix[cancer], "_PHO_formatted_normalized_replicate_averaged_imputed.txt"),
               data.table = F)
  background <- unique(Pho$Gene)
  background_ids <- bitr(background, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
  
  for (diffexp_type in c("up", "down")) {
    ## get gene list to test
    genes_test <- unique(Pho.diffexp_can$SUBSTRATE[Pho.diffexp_can$diffexp_type == diffexp_type])
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

# KEGG Gene Set Enrichment Analysis ---------------------------------------

# KEGG Module Gene Set Enrichment Analysis --------------------------------

# gmtfile <- system.file("extdata", package="clusterProfiler")
# c5 <- read.gmt(gmtfile)
# 
# barplot(kk)
# 
# kk <- enrichKEGG_results[['BRCA']][['up']]
# 
# browseKEGG(x = enrichKEGG_results[['BRCA']][['up']], pathID = 'hsa01200')
# browseKEGG(x = enrichKEGG_results[['BRCA']][['up']], pathID = 'hsa00030')
# browseKEGG(x = enrichKEGG_results[['BRCA']][['up']], pathID = 'hsa04610')
# 
# 
# gmtfile <- system.file("extdata", "c5.cc.v5.0.entrez.gmt", package="clusterProfiler")
# c5 <- read.gmt(gmtfile)
# 
# egmt <- enricher(gene, TERM2GENE=c5)
# head(egmt)
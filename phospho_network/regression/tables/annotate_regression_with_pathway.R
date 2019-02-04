# Yige Wu @ WashU 2018 July
## annotate the genes in the regression table to pathways
## allow for duplicate annotation for one gene

# source ------------------------------------------------------------------
baseD = "/Users/yigewu/Box\ Sync/"
wd <- getwd()
if (wd != baseD) {
  setwd(baseD)
}
source("./cptac2p_analysis/phospho_network/phospho_network_shared.R")
source("./cptac2p_analysis/dependencies/tables/load_TCGA_pathways.R")

# input regression table --------------------------------------------------
if (!exists("table2pathway")) {
  stop("input regression table first and named after table2pathway!")
}
table2pathway$row_id <- paste0(table2pathway$pair, ":", table2pathway$Cancer)

# set variables -----------------------------------------------------------
pathway_order <- c("TP53", "Cell Cycle","WNT", "PI3K","RTK RAS", "MAPK", "SWI/SNF", "NOTCH", "MMR", "HIPPO", "TGF-Beta", "other")

# annotate the substrates to oncogenic pathways ---------------------------
sub_genes2pathways <- map2TCGApathwaways(gene_list = unique(table2pathway$SUB_GENE), pathway_list = tcga_pathways_pluskegg_and_pathway)
sub_genes2pathways <- data.frame(SUB_GENE = rep(x = names(sub_genes2pathways), sapply(X = sub_genes2pathways, FUN = function(x) length(x))), 
                                 SUB_GENE.path = unlist(sub_genes2pathways, use.names = F))
tmp <- as.vector(sub_genes2pathways$SUB_GENE.path)
tmp[is.na(tmp)] <- "other"
sub_genes2pathways$SUB_GENE.path <- tmp
sub_genes2pathways <- unique(sub_genes2pathways)

table2pathway$SUB_GENE.path <- NULL
table2pathway <- unique(table2pathway)
table2pathway <- merge(table2pathway, sub_genes2pathways, all.x = T)

tmp <- as.vector(table2pathway$SUB_GENE.path)
tmp[is.na(tmp)] <- "other"
table2pathway$SUB_GENE.path <- tmp
table2pathway <- unique(table2pathway)

# annotate the substrates to oncogenic pathways ---------------------------
GENEs2pathways <- map2TCGApathwaways(gene_list = unique(table2pathway$GENE), pathway_list = tcga_pathways_pluskegg_and_pathway)
GENEs2pathways <- data.frame(GENE = rep(x = names(GENEs2pathways), sapply(X = GENEs2pathways, FUN = function(x) length(x))), 
                             GENE.path = unlist(GENEs2pathways, use.names = F))
tmp <- as.vector(GENEs2pathways$GENE.path)
tmp[is.na(tmp)] <- "other"
GENEs2pathways$GENE.path <- tmp
GENEs2pathways <- unique(GENEs2pathways)

table2pathway$GENE.path <- NULL
table2pathway <- unique(table2pathway)
table2pathway <- merge(table2pathway, GENEs2pathways, all.x = T)

tmp <- as.vector(table2pathway$GENE.path)
tmp[is.na(tmp)] <- "other"
table2pathway$GENE.path <- tmp
table2pathway <- unique(table2pathway)


## correct for duplicate annotation
# tmp[table2pathway$SUB_GENE %in% c("STAT3", "JAK1")] <- "JAK-STAT"
# tmp[table2pathway$SUB_GENE %in% c("ARID1A", "SMARCA4", "PBRM1", "ARID1B")] <- "SWI/SNF"
# tmp[table2pathway$SUB_GENE %in% c("KRAS", "IRS1", "BRAF", "ARAF", "RAF1")] <- "RTK RAS"
# tmp[table2pathway$SUB_GENE %in% c("AXIN1", "APC", "TCF7L2")] <- "WNT"
# tmp[table2pathway$SUB_GENE %in% c("BCL2", "BCL2L1", "TP53")] <- "TP53"
# tmp[table2pathway$SUB_GENE %in% c("MAP3K8", "MAP3K5", "MAP3K4",
#                                  "MAP2K1", "MAP2K2", "MAP2K4", "MAP3K1", "MAPK1", "MAPK3", "MAPK9", "MAPK14", "MAPK13")] <- "MAPK"
# tmp[table2pathway$SUB_GENE %in% c("GSK3B")] <- "PI3K"
# table2pathway$SUB_GENE.path <- tmp
# table2pathway <- unique(table2pathway)


# 
# table2pathway2correct <- table2pathway[table2pathway$row_id %in% table2pathway$row_id[duplicated(table2pathway$row_id)],]
# table2pathway2correct <- table2pathway2correct[order(table2pathway2correct$row_id),]
# 
# genes2resolve <- table(unique(table2pathway2correct[, c("SUB_GENE", "SUB_GENE.path")]))

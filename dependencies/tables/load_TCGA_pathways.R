# Yige Wu @ WashU 2018 Aug
# load TCGA pathway templates
# ref: https://ars.els-cdn.com/content/image/1-s2.0-S0092867418303593-mmc3.xlsx


source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')

pathway_names <- c("Cell Cycle", "HIPPO", "MYC", "NOTCH", "NRF2", "PI3K", "TGF-Beta", "RTK RAS", "TP53", "WNT")

tcga_pathways <- list()
gene_count <- 0
for (pathway_name in pathway_names) {
  if (pathway_name == "Cell Cycle") {
    curated_gene_tab <- read_excel("./Ding_Lab/Projects_Current/TCGA_data/gene_lists/TCGA_pathways_templates.xlsx", 
                                   sheet = pathway_name, skip = 2)
  } else {
    curated_gene_tab <- read_excel("./Ding_Lab/Projects_Current/TCGA_data/gene_lists/TCGA_pathways_templates.xlsx", 
                                   sheet = pathway_name)
  }
  curated_gene_tab <- data.frame(curated_gene_tab)
  tcga_pathways[[pathway_name]] <- curated_gene_tab
  gene_count <- gene_count + nrow(curated_gene_tab)
}
print(gene_count)
## 335 gene

map2TCGApathwaways <- function(gene_list, pathway_list) {
  mapped_path <- list()
  for (gene in gene_list) {
    tmp <- NULL
    for (pathway_name in names(pathway_list)) {
      tmp <- c(tmp, ifelse(gene %in% pathway_list[[pathway_name]]$Gene, pathway_name, NA))
    }
    tmp <- tmp[!is.na(tmp)]
    if (!is.null(tmp)) {
      mapped_path[[gene]] <- tmp
    }
  }
  return(mapped_path)
}


# combine KEGG pathway genes with TCGA pathway genes --------
pathway_names4kegg <- pathway_names
pathway_names4kegg[pathway_names4kegg == "RTK RAS"] <- "RAS"
pathway_names4kegg[pathway_names4kegg == "TP53"] <- "p53"
pathway_names4kegg[pathway_names4kegg == "MYC"] <- "XXXX"

load(paste0(pan3can_shared_dataD, "Gene_family/2015-08-01_Gene_Set.RData"))
all_keggs <- names(KEGG)
keggs <- NULL
key_names <- NULL
for (key in pathway_names4kegg) {
  tmp <- all_keggs[grepl(key, all_keggs, ignore.case = T)]
  if (length(tmp) > 0){
    if (key == "RAS") {
      tmp <- tmp[tmp != "hsa03020\tRNA polymerase"]
    }
    keggs <- c(keggs, tmp)
    key_names <- c(key_names, pathway_names[pathway_names4kegg == key])
  }
}
names(keggs) <- key_names
keggs

tcga_pathways_pluskegg <- tcga_pathways
gene_count <- 0
for (pathway_name in names(keggs)) {
  genes_orig <- tcga_pathways[[pathway_name]]$Gene
  genes_new <- unique(c(as.vector(genes_orig), KEGG[[keggs[pathway_name]]]))
  if ("CTNNB1" %in% genes_new & pathway_name != "WNT") {
    genes_new <- genes_new[genes_new != "CTNNB1"]
  }
  tcga_pathways_pluskegg[[pathway_name]] <- data.frame(Gene = genes_new)
  gene_count <- gene_count + length(genes_new)
}
gene_count

## adding DNA repair pathway, spliceosome, ubiquitination, metabolic pathways
keggs2add <- NULL
keys2add <- c("Mismatch repair", "Spliceosome", "MAPK")
for (key in keys2add) {
  tmp <- all_keggs[grepl(key, all_keggs, ignore.case = T)]
  keggs2add <- c(keggs2add, tmp)
}
names(keggs2add) <- keys2add
tcga_pathways_pluskegg_and_pathway <- tcga_pathways_pluskegg
for (pathway_name in keys2add) {
  tcga_pathways_pluskegg_and_pathway[[pathway_name]] <- data.frame(Gene = KEGG[[keggs2add[pathway_name]]])
}
tcga_pathways_pluskegg_and_pathway[["MAPK"]] <- data.frame(Gene = c("MAP3K5", "MAP3K8", "BRAF", "ARAF", "RAF1", "MAP2K1", "MAP2K2", "MAP2K4", "MAP3K1", "MAPK1", "MAPK3", "MAPK9", "MAPK14", "MAPK13"))
# set color for each pathway ---------------------------------------------
# tcga_path_colors <- c(brewer.pal(n = length(pathway_names), name = "Dark2"), brewer.pal(n = length(pathway_names), name = "Set1"))
# tcga_path_colors <- tcga_path_colors[1:length(pathway_names)]
tcga_path_colors <- c("#fb9a99", "#1f78b4", "#b2df8a", "#33a02c", "#a6cee3", "#e31a1c", "#fdbf6f", "#ff7f00", "#b15928", "#6a3d9a", "black", "grey", "#c51b7d", "#8c510a", "#01665e", "#c51b7d")
names(tcga_path_colors) <- c(pathway_names, "NA", "NAgrey", keys2add, "driver")

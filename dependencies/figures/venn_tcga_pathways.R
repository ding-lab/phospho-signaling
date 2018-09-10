# Yige Wu @ WashU 2018 Aug
# check the overlap of TCGA pathway templates and with 
# ref: https://ars.els-cdn.com/content/image/1-s2.0-S0092867418303593-mmc3.xlsx


source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')

pathway_names <- c("Cell Cycle", "HIPPO", "MYC", "NOTCH", "NRF2", "PI3K", "TGF-Beta", "RTK RAS", "TP53", "WNT")

tcga_pathways <- list()
genes_all <- NULL
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
  genes_all <- unique(c(genes_all, as.vector(curated_gene_tab$Gene)))
}
length(genes_all)

## check how each pathway members are overlapping with 
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
  genes_all <- unique(c(genes_all, as.vector(curated_gene_tab$Gene)))
}

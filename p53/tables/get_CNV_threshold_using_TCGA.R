# Yige Wu @ WashU Oct 2018
# compare BICSEQ2 with TCGA result

# source ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/Ding_Lab/Projects_Current/CPTAC/CPTACIII/cptac3_cnv/cptac3_cnv_analysis/cptac3_cnv_shared.R')
# source('/Users/yigewu/Box Sync/Ding_Lab/Projects_Current/CPTAC/CPTACIII/cptac3_cnv/cptac3_cnv_analysis/qc/figures/plot_cnv_comparision.R')
library(biomaRt)
library(dplyr)
library(tidyverse)


# inputs ------------------------------------------------------------------
## input ensembl
# ensembl = useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl",host="www.ensembl.org")
# ensembl_attr <- listAttributes(ensembl)
# ensembl_filter <- listFilters(ensembl)

## get function for calculating TCGA
tcga_pct <- function(tcga_tbl, gene_name, cnv_type){
  gene_values <- tcga_tbl %>% filter(Hugo_Symbol == gene_name) %>% dplyr::select(-Hugo_Symbol) %>% dplyr::select(-Entrez_Gene_Id)
  if (cnv_type == "Amplification"){
    return_pct <- mean(gene_values >= 1, na.rm = T)*100
  } else if (cnv_type == "Deletion") {
    return_pct <- mean(gene_values <= -1, na.rm = T)*100
  } else {
    return_pct <- mean(gene_values == 0, na.rm = T)*100
  }
  return(return_pct)
}

tcga_pct_deep <- function(tcga_tbl, gene_name, cnv_type){
  gene_values <- tcga_tbl %>% filter(Hugo_Symbol == gene_name) %>% dplyr::select(-Hugo_Symbol) %>% dplyr::select(-Entrez_Gene_Id)
  if (cnv_type == "Amplification"){
    return_pct <- mean(gene_values >= 2, na.rm = T)*100
  } else if (cnv_type == "Deletion") {
    return_pct <- mean(gene_values <= -2, na.rm = T)*100
  } else {
    return_pct <- mean(gene_values == 0, na.rm = T)*100
  }
  return(return_pct)
}

## get function for calculating TCGA
cptac_pct <- function(cptac_tbl, gene_name, cnv_type,
                      amp_threshold, del_threshold){
  gene_values <- cptac_tbl %>% filter(gene == gene_name) %>% dplyr::select(-gene)
  if (cnv_type == "Amplification") {
    return_pct <- mean(gene_values > amp_threshold, na.rm = T)*100
  } else if (cnv_type == "Deletion") {
    return_pct <- mean(gene_values < del_threshold, na.rm = T)*100
  } else {
    return_pct <- mean(gene_values >= del_threshold & gene_values <= amp_threshold, na.rm = T)*100
  }
  return(return_pct)
}




# get ready for transforming TCGA OV CNA table entrezgene ids -------------
# cancer = "OV"
# tcga_tbl <- fread(file = paste0(cbioportalD, cancer, "/", tolower(cancer), "_tcga_pub/data_CNA.txt"), data.table = F)
# mapTab = getBM(attributes = c("hgnc_symbol", "entrezgene"), filters = "entrezgene", 
#                values = as.vector(tcga_tbl$Entrez_Gene_Id), mart = ensembl, uniqueRows=FALSE)
# mapTab <- unique(mapTab)
# mapTab <- mapTab[!duplicated(mapTab$entrezgene),]
# rownames(mapTab) <- as.vector(mapTab$entrezgene)

# set variables -----------------------------------------------------------


del_threshold <- log2(0.9)
amp_threshold <- log2(1.1)

# business ----------------------------------------------------------------
del_pct_df_sup <- NULL
# for (cancer in c("UCEC", "CCRCC", "LIHC", "CO", "BRCA", "OV")) {
for (cancer in c("BRCA")) {
  ## input CPTAC CNVs
  if (cancer %in% c("UCEC", "CCRCC", "LIHC")) {
    cptac_tbl <- fread(input = paste0("../../../TP53_shared_data/resources/Somatic_Copy_Number/somatic_CNA." , cancer, ".partID.txt"), data.table = F)
  }
  if (cancer %in% c("CO", "BRCA", "OV")) {
    cptac_tbl <- fread(input = paste0("../../CPTAC_Prospective_Samples/copy_number/gatk/v1.3.swap_contamination_fixed/prospective_somatic/gene_level/v1.3.CPTAC2_prospective.2018-03-19/", toupper(substr(cancer, start = 1, stop = 2)), "/gene_level_CNV.", substr(cancer, start = 1, stop = 2), ".v1.3.2018-03-19.tsv"), data.table = F)
    colnames(cptac_tbl)[1] <- "gene"
  }
  
  ## input cBioPortal gene-level CNV
  if (cancer %in% c("BRCA", "OV", "UCEC")) {
    tcga_tbl <- fread(file = paste0(cbioportalD, cancer, "/", tolower(cancer), "_tcga_pub/data_CNA.txt"), data.table = F)
    # tcga_tbl <- read_delim( paste0(cbioportalD, cancer, "/", tolower(cancer), "_tcga_pub/data_CNA.txt"), 
    #                        "\t", escape_double = FALSE, trim_ws = TRUE)
  }
  if (cancer %in% c("CCRCC")) {
    tcga_tbl <- fread(file = paste0(cbioportalD, cancer, "/", "kirc", "_tcga_pub/data_CNA.txt"), data.table = F)
  }
  if (cancer %in% c("LIHC")) {
    tcga_tbl <- fread(file = paste0(cbioportalD, cancer, "/", tolower(cancer), "_tcga/data_CNA.txt"), data.table = F)
  }
  if (cancer == "CO") {
    tcga_tbl <- fread(file = paste0(cbioportalD, cancer, "/", "coadread", "_tcga_pub/data_CNA.txt"), data.table = F)
  }
  if (cancer == "OV") {
    tcga_tbl <- fread(file = paste0(cbioportalD, cancer, "/", "ov_tcga_pan_can_atlas_2018", "/data_CNA.txt"), data.table = F)
    tcga_tbl <- tcga_tbl[!is.na(tcga_tbl$Hugo_Symbol),]
  }
  
  breaks <- seq(-1.5, 0, 0.05)
  del_pct <- NULL
  for (del_threshold in breaks) {
    del_pct <-   c(del_pct, cptac_pct(cptac_tbl = cptac_tbl, gene_name = "TP53", cnv_type = "Deletion", del_threshold = del_threshold))
  }
  del_pct_df <- data.frame(del_threshold = breaks, del_pct = del_pct, 
                           del_pct_diff_tcga = (del_pct - tcga_pct(tcga_tbl = tcga_tbl, "TP53", "Deletion")),
                           del_pct_diff_tcga_deep = (del_pct - tcga_pct_deep(tcga_tbl = tcga_tbl, "TP53", "Deletion")))
  del_pct_df$cancer <- cancer
  del_pct_df_sup <- rbind(del_pct_df, del_pct_df_sup)
  
}

tcga_pct(tcga_tbl = tcga_tbl, gene_name = "MET", cnv_type = "Amplification")

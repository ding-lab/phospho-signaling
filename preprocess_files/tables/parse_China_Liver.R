# Yige Wu @ WashU 2018 Aug
## For formating CCRCC data in data freeze to be in concord with CDAP proteomics data and genomics data


# source ------------------------------------------------------------------
wd <- getwd()
if (wd != "/Users/yigewu/Box Sync") {
  setwd("/Users/yigewu/Box Sync")
}
source("./cptac2p_analysis/phospho_network/phospho_network_shared.R")
library(readxl)
# set variables -----------------------------------------------------------
cancer <- "LIHC"
pipeline_type <- "PGDAC"
gene_col_names <- NULL
gene_col_names["PRO"] <- "Gene names"
gene_col_names["PHO"] <- "Gene names"
gene_col_names["collapsed_PHO"] <- "protein"

# input clinical data -----------------------------------------------------
meta_tab <- read_excel("./Ding_Lab/Projects_Current/TP53_shared_data/resources/TP53_samples_classification/HCC Clinical information and follow-up information.xlsx", 
                        sheet = "Table S1 Clinical and follow-up")
meta_tab <- data.frame(meta_tab)
meta_tab$partID <- paste0("LIHC", meta_tab$Case.ID)
rownames(meta_tab) <- as.vector(meta_tab$partID)
meta_tab$sampID.tumor <- str_split_fixed(string = meta_tab$Tumor..T..sample.ID, pattern = "T", 2)[,2]
meta_tab$sampID.normal <- str_split_fixed(string = meta_tab$Adjacent.Non.tumor.liver.tissue..N..sample.ID, pattern = "N", 2)[,2]

exclude_sampIDs <- c("")
# rewrite proteomics/phosphoproteomics data split by tumor and normal ----------------------------------
for (sample_type in c("tumor")) {
  ## write protein data
  expression_type <- "PRO"
  pro_data <- read_xlsx(path = "Ding_Lab/Projects_Current/ICPC_liver_proteogenomic/not_filtered_protein_phospho/Proteome_HCC_proteinGroups_33sets.xlsx", sheet = "10783 protein groups")
  pro_data <- pro_data[!is.na(pro_data$`Gene names`),]
  sampIDs <- intersect(colnames(pro_data), c(as.vector(meta_tab$sampID.tumor)))
  sampIDs <- sampIDs[!(sampIDs %in% paste0("T", exclude_sampIDs))]
  partIDs <- sapply(X = sampIDs, FUN = function(sampID, meta_tab) {
    partID <- unlist(meta_tab$partID[meta_tab$sampID.normal == sampID | meta_tab$sampID.tumor == sampID])[1]
    return(partID)
  }, meta_tab = meta_tab)
  
  partIDs <- partIDs[!is.na(partIDs)]
  sampIDs <- names(partIDs)

  ## log2 transform and median cencerting
  exp_mat <- pro_data[, sampIDs]
  exp_mat <- log2(exp_mat)
  head(exp_mat)
  per_sample_median <- sapply(X = colnames(exp_mat), FUN = function(i, mat) median(unlist(mat[,i]), na.rm = T), mat = exp_mat)
  exp_mat_MD <- exp_mat - matrix(rep(per_sample_median, nrow(exp_mat)), byrow = T, nrow = nrow(exp_mat))
  colnames(exp_mat_MD) <- partIDs
  exp_mat.partID <- data.frame(Gene = unlist(pro_data[, gene_col_names[expression_type]]))
  exp_mat.partID <- cbind(exp_mat.partID, exp_mat_MD)
  head(exp_mat.partID)
  write.table(x = exp_mat.partID, file = paste0(makeOutDir(resultD = resultD), cancer, "_", expression_type, "_", sample_type, "_", pipeline_type, "_", "MD", "_", "partID", ".txt"), quote = F, row.names = F, sep = "\t")
  rm(pro_data)
  rm(exp_mat.partID)
  
  ## write collapsed_PHOtein data
  expression_type <- "collapsed_PHO"
  collapsed_PHO_data <- fread("Ding_Lab/Projects_Current/TP53_shared_data/resources/Phospho/gene_level/LIHC_tumor.log2_NA50_collapse_medianHCC_phosproteomics_165.tsv", data.table = F)
  if (any(is.na(collapsed_PHO_data$protein))) {
    collapsed_PHO_data <- collapsed_PHO_data[!is.na(collapsed_PHO_data$protein),]
  }
  sampleIDs_old <- colnames(collapsed_PHO_data)
  partIDs <- sapply(X = sampleIDs_old, FUN = function(sampID_old, meta_tab) {
    sampID <- str_split(string = sampID_old, pattern = "T")[[1]][2]
    partID <- unique(unlist(meta_tab$partID[meta_tab$sampID.normal == sampID | meta_tab$sampID.tumor == sampID])[1])
    return(partID)
  }, meta_tab = meta_tab)
  
  partIDs <- partIDs[!is.na(partIDs)]
  sampleIDs_old <- names(partIDs)
  
  ## log2 transform and median cencerting
  exp_mat.partID <- collapsed_PHO_data[, c(gene_col_names[expression_type], sampleIDs_old)]
  colnames(exp_mat.partID) <- c("Gene", partIDs)
  head(exp_mat.partID)
  write.table(x = exp_mat.partID, file = paste0(makeOutDir(resultD = resultD), cancer, "_", expression_type, "_", sample_type, "_", pipeline_type, "_", "MD", "_", "partID", ".txt"), quote = F, row.names = F, sep = "\t")
  rm(exp_mat.partID)
  
  
  # ## write phosphorylatin data
  # # pho_data <- read_xlsx(path = "Ding_Lab/Projects_Current/ICPC_liver_proteogenomic/not_filtered_protein_phospho/Phosphoproteome_HCC_Phospho (STY)Sites_33sets.xlsx", sheet = "Phospho (STY)Sites_")
  # 
  # pho_data <- fread(input = "Ding_Lab/Projects_Current/TP53_shared_data/resources/Phospho/site_level/LIHC_tumor.log2_noimputation_NA50_median_HCC_phosphosite_159.tsv", data.table = F)
  # head(pho_data)
  # sampIDs <- colnames(pho_data)[!(colnames(pho_data) %in% c("protein"))]
  # sampIDs <- sampIDs[!(sampIDs %in% paste0("T", exclude_sampIDs))]
  # partIDs <- sapply(X = sampIDs, FUN = function(sampID, meta_tab) {
  #   tumorID <- str_split_fixed(string = sampID, pattern = "T", n = 2)[2]
  #   partID <- unlist(meta_tab$partID[meta_tab$Tumor.ID.in.Proteome.experiment == tumorID])[1]
  #   return(partID)
  # }, meta_tab = meta_tab)
  # tmp <- pho_data$protein
  # tmp <- str_split_fixed(string = tmp, pattern = "\\:", 2)
  # pho_head <- data.frame(Gene = tmp[,1], Phosphosite = tmp[,2])
  # head(pho_head)
  # partIDs <- partIDs[!is.na(partIDs)]
  # sampIDs <- names(partIDs)
  # pho2w <- pho_data[, sampIDs]
  # pho2w.partID <- cbind(pho_head, pho2w)
  # colnames(pho2w.partID) <- c("Gene", "Phosphosite", partIDs)
  # expression_type <- "PHO"
  # write.table(x = pho2w.partID, file = paste0(makeOutDir(resultD = resultD), cancer, "_", expression_type, "_", sample_type, "_", pipeline_type, "_", "MD", "_", "partID", ".txt"), quote = F, row.names = F, sep = "\t")
  # 
  # ## write RNA data
  # rna_tab <- read_delim("Ding_Lab/Projects_Current/TP53_shared_data/resources/rna/LIHC_UQ_RSEM_hugo_tumor.tsv", 
  #                       "\t", escape_double = FALSE, trim_ws = TRUE)
  # sampIDs <- colnames(rna_tab)[!(colnames(rna_tab) %in% c("X1"))]
  # sampIDs <- sampIDs[!(sampIDs %in% exclude_sampIDs)]
  # 
  # partIDs <- sapply(X = sampIDs, FUN = function(sampID, meta_tab) {
  #   tumorID <- sampID
  #   partID <- unlist(meta_tab$partID[meta_tab$Tumor.ID.in.Proteome.experiment == tumorID])[1]
  #   return(partID)
  # }, meta_tab = meta_tab)
  # partIDs <- partIDs[!is.na(partIDs)]
  # sampIDs <- names(partIDs)
  # rna2w <- rna_tab[, c("X1", sampIDs)]
  # head(rna2w)
  # colnames(rna2w) <- c("gene", partIDs)
  # head(rna2w)
  # rna2w <- data.frame(rna2w)
  # head(rna2w)
  # expression_type <- "RNA"
  # write.table(x = rna2w, file = paste0(makeOutDir(resultD = resultD), cancer, "_", expression_type, "_", sample_type, "_", pipeline_type, "_", "RSEM", "_", "partID", ".txt"), quote = F, row.names = F, sep = "\t")
  # rm(rna_tab)
  # 
  ## reformat the maf file
  maf <- fread(input = "./Ding_Lab/Projects_Current/TP53_shared_data/resources/Somatic_Mutation/MAF/LIHC_AllSamples.filtered.mutect2.maf", data.table = F)
  sampIDs <- maf$Tumor_Sample_Barcode
  partIDs <- sapply(X = sampIDs, FUN = function(sampID, meta_tab) {
    partID <- unlist(meta_tab$partID[meta_tab$sampID.normal == sampID | meta_tab$sampID.tumor == sampID])[1]
    return(partID)
  }, meta_tab = meta_tab)
  nrow(maf)
  maf <- maf[!is.na(partIDs),]
  nrow(maf)
  maf$partID <- partIDs[!is.na(partIDs)]
  maf <- maf[!(maf$Tumor_Sample_Barcode %in% exclude_sampIDs),]
  maf$Tumor_Sample_Barcode <- paste0(maf$partID, "_T")
  write.table(x = maf, file = paste0(makeOutDir(resultD = resultD),  "LIHC_AllSamples.filtered.mutect2.partID.maf"), quote = F, row.names = F, sep = "\t")

  ## reformat the CNV file
  cna <- read_delim("Ding_Lab/Projects_Current/TP53_shared_data/resources/Somatic_Copy_Number/gene_level_CNV.just4anno_LIHC_cnv.v1.2.tsv",
                    "\t", escape_double = FALSE, trim_ws = TRUE)
  head(cna)
  sampIDs <- colnames(cna)[!(colnames(cna) %in% c("gene"))]
  sampIDs <- sampIDs[!(sampIDs %in% exclude_sampIDs)]
  partIDs <- sapply(X = sampIDs, FUN = function(sampID, meta_tab) {
    partID <- unlist(meta_tab$partID[meta_tab$sampID.normal == sampID | meta_tab$sampID.tumor == sampID])[1]
    return(partID)
  }, meta_tab = meta_tab)
  partIDs <- partIDs[!is.na(partIDs)]
  sampIDs <- names(partIDs)
  cna2w <- cna[, c("gene", sampIDs)]
  colnames(cna2w) <- c("gene", partIDs)
  cna2w <- data.frame(cna2w)
  head(cna2w)
  write.table(x = cna2w, file = paste0(makeOutDir(resultD = resultD), "somatic_CNA.", cancer, ".partID.txt"), quote = F, row.names = F, sep = "\t")

}

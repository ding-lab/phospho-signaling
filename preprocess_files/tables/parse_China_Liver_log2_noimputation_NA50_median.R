# Yige Wu @ WashU 2018 Aug
## For formating CCRCC data in data freeze to be in concord with CDAP proteomics data and genomics data


# source ------------------------------------------------------------------
wd <- getwd()
if (wd != "/Users/yigewu/Box Sync") {
  setwd("/Users/yigewu/Box Sync")
}
source("./cptac2p_analysis/phospho_network/phospho_network_shared.R")

# set variables -----------------------------------------------------------
cancer <- "LIHC"
pipeline_type <- "PGDAC"


# input clinical data -----------------------------------------------------
library(readxl)
meta_tab <- read_excel("./Ding_Lab/Projects_Current/TP53_shared_data/resources/TP53_samples_classification/HCC_Clinical_information_follow-up_information.xlsx", 
                        sheet = "Table 1. Clinical Information")
meta_tab <- data.frame(meta_tab)
meta_tab$partID <- paste0("LIHC", meta_tab$Case.ID)
rownames(meta_tab) <- as.vector(meta_tab$partID)

exclude_sampIDs <- c("")
# rewrite proteomics/phosphoproteomics data split by tumor and normal ----------------------------------
for (sample_type in c("tumor")) {
  ## write protein data
  pro_data <- fread(input = "Ding_Lab/Projects_Current/TP53_shared_data/resources/Protein/LIHC_tumor.log2_filter_N50_median_corrected_HCC_proteomics_165.tsv", data.table = F)
  sampIDs <- colnames(pro_data)[!(colnames(pro_data) %in% c("protein"))]
  sampIDs <- sampIDs[!(sampIDs %in% paste0("T", exclude_sampIDs))]
  partIDs <- sapply(X = sampIDs, FUN = function(sampID, meta_tab) {
    tumorID <- str_split_fixed(string = sampID, pattern = "T", n = 2)[2]
    partID <- unlist(meta_tab$partID[meta_tab$Tumor.ID.in.Proteome.experiment == tumorID])[1]
    return(partID)
  }, meta_tab = meta_tab)
  partIDs <- partIDs[!is.na(partIDs)]
  sampIDs <- names(partIDs)
  pro2w <- pro_data[, c("protein", sampIDs)]
  pro2w.partID <- pro2w
  colnames(pro2w.partID) <- c("Gene", partIDs)
  head(pro2w.partID)
  expresson_type <- "PRO"
  write.table(x = pro2w.partID, file = paste0(makeOutDir(resultD = resultD), cancer, "_", expresson_type, "_", sample_type, "_", pipeline_type, "_", "MD", "_", "partID", ".txt"), quote = F, row.names = F, sep = "\t")
  
  ## write phosphorylatin data
  pho_data <- fread(input = "Ding_Lab/Projects_Current/TP53_shared_data/resources/Phospho/site_level/LIHC_tumor.log2_noimputation_NA50_median_HCC_phosphosite_159.tsv", data.table = F)
  head(pho_data)
  sampIDs <- colnames(pho_data)[!(colnames(pho_data) %in% c("protein"))]
  sampIDs <- sampIDs[!(sampIDs %in% paste0("T", exclude_sampIDs))]
  partIDs <- sapply(X = sampIDs, FUN = function(sampID, meta_tab) {
    tumorID <- str_split_fixed(string = sampID, pattern = "T", n = 2)[2]
    partID <- unlist(meta_tab$partID[meta_tab$Tumor.ID.in.Proteome.experiment == tumorID])[1]
    return(partID)
  }, meta_tab = meta_tab)
  tmp <- pho_data$protein
  tmp <- str_split_fixed(string = tmp, pattern = "\\:", 2)
  pho_head <- data.frame(Gene = tmp[,1], Phosphosite = tmp[,2])
  head(pho_head)
  partIDs <- partIDs[!is.na(partIDs)]
  sampIDs <- names(partIDs)
  pho2w <- pho_data[, sampIDs]
  pho2w.partID <- cbind(pho_head, pho2w)
  colnames(pho2w.partID) <- c("Gene", "Phosphosite", partIDs)
  expresson_type <- "PHO"
  write.table(x = pho2w.partID, file = paste0(makeOutDir(resultD = resultD), cancer, "_", expresson_type, "_", sample_type, "_", pipeline_type, "_", "MD", "_", "partID", ".txt"), quote = F, row.names = F, sep = "\t")
  
  ## write RNA data
  rna_tab <- read_delim("Ding_Lab/Projects_Current/TP53_shared_data/resources/rna/LIHC_UQ_RSEM_hugo_tumor.tsv", 
                        "\t", escape_double = FALSE, trim_ws = TRUE)
  sampIDs <- colnames(rna_tab)[!(colnames(rna_tab) %in% c("X1"))]
  sampIDs <- sampIDs[!(sampIDs %in% exclude_sampIDs)]
  
  partIDs <- sapply(X = sampIDs, FUN = function(sampID, meta_tab) {
    tumorID <- sampID
    partID <- unlist(meta_tab$partID[meta_tab$Tumor.ID.in.Proteome.experiment == tumorID])[1]
    return(partID)
  }, meta_tab = meta_tab)
  partIDs <- partIDs[!is.na(partIDs)]
  sampIDs <- names(partIDs)
  rna2w <- rna_tab[, c("X1", sampIDs)]
  head(rna2w)
  colnames(rna2w) <- c("gene", partIDs)
  head(rna2w)
  rna2w <- data.frame(rna2w)
  head(rna2w)
  expresson_type <- "RNA"
  write.table(x = rna2w, file = paste0(makeOutDir(resultD = resultD), cancer, "_", expresson_type, "_", sample_type, "_", pipeline_type, "_", "RSEM", "_", "partID", ".txt"), quote = F, row.names = F, sep = "\t")
  rm(rna_tab)
  
  ## reformat the maf file
  maf <- fread(input = "./Ding_Lab/Projects_Current/TP53_shared_data/resources/Somatic_Mutation/MAF/LIHC_AllSamples.filtered.mutect2.maf", data.table = F) 
  sampIDs <- maf$Tumor_Sample_Barcode
  partIDs <- sapply(X = sampIDs, FUN = function(sampID, meta_tab) {
    paratumorID <- sampID
    partID <- unlist(meta_tab$partID[meta_tab$Paratumor.ID.in.Proteome.experiment == paratumorID])[1]
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
    tumorID <- sampID
    partID <- unlist(meta_tab$partID[meta_tab$Tumor.ID.in.Proteome.experiment == tumorID])[1]
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

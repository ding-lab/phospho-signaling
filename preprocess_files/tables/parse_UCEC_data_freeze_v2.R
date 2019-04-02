# Yige Wu @ WashU 2018 Aug
## For formating UCEC data in data freeze to be in concord with CDAP proteomics data and genomics data


# source ------------------------------------------------------------------
wd <- getwd()
if (wd != "/Users/yigewu/Box Sync") {
  setwd("/Users/yigewu/Box Sync")
}
source("./cptac2p_analysis/phospho_network/phospho_network_shared.R")

# set variables -----------------------------------------------------------
cancer <- "UCEC"
pipeline_type <- "PGDAC"

## input UCEC meta data
meta_tab <- fread("./Ding_Lab/Projects_Current/CPTAC/PGDAC/Endometrium_CPTAC3/01_Data_tables/Baylor_DataFreeze_V2/UCEC_CPTAC3_meta_table_v2.0.txt", data.table = F)
rownames(meta_tab) <- meta_tab$idx


# input bobo’s PTMcosmos to add peptide ID --------------------------------
uniprot_ensembl_map <- fread("./cptac2p/resources/PTMcosmos/Data dump/2019-01-11/all_uniprot_ensembl_id_mappings.tsv", data.table = F)
ptm_sites_map <- fread("./cptac2p/resources/PTMcosmos/Data dump/2019-01-11/all_ptm_sites.all_txs.tsv", data.table = F)


# rewrite proteomics/phosphoproteomics data split by tumor and normal ----------------------------------
for (sample_type in c("tumor")) {
  pro_data <- fread(input = "Ding_Lab/Projects_Current/CPTAC/PGDAC/Endometrium_CPTAC3/01_Data_tables/Baylor_DataFreeze_V2/UCEC_proteomics_PNNL_ratio_median_polishing_log2_V2.0.cct", data.table = F)
  head(pro_data)
  col_names_origin <- colnames(pho_data)[!(colnames(pho_data) %in% c("idx"))]
  col_names_keep <- meta_tab$idx[grepl(x = meta_tab$Proteomics_Tumor_Normal, pattern = sample_type, ignore.case = T)]
  data2w <- pro_data[, c("idx", col_names_keep)]
  col_names_sample_type2w <- meta_tab[col_names_keep, "Proteomics_Participant_ID"]
  colnames(data2w) <- c("Gene", col_names_sample_type2w)
  col_names_sample_type2w <- col_names_sample_type2w[!grepl(pattern = "rep", x = col_names_sample_type2w) & !(col_names_sample_type2w %in% c("C3L-00084", "C3L-01284", "C3N-01001", "C3L-00356", "C3L-00938", "C3L-01247","C3L-01253"))]
  data2w <- data2w[, c("Gene", col_names_sample_type2w)]
  expresson_type <- "PRO"
  # write.table(x = data2w, file = paste0(makeOutDir(resultD = resultD), cancer, "_", expresson_type, "_", sample_type, "_", pipeline_type, "_", "median_polishing", "_", "partID", ".txt"), quote = F, row.names = F, sep = "\t")
  
  pho_data <- fread(input = "Ding_Lab/Projects_Current/CPTAC/PGDAC/Endometrium_CPTAC3/01_Data_tables/Baylor_DataFreeze_V2/UCEC_phosphoproteomics_PNNL_ratio_median_polishing_site_level_log2_V2.0.cct", data.table = F)
  head(pho_data)
  col_names_origin <- colnames(pho_data)[!(colnames(pho_data) %in% c("idx"))]
  col_names_keep <- meta_tab$idx[grepl(x = meta_tab$Proteomics_Tumor_Normal, pattern = sample_type, ignore.case = T)]
  data2w <- pho_data[, c("idx", col_names_keep)]
  col_names_sample_type2w <- meta_tab[col_names_keep, "Proteomics_Participant_ID"]
  colnames(data2w) <- c("idx", col_names_sample_type2w)
  col_names_sample_type2w <- col_names_sample_type2w[!grepl(pattern = "rep", x = col_names_sample_type2w) & !(col_names_sample_type2w %in% c("C3L-00084", "C3L-01284", "C3N-01001", "C3L-00356", "C3L-00938", "C3L-01247","C3L-01253"))]
  data2w <- data2w[, c("idx", col_names_sample_type2w)]
  phosphosite_head <- str_split_fixed(string = data2w$idx, pattern = "-", 2)
  phosphosite_head %>% head()
  data2w$Gene <- phosphosite_head[,1]
  data2w$Phosphosite <- phosphosite_head[,2]
  data2w <- data2w[, c("Gene", "Phosphosite", col_names_sample_type2w)]
  expresson_type <- "PHO"
  stop("")
  # write.table(x = data2w, file = paste0(makeOutDir(resultD = resultD), cancer, "_", expresson_type, "_", sample_type, "_", pipeline_type, "_", "median_polishing", "_", "partID", ".txt"), quote = F, row.names = F, sep = "\t")
  
  phog_data <- fread(input = "Ding_Lab/Projects_Current/CPTAC/PGDAC/Endometrium_CPTAC3/01_Data_tables/Baylor_DataFreeze_V2/UCEC_phosphoproteomics_PNNL_ratio_median_polishing_gene_level_log2_V2.0.cct", data.table = F)
  head(phog_data)
  col_names_origin <- colnames(pho_data)[!(colnames(pho_data) %in% c("idx"))]
  col_names_keep <- meta_tab$idx[grepl(x = meta_tab$Proteomics_Tumor_Normal, pattern = sample_type, ignore.case = T)]
  data2w <- phog_data[, c("idx", col_names_keep)]
  col_names_sample_type2w <- meta_tab[col_names_keep, "Proteomics_Participant_ID"]
  colnames(data2w) <- c("Gene", col_names_sample_type2w)
  col_names_sample_type2w <- col_names_sample_type2w[!grepl(pattern = "rep", x = col_names_sample_type2w) & !(col_names_sample_type2w %in% c("C3L-00084", "C3L-01284", "C3N-01001", "C3L-00356", "C3L-00938", "C3L-01247","C3L-01253"))]
  data2w <- data2w[, c("Gene", col_names_sample_type2w)]
  expresson_type <- "collapsed_PHO"
  # write.table(x = data2w, file = paste0(makeOutDir(resultD = resultD), cancer, "_", expresson_type, "_", sample_type, "_", pipeline_type, "_", "median_polishing", "_", "partID", ".txt"), quote = F, row.names = F, sep = "\t")
}

# Parse Baylor’s CNV ------------------------------------------------------
cna <- fread(input = "./Ding_Lab/Projects_Current/CPTAC/PGDAC/Endometrium_CPTAC3/01_Data_tables/Baylor_DataFreeze_V2/UCEC_WGS_CNA_gene_level_V2.0.cct", data.table = F)
colnames(cna) <- c("gene", meta_tab[colnames(cna)[-1], "Proteomics_Participant_ID"])
head(cna)
# write.table(x = cna, file = paste0(makeOutDir(resultD = resultD), "somatic_CNA", ".", cancer, ".partID.txt"), quote = F, row.names = F, sep = "\t")


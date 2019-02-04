# Yige Wu @ WashU 2018 Aug
## For formating UCEC CDAP proteomics data and genomics data
## output file name format: {cancer}_{expression_type}_{sample_type}_{pipeline_type}_{unnormalized/normalization_method}_{aliquot_id/participant_id_labeled}.txt

# source ------------------------------------------------------------------
wd <- getwd()
if (wd != "/Users/yigewu/Box Sync") {
  setwd("/Users/yigewu/Box Sync")
}
source("./cptac2p_analysis/phospho_network/phospho_network_shared.R")
source("./cptac2p_analysis/preprocess_files/preprocess_files_shared.R")

# set variables -----------------------------------------------------------
cancer <- "UCEC"
pipeline_type <- "CDAP"

## input UCEC meta data
meta_tab <- fread("./Ding_Lab/Projects_Current/CPTAC/PGDAC/Endometrium_CPTAC3/01_Data_tables/Baylor_DataFreeze_V2/UCEC_CPTAC3_meta_table_v2.0.txt", data.table = F)
rownames(meta_tab) <- meta_tab$idx
meta_tab %>% head()

# process protein data ------------------------------------------------------
expresson_type <- "PRO"
pro_data <- fread(input = "Ding_Lab/Projects_Current/CPTAC/cptac_shared/7_CPTAC3_Uterine_Corpus_Endometrial_Carcinoma/A_CPTAC3_hg38_CDAP/UCEC_SummaryReports/CPTAC3_Uterine_Corpus_Endometrial_Carcinoma_Proteome.tmt10.tsv", data.table = F)
## keep only the unshared log ratios
col_names_origin <- colnames(pro_data)
col_names_keep <- col_names_origin[grepl(x = col_names_origin, pattern = "Unshared")]
col_names_keep
col_names_keep_sim <- str_split_fixed(string = col_names_keep, pattern = " ", n = 4)[,1]
col_names_keep_sim %>%
  length()
pro_data <- pro_data[,c("Gene", col_names_keep)]
colnames(pro_data) <- c("Gene", col_names_keep_sim)

## remove the first 3 rows
rownames(pro_data) <- pro_data$Gene
pro_data <- pro_data[rownames(pro_data)[!(rownames(pro_data) %in% c("Mean", "Median", "StdDev"))],]

## check whether this file contains all the samples
meta_tab[!(meta_tab$Proteomics_Aliquot_ID %in% col_names_keep_sim),] %>% head()
meta_tab[(meta_tab$Proteomics_Aliquot_ID %in% col_names_keep_sim) & meta_tab$Proteomics_Tumor_Normal == "Tumor",] %>% nrow()
meta_tab$Proteomics_Tumor_Normal %>%
  unique()
### 104 tumor samples, should be good
### save tumor only
for (sample_type in c("tumor")) {
  col_names_sample_type <- meta_tab$Proteomics_Aliquot_ID[(meta_tab$Proteomics_Aliquot_ID %in% col_names_keep_sim) & grepl(pattern = sample_type, x = meta_tab$Proteomics_Tumor_Normal, ignore.case = T)]
  col_names_sample_type2w <- meta_tab$Proteomics_Participant_ID[(meta_tab$Proteomics_Aliquot_ID %in% col_names_keep_sim) & grepl(pattern = sample_type, x = meta_tab$Proteomics_Tumor_Normal, ignore.case = T)]
  data2w <- pro_data[, c("Gene", col_names_sample_type)]
  colnames(data2w) <- c("Gene", col_names_sample_type2w)
  col_names_sample_type2w <- col_names_sample_type2w[!grepl(pattern = "rep", x = col_names_sample_type2w)]
  data2w <- data2w[, c("Gene", col_names_sample_type2w)]
  write.table(x = data2w, file = paste0(makeOutDir(resultD = resultD), cancer, "_", expresson_type, "_", sample_type, "_", pipeline_type, "_", "unnormalized", "_", "partID", ".txt"), quote = F, row.names = F, sep = "\t")
  
  ## scale by sd and mean
  pro_data_mat <- data2w[, col_names_sample_type2w]
  pro_data_mat.n <- scale_by_sample(m = pro_data_mat)
  data2w <- data.frame(Gene = data2w$Gene)
  data2w <- cbind(data2w, pro_data_mat.n)
  write.table(x = data2w, file = paste0(makeOutDir(resultD = resultD), cancer, "_", expresson_type, "_", sample_type, "_", pipeline_type, "_", "scaled", "_", "partID", ".txt"), quote = F, row.names = F, sep = "\t")
}
## next do QC plots


# rewrite phosphorylation data ----------------------------------------------
pho_data <- fread(input = "Ding_Lab/Projects_Current/CPTAC/cptac_shared/7_CPTAC3_Uterine_Corpus_Endometrial_Carcinoma/A_CPTAC3_hg38_CDAP/UCEC_SummaryReports/CPTAC3_Uterine_Corpus_Endometrial_Carcinoma_Phosphoproteome.phosphosite.tmt10.tsv", data.table = F)
pho_data %>% head()
## phosphosite data didn't distinguish between shared and unshared log ratios
col_names_origin <- colnames(pho_data)
col_names_keep <- col_names_origin[!(col_names_origin %in% c("Peptide", "Gene", "Organism", "Phosphosite"))]
col_names_keep
col_names_keep_sim <- str_split_fixed(string = col_names_keep, pattern = " ", n = 4)[,1]
col_names_keep_sim %>%
  length()
pho_data <- pho_data[,c("Gene", "Phosphosite", col_names_keep)]
colnames(pho_data) <- c("Gene", "Phosphosite", col_names_keep_sim)

## edit the phosphosite column
phosphosite_tmp <- str_split_fixed(string = pho_data$Phosphosite, pattern = ":", 2)
phosphosite_tmp %>% head()
pho_data$Phosphosite <- toupper(phosphosite_tmp[,2])
pho_data$Peptide_ID <- phosphosite_tmp[,1]
pho_data <- pho_data[,c("Gene", "Phosphosite", "Peptide_ID", col_names_keep_sim)]
head(pho_data)
expresson_type <- "PHO"
for (sample_type in c("tumor")) {
  col_names_sample_type <- meta_tab$Proteomics_Aliquot_ID[(meta_tab$Proteomics_Aliquot_ID %in% col_names_keep_sim) & grepl(pattern = sample_type, x = meta_tab$Proteomics_Tumor_Normal, ignore.case = T)]
  col_names_sample_type2w <- meta_tab$Proteomics_Participant_ID[(meta_tab$Proteomics_Aliquot_ID %in% col_names_keep_sim) & grepl(pattern = sample_type, x = meta_tab$Proteomics_Tumor_Normal, ignore.case = T)]
  data2w <- pho_data[, c("Gene", "Phosphosite", "Peptide_ID", col_names_sample_type)]
  colnames(data2w) <- c("Gene", "Phosphosite", "Peptide_ID", col_names_sample_type2w)
  col_names_sample_type2w <- col_names_sample_type2w[!grepl(pattern = "rep", x = col_names_sample_type2w)]
  data2w <- data2w[, c("Gene", "Phosphosite", "Peptide_ID", col_names_sample_type2w)]
  write.table(x = data2w, file = paste0(makeOutDir(resultD = resultD), cancer, "_", expresson_type, "_", sample_type, "_", pipeline_type, "_", "unnormalized", "_", "partID", ".txt"), quote = F, row.names = F, sep = "\t")
  
  ## scale by sd and mean
  pho_data_mat <- data2w[, col_names_sample_type2w]
  pho_data_mat.n <- scale_by_sample(m = pho_data_mat)
  data2w <- data2w[, c("Gene", "Phosphosite", "Peptide_ID")]
  data2w <- cbind(data2w, pho_data_mat.n)
  write.table(x = data2w, file = paste0(makeOutDir(resultD = resultD), cancer, "_", expresson_type, "_", sample_type, "_", pipeline_type, "_", "scaled", "_", "partID", ".txt"), quote = F, row.names = F, sep = "\t")
}

# Yige Wu @ WashU 2018 Aug
## For formating CCRCC data in data freeze to be in concord with CDAP proteomics data and genomics data


# source ------------------------------------------------------------------
wd <- getwd()
if (wd != "/Users/yigewu/Box Sync") {
  setwd("/Users/yigewu/Box Sync")
}
source("./cptac2p_analysis/phospho_network/phospho_network_shared.R")

# set variables -----------------------------------------------------------
cancer <- "CCRCC"
pipeline_type <- "PGDAC"

# rewrite protein data (Umich) ------------------------------------------------------
pro_data <- fread(input = "Ding_Lab/Projects_Current/CPTAC/PGDAC/JHU_PCC_shared_data/CCRCC_analysis/proteomics/global/Michigan/W_Batch1-5/6_CPTAC3_Clear_Cell_Renal_Cell_Carcinoma_W_abundance_groupby=0_protNorm=2.tsv", data.table = F)
pro_data <- data.frame(pro_data)
ncol(pro_data)
samples <- colnames(pro_data)[!(colnames(pro_data) %in% c("Index", "NumberPSM", "Proteins", "ReferenceIntensity"))]
pro_values.a <- pro_data[, samples]
pro_values.r <- pro_data[, "ReferenceIntensity"]
pro_values.a_r <- pro_values.a - pro_values.r
pro_data <- cbind(data.frame(Gene = pro_data$Index), pro_values.a_r)
write.table(x = pro_data, file = paste0(makeOutDir(resultD = resultD), "6_CPTAC3_Clear_Cell_Renal_Cell_Carcinoma_W_abundance_groupby=0_protNorm=2.tsv.formatted.txt"), quote = F, row.names = F, sep = "\t")


# Rewrite phosphorylation gene level --------------------------------------
phog_data <- fread(input = "Ding_Lab/Projects_Current/CPTAC/PGDAC/JHU_PCC_shared_data/CCRCC_analysis/proteomics/phospho/Michigan/P_Batch1-5/6_CPTAC3_Clear_Cell_Renal_Cell_Carcinoma_P_abundance_groupby=0_protNorm=2.tsv", data.table = F)
phog_data <- data.frame(phog_data)
ncol(phog_data)
head(phog_data)
samples <- colnames(phog_data)[!(colnames(phog_data) %in% c("Index", "NumberPSM", "Proteins", "ReferenceIntensity"))]
pro_values.a <- phog_data[, samples]
pro_values.r <- phog_data[, "ReferenceIntensity"]
pro_values.a_r <- pro_values.a - pro_values.r
phog_data <- cbind(data.frame(Gene = phog_data$Index), pro_values.a_r)
write.table(x = phog_data, file = paste0(makeOutDir(resultD = resultD), "6_CPTAC3_Clear_Cell_Renal_Cell_Carcinoma_P_abundance_groupby=0_protNorm=2.tsv.formatted.txt"), quote = F, row.names = F, sep = "\t")

# rewrite phosphorylation data ----------------------------------------------
pho_data <- fread(input = "Ding_Lab/Projects_Current/CPTAC/PGDAC/JHU_PCC_shared_data/CCRCC_analysis/proteomics/phospho/Michigan/P_Batch1-5/6_CPTAC3_Clear_Cell_Renal_Cell_Carcinoma_P_abundance_groupby=3_protNorm=2.tsv", data.table = F)
pho_data <- data.frame(pho_data)
pho_values.a <- pho_data[, samples]
pho_values.r <- pho_data[, "ReferenceIntensity"]
pho_values.a_r <- pho_values.a - pho_values.r
pho_data <- cbind(data.frame(Gene = pho_data$Gene, Phosphosite = str_split_fixed(string = pho_data$Index, pattern = "_", n = 7)[,7]), pho_values.a_r)
pho_data <- pho_data[pho_data$Phosphosite != "",]
write.table(x = pho_data, file = paste0(makeOutDir(resultD = resultD), "6_CPTAC3_Clear_Cell_Renal_Cell_Carcinoma_P_abundance_groupby=3_protNorm=2.tsv.formatted.txt"), quote = F, row.names = F, sep = "\t")

# input sample mapping file -----------------------------------------------
sample_map <- fread(input = "./Ding_Lab/Projects_Current/CPTAC/PGDAC/A6_PGDAC_Data_Clear_Cell_Renal_Cell_Carcinoma/MountSinai/CPTAC3_CCRCC/Batch1-4/CPTAC-3 CCRCC samples tumor-normal info.csv", data.table = F)
sample_map <- unique(sample_map[, c("Aliquot ID", "Case ID", "Tissue Type")])
rownames(sample_map) <- sample_map$`Aliquot ID`
sample_map <- sample_map[!(sample_map$`Case ID` %in% c("C3N-01175", "C3N-01180", "C3N-00313", "C3N-00435", "C3N-00832")),]

# rewrite proteomics/phosphoproteomics data split by tumor and normal ----------------------------------
for (sample_type in c("tumor", "normal")) {
  sampIDs <- sample_map$`Aliquot ID`[grepl(x = sample_map$`Tissue Type`, pattern= sample_type, ignore.case = T)]
  sampIDs <- intersect(sampIDs, samples)
  partIDs <- sample_map[sampIDs, "Case ID"]
  
  ## write protein data
  pro2w <- pro_data[, c("Gene", sampIDs)]
  pro2w.partID <- pro2w
  colnames(pro2w.partID) <- c("Gene", partIDs)
  expresson_type <- "PRO"
  write.table(x = pro2w.partID, file = paste0(makeOutDir(resultD = resultD), cancer, "_", expresson_type, "_", sample_type, "_", pipeline_type, "_", "MD_MAD", "_", "partID", ".txt"), quote = F, row.names = F, sep = "\t")
  
  ## write gene-level phosphorylation data
  phog2w <- phog_data[, c("Gene", sampIDs)]
  phog2w.partID <- phog2w
  colnames(phog2w.partID) <- c("Gene", partIDs)
  expresson_type <- "collapsed_PHO"
  write.table(x = phog2w.partID, file = paste0(makeOutDir(resultD = resultD), cancer, "_", expresson_type, "_", sample_type, "_", pipeline_type, "_", "MD_MAD", "_", "partID", ".txt"), quote = F, row.names = F, sep = "\t")
  
  
  ## write phosphorylatin data
  pho2w <- pho_data[, c("Gene", "Phosphosite", sampIDs)]
  pho2w.partID <- pho2w
  colnames(pho2w.partID) <- c("Gene", "Phosphosite", partIDs)
  expresson_type <- "PHO"
  write.table(x = pho2w.partID, file = paste0(makeOutDir(resultD = resultD), cancer, "_", expresson_type, "_", sample_type, "_", pipeline_type, "_", "MD_MAD", "_", "partID", ".txt"), quote = F, row.names = F, sep = "\t")
}

# Parse Michigan’s CNV ------------------------------------------------------
cna <- fread(input = "./Ding_Lab/Projects_Current/CPTAC/PGDAC/JHU_PCC_shared_data/CCRCC_analysis/genomics/WXS_CNV/michigan/kirc_cnv_genelr.csv", data.table = F)
cna <- cna[, c("gene_name", colnames(cna)[grepl(pattern = "lr_", x = colnames(cna))])]
colnames(cna) <- str_split_fixed(string = colnames(cna), pattern = "lr_", 2)[,2]
colnames(cna)[1] <- "gene"
cna.val <- as.matrix(cna[,-1])
cna.val[is.infinite(cna.val)] <- NA
colnames_tmp <- colnames(cna)
cna <- cbind(data.frame(gene = cna$gene), data.frame(cna.val))
colnames(cna) <- colnames_tmp
cancer <- "CCRCC"
write.table(x = cna, file = paste0(makeOutDir(resultD = resultD), "somatic_CNA", ".", cancer, ".partID.txt"), quote = F, row.names = F, sep = "\t")



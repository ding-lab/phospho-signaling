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

## input UCEC meta data
meta_tab <- fread("./Ding_Lab/Projects_Current/CPTAC/PGDAC/Endometrium_CPTAC3/01_Data_tables/UCEC_V1.1/UCEC_CPTAC3_meta_table.txt", data.table = F)
rownames(meta_tab) <- meta_tab$idx

# rewrite protein data ------------------------------------------------------
pro_data <- fread(input = "Ding_Lab/Projects_Current/CPTAC/PGDAC/Endometrium_CPTAC3/01_Data_tables/UCEC_V1/UCEC_proteomics_V1.cct", data.table = F)
pro_data <- data.frame(pro_data)
colnames(pro_data)[1] <- "Gene"
write.table(x = pro_data, file = paste0(makeOutDir(resultD = resultD), "UCEC_proteomics_V1.cct.formatted.txt"), quote = F, row.names = F, sep = "\t")

# rewrite phosphorylation data ----------------------------------------------
pho_data <- fread(input = "Ding_Lab/Projects_Current/CPTAC/PGDAC/Endometrium_CPTAC3/01_Data_tables/UCEC_V1/UCEC_phosphoproteomics_site_level_V1.cct", data.table = F)
pho_data <- data.frame(pho_data)
pho_data$Gene <- str_split_fixed(string = pho_data$idx, pattern = "-", 2)[,1]
pho_data$Phosphosite <- toupper(str_split_fixed(string = pho_data$idx, pattern = "-", 2)[,2])
pho_data <- pho_data[, c("Gene", "Phosphosite", colnames(pho_data)[!(colnames(pho_data) %in% c("idx", "Gene", "Phosphosite"))])]
write.table(x = pho_data, file = paste0(makeOutDir(resultD = resultD), "UCEC_phosphoproteomics_site_level_V1.cct.formatted.txt"), quote = F, row.names = F, sep = "\t")

# rewrite proteomics/phosphoproteomics data split by tumor and normal ----------------------------------
for (sample_type in c("tumor", "normal")) {
  samples <- colnames(pro_data)[!(colnames(pro_data) %in% c("idx")) & !grepl(pattern = "Mix", x = (colnames(pro_data)))  & !grepl(pattern = "rep", x = (colnames(pro_data)))]
  samples <- samples[grepl(pattern = paste0('\\.', toupper(substr(x = sample_type, start = 1, stop = 1))), x = samples)]
  sampIDs <- samples
  tmp <- str_split_fixed(string = samples, pattern = "\\.", n = 3)
  partIDs <- paste0(tmp[,1], "-", tmp[,2])
  
  ## write protein data
  pro2w <- pro_data[, c("Gene", sampIDs)]
  write.table(x = pro2w, file = paste0(makeOutDir(resultD = resultD), "UCEC_proteomics_V1.cct.formatted", ".", sample_type, ".txt"), quote = F, row.names = F, sep = "\t")
  
  pro2w.partID <- pro2w
  colnames(pro2w.partID) <- c("Gene", partIDs)
  write.table(x = pro2w.partID, file = paste0(makeOutDir(resultD = resultD), "PRO", ".", cancer, ".", sample_type, ".partID.txt"), quote = F, row.names = F, sep = "\t")
  
  ## write phosphorylatin data
  pho2w <- pho_data[, c("Gene", "Phosphosite", sampIDs)]
  write.table(x = pho2w, file = paste0(makeOutDir(resultD = resultD), "UCEC_phosphoproteomics_site_level_V1.cct.formatted", ".", sample_type, ".txt"), quote = F, row.names = F, sep = "\t")
  
  pho2w.partID <- pho2w
  colnames(pho2w.partID) <- c("Gene", "Phosphosite", partIDs)
  write.table(x = pho2w.partID, file = paste0(makeOutDir(resultD = resultD), "PHO", ".", cancer, ".", sample_type, ".partID.txt"), quote = F, row.names = F, sep = "\t")
}


# Parse Baylor’s CNV ------------------------------------------------------
cna <- fread(input = "./Ding_Lab/Projects_Current/CPTAC/PGDAC/Endometrium_CPTAC3/01_Data_tables/UCEC_V1.1/UCEC_CNA_gene_level_hg19_V1.1.cct", data.table = F)
colnames(cna) <- c("gene", meta_tab[colnames(cna)[-1], "Proteomics_Participant_ID"])
write.table(x = cna, file = paste0(makeOutDir(resultD = resultD), "somatic_CNA", ".", cancer, ".partID.txt"), quote = F, row.names = F, sep = "\t")


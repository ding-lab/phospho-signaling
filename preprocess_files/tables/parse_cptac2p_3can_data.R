# Yige Wu @ WashU 2018 Aug
## For formating CCRCC data in data freeze to be in concord with CDAP proteomics data and genomics data


# source ------------------------------------------------------------------
wd <- getwd()
if (wd != "/Users/yigewu/Box Sync") {
  setwd("/Users/yigewu/Box Sync")
}
source("./cptac2p_analysis/phospho_network/phospho_network_shared.R")
pipeline_type <- "CDAP"

# rewrite proteomics/phosphoproteomics data split by tumor and normal ----------------------------------
for (sample_type in c("tumor", "normal")) {
  for (cancer in cancers_sort) {
    pro_data <- fread(input = paste0("./Ding_Lab/Projects_Current/CPTAC/cptac_shared/", cancer, "/", prefix[cancer], "_PRO_formatted_normalized_", ifelse(sample_type == "tumor", "noControl", "Control"),".txt"), data.table = F)
    sampIDs <- colnames(pro_data)[!(colnames(pro_data) %in% c("Gene"))]
    partIDs <- sampID2partID(sampleID_vector = sampIDs, sample_map = loadSampMap())
    
    pro2w.partID <- pro_data
    colnames(pro2w.partID) <- c("Gene", partIDs)
    expresson_type <- "PRO"
    write.table(x = pro2w.partID, file = paste0(makeOutDir(resultD = resultD), cancer, "_", expresson_type, "_", sample_type, "_", pipeline_type, "_", "scaled", "_", "partID", ".txt"), quote = F, row.names = F, sep = "\t")
    
    pho_data <- fread(input = paste0("./Ding_Lab/Projects_Current/CPTAC/cptac_shared/", cancer, "/", prefix[cancer], "_PHO_formatted_normalized_", ifelse(sample_type == "tumor", "noControl", "Control"),".txt"), data.table = F)
    pho_data$Peptide_ID <- toupper(str_split_fixed(string = pho_data$Phosphosite, pattern = "\\:", 2)[,1])
    pho_data$Phosphosite <- toupper(str_split_fixed(string = pho_data$Phosphosite, pattern = "\\:", 2)[,2])
    
    pho2w.partID <- pho_data
    pho2w.partID <- pho2w.partID[,c("Gene", "Phosphosite", "Peptide_ID", sampIDs)]
    colnames(pho2w.partID) <- c("Gene", "Phosphosite", "Peptide_ID", partIDs)
    expresson_type <- "PHO"
    write.table(x = pho2w.partID, file = paste0(makeOutDir(resultD = resultD), cancer, "_", expresson_type, "_", sample_type, "_", pipeline_type, "_", "scaled", "_", "partID", ".txt"), quote = F, row.names = F, sep = "\t")
    
    pho_collapsed_data <- fread(input = paste0("./Ding_Lab/Projects_Current/CPTAC/cptac_shared/", cancer, "/", prefix[cancer], "_collapsed_PHO_formatted_normalized_replicate_averaged_", ifelse(sample_type == "tumor", "Tumor", "Normal"),".txt"), data.table = F)
    pho_collapsed2w.partID <- pho_collapsed_data
    colnames(pho_collapsed2w.partID) <- c("Gene", partIDs)
    expresson_type <- "collapsed_PHO"
    write.table(x = pho_collapsed2w.partID, file = paste0(makeOutDir(resultD = resultD), cancer, "_", expresson_type, "_", sample_type, "_", pipeline_type, "_", "scaled", "_", "partID", ".txt"), quote = F, row.names = F, sep = "\t")
  }
}

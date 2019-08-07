# Yige Wu @ WashU 2019 Apr
## mark patient with TP53 protein expression status and only deep deletion
### For samples without TP53 protein level the protein status will be NA

# source ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/p53/TP53_shared.R')


# set variables -----------------------------------------------------------
# cancers2process <- "BRCA"
# cancers2process <- "OV"
# cancers2process <- "CO"
# cancers2process <- "UCEC"
# cancers2process <- "CCRCC"
cancers2process <- c("LIHC", "CCRCC", "UCEC", "CO", "OV", "BRCA")

gene_altered <- "TP53"

meta_tab <- read_excel("./Ding_Lab/Projects_Current/TP53_shared_data/resources/TP53_samples_classification/HCC Clinical information and follow-up information.xlsx", 
                       sheet = "Table S1 Clinical and follow-up")
meta_tab <- data.frame(meta_tab)
meta_tab$partID <- paste0("LIHC", meta_tab$Case.ID)
rownames(meta_tab) <- as.vector(meta_tab$partID)
meta_tab$sampID.tumor <- str_split_fixed(string = meta_tab$Tumor..T..sample.ID, pattern = "T", 2)[,2]
meta_tab$sampID.normal <- str_split_fixed(string = meta_tab$Adjacent.Non.tumor.liver.tissue..N..sample.ID, pattern = "N", 2)[,2]

getChinaLiverSampID2PartID <- function(sampIDs) {
  partIDs <- sapply(sampIDs, FUN = function(sampID, meta_tab) {
    new_partID <- meta_tab$partID[meta_tab$sampID.normal == sampID | meta_tab$sampID.tumor == sampID]
    return(new_partID)
  }, meta_tab = meta_tab)
  return(partIDs)
}


for (cancer_tmp in cancers2process) {
  # input previous patient table --------------------------------------------
  part_tab <- fread(input = paste0("Ding_Lab/Projects_Current/TP53_shared_data/resources/TP53_samples_classification/", cancer_tmp, "_tp53_mut_classification_table.txt"), data.table = F)
  if (cancer_tmp == "LIHC") {
    part_tab$V1 <- getChinaLiverSampID2PartID(part_tab$V1)
  }
  # input protein data ------------------------------------------------------
  pro_tab <- loadTP53proteomics(cancer = cancer_tmp, expression_type = "PRO")
  pro_tab <- pro_tab %>%
    filter(Gene == gene_altered)
  pro_tab.m <- melt(data = pro_tab)
 
  # add the column for high p53 protein -------------------------------------
  ## high expressor is defined as those higher than 90% quantile
  exp_high <- quantile(x = pro_tab.m$value, probs = 0.9, na.rm = T)
  exp_partIDs <- unique(pro_tab.m$variable)
  exp_high_partIDs <- pro_tab.m$variable[pro_tab.m$value >= exp_high]
  part_tab$p53_protein_status <- "NA"
  part_tab$p53_protein_status[part_tab$V1 %in% exp_partIDs] <- "< 90% percentile"
  part_tab$p53_protein_status[part_tab$V1 %in% exp_high_partIDs] <- ">= 90% percentile"
  
  # add the column for deep deletion ----------------------------------------
  cna_tab <- loadTP53DeepDeletion(cancer = cancer_tmp)
  cna_tab.m <- melt(data = cna_tab, id.vars = "gene")
  deep_deletion_partIDs <- cna_tab.m$variable[cna_tab.m$value == "deletion"]
  part_tab$Deep_Deletion <- FALSE
  part_tab$Deep_Deletion[part_tab$V1 %in% deep_deletion_partIDs] <- TRUE
  
  # write patient table -----------------------------------------------------
  part_tab$Classification_complex <- "WT"
  part_tab$Classification_complex[part_tab$Missense == T & part_tab$p53_protein_status == ">= 90% percentile"] <- "Missense_protein_level_high"
  part_tab$Classification_complex[part_tab$Missense == T & part_tab$p53_protein_status == "< 90% percentile"] <- "Missense_protein_level_average"
  part_tab$Classification_complex[part_tab$Missense == T & part_tab$p53_protein_status == "NA"] <- "Missense_protein_level_NA"
  part_tab$Classification_complex[(part_tab$Missense == F & part_tab$Truncation == T) | (part_tab$Missense == F & part_tab$Other_mutations & grepl(pattern = "S", x = part_tab$TP53)) ] <- "Truncation"
  part_tab$Classification_complex[(part_tab$Missense == F & part_tab$Truncation == F & part_tab$Deep_Deletion == T) | (part_tab$Missense == F & part_tab$Truncation == F & part_tab$Other_mutations & grepl(pattern = "I", x = part_tab$TP53))] <- "Deletion"
  part_tab$Deletion <- NULL
  
  table(part_tab$Classification_complex) %>%
    print()
  
  write.table(x = part_tab, file = paste0(makeOutDir(resultD = resultD), cancer_tmp, "_tp53_mut_protein_classification_table.txt"), quote = F, sep = "\t", row.names = F)
}

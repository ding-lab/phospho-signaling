# Yige Wu @ WashU Mar. 2019
## get the number of unique phosphosites that were sufficiently detected and variably expressed

# source ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_plotting.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')

# set cancers to process --------------------------------------------------
cancers2process <- c("BRCA", "OV", "CO", "UCEC", "CCRCC")


# set the least number of samples detected --------------------------------
num_nonNA <- 20

# get the phosphosites that are sufficiently detected across cohorts --------
pho_head_cancers <- NULL
for (cancer in cancers2process) {
  if (cancer %in% c("BRCA", "OV", "CO")) {
    pho_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "PHO", sample_type = "tumor", pipeline_type = "CDAP", norm_type = "scaled")
  } else if (cancer == "UCEC") {
    pho_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "PHO", sample_type = "tumor", pipeline_type = "PGDAC", norm_type = "median_polishing")
  } else if (cancer == "CCRCC") {
    pho_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "PHO", sample_type = "tumor", pipeline_type = "PGDAC", norm_type = "MD_MAD")
  } else if (cancer == "LIHC") {
    pho_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "PHO", sample_type = "tumor", pipeline_type = "PGDAC", norm_type = "MD")
  }
  partIDs <- colnames(pho_tab)
  partIDs <- partIDs[!(partIDs %in% c("Gene", "Phosphosite", "Peptide_ID"))]
  
  ## detected in at least 20 samples
  pho_tab <- pho_tab[rowSums(!is.na(pho_tab[,partIDs])) >= num_nonNA,]
  pho_head <- data.frame(Gene = pho_tab$Gene, Phosphosite = pho_tab$Phosphosite)
  pho_head_cancers <- rbind(pho_head_cancers, pho_head)
}

tmp <- str_split_fixed(string = pho_head_cancers$Phosphosite, pattern = "[STY]", 3)
pho_head_cancers$is.solo <- (tmp[,3] == "")
pho_head_cancers$p_coord <- tmp[,2]
pho_head_cancers$residue <- sapply(X = pho_head_cancers$Phosphosite, FUN = function(phosphosite) substr(x = phosphosite, start = 1, stop = 1))
pho_head_cancers <- pho_head_cancers %>%
  mutate(id = paste0(Gene, "_", Phosphosite)) %>%
  filter(is.solo == T) %>%
  unique

pho_head_cancers %>%
  head

pho_head_cancers %>%
  select(id) %>%
  unique %>%
  nrow()


# get the phosphosites that are sufficiently detected and variabley expressed across cohorts --------
pho_head_cancers <- NULL
for (cancer in cancers2process) {
  if (cancer %in% c("BRCA", "OV", "CO")) {
    pho_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "PHO", sample_type = "tumor", pipeline_type = "CDAP", norm_type = "scaled")
    pro_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "PRO", sample_type = "tumor", pipeline_type = "CDAP", norm_type = "scaled")
    
  } else if (cancer == "UCEC") {
    pho_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "PHO", sample_type = "tumor", pipeline_type = "PGDAC", norm_type = "median_polishing")
    pro_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "PRO", sample_type = "tumor", pipeline_type = "PGDAC", norm_type = "median_polishing")
    
  } else if (cancer == "CCRCC") {
    pho_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "PHO", sample_type = "tumor", pipeline_type = "PGDAC", norm_type = "MD_MAD")
    pro_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "PRO", sample_type = "tumor", pipeline_type = "PGDAC", norm_type = "MD_MAD")
    
  } else if (cancer == "LIHC") {
    pho_tab <- loadParseProteomicsData(cancer = cancer, expression_type  = "PHO", sample_type = "tumor", pipeline_type = "PGDAC", norm_type = "MD")
  }
  partIDs <- colnames(pho_tab)
  partIDs <- partIDs[!(partIDs %in% c("Gene", "Phosphosite", "Peptide_ID"))]
  
  ## detected in at least 20 samples
  pho_tab <- pho_tab[rowSums(!is.na(pho_tab[,partIDs])) >= num_nonNA,]
  
  ## get the standard deviations per site
  pho_sds <- sapply(1:nrow(pho_tab), FUN = function(n, mat) sd(mat[n,], na.rm = T), mat = as.matrix(pho_tab[, partIDs]))
  boxplot(pho_sds)
  pho_sds_iqr <- IQR(pho_sds)
  
  # pho_sds_outlier_low <- quantile(x = pho_sds, probs = 0.25, na.rm = T) - 1*pho_sds_iqr
  pho_sds_outlier_low <- quantile(x = pho_sds, probs = 0.5, na.rm = T) - 1*pho_sds_iqr
  # print(pho_sds_outlier_low)
  # print(length(which(pho_sds < pho_sds_outlier_low)))
  # test <- pho_tab[pho_sds <= pho_sds_outlier_low,]
  # table(test$Gene)
  # stop("")
  
  pho_tab <- pho_tab[pho_sds >= pho_sds_outlier_low,]
  # print(nrow(pho_tab))
  # print(length(unique(pho_tab$Gene)))
  print(length(unique(pro_tab$Gene)))
  # pho_head <- data.frame(Gene = pho_tab$Gene, Phosphosite = pho_tab$Phosphosite)
  # pho_head_cancers <- rbind(pho_head_cancers, pho_head)
}

tmp <- str_split_fixed(string = pho_head_cancers$Phosphosite, pattern = "[STY]", 3)
pho_head_cancers$is.solo <- (tmp[,3] == "")
pho_head_cancers$p_coord <- tmp[,2]
pho_head_cancers$residue <- sapply(X = pho_head_cancers$Phosphosite, FUN = function(phosphosite) substr(x = phosphosite, start = 1, stop = 1))
pho_head_cancers <- pho_head_cancers %>%
  mutate(id = paste0(Gene, "_", Phosphosite)) %>%
  filter(is.solo == T) %>%
  unique

pho_head_cancers %>%
  head

pho_head_cancers %>%
  select(id) %>%
  unique %>%
  nrow()

pho_head_cancers %>%
  select(Gene) %>%
  unique %>%
  nrow()

write.table(pho_head_cancers, file = paste0(makeOutDir(resultD = resultD), "pho_head_cancers_nonNA", num_nonNA, "1IQR_below_Median.txt"), quote = F, row.names = F, sep = "\t")


# # input genomically mapped phosphosites -----------------------------------
# pho_head_mapped_previous <- fread(input = "./cptac2p/analysis_results/preprocess_files/tables/convert_cptac_pho2genomic_coord/pho_head_cancers_single_site_wpid_hg38.txt", data.table = F)
# pho_head_mapped_previous %>%
#   head()
# pho_head_mapped_previous <- pho_head_mapped_previous %>%
#   mutate(coordinate_hg38 = paste0("chr", seq_region_name, ":", start, "-", end))
# pho_head_mapped_previous %>%
#   head()
# 
# pho_head_mapped_previous2merge <- pho_head_mapped_previous %>%
#   select(id, coordinate_hg38) %>%
#   unique()
# # get phosphosites with unique genomic coordination -----------------------
# pho_head_cancers <- merge(pho_head_cancers, pho_head_mapped_previous2merge, all.x = T)
# 
# pho_head_cancers %>%
#   select(id) %>%
#   unique %>%
#   nrow()
# 
# pho_head_cancers_wcoord <- pho_head_cancers %>%
#   filter(coordinate_hg38 != "chrNA:NA-NA")
# 
# pho_head_cancers_wcoord %>%
#   select(coordinate_hg38) %>%
#   unique() %>%
#   nrow()

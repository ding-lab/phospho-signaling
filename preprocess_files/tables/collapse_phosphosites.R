# Yige Wu @ WashU Dec 2017 
# collapse the phosphosite-level phosphorylation data to protein-level phosphorylation


# source ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/preprocess_files/preprocess_files_shared.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')

# set variables -------------------------------------------------------------------
cancers2process <- c("UCEC")

# input -------------------------------------------------------------------
clinical <- loadSampMap()

# loop around 3 cancers ---------------------------------------------------
for (cancer in cancers2process) { 
  sink(paste0(cptac_sharedD, cancer, "/", 'collapsed_PHO_preprocessing_log.txt'), append=FALSE, split=FALSE)
  
  ## process proteome
  print(paste0('-Start processing ', cancer, ' phosphosite data'))
  
  path_tmp <- paste0(cptac_sharedD, cancer, "/", prefix[cancer], "_PHO_formatted_normalized_replicate_averaged.txt")
  pho_raw <- fread(path_tmp, data.table = F)
  
  pho_raw.f <- pho_raw[,!(colnames(pho_raw) %in% c("Gene", "Phosphosite"))]
  rownames(pho_raw.f) <- pho_raw$Phosphosite
  
  ## log the number of row phosphosites and phosphoprotein
  print(paste0('-- ', length(unique(pho_raw$Phosphosite)), ' phosphosites'))
  
  pho_graw = collapseRows(pho_raw.f, rowGroup=pho_raw$Gene, rowID=pho_raw$Phosphosite, connectivityBasedCollapsing = T)$datETcollapsed
  pho_graw <- cbind(rownames(pho_graw), pho_graw)
  colnames(pho_graw)[1] <- "Gene"
  pho_graw <- as.data.frame(pho_graw)
  
  path_tmp <- paste0(cptac_sharedD, cancer,"/",prefix[cancer], "_collapsed_PHO_formatted_normalized_replicate_averaged.txt")
  write.table(pho_graw, row.names = F, quote=F, sep = '\t', file = path_tmp)
  
  ## log the number of collapsed phosphoprotein
  
  
  phog_t <- get_tumor(expression = pho_graw,clinical.m = clinical)
  tumorSampIDs <- colnames(phog_t)
  phog_t <- pho_graw[, c('Gene', tumorSampIDs)]
  path_tmp <- paste0(cptac_sharedD, cancer,"/",prefix[cancer], "_collapsed_PHO_formatted_normalized_replicate_averaged_Tumor.txt")
  write.table(phog_t, row.names = F, quote=F, sep = '\t', file = path_tmp)

  phog_n <- get_normal(expression = pho_graw,clinical.m = clinical)
  normalSampIDs <- colnames(phog_n)
  phog_n <- pho_graw[, c('Gene', normalSampIDs)]
  path_tmp <- paste0(cptac_sharedD, cancer,"/",prefix[cancer], "_collapsed_PHO_formatted_normalized_replicate_averaged_Normal.txt")
  write.table(phog_n, row.names = F, quote=F, sep = '\t', file = path_tmp)  
  
  sink()
  closeAllConnections()
}

# for (cancer in c("UCEC")) { 
#   path_tmp <- paste0(cptac_sharedD, cancer, "/", prefix[cancer], "_PHO_formatted_normalized_replicate_averaged.txt")
#   pho_raw <- fread(path_tmp, data.table = F)
#   
#   pho_raw.f <- pho_raw[,!(colnames(pho_raw) %in% c("Gene", "Phosphosite"))]
#   rownames(pho_raw.f) <- pho_raw$Phosphosite
#   
#   pho_graw = collapseRows(pho_raw.f, rowGroup=pho_raw$Gene, rowID=pho_raw$Phosphosite, connectivityBasedCollapsing = T)$datETcollapsed
#   pho_graw <- cbind(rownames(pho_graw), pho_graw)
#   colnames(pho_graw)[1] <- "Gene"
#   pho_graw <- as.data.frame(pho_graw)
#   
#   path_tmp <- paste0(cptac_sharedD, cancer,"/",prefix[cancer], "_collapsed_PHO_formatted_normalized_replicate_averaged.txt")
#   write.table(pho_graw, row.names = F, quote=F, sep = '\t', file = path_tmp)
#   
#   phog_t <- get_tumor(expression = pho_graw,clinical.m = clinical)
#   tumorSampIDs <- colnames(phog_t)
#   phog_t <- pho_graw[, c('Gene', tumorSampIDs)]
#   path_tmp <- paste0(cptac_sharedD, cancer,"/",prefix[cancer], "_collapsed_PHO_formatted_normalized_replicate_averaged_Tumor.txt")
#   write.table(phog_t, row.names = F, quote=F, sep = '\t', file = path_tmp)
#   
#   phog_n <- get_normal(expression = pho_graw,clinical.m = clinical)
#   normalSampIDs <- colnames(phog_n)
#   phog_n <- pho_graw[, c('Gene', normalSampIDs)]
#   path_tmp <- paste0(cptac_sharedD, cancer,"/",prefix[cancer], "_collapsed_PHO_formatted_normalized_replicate_averaged_Normal.txt")
#   write.table(phog_n, row.names = F, quote=F, sep = '\t', file = path_tmp)  
# 
# }

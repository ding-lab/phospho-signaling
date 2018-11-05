# Yige Wu @ WashU 2018 Aug
## generate a matrix of gene*sample for focal_data_by_genes.txt from GISTIC2.0

# source ------------------------------------------------------------------
source("./cptac2p_analysis/phospho_network/phospho_network_shared.R")

# Set variables -----------------------------------------------------------
num_genoalt_thres <- 4
amp_thres <- 0.1
del_thres <- -0.1

# inputs ------------------------------------------------------------------
## input enzyme-substrate table
ptms_site_pairs_sup <- read_csv(paste0(ppnD, "compile_enzyme_substrate/tables/compile_omnipath_networkin_depod_signor_manual/Omnipath_NetworKin_DEPOD_SignorNotSiteMapped_manualAdded.csv"))
ptms_site_pairs_sup <- ptms_site_pairs_sup[ptms_site_pairs_sup$Source != "NetKIN" | (ptms_site_pairs_sup$Source == "NetKIN" & ptms_site_pairs_sup$networkin_score >= 5),]
ptms_site_pairs_sup <- data.frame(ptms_site_pairs_sup)
kinases <- unique(ptms_site_pairs_sup$GENE[ptms_site_pairs_sup$enzyme_type == "kinase"])
phosphatases <- unique(ptms_site_pairs_sup$GENE[ptms_site_pairs_sup$enzyme_type == "phosphatase"])
substrates <- unique(ptms_site_pairs_sup$SUB_GENE)

# get the gene list covering all the enzymes and substrates ---------------
genes4mat <- unique(c(kinases, phosphatases, substrates))
length(genes4mat)

# loop by cancer ----------------------------------------------------------
for (cancer in cancers_sort) {
  ## input focal CNA values
  ### if it is breast cancer data, irregular columns names actually won't overlap with proteomics data
  cna <- fread(input = paste0(cptac2pD, "copy_number/gistic/outputs/cptac2p.", tolower(substr(x = cancer, start = 1, stop = 2)), "/somatic/run_gistic2/0.25/0.99/armpeel1/somatic_clean/focal_data_by_genes.txt"),
               data.table = F)
  nrow(cna)
  cna <- cna[cna$`Gene Symbol` %in% genes4mat,]
  nrow(cna)
  partIDs <- colnames(cna)
  partIDs <- partIDs[!(partIDs %in% c("Gene Symbol","Gene ID","Cytoband"))]
  
  ### input by-case threshold for high-level CNA
  cna_thres <- fread(input = paste0(cptac2pD, "copy_number/gistic/outputs/cptac2p.", tolower(substr(x = cancer, start = 1, stop = 2)), "/somatic/run_gistic2/0.25/0.99/armpeel1/somatic_clean/sample_cutoffs.txt"),
                     data.table = F)
  deep_del_thres_mat <- matrix(data = rep(cna_thres$Low, nrow(cna)), nrow = nrow(cna), byrow = T)
  deep_amp_thres_mat <- matrix(data = rep(cna_thres$High, nrow(cna)), nrow = nrow(cna), byrow = T)
  
  if (any(cna_thres$Low > del_thres) | any(cna_thres$High < amp_thres)) {
    stop("meow")
  }
  
  ## generate low-level CNA thresholded matrix
  shallow_del_thres_mat <- matrix(data = rep(del_thres, nrow(cna)*length(partIDs)), nrow = nrow(cna), byrow = T)
  shallow_amp_thres_mat <- matrix(data = rep(amp_thres,  nrow(cna)*length(partIDs)), nrow = nrow(cna), byrow = T)
  
  shallow_del_thresholded_mat <- (cna[, partIDs] <= shallow_del_thres_mat & cna[, partIDs] > deep_del_thres_mat)
  shallow_amp_thresholded_mat <- (cna[, partIDs] >= shallow_amp_thres_mat & cna[, partIDs] < deep_amp_thres_mat)

  Hugo_Symbol <- data.frame(Hugo_Symbol = as.vector(cna$`Gene Symbol`))
  ## generate high-level CNA thresholded matrix
  deep_del_thresholded_mat <- (cna[, partIDs] <= deep_del_thres_mat)
  deep_amp_thresholded_mat <- (cna[, partIDs] >= deep_amp_thres_mat)
  
  fn <- paste0(makeOutDir(resultD = resultD), cancer, "_deep_deletions_in_enzyme_substrate.txt")
  write.table(x = cbind(Hugo_Symbol, deep_del_thresholded_mat), file = fn, row.names = F, quote = F, sep = '\t')
  
  fn <- paste0(makeOutDir(resultD = resultD), cancer, "_deep_amplifications_in_enzyme_substrate.txt")
  write.table(x = cbind(Hugo_Symbol, deep_amp_thresholded_mat), file = fn, row.names = F, quote = F, sep = '\t')
  
  fn <- paste0(makeOutDir(resultD = resultD), cancer, "_shallow_deletions_", del_thres, "_in_enzyme_substrate.txt")
  write.table(x = cbind(Hugo_Symbol, shallow_del_thresholded_mat), file = fn, row.names = F, quote = F, sep = '\t')
  
  fn <- paste0(makeOutDir(resultD = resultD), cancer, "_shallow_amplifications_", amp_thres, "_in_enzyme_substrate.txt")
  write.table(x = cbind(Hugo_Symbol, shallow_amp_thresholded_mat), file = fn, row.names = F, quote = F, sep = '\t')
  
  shallow_del_count <- rowSums(shallow_del_thresholded_mat); names(shallow_del_count) <- cna$`Gene Symbol`
  shallow_amp_count <- rowSums(shallow_amp_thresholded_mat); names(shallow_amp_count) <- cna$`Gene Symbol`
  deep_del_count <- rowSums(deep_del_thresholded_mat); names(deep_del_count) <- cna$`Gene Symbol`
  deep_amp_count <- rowSums(deep_amp_thresholded_mat); names(deep_amp_count) <- cna$`Gene Symbol`
  
  ## number of genes with at least 4 deep deletion
  deep_del_count_thresholded <- deep_del_count[deep_del_count >= num_genoalt_thres]
  deep_amp_count_thresholded <- deep_amp_count[deep_amp_count >= num_genoalt_thres]
  
  print(length(deep_del_count_thresholded))
  print(length(deep_amp_count_thresholded))
  
  ## number of kinases genes with at least 4 deep deletions/amplification
  kinases_deep_amp_count <- deep_amp_count_thresholded[names(deep_amp_count_thresholded) %in% kinases]
  print(kinases_deep_amp_count)
  kinases_deep_del_count <- deep_del_count_thresholded[names(deep_del_count_thresholded) %in% kinases]
  print(kinases_deep_del_count)
  
  ## number of phosphatases genes with at least 4 deep deletions/amplification
  phosphatases_deep_amp_count <- deep_amp_count_thresholded[names(deep_amp_count_thresholded) %in% phosphatases]
  print(phosphatases_deep_amp_count)
  phosphatases_deep_del_count <- deep_del_count_thresholded[names(deep_del_count_thresholded) %in% phosphatases]
  print(phosphatases_deep_del_count)
}


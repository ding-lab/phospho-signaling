# Yige Wu @ WashU 2018 Aug
## generate a matrix of gene*sample for focal_data_by_genes.txt from GISTIC2.0

# source ------------------------------------------------------------------
source("Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R")

# Set variables -----------------------------------------------------------
num_genoalt_thres <- 4
# amp_thres <- 0.1
# del_thres <- -0.1
amp_thres <- log2(1.1)
del_thres <- log2(0.9)

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
  ## input CNA values
  ### if it is breast cancer data, irregular columns names actually won't overlap with proteomics data
  cna <- fread(input = paste0(cptac2p_genomicD, "copy_number/gatk/v1.3.swap_contamination_fixed/prospective_somatic/gene_level/v1.3.CPTAC2_prospective.2018-03-19/", toupper(substr(cancer, start = 1, stop = 2)), "/gene_level_CNV.", substr(cancer, start = 1, stop = 2), ".v1.3.2018-03-19.tsv"), data.table = F)
  nrow(cna)
  cna <- cna[cna$gene %in% genes4mat,]
  nrow(cna)
  partIDs <- colnames(cna)
  partIDs <- partIDs[!(partIDs %in% c("gene","Gene ID","Cytoband"))]
  
  ## generate low-level CNA thresholded matrix
  shallow_del_thres_mat <- matrix(data = rep(del_thres, nrow(cna)*length(partIDs)), nrow = nrow(cna), byrow = T)
  shallow_amp_thres_mat <- matrix(data = rep(amp_thres,  nrow(cna)*length(partIDs)), nrow = nrow(cna), byrow = T)
  
  shallow_del_thresholded_mat <- (cna[, partIDs] <= shallow_del_thres_mat)
  shallow_amp_thresholded_mat <- (cna[, partIDs] >= shallow_amp_thres_mat)

  Hugo_Symbol <- data.frame(Hugo_Symbol = as.vector(cna$`gene`))

  fn <- paste0(makeOutDir(resultD = resultD), cancer, "_shallow_deletions_", del_thres, "_in_enzyme_substrate.txt")
  write.table(x = cbind(Hugo_Symbol, shallow_del_thresholded_mat), file = fn, row.names = F, quote = F, sep = '\t')
  
  fn <- paste0(makeOutDir(resultD = resultD), cancer, "_shallow_amplifications_", amp_thres, "_in_enzyme_substrate.txt")
  write.table(x = cbind(Hugo_Symbol, shallow_amp_thresholded_mat), file = fn, row.names = F, quote = F, sep = '\t')
  
  shallow_del_count <- rowSums(shallow_del_thresholded_mat); names(shallow_del_count) <- cna$`gene`
  shallow_amp_count <- rowSums(shallow_amp_thresholded_mat); names(shallow_amp_count) <- cna$`gene`
}


# Yige Wu @ WashU 2018 Aug
## generate a matrix of gene*sample for focal_data_by_genes.txt from GISTIC2.0

# source ------------------------------------------------------------------
source("Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R")

# Set variables -----------------------------------------------------------
num_genoalt_thres <- 4
# amp_thres <- 0.1
# del_thres <- -0.1
amp_thres_cans <- c(log2(1.1), log2(1.1), log2(1.1), 0.1, 0.1); names(amp_thres_cans) <- c(cancers_sort, "UCEC", "CCRCC")
del_thres_cans <- c(log2(0.9), log2(0.9), log2(0.9), -0.1, -0.1); names(del_thres_cans) <- c(cancers_sort, "UCEC", "CCRCC")

# inputs ------------------------------------------------------------------
## input enzyme-substrate table
ptms_site_pairs_sup <- read_csv(paste0(ppnD, "compile_enzyme_substrate/tables/compile_omnipath_networkin_depod_signor_manual/Omnipath_NetworKin_DEPOD_SignorNotSiteMapped_manualAdded.csv"))
ptms_site_pairs_sup <- ptms_site_pairs_sup[ptms_site_pairs_sup$Source != "NetKIN" | (ptms_site_pairs_sup$Source == "NetKIN" & ptms_site_pairs_sup$networkin_score >= 5),]
ptms_site_pairs_sup <- data.frame(ptms_site_pairs_sup)
kinases <- unique(ptms_site_pairs_sup$GENE[ptms_site_pairs_sup$enzyme_type == "kinase"])
phosphatases <- unique(ptms_site_pairs_sup$GENE[ptms_site_pairs_sup$enzyme_type == "phosphatase"])
substrates <- unique(ptms_site_pairs_sup$SUB_GENE)

## input complex table
complex_pair_tab <- fread(input = "./cptac2p/analysis_results/phospho_network/compile_enzyme_substrate/tables/parse_corum_signor_reactome/sup_complex_pair_uniq.txt", data.table = F)

# parse UCEC Baylor copy number -------------------------------------------
## input UCEC meta data
meta_tab <- fread("./Ding_Lab/Projects_Current/CPTAC/PGDAC/Endometrium_CPTAC3/01_Data_tables/UCEC_V1.1/UCEC_CPTAC3_meta_table.txt", data.table = F)
rownames(meta_tab) <- meta_tab$idx

# get the gene list covering all the enzymes and substrates ---------------
genes4mat <- unique(c(kinases, phosphatases, substrates, as.vector(complex_pair_tab$geneA), as.vector(complex_pair_tab$geneB)))
length(genes4mat)

# process prospective cancers ----------------------------------------------------------
for (cancer in c("UCEC", "CCRCC", cancers_sort)) {
  ## get CNA thresholds
  amp_thres <- amp_thres_cans[cancer]
  del_thres <- del_thres_cans[cancer]
  
  ## input CNA values
  if (cancer %in% cancers_sort) {
    ### if it is breast cancer data, irregular columns names actually won't overlap with proteomics data
    cna <- fread(input = paste0(cptac2p_genomicD, "copy_number/gatk/v1.3.swap_contamination_fixed/prospective_somatic/gene_level/v1.3.CPTAC2_prospective.2018-03-19/", toupper(substr(cancer, start = 1, stop = 2)), "/gene_level_CNV.", substr(cancer, start = 1, stop = 2), ".v1.3.2018-03-19.tsv"), data.table = F)
    partIDs <- colnames(cna)
    partIDs <- partIDs[!(partIDs %in% c("gene","Gene ID","Cytoband"))]
  }
  if (cancer %in% c("UCEC", "CCRCC")) {
    cna <- fread(input = paste0("./cptac2p/analysis_results/preprocess_files/tables/parse_", cancer, "_data_freeze/somatic_CNA.", cancer, ".partID.txt"), data.table = F)
    partIDs <- colnames(cna)[-1]
  }
  cna <- cna[cna$gene %in% genes4mat,]

  ## generate low-level CNA thresholded matrix
  shallow_del_thres_mat <- matrix(data = rep(del_thres, nrow(cna)*length(partIDs)), nrow = nrow(cna), byrow = T)
  shallow_amp_thres_mat <- matrix(data = rep(amp_thres,  nrow(cna)*length(partIDs)), nrow = nrow(cna), byrow = T)
  
  shallow_del_thresholded_mat <- (cna[, partIDs] <= shallow_del_thres_mat)
  shallow_amp_thresholded_mat <- (cna[, partIDs] >= shallow_amp_thres_mat)

  Hugo_Symbol <- data.frame(Hugo_Symbol = as.vector(cna$`gene`))

  tab2w <- cbind(Hugo_Symbol, shallow_del_thresholded_mat)
  fn <- paste0(makeOutDir(resultD = resultD), cancer, "_shallow_deletions_", del_thres, ".txt")
  write.table(x = tab2w, file = fn, row.names = F, quote = F, sep = '\t')
  
  tab2w <- cbind(Hugo_Symbol, shallow_amp_thresholded_mat)
  fn <- paste0(makeOutDir(resultD = resultD), cancer, "_shallow_amplifications_", amp_thres, ".txt")
  write.table(x = tab2w, file = fn, row.names = F, quote = F, sep = '\t')
}


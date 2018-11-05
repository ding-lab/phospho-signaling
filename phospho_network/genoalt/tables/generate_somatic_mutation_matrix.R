# Yige Wu @ WashU 2018 Aug
## generate a matrix of gene*sample for somatic mutations


# source ------------------------------------------------------------------
source("Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R")


# Set variables -----------------------------------------------------------
num_genoalt_thres <- 4

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

# get the gene list covering all the enzymes and substrates ---------------
genes4mat <- unique(c(kinases, phosphatases, substrates, as.vector(complex_pair_tab$geneA), as.vector(complex_pair_tab$geneB)))
length(genes4mat)

# loop by cancer ----------------------------------------------------------
for (cancer in c(cancers_sort, "UCEC")) {
  maf <- loadMaf(cancer = cancer, maf_files = maf_files)
  nrow(maf)
  maf <- maf[maf$Hugo_Symbol %in% genes4mat,]
  nrow(maf)
  maf$sampID <- str_split_fixed(string = maf$Tumor_Sample_Barcode, pattern = "_", 2)[,1]
  mut_mat <- dcast(data = maf, Hugo_Symbol ~ sampID, fun =  function(x) {
    variant_class <- paste0(unique(x), collapse = ",")
    return(variant_class)
  },value.var = "Variant_Classification", drop=TRUE)
  rownames(mut_mat) <- as.vector(mut_mat$Hugo_Symbol)
  print(nrow(mut_mat))
  
  mut_count <- rowSums(x = ((mut_mat[, -1] != "") & (mut_mat[, -1] != "Silent")))
  mut_count_thres <- mut_count[mut_count >= num_genoalt_thres]
  
  ## number of genes with at least 4 nonsynonymous mutations
  print(length(mut_count_thres))
  
  ## number of kinases genes with at least 4 nonsynonymous mutations
  kinases_mut_count <- mut_count_thres[names(mut_count_thres) %in% kinases]
  print(kinases_mut_count)
  
  ## number of phosphatases genes with at least 4 nonsynonymous mutations
  phosphatases_mut_count <- mut_count_thres[names(mut_count_thres) %in% phosphatases]
  print(phosphatases_mut_count)
  
  fn <- paste0(makeOutDir(resultD = resultD), cancer, "_somatic_mutation_in_enzyme_substrate.txt")
  write.table(x = mut_mat, file = fn, row.names = F, quote = F, sep = '\t')
}


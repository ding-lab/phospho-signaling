# Yige Wu @ WashU 2018 Aug
## generate a matrix of gene*sample for somatic mutations


# source ------------------------------------------------------------------
source("./cptac2p_analysis/phospho_network/phospho_network_shared.R")

# Set variables -----------------------------------------------------------
num_genoalt_thres <- 4

# inputs ------------------------------------------------------------------
complex_pair_tab <- fread(input = "./cptac2p/analysis_results/phospho_network/compile_enzyme_substrate/tables/parse_corum_signor_reactome/sup_complex_pair_uniq.txt", data.table = F)

# get the gene list covering all the enzymes and substrates ---------------
genes4mat <- unique(c(as.vector(complex_pair_tab$geneA), as.vector(complex_pair_tab$geneB)))
length(genes4mat)

# loop by cancer ----------------------------------------------------------
for (cancer in "UCEC") {
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
  
  fn <- paste0(makeOutDir(resultD = resultD), cancer, "_somatic_mutation_in_complex.txt")
  write.table(x = mut_mat, file = fn, row.names = F, quote = F, sep = '\t')
}


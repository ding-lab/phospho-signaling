# Yige Wu @ WashU 2018 Aug
## generate a matrix of gene*sample for somatic mutations

# get the gene list covering all the enzymes and substrates ---------------
genes4mat <- unique(unlist(pair_tab))
length(genes4mat)

# loop by cancer ----------------------------------------------------------
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

fn <- paste0(makeOutDir(resultD = resultD), cancer, "_somatic_mutation.txt")
write.table(x = mut_mat, file = fn, row.names = F, quote = F, sep = '\t')



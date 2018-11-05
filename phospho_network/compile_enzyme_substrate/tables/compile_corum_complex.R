# Yige Wu @ WashU 2018 Apr
# compile complex info from Signor


# source ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
library(biomaRt)

# set variables -----------------------------------------------------------

# inputs ------------------------------------------------------------------
## input ensembl
ensembl = useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl",host="www.ensembl.org")
ensembl_attr <- listAttributes(ensembl)
ensembl_filter <- listFilters(ensembl)

## input signor complex information
corum_complexes <- fread(input = paste0("./Ding_Lab/Projects_Current/PPI/PPI_shared_data/protein_complex/CPRUM/allComplexes_CORUM.txt"), data.table = F)
corum_complexes <- corum_complexes[corum_complexes$Organism == "Human",]

# extract the list -----------------------------
corum_complex_tab <- corum_complexes[, c("ComplexID", "subunits(Gene name)")]

## extract the list
corum_complex <- lapply(X = corum_complex_tab$ComplexID, FUN = function(id, tab) unlist(strsplit(x =  as.vector(tab$`subunits(Gene name)`[tab$`ComplexID` == id]), split = "\\;")), tab = corum_complex_tab)
names(corum_complex) <- as.vector(corum_complex_tab$`ComplexID`)

corum_complex_geneSymbols <- corum_complex

# compile a data frame get protein pairs within a complex----------------------------------------------------
## get the number of possible pairs
num_pairs <- sum(sapply(X = corum_complex_geneSymbols, function(c) length(c)*(length(c)-1)))
geneA = vector(mode = "character", length = num_pairs)
geneB = vector(mode = "character", length = num_pairs)
corum_ID = vector(mode = "character", length = num_pairs)

count <- 0
for (i in 1:length(corum_complex_geneSymbols)) {
  corum_id <- names(corum_complex_geneSymbols)[i]
  for (genea in corum_complex_geneSymbols[[corum_id]]) {
    for (geneb in corum_complex_geneSymbols[[corum_id]][corum_complex_geneSymbols[[corum_id]] != genea]) {
      count <- count + 1
      geneA[count] <- genea
      geneB[count] <- geneb
      corum_ID[count] <- corum_id
    }
  }
}

corum_complex_pair_tab <- data.frame(geneA = geneA, 
                                      geneB = geneB,
                                      corum_ID = corum_ID)

write.table(x = corum_complex_pair_tab, file = paste0(makeOutDir(resultD = resultD), "corum_complex_pair_tab.txt"), quote = F, sep = "\t", row.names = F)


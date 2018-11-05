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
signor_complexes <- fread(input = paste0(pan3can_shared_dataD, "Phospho_databases/Signor/SIGNOR_complexes.csv"), data.table = F)

# combine the complex and protein family list -----------------------------
## combine the data frame first
signor_complex_tab <- signor_complexes[, c("SIGNOR ID", "LIST OF ENTITIES")]

## extract the list
signor_complex <- lapply(X = signor_complex_tab$`SIGNOR ID`, FUN = function(id, tab) unlist(strsplit(x =  as.vector(tab$`LIST OF ENTITIES`[tab$`SIGNOR ID` == id]), split = ",| ")), tab = signor_complex_tab)
names(signor_complex) <- as.vector(signor_complex_tab$`SIGNOR ID`)

## clean out the empty strings and small molecules
signor_complex <- lapply(X = signor_complex, FUN = function(listofentity) unique(listofentity[listofentity != "" & !grepl(pattern = "CID", x = listofentity)]))

## expand the embedded the family and complex
allentities <- unlist(signor_complex)
while (any(grepl(pattern = "signor", x = allentities))) {
  for (i in names(signor_complex)) {
    listofentity <- signor_complex[[i]]
    signor_ids <- listofentity[grepl(pattern = "signor", x = listofentity, ignore.case = T)]
    if (length(signor_ids) > 0){
      for (signor_id in signor_ids) {
        listofentity <- listofentity[listofentity != signor_id]
        listofentity <- unique(c(listofentity, signor_complex[[signor_id]]))
        signor_complex[[i]] <- listofentity
      }
    }
  }
  allentities <- unlist(signor_complex)
}

## map above list with gene names
allentities_uniprots <- unique(allentities)
length(allentities_uniprots)
mapTab = getBM(attributes = c("hgnc_symbol", "uniprotswissprot"), filters = "uniprotswissprot", 
               values = allentities_uniprots, mart = ensembl, uniqueRows=FALSE)
mapTab <- unique(mapTab)
nrow(mapTab)
allentities_uniprots_unmapped <- allentities_uniprots[!(allentities_uniprots %in% mapTab$uniprotswissprot)]
allentities_uniprots_unmapped <- data.frame(hgnc_symbol = c("GUCY1B2", "TRBC1", "TRAC"),
                                            uniprotswissprot = c("O75343", "P01850", "P01848"))
## unmapped could be a product of a pseudogene
mapTab <- rbind(mapTab, allentities_uniprots_unmapped)
nrow(mapTab)
mapTab[(mapTab$uniprotswissprot %in% mapTab$uniprotswissprot[duplicated(mapTab$uniprotswissprot)]),]

## get uniprots named by gene name
allentities_uniprots <- as.vector(mapTab$uniprotswissprot)
names(allentities_uniprots) <- as.vector(mapTab$hgnc_symbol)

## add all gene symbols into a new list, even the multiple genes to one uniprot ID
signor_complex_geneSymbols <- sapply(X = signor_complex,
                                        FUN = function(listofentity, tab) unique(as.vector(tab$hgnc_symbol[tab$uniprotswissprot %in% listofentity])), tab = mapTab)



# compile a data frame get protein pairs within a complex----------------------------------------------------
## get the number of possible pairs
num_pairs <- sum(sapply(X = signor_complex_geneSymbols, function(c) length(c)*(length(c)-1)))
geneA = vector(mode = "character", length = num_pairs)
geneB = vector(mode = "character", length = num_pairs)
SIGNOR_ID = vector(mode = "character", length = num_pairs)

count <- 0
for (i in 1:length(signor_complex_geneSymbols)) {
  signor_id <- names(signor_complex_geneSymbols)[i]
  for (genea in signor_complex_geneSymbols[[signor_id]]) {
    for (geneb in signor_complex_geneSymbols[[signor_id]][signor_complex_geneSymbols[[signor_id]] != genea]) {
      count <- count + 1
      geneA[count] <- genea
      geneB[count] <- geneb
      SIGNOR_ID[count] <- signor_id
    }
  }
}

signor_complex_pair_tab <- data.frame(geneA = geneA, 
                                      geneB = geneB,
                                      SIGNOR_ID = SIGNOR_ID)
write.table(x = signor_complex_pair_tab, file = paste0(makeOutDir(resultD = resultD), "signor_complex_pair_tab.txt"), quote = F, sep = "\t", row.names = F)

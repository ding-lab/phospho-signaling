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
reactome_complexes <- fread(input = paste0("./Ding_Lab/Projects_Current/PPI/PPI_shared_data/protein_complex/Reactome/ComplexParticipantsPubMedIdentifiers_Human_Reactome.txt"), data.table = F)

# combine the complex and protein family list -----------------------------
## combine the data frame first
reactome_complex_tab <- reactome_complexes[, c("identifier", "participants")]

## extract the list
reactome_complex <- lapply(X = reactome_complex_tab$`identifier`, FUN = function(id, tab) unlist(strsplit(x =  as.vector(tab$`participants`[tab$`identifier` == id]), split = "\\|")), tab = reactome_complex_tab)
names(reactome_complex) <- as.vector(reactome_complex_tab$`identifier`)

## clean out the empty strings and small molecules
reactome_complex <- lapply(X = reactome_complex, FUN = function(listofentity) unique(listofentity[listofentity != "" & !grepl(pattern = "CID", x = listofentity)]))

allentities <- unlist(reactome_complex)
while (any(!grepl(pattern = "uniprot", x = allentities))) {
  for (i in names(reactome_complex)) {
    listofentity <- reactome_complex[[i]]
    if (!any(!grepl(pattern = "uniprot", x = allentities)) & length(allentities) > 1) {
      reactome_ids <- listofentity[grepl(pattern = "uniprot", x = listofentity, ignore.case = T)]
      if (length(reactome_ids) >= 2){
        for (reactome_id in reactome_ids) {
          listofentity <- listofentity[listofentity != reactome_id]
          listofentity <- unique(c(listofentity, reactome_complex[[reactome_id]]))
          reactome_complex[[i]] <- listofentity
        }
      }
    } else {
      reactome_complex[[i]] <- NULL
    }
  }
  allentities <- unlist(reactome_complex)
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
reactome_complex_geneSymbols <- sapply(X = reactome_complex,
                                        FUN = function(listofentity, tab) unique(as.vector(tab$hgnc_symbol[tab$uniprotswissprot %in% listofentity])), tab = mapTab)



# compile a data frame get protein pairs within a complex----------------------------------------------------
## get the number of possible pairs
num_pairs <- sum(sapply(X = reactome_complex_geneSymbols, function(c) length(c)*(length(c)-1)))
geneA = vector(mode = "character", length = num_pairs)
geneB = vector(mode = "character", length = num_pairs)
reactome_ID = vector(mode = "character", length = num_pairs)

count <- 0
for (i in 1:length(reactome_complex_geneSymbols)) {
  reactome_id <- names(reactome_complex_geneSymbols)[i]
  for (genea in reactome_complex_geneSymbols[[reactome_id]]) {
    for (geneb in reactome_complex_geneSymbols[[reactome_id]][reactome_complex_geneSymbols[[reactome_id]] != genea]) {
      count <- count + 1
      geneA[count] <- genea
      geneB[count] <- geneb
      reactome_ID[count] <- reactome_id
    }
  }
}

reactome_complex_pair_tab <- data.frame(geneA = geneA, 
                                      geneB = geneB,
                                      reactome_ID = reactome_ID)

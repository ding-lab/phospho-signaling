# Yige Wu @ WashU 2018 Apr
# compile protein-level and site-level ks pairs from Omnipath + networkin + DEPOD


# source ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
library(biomaRt)

# set variables -----------------------------------------------------------
entitytypes2process <- c("protein", "proteinfamily", "complex")

# inputs ------------------------------------------------------------------
## input ensembl
ensembl = useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl",host="www.ensembl.org")
ensembl_attr <- listAttributes(ensembl)
ensembl_filter <- listFilters(ensembl)

## inputs omnipath + networkin + depod
ptms_pairs <- fread(input = paste0(ppnD, "compile_enzyme_substrate/tables/compile_omnipath_networkin_depod/omnipath_networkin_enzyme_substrate_site_level_union_source_cleaned_addDEPOD_extended.csv"),
                    data.table = F)
## input signor
signor <- fread(input = paste0(pan3can_shared_dataD, "Phospho_databases/Signor/human_phosphorylations_18_08_18.tsv"), data.table = F)

## input signor complex information
signor_complexes <- fread(input = paste0(pan3can_shared_dataD, "Phospho_databases/Signor/SIGNOR_complexes.csv"), data.table = F)

## input signor protein family information
signor_pf <- fread(input = paste0(pan3can_shared_dataD, "Phospho_databases/Signor/SIGNOR_PF.csv"), data.table = F)

# combine the complex and protein family list -----------------------------
## combine the data frame first
signor_complex_pf_tab <- rbind(signor_complexes[, c("SIGNOR ID", "LIST OF ENTITIES")], signor_pf[, c("SIGNOR ID", "LIST OF ENTITIES")])

## extract the list
signor_complex_pf <- lapply(X = signor_complex_pf_tab$`SIGNOR ID`, FUN = function(id, tab) unlist(strsplit(x =  as.vector(tab$`LIST OF ENTITIES`[tab$`SIGNOR ID` == id]), split = ",| ")), tab = signor_complex_pf_tab)
names(signor_complex_pf) <- as.vector(signor_complex_pf_tab$`SIGNOR ID`)

## clean out the empty strings and small molecules
signor_complex_pf <- lapply(X = signor_complex_pf, FUN = function(listofentity) unique(listofentity[listofentity != "" & !grepl(pattern = "CID", x = listofentity)]))

## expand the embedded the family and complex
allentities <- unlist(signor_complex_pf)
while (any(grepl(pattern = "signor", x = allentities))) {
  for (i in names(signor_complex_pf)) {
    listofentity <- signor_complex_pf[[i]]
    signor_ids <- listofentity[grepl(pattern = "signor", x = listofentity, ignore.case = T)]
    if (length(signor_ids) > 0){
      for (signor_id in signor_ids) {
        listofentity <- listofentity[listofentity != signor_id]
        listofentity <- unique(c(listofentity, signor_complex_pf[[signor_id]]))
        signor_complex_pf[[i]] <- listofentity
      }
    }
  }
  allentities <- unlist(signor_complex_pf)
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
signor_complex_pf_geneSymbols <- sapply(X = signor_complex_pf,
                                        FUN = function(listofentity, tab) unique(as.vector(tab$hgnc_symbol[tab$uniprotswissprot %in% listofentity])), tab = mapTab)

# process protein complex --------------------------------------
## process entityA with protein complexes
signor_acomplex <- signor[signor$TYPEA == "complex",]
signor_anotcomplex <- signor[signor$TYPEA != "complex",]

## fill in the members of the complex into entity A
signor_acomplex_new <- NULL
for (i in 1:nrow(signor_acomplex)) {
  ## get each line
  signor_tmp <- signor_acomplex[i,]
  
  ## get signor id
  id_tmp <- as.character(signor_tmp$IDA)
  
  ## get the vector of complex members
  entitiesA <- signor_complex_pf_geneSymbols[[id_tmp]]
  
  ## create the new table with the line
  signor_tmp_new <- signor_tmp[rep(c(1), length(entitiesA)),]
  signor_tmp_new$ENTITYA <- entitiesA
  colnames_signor <- colnames(signor_tmp_new)
  signor_tmp_new$IDA <- allentities_uniprots[entitiesA]
  
  ## combine to the cumulating new table
  signor_acomplex_new <- rbind(signor_acomplex_new, signor_tmp_new)
}

signor <- rbind(signor_acomplex_new, signor_anotcomplex)
## process entityB with protein complexes
signor_bcomplex <- signor[signor$TYPEB == "complex",]
signor_bnotcomplex <- signor[signor$TYPEB != "complex",]

## fill in the members of the complex into entity B
signor_bcomplex_new <- NULL
for (i in 1:nrow(signor_bcomplex)) {
  ## get each line
  signor_tmp <- signor_bcomplex[i,]
  
  ## get signor id
  id_tmp <- as.character(signor_tmp$IDB)
  
  ## get the vector of complex members
  entitiesB <- signor_complex_pf_geneSymbols[[id_tmp]]
  
  ## create the new table with the line
  signor_tmp_new <- signor_tmp[rep(c(1), length(entitiesB)),]
  signor_tmp_new$ENTITYB <- entitiesB
  colnames_signor <- colnames(signor_tmp_new)
  signor_tmp_new$IDB <- allentities_uniprots[entitiesB]
  
  ## combine to the cumulating new table
  signor_bcomplex_new <- rbind(signor_bcomplex_new, signor_tmp_new)
}
signor <- rbind(signor_bcomplex_new, signor_bnotcomplex)


# process protein proteinfamily --------------------------------------
## process entityA with protein proteinfamilyes
signor_aproteinfamily <- signor[signor$TYPEA == "proteinfamily",]
signor_anotproteinfamily <- signor[signor$TYPEA != "proteinfamily",]

if (nrow(signor_aproteinfamily) > 0) {
  ## fill in the members of the proteinfamily into entity A
  signor_aproteinfamily_new <- NULL
  for (i in 1:nrow(signor_aproteinfamily)) {
    ## get each line
    signor_tmp <- signor_aproteinfamily[i,]
    
    ## get signor id
    id_tmp <- as.character(signor_tmp$IDA)
    
    ## get the vector of proteinfamily members
    entitiesA <- signor_complex_pf_geneSymbols[[id_tmp]]
    
    ## create the new table with the line
    signor_tmp_new <- signor_tmp[rep(c(1), length(entitiesA)),]
    signor_tmp_new$ENTITYA <- entitiesA
    colnames_signor <- colnames(signor_tmp_new)
    signor_tmp_new$IDA <- allentities_uniprots[entitiesA]
    
    ## combine to the cumulating new table
    signor_aproteinfamily_new <- rbind(signor_aproteinfamily_new, signor_tmp_new)
  }
  
  signor <- rbind(signor_aproteinfamily_new, signor_anotproteinfamily)
}

## process entityB with protein proteinfamilyes
signor_bproteinfamily <- signor[signor$TYPEB == "proteinfamily",]
signor_bnotproteinfamily <- signor[signor$TYPEB != "proteinfamily",]

## fill in the members of the proteinfamily into entity B
if (nrow(signor_bproteinfamily) > 0) {
  signor_bproteinfamily_new <- NULL
  for (i in 1:nrow(signor_bproteinfamily)) {
    ## get each line
    signor_tmp <- signor_bproteinfamily[i,]
    
    ## get signor id
    id_tmp <- as.character(signor_tmp$IDB)
    
    ## get the vector of proteinfamily members
    entitiesB <- signor_complex_pf_geneSymbols[[id_tmp]]
    
    ## create the new table with the line
    signor_tmp_new <- signor_tmp[rep(c(1), length(entitiesB)),]
    signor_tmp_new$ENTITYB <- entitiesB
    colnames_signor <- colnames(signor_tmp_new)
    signor_tmp_new$IDB <- allentities_uniprots[entitiesB]
    
    ## combine to the cumulating new table
    signor_bproteinfamily_new <- rbind(signor_bproteinfamily_new, signor_tmp_new)
  }
  signor <- rbind(signor_bproteinfamily_new, signor_bnotproteinfamily)
}
write.table(x = signor, file = paste0(pan3can_shared_dataD, "Phospho_databases/Signor/human_phosphorylations_18_08_18_proteinfamily_complex_expanded.txt"), quote = F, row.names = F, sep = "\t")

# combine tables processed from signor and format as previous k-s table ------------------------------------
signor_tab <- data.frame(GENE = as.vector(signor$ENTITYA), SUB_GENE = as.vector(signor$ENTITYB), SUB_MOD_RSD = "NULL", pair = "NULL",
                         KINASE = as.vector(signor$ENTITYA), KIN_ACC_ID = as.vector(signor$IDA), KIN_ORGANISM = "human",
                         SUBSTRATE = as.vector(signor$ENTITYB), SUB_GENE_ID = "NULL", SUB_ACC_ID = as.vector(signor$IDB), SUB_ORGANISM = "human", 
                         SITE_GRP_ID = "NULL", SITE_...7_AA = "NULL", Source = "Signor", networkin_score = Inf, enzyme_type = ifelse(signor$MECHANISM == "phosphorylation", "kinase", "phosphatase") )
signor_tab <- unique(signor_tab)
nrow(signor_tab)
# check the proportion that overlaps with previous table ------------------
ptms_pairs$pair_pro <- paste0(ptms_pairs$GENE, ":", ptms_pairs$SUB_GENE)
signor_tab$pair_pro <- paste0(signor_tab$GENE, ":", signor_tab$SUB_GENE)

## only comparing the protein pair part
signor_tab2add <- signor_tab[!(signor_tab$pair_pro %in% ptms_pairs$pair_pro),]
nrow(signor_tab2add)

## add protein pairs not included
ptms_pairs <- rbind(signor_tab2add, ptms_pairs)
write.table(x = ptms_pairs, file = paste0(makeOutDir(resultD = resultD), "omnipath_networkin_DEPOD_SignorNotSiteMapped.csv"), quote = F, row.names = F, sep = ",")



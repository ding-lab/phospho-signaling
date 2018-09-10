# Yige Wu @ WashU 2018 Jan
# convert site-level ks pairs from Omnipath to genomic positions
# OmniPath output only has uniprot Ids and residue offset, the strategy is to first convert to ensembl peptide id then query for genomic coordication

# source -------------------------------------------------------------------
## for convert uniprot ID to ensembl peptide ID
source("https://bioconductor.org/biocLite.R")
# biocLite("devtools")
# BiocInstaller::biocLite('grimbough/biomaRt')
library(biomaRt)

## for convert ensembl peptide ID + residue to genomic coord
library(httr)
library(jsonlite)
library(xml2)

## others
library(data.table)
library(readr)
source('/Users/yigewu/Box Sync/cptac2p_analysis/preprocess_files/preprocess_files.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')



# shared inputs -----------------------------------------------------------

# convert uniprot ID to ensembl peptide ID ------------------------------------------------------------------
ensembl = useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl",host="www.ensembl.org")
filters = listFilters(ensembl)

## input omnipath
op <- read_delim(paste0(pan3can_shared_dataD, "Phospho_databases/OmniPath/ptms.txt"),"\t", escape_double = FALSE, trim_ws = TRUE)

## input manually curated gene symbol mapping
omnipath_all_gene_uniprot_mapped <- read_delim(paste0(ppnD, "compile_enzyme_substrate/compile_omnipath/omnipath_all_gene_uniprot_mapped.tab"),
                                               "\t", escape_double = FALSE, trim_ws = TRUE)
omnipath_all_gene_uniprot_mapped <- rbind(omnipath_all_gene_uniprot_mapped,
                                          data.frame(From = "P62158", To = "CALM1"))

## add gene symbols for enzymes and substrates
op <- merge(op, omnipath_all_gene_uniprot_mapped, by.x = c("enzyme"), by.y = c("From"), all.x = T)
which(is.na(op$To))
colnames(op)[ncol(op)] <- "GENE"
op <- merge(op, omnipath_all_gene_uniprot_mapped, by.x = c("substrate"), by.y = c("From"), all.x = T)
which(is.na(op$To))
colnames(op)[ncol(op)] <- "SUB_GENE"


# examine which phosphosites are detected by global phosphoproteom --------
# for (cancer in c("BRCA", "OV", "CO")) {
# # for (cancer in c("OV", "CO")) {
#   ## input CPTAC prospective data
#   pho <- fread(paste0(inputD, cancer, "/", prefix[cancer], "_", "PHO", "_formatted_normalized_replicate_averaged.txt"),
#                data.table = F)
#   pho_head <- formatPhosphosite(phosphosite_vector = as.vector(pho$Phosphosite), gene_vector = as.vector(pho$Gene))
#   op[, paste0(cancer, "_sub_gene")] <- (op$SUB_GENE %in% pho_head$SUBSTRATE)
#   op[, paste0(cancer, "_sub_site")] <- FALSE
#   indata <- sapply(1:nrow(op), FUN = function(n, tab, data) {
#     residue <- tab$residue_type[n]
#     pos <- tab$residue_offset[n]
#     rsd <- paste0(residue, pos)
#     gene <- tab$SUB_GENE[n]
#     row_rsd <- which(data$SUBSTRATE == gene & data$SUB_MOD_RSD == rsd)
#     if (length(row_rsd) > 0) {
#       tmp <- TRUE
#     } else {
#       tmp <- FALSE
#     }
#     return(tmp)
#   }, tab = op, data = pho_head)
#   op[, paste0(cancer, "_sub_site")] <- indata
# }
# write.table(x = op,
#             file = paste0(makeOutDir(), "op_w.phosphosite.availability.txt"),
#             row.names = F, col.names = T, quote = F, sep = ",")
op <- fread(input = "cptac2p/cptac_shared/analysis_results/preprocess_files/tables/convert_omnipath2genomic_coord/op_w.phosphosite.availability.txt", data.table = F)


# take those substrates at least one phosphosite has been detected but not all are detected and convert to genomic positions --------
# ## convert substrate uniprot IDs to ensembl peptide ID
# op_substrate_uniprots <- unique(c(as.vector(op$substrate[(op$BRCA_sub_gene & !op$BRCA_sub_site) | (op$OV_sub_gene & !op$OV_sub_site) | (op$CO_sub_gene & !op$CO_sub_site)])))
# op_substrate_uniprots_mapTab = getBM(attributes = c("ensembl_peptide_id", "uniprotswissprot", "hgnc_symbol"), 
#                filters = "uniprotswissprot", 
#                values = op_substrate_uniprots, mart = ensembl, uniqueRows=T)
# 
# ## add in the un-converted uniprot IDs
# uniprots_unconverted <- op_substrate_uniprots[!(op_substrate_uniprots %in% op_substrate_uniprots_mapTab$uniprotswissprot)]
# uniprots_unconverted_genes <- unique(op[op$substrate %in% uniprots_unconverted, c("SUB_GENE", "substrate")])
# op_substrate_uniprots_mapTab_patch <- getBM(attributes = c("ensembl_peptide_id", "uniprotswissprot", "hgnc_symbol"), 
#                                             filters = "hgnc_symbol", 
#                                             values = unique(uniprots_unconverted_genes$SUB_GENE), mart = ensembl, uniqueRows=T)
# op_substrate_uniprots_mapTab_patch$uniprotswissprot <- NULL
# op_substrate_uniprots_mapTab_patch <- merge(uniprots_unconverted_genes, op_substrate_uniprots_mapTab_patch, by.x = c("SUB_GENE"), by.y = c("hgnc_symbol"), all.x = T)
# colnames(op_substrate_uniprots_mapTab_patch) <- c("hgnc_symbol", "uniprotswissprot", "ensembl_peptide_id")
# op_substrate_uniprots_mapTab <- rbind(op_substrate_uniprots_mapTab_patch[,colnames(op_substrate_uniprots_mapTab)],
#                                       op_substrate_uniprots_mapTab)
# 
# ## get the peptide sequence and length for each ensembl peptide IDs, sort by length
# peptide_seq <- sapply(as.vector(op_substrate_uniprots_mapTab$ensembl_peptide_id), FUN = function(pid) {
#   r <- GET(paste("http://rest.ensembl.org", paste0("/sequence/id/", pid, "?"), sep = ""), content_type("application/json"))
#   stop_for_status(r)
#   seq <- fromJSON(toJSON(content(r)))$seq
#   print(pid)
#   return(seq)
# })
# op_substrate_uniprots_mapTab$peptide_seq <- peptide_seq
# op_substrate_uniprots_mapTab$peptide_len <- sapply(as.vector(op_substrate_uniprots_mapTab$peptide_seq), FUN = function(s) {
#   len <- nchar(s)
#   return(len)
# })
# op_substrate_uniprots_mapTab <- op_substrate_uniprots_mapTab[order(op_substrate_uniprots_mapTab$hgnc_symbol, op_substrate_uniprots_mapTab$peptide_len, decreasing = T),]
# write.table(x = op_substrate_uniprots_mapTab,
#             file = paste0(makeOutDir(), "op_substrate_uniprots_mapTab.txt"),
#             row.names = F, col.names = T, quote = F, sep = ",")
op_substrate_uniprots_mapTab <- fread(input = "cptac2p/cptac_shared/analysis_results/preprocess_files/tables/convert_omnipath2genomic_coord/op_substrate_uniprots_mapTab.txt", data.table = F)
op_substrate_uniprots_mapTab <- op_substrate_uniprots_mapTab[order(op_substrate_uniprots_mapTab$hgnc_symbol, op_substrate_uniprots_mapTab$peptide_len, decreasing = T),]
peptide_seq <- as.vector(op_substrate_uniprots_mapTab$peptide_seq)
peptide_len <- as.vector(op_substrate_uniprots_mapTab$peptide_len)
peptide_gene <- as.vector(op_substrate_uniprots_mapTab$hgnc_symbol)
peptide_id <- as.vector(op_substrate_uniprots_mapTab$ensembl_peptide_id)


# filter multiple ensembl peptide IDs mapping to one uniprot IDs by matching the residue on the peptide sequence -----------------------------------------------------------------------
### first try to find at least one sensembl peptide IDs whose sequence matches all the phosphosites
### if none found, use the longest to shortest peptide ID to convert to genomic coordination, use longer peptide when possible
# op_substrate_siteTab <- unique(op[(op$BRCA_sub_gene & !op$BRCA_sub_site) | (op$OV_sub_gene & !op$OV_sub_site) | (op$CO_sub_gene & !op$CO_sub_site), c("SUB_GENE", "residue_type", "residue_offset")])
# matched_ens_pid <- vector(mode = "character", length = nrow(op_substrate_siteTab))
# for (i in 1:nrow(op_substrate_siteTab)) {
#   gene <- op_substrate_siteTab$SUB_GENE[i]
#   rsd <- op_substrate_siteTab$residue_type[i]
#   offset <- op_substrate_siteTab$residue_offset[i]
#   ens_pids_rows <- which(peptide_gene == gene)
#   for (j in ens_pids_rows) {
#     rsd_j <- substr(x = peptide_seq[j], start = offset, stop = offset)
#     if (!is.na(rsd_j) & rsd_j == rsd) {
#       matched_ens_pid[i] <- peptide_id[j]
#       break
#     }
#   }
# }
# op_substrate_siteTab$matched_ens_pid <- matched_ens_pid
# op_substrate_siteTab_match1 <- op_substrate_siteTab[op_substrate_siteTab$matched_ens_pid != "",]
# op_substrate_siteTab_match2 <- op_substrate_siteTab[op_substrate_siteTab$matched_ens_pid == "",]
# 
# ## take the unmatched genes and use gene symbol (rather than uniprot Ids) to get all ensembl peptide ids
# genes_unmatched <- unique(op_substrate_siteTab_match2$SUB_GENE)
# genes_unmatched_mapTab <- getBM(attributes = c("ensembl_peptide_id", "hgnc_symbol"), 
#                                             filters = "hgnc_symbol", 
#                                             values = genes_unmatched, mart = ensembl, uniqueRows=T)
# genes_unmatched_mapTab <- genes_unmatched_mapTab[!(genes_unmatched_mapTab$ensembl_peptide_id == "") & !(genes_unmatched_mapTab$ensembl_peptide_id %in% op_substrate_uniprots_mapTab$ensembl_peptide_id),]
# 
# ## get peptide sequences for these new ensemble peptides
# genes_unmatched_mapTab$peptide_seq <- sapply(as.vector(genes_unmatched_mapTab$ensembl_peptide_id), FUN = function(pid) {
#   r <- GET(paste("http://rest.ensembl.org", paste0("/sequence/id/", pid, "?"), sep = ""), content_type("application/json"))
#   stop_for_status(r)
#   seq <- fromJSON(toJSON(content(r)))$seq
#   print(pid)
#   return(seq)
# })
# genes_unmatched_mapTab$peptide_len <- sapply(as.vector(genes_unmatched_mapTab$peptide_seq), FUN = function(s) {
#   len <- nchar(s)
#   return(len)
# })
# genes_unmatched_mapTab <- genes_unmatched_mapTab[order(genes_unmatched_mapTab$hgnc_symbol, genes_unmatched_mapTab$peptide_len, decreasing = T),]
# 
# ## do another round of matching
# op_substrate_siteTab_match2$matched_ens_pid <- NULL
# matched_ens_pid <- vector(mode = "character", length = nrow(op_substrate_siteTab_match2))
# for (i in 1:nrow(op_substrate_siteTab_match2)) {
#   gene <- op_substrate_siteTab_match2$SUB_GENE[i]
#   rsd <- op_substrate_siteTab_match2$residue_type[i]
#   offset <- op_substrate_siteTab_match2$residue_offset[i]
#   ens_pids_rows <- which(as.vector(genes_unmatched_mapTab$hgnc_symbol) == gene)
#   for (j in ens_pids_rows) {
#     rsd_j <- substr(x = as.vector(genes_unmatched_mapTab$peptide_seq)[j], start = offset, stop = offset)
#     if (!is.na(rsd_j) & rsd_j == rsd) {
#       matched_ens_pid[i] <- as.vector(genes_unmatched_mapTab$ensembl_peptide_id)[j]
#       break
#     }
#   }
# }
# op_substrate_siteTab_match2$matched_ens_pid <- matched_ens_pid
# 
# ## merge
# op_substrate_siteTab <- rbind(op_substrate_siteTab_match1, op_substrate_siteTab_match2)
# write.table(x = op_substrate_siteTab,
#             file = paste0(makeOutDir(), "op_substrate_siteTab.txt"),
#             row.names = F, col.names = T, quote = F, sep = ",")
op_substrate_siteTab <- fread(input = "cptac2p/cptac_shared/analysis_results/preprocess_files/tables/convert_omnipath2genomic_coord/op_substrate_siteTab.txt", data.table = F)

# convert ensembl peptide ID + residue to genomic coord ------------------------------------------------------------------
op_substrate_siteTab_wpid <- op_substrate_siteTab[op_substrate_siteTab$matched_ens_pid != "",]

# test <- lapply(1:10, FUN = function(n, pids, offsets) {
coords <- lapply(1:nrow(op_substrate_siteTab_wpid), FUN = function(n, pids, offsets) {
  r <- GET(paste("http://rest.ensembl.org", paste0("/map/translation/", pids[n], "/", offsets[n], "..", offsets[n], "?"), sep = ""), content_type("application/json"))
  stop_for_status(r)
  tmp <- fromJSON(toJSON(content(r)))$mappings
  print(pids[n])
  return(tmp)
}, pids = as.vector(op_substrate_siteTab_wpid$matched_ens_pid), offsets = as.vector(op_substrate_siteTab_wpid$residue_offset))

## TODO: resolve incorrect coordinates, e.g. i = 19
df <- matrix(NA, nrow=length(coords), ncol = ncol(coords[[1]]))
for (i in 1:length(coords)) {
  if (nrow(coords[[i]]) == 1) {
    df[i,] <- matrix(unlist(coords[[i]]), nrow = 1, byrow = T)
  }
}
df <- data.frame(df,stringsAsFactors=FALSE)
colnames(df) <- colnames(coords[[1]])
op_substrate_siteTab_wpid <- cbind(op_substrate_siteTab_wpid, df)
write.table(x = op_substrate_siteTab_wpid,
            file = paste0(makeOutDir(), "op_substrate_siteTab_wpid.txt"),
            row.names = F, col.names = T, quote = F, sep = ",")




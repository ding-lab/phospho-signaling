# Yige Wu @ WashU 2018 Jan
# convert site-level ks pairs from NetWorKin to genomic positions
# use ensembl peptide ids to convert to genomic positions

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



# inputs ------------------------------------------------------------------
## input CPTAC phosphosite table
pho_head_cancers <- NULL
for (cancer in c("BRCA", "OV", "CO", "UCEC", "CCRCC")) {
  pho <- fread(paste0(inputD, cancer, "/", prefix[cancer], "_", "PHO", "_formatted.txt"),
               data.table = F)
  pho_head <- formatPhosphosite(phosphosite_vector = as.vector(pho$Phosphosite), gene_vector = as.vector(pho$Gene))
  pho_head_cancers <- rbind(pho_head_cancers, pho_head)
  pho_head_cancers <- unique(pho_head_cancers)
}
pho_head_cancers$offset <- str_split_fixed(string = pho_head_cancers$SUB_MOD_RSD, pattern = '[STY]', 3)[,2]
pho_head_cancers$offset.2 <- str_split_fixed(string = pho_head_cancers$SUB_MOD_RSD, pattern = '[STY]', 3)[,3]
pho_head_cancers_sim <- pho_head_cancers[pho_head_cancers$offset.2 == "",]

## input NetWorKin table
networkin <- read_delim(paste0(pan3can_shared_dataD, "Phospho_databases/NetworKIN/networkin_human_predictions_3.1.tsv"), "\t", escape_double = FALSE, trim_ws = TRUE)
### filter networkin table with prediction score
networkin_cpatc <- networkin[networkin$networkin_score >= 2,]
networkin_cpatc <- data.frame(networkin_cpatc)
### filter networkin table with cptac gene names
networkin_cpatc <- networkin_cpatc[networkin_cpatc$substrate_name %in% pho_head_cancers_sim$SUBSTRATE,]
networkin_cpatc$substrate_ens_pid <- str_split_fixed(string = networkin_cpatc$X.substrate, pattern = '[ ()]', n = 4)[,3]
networkin_cpatc_sites <- unique(networkin_cpatc[, c("substrate_name", "position", "substrate_ens_pid")])

# convert networkin phosphosites to genomic coordinates ----------------------------
# coords <- list()
# pids = as.vector(networkin_cpatc_sites$substrate_ens_pid)
# offsets = as.vector(networkin_cpatc_sites$position)
# for (n in 1:nrow(networkin_cpatc_sites)) {
#   if (!is.na(pids[n])) {
#     r <- GET(paste0("http://rest.ensembl.org", paste0("/map/translation/", pids[n], "/", offsets[n], "..", offsets[n], "?"), sep = ""), content_type("application/json"))
#     msg <- try(stop_for_status(r))
#     if (!grepl(x = msg[1], pattern = "Error")) {
#       tmp <- fromJSON(toJSON(content(r)))$mappings
#       coords[[n]] <- tmp
#     }
#     print(n)
#   } 
# }
# saveRDS(object = coords, file = paste0(makeOutDir(), "coords.RDS"))

df <- matrix(data = NA, nrow = nrow(networkin_cpatc_sites), ncol = ncol(coords[[1]]))
for (i in 1:nrow(networkin_cpatc_sites)) {
  if (!is.null(coords[[i]])) {
    if (nrow(coords[[i]]) == 1) {
      df[i,] <- matrix(unlist(coords[[i]]), nrow = 1, byrow = T)
    }
  }
}
df <- data.frame(df,stringsAsFactors=FALSE)
colnames(df) <- colnames(coords[[1]])
networkin_cpatc_sites <- cbind(networkin_cpatc_sites, df)
write.table(x = networkin_cpatc_sites,
            file = paste0(makeOutDir(), "networkin_cpatc_sites_wpid.txt"),
            row.names = F, col.names = T, quote = F, sep = ",")

rm(coords)

  

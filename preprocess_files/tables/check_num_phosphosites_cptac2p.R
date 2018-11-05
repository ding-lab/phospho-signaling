# Yige Wu @ WashU 2018 Jan
# convert CPTAC global proteomids phosphosites to genomic positions
# refseq IDs to ensembl peptide ids to genomic coordinates

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
ensembl = useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl",host="www.ensembl.org") ## 38
filters = listFilters(ensembl)
ensembl37 = useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl",host="grch37.ensembl.org")

# gather all phosphosites -------------------------------------------------
pho_head_cancers <- NULL
for (cancer in c("BRCA", "OV", "CO")) {
  pho <- fread(paste0(cptac_sharedD, cancer, "/", prefix[cancer], "_", "PHO", "_formatted.txt"),
               data.table = F)
  pho_head <- formatPhosphosite(phosphosite_vector = as.vector(pho$Phosphosite), gene_vector = as.vector(pho$Gene))
  pho_head_cancers <- rbind(pho_head_cancers, pho_head)
  pho_head_cancers <- unique(pho_head_cancers)
}
pho_head_cancers$offset <- str_split_fixed(string = pho_head_cancers$SUB_MOD_RSD, pattern = '[STY]', 3)[,2]
pho_head_cancers$offset.2 <- str_split_fixed(string = pho_head_cancers$SUB_MOD_RSD, pattern = '[STY]', 3)[,3]
pho_head_cancers_sim <- pho_head_cancers[pho_head_cancers$offset.2 == "",]


# gather filtered phosphosites -------------------------------------------------
num_nonNA <- 15
pho_head_cancers <- NULL
for (cancer in c("BRCA", "OV", "CO")) {
  pho <- fread(paste0(cptac_sharedD, cancer, "/", prefix[cancer], "_", "PHO", "_formatted.txt"),
               data.table = F)
  pho <- pho[rowSums(!is.na(pho)) >= (2+num_nonNA),]
  pho_head <- formatPhosphosite(phosphosite_vector = as.vector(pho$Phosphosite), gene_vector = as.vector(pho$Gene))
  pho_head_cancers <- rbind(pho_head_cancers, pho_head)
  pho_head_cancers <- unique(pho_head_cancers)
}
pho_head_cancers$offset <- str_split_fixed(string = pho_head_cancers$SUB_MOD_RSD, pattern = '[STY]', 3)[,2]
pho_head_cancers$offset.2 <- str_split_fixed(string = pho_head_cancers$SUB_MOD_RSD, pattern = '[STY]', 3)[,3]
pho_head_cancers_sim <- pho_head_cancers[pho_head_cancers$offset.2 == "",]
nrow(pho_head_cancers_sim)



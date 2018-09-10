# Yige Wu @ WashU March 2018
# plot heatmap showing mutations, gene expression and protein expression for MMR genes

# source ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/preprocess_files/preprocess_files.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')

## input MMR genes
mmr_gene <- read_delim("~/Box Sync/MSI_CPTAC/Data/mmr_gene.txt","\t", escape_double = FALSE, col_names = FALSE,trim_ws = TRUE)
mmr_gene <- mmr_gene$X1

## input merged CPTAC somatic maf
test <- fread(input = "./MSI_CPTAC/Data/mergerdMMR_CPTAC.somatic.maf", data.table = F, nrows = 10)
somatic <- fread(input = "./MSI_CPTAC/Data/MAF/mergerdCPTAC.somatic.maf", data.table = F, header = F, col.names = colnames(test))
somatic$Cancer_Type

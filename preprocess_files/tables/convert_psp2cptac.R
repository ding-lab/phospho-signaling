# Yige Wu @ WashU 2018 Apr
# convert site-level ks pairs from PhosphositePlus by matching genomic coordinations with CPTAC data

# source -------------------------------------------------------------------
library(data.table)
library(readr)
source('/Users/yigewu/Box Sync/cptac2p_analysis/preprocess_files/preprocess_files.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')

# inputs ------------------------------------------------------------------
op_cptac <- fread(input = "cptac2p/cptac_shared/analysis_results/preprocess_files/tables/convert_omnipath_site2cptac/op_cptac_unconverted.txt", data.table = F, sep = ",")
op_cptac$SUB_MOD_RSD <- paste0(op_cptac$residue_type, op_cptac$residue_offset)
op_cptac$SUB_MOD_RSD_cptac <- paste0(op_cptac$residue_type, op_cptac$residue_offset_cptac)

for( enzyme_type in c("kinase")) {
  es_tab <- load_ks_table(protein = enzyme_type)
  es_tab <- merge(es_tab, op_cptac, by.x = c("KIN_ACC_ID", "SUB_ACC_ID", "SUB_MOD_RSD"), by.y = c("enzyme", "substrate", "SUB_MOD_RSD"), all.x = T)
  
}
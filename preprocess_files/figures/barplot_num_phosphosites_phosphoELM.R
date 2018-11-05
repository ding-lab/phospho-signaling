# Yige Wu @ WashU 2018 July
## show the number of 

# source ------------------------------------------------------------------
source("Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R")
library(readr)

# set variables -----------------------------------------------------------
cancers_sort <- c("BRCA", "OV", "CO")
# input phosphosite genomic coordinates (hg38) -----------------------------------
pho_head_geno_pos <- fread(input = paste0("./cptac2p/analysis_results/preprocess_files/tables/convert_cptac_pho2genomic_coord/pho_head_cancers_single_site_wpid.txt"), data.table = F)
# pho_head_geno_pos <- fread(input = paste0("./cptac2p/analysis_results/preprocess_files/tables/convert_cptac_pho2genomic_coord/pho_head_cancers_single_site_wpid_hg38_hg19.txt"), data.table = F)
pho_head_geno_pos$id <- paste0(pho_head_geno_pos$SUBSTRATE, ":", pho_head_geno_pos$SUB_MOD_RSD)
pho_head_geno_pos <- pho_head_geno_pos[!duplicated(pho_head_geno_pos$id),]
rownames(pho_head_geno_pos) <- pho_head_geno_pos$id
pho_head_geno_pos$geno_pos <- paste0(pho_head_geno_pos$seq_region_name, ":", pho_head_geno_pos$start, "-", pho_head_geno_pos$end)

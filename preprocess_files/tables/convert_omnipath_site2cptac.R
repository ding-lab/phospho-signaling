# Yige Wu @ WashU 2018 Jan
# convert site-level ks pairs from Omnipath by matching genomic coordinations with CPTAC data

# source -------------------------------------------------------------------
library(data.table)
library(readr)
source('/Users/yigewu/Box Sync/cptac2p_analysis/preprocess_files/preprocess_files.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')

# inputs ------------------------------------------------------------------
op <- fread(input = "cptac2p/analysis_results/preprocess_files/tables/convert_omnipath2genomic_coord/op_w.phosphosite.availability.txt", data.table = F)
op <- data.frame(op)
op2add <- op[(op$BRCA_sub_gene & !op$BRCA_sub_site) | (op$OV_sub_gene & !op$OV_sub_site) | (op$CO_sub_gene & !op$CO_sub_site),]
op_matched <- op[!((op$BRCA_sub_gene & !op$BRCA_sub_site) | (op$OV_sub_gene & !op$OV_sub_site) | (op$CO_sub_gene & !op$CO_sub_site)) & (op$BRCA_sub_site | op$OV_sub_site | op$CO_sub_site),]

## input omnipath+networkin table with genomic coordinations
op_substrate_siteTab_wpid <- fread(input = "cptac2p/analysis_results/preprocess_files/tables/convert_omnipath2genomic_coord/op_substrate_siteTab_wpid.txt", data.table = F)

## input CPTAC phosphosites with genomic coordinations
pho_head_cancers_sim <- fread(input = "cptac2p/analysis_results/preprocess_files/tables/convert_cptac_pho2genomic_coord/pho_head_cancers_single_site_wpid.txt", data.table = F)

# edit the residue offset according to the genomic-matched cptac phosphosites ----------------
residue_offset_cptac <- vector(mode = "numeric", length = nrow(op_substrate_siteTab_wpid)) + NA
for (i in 1:nrow(op_substrate_siteTab_wpid)) {
  matched_row <- which(pho_head_cancers_sim$seq_region_name == op_substrate_siteTab_wpid$seq_region_name[i] & pho_head_cancers_sim$start == op_substrate_siteTab_wpid$start[i])
  if (length(matched_row) > 0) {
    residue_offset_cptac[i] <- pho_head_cancers_sim$offset[matched_row[1]]
  }
}
op_substrate_siteTab_wpid$residue_offset_cptac <- residue_offset_cptac
op_substrate_siteTab_wpid <- op_substrate_siteTab_wpid[op_substrate_siteTab_wpid$residue_offset_cptac > 0,]
op2add <- merge(op2add, op_substrate_siteTab_wpid[, c("SUB_GENE", "residue_type", "residue_offset", "residue_offset_cptac")], by = c("SUB_GENE", "residue_type", "residue_offset"))
op2add <- data.frame(op2add)

tmp1 <- op2add
tmp2 <- op_matched
tmp2$residue_offset_cptac <- tmp2$residue_offset
tmp <- rbind(tmp1, tmp2)
write.table(x = tmp,
            file = paste0(makeOutDir(), "op_cptac_unconverted.txt"),
            row.names = F, col.names = T, quote = F, sep = ",")

## merge
op2add$residue_offset <- op2add$residue_offset_cptac
op2add$residue_offset_cptac <- NULL
op_cptac <- rbind(op_matched, op2add)
write.table(x = op_cptac,
            file = paste0(makeOutDir(), "op_cptac_converted.txt"),
            row.names = F, col.names = T, quote = F, sep = ",")

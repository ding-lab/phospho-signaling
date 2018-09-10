# Yige Wu @ WashU 2018 Jan
# convert site-level ks pairs from NetWorKin by matching genomic coordinations with CPTAC data

# source -------------------------------------------------------------------
library(data.table)
library(readr)
source('/Users/yigewu/Box Sync/cptac2p_analysis/preprocess_files/preprocess_files.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')

# inputs ------------------------------------------------------------------
## input CPTAC phosphosites with genomic coordinations
pho_head_cancers_sim <- fread(input = "cptac2p/cptac_shared/analysis_results/preprocess_files/tables/convert_cptac_pho2genomic_coord/pho_head_cancers_single_site_wpid.txt", data.table = F)

## input networkin table
networkin <- read_delim(paste0("./pan3can_shared_data/Phospho_databases/NetworKIN/networkin_human_predictions_3.1.tsv"), "\t", escape_double = FALSE, trim_ws = TRUE)
### filter networkin table with prediction score
networkin_cpatc <- networkin[networkin$networkin_score >= 2,]
networkin_cpatc <- data.frame(networkin_cpatc)
networkin_cpatc <- networkin_cpatc[networkin_cpatc$substrate_name %in% pho_head_cancers_sim$SUBSTRATE,]

## add residue type
residue_type <- vector("character", length = nrow(networkin_cpatc))
for (residue in c("s", "t", "y")) {
  residue_type[grepl(x = as.vector(networkin_cpatc$sequence), pattern = residue, ignore.case = F)] <- toupper(residue)
}
networkin_cpatc$residue_type <- residue_type

## input networkin phosphosites with genomic coordinates
networkin_cpatc_sites <- fread(input = "cptac2p/cptac_shared/analysis_results/preprocess_files/tables/convert_networkin2genomic_coord/networkin_cpatc_sites_wpid.txt", data.table = F)

# annotate netwokin table to whether it's in the dataset ----------------
# for (cancer in c("BRCA", "OV", "CO")) {
#   ## input CPTAC prospective data
#   pho <- fread(paste0(inputD, cancer, "/", prefix[cancer], "_", "PHO", "_formatted_normalized_replicate_averaged.txt"),
#                data.table = F)
#   pho_head <- formatPhosphosite(phosphosite_vector = as.vector(pho$Phosphosite), gene_vector = as.vector(pho$Gene))
#   networkin_cpatc[, paste0(cancer, "_sub_gene")] <- (networkin_cpatc$substrate_name %in% pho_head$SUBSTRATE)
#   networkin_cpatc[, paste0(cancer, "_sub_site")] <- FALSE
#   indata <- sapply(1:nrow(networkin_cpatc), FUN = function(n, tab, data) {
#     residue <- tab$residue_type[n]
#     pos <- tab$position[n]
#     rsd <- paste0(residue, pos)
#     gene <- tab$substrate_name[n]
#     row_rsd <- which(data$SUBSTRATE == gene & data$SUB_MOD_RSD == rsd)
#     if (length(row_rsd) > 0) {
#       tmp <- TRUE
#     } else {
#       tmp <- FALSE
#     }
#     return(tmp)
#   }, tab = networkin_cpatc, data = pho_head)
#   networkin_cpatc[, paste0(cancer, "_sub_site")] <- indata
# }
# write.table(x = networkin_cpatc,
#             file = paste0(makeOutDir(), "nk_w.phosphosite.availability.txt"),
#             row.names = F, col.names = T, quote = F, sep = ",")

# divide netwokin table into already matched and unmatched ----------------
networkin2add <- networkin_cpatc[(networkin_cpatc$BRCA_sub_gene & !networkin_cpatc$BRCA_sub_site) | (networkin_cpatc$OV_sub_gene & !networkin_cpatc$OV_sub_site) | (networkin_cpatc$CO_sub_gene & !networkin_cpatc$CO_sub_site),]
networkin_matched <- networkin_cpatc[!((networkin_cpatc$BRCA_sub_gene & !networkin_cpatc$BRCA_sub_site) | (networkin_cpatc$OV_sub_gene & !networkin_cpatc$OV_sub_site) | (networkin_cpatc$CO_sub_gene & !networkin_cpatc$CO_sub_site)) & (networkin_cpatc$BRCA_sub_site | networkin_cpatc$OV_sub_site | networkin_cpatc$CO_sub_site),]

## filter the sites by unmatched sites
networkin_cpatc_sites <- networkin_cpatc_sites[(networkin_cpatc_sites$substrate_name %in% networkin2add$substrate_name) & !is.na(networkin_cpatc_sites$end),]

# edit the residue offset according to the genomic-matched cptac phosphosites ----------------
residue_offset_cptac <- vector(mode = "numeric", length = nrow(networkin_cpatc_sites)) + NA
for (i in 1:nrow(networkin_cpatc_sites)) {
  matched_row <- which(pho_head_cancers_sim$seq_region_name == networkin_cpatc_sites$seq_region_name[i] & pho_head_cancers_sim$start == networkin_cpatc_sites$start[i])
  if (length(matched_row) > 0) {
    residue_offset_cptac[i] <- pho_head_cancers_sim$offset[matched_row[1]]
  }
}
networkin_cpatc_sites$residue_offset_cptac <- residue_offset_cptac
networkin_cpatc_sites <- networkin_cpatc_sites[!is.na(networkin_cpatc_sites$residue_offset_cptac),]
networkin2add <- merge(networkin2add, networkin_cpatc_sites[, c("substrate_name", "position", "residue_offset_cptac")], by = c("substrate_name", "position"))
networkin2add$position <- networkin2add$residue_offset_cptac
networkin2add$residue_offset_cptac <- NULL

## merge
networkin_cpatc <- rbind(networkin_matched, networkin2add)
write.table(x = networkin_cpatc,
            file = paste0(makeOutDir(), "networkin_cptac_converted.txt"),
            row.names = F, col.names = T, quote = F, sep = ",")

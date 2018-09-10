# Yige Wu @ WashU 2018 Jan
# compile protein-level and site-level ks pairs from Omnipath

# souce ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/cptac2p_analysis_shared.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
source("https://bioconductor.org/biocLite.R")
biocLite("devtools")
BiocInstaller::biocLite('grimbough/biomaRt')
library(biomaRt)

# inputs ------------------------------------------------------------------
makeOutDir()
# compile protein level ---------------------------------------------------
## input Omnipath PTM sites
ptms <- read_delim("~/Box Sync/pan3can_shared_data/Phospho_databases/OmniPath/ptms.txt","\t", escape_double = FALSE, trim_ws = TRUE)

## input NetKIN predicted enzyme-substrate pairs

## input DEPOD
ps_orig <- load_ks_table(protein = "phosphotase")

## transform omnipath uniprot ids to gene names
ptms_pro_pairs <- unique(ptms[, c("enzyme", "substrate", "modification")])
length(unique(c(ptms_pro_pairs$enzyme, ptms_pro_pairs$substrate)))
nrow(ptms_pro_pairs[ptms_pro_pairs$modification == "phosphorylation",])
nrow(ptms_pro_pairs[ptms_pro_pairs$modification == "dephosphorylation",])
# write.table(x = unique(c(ptms_pro_pairs$enzyme, ptms_pro_pairs$substrate)), 
#             file = paste0(makeOutDir(), "omnipath_all_gene_uniprot.txt"),
#             row.names = F, col.names = F, quote = F)
omnipath_all_gene_uniprot_mapped <- read_delim("~/Box Sync/cptac2p/cptac_shared/analysis_results/phospho_network/compile_enzyme_substrate/compile_omnipath/omnipath_all_gene_uniprot_mapped.tab",
                                               "\t", escape_double = FALSE, trim_ws = TRUE)
omnipath_all_gene_uniprot_mapped <- rbind(omnipath_all_gene_uniprot_mapped,
                                          data.frame(From = "P62158", To = "CALM1"))
ptms_pro_pairs <- merge(ptms_pro_pairs, omnipath_all_gene_uniprot_mapped, by.x = c("enzyme"), by.y = c("From"), all.x = T)
colnames(ptms_pro_pairs)[ncol(ptms_pro_pairs)] <- "GENE"
ptms_pro_pairs <- merge(ptms_pro_pairs, omnipath_all_gene_uniprot_mapped, by.x = c("substrate"), by.y = c("From"), all.x = T)
colnames(ptms_pro_pairs)[ncol(ptms_pro_pairs)] <- "SUB_GENE"
ptms_pro_pairs$source <- "OmniPath"

## modify in kinase-substrate table from Phosphositeplus
ks_orig <- load_ks_table(protein = "kinase")
ks_orig_pro_pairs <- unique(ks_orig[, c("KIN_ACC_ID", "SUB_ACC_ID", "GENE", "SUB_GENE")])
ks_orig_pro_pairs$modification <- "phosphorylation"
colnames(ks_orig_pro_pairs) <- c("enzyme", "substrate", "GENE", "SUB_GENE", "modification")
ks_orig_pro_pairs$source <- "PhosphositePlus"

## modify in phosphotase-substrate table from DEPOD
ps_orig_pro_pairs <- unique(ps_orig[, c("Phosphatase_UniProtAC_human", "Substrate_UniProtAC_ref")])
colnames(ps_orig_pro_pairs) <- c("enzyme", "substrate")
ps_orig_pro_pairs <- merge(ps_orig_pro_pairs, omnipath_all_gene_uniprot_mapped, by.x = c("enzyme"), by.y = c("From"), all.x = T)
colnames(ps_orig_pro_pairs)[ncol(ps_orig_pro_pairs)] <- "GENE"
ps_orig_pro_pairs <- merge(ps_orig_pro_pairs, omnipath_all_gene_uniprot_mapped, by.x = c("substrate"), by.y = c("From"), all.x = T)
colnames(ps_orig_pro_pairs)[ncol(ps_orig_pro_pairs)] <- "SUB_GENE"
ps_orig_pro_pairs$modification <- "dephosphorylation"
ps_orig_pro_pairs$source <- "DEPOD"

## merge enzyme-substrate table
ptms_pro_pairs_sup <- merge(unique(ptms_pro_pairs[, c("GENE", "SUB_GENE", "modification", "source")]), 
              unique(ks_orig_pro_pairs[, c("GENE", "SUB_GENE", "modification", "source")]), 
              by = c("GENE", "SUB_GENE", "modification"), all = T, suffixes = c(".omni", ".pps"))
ptms_pro_pairs_sup <- merge(ptms_pro_pairs_sup, 
                            unique(ps_orig_pro_pairs[, c("GENE", "SUB_GENE", "modification", "source")]), 
                            by = c("GENE", "SUB_GENE", "modification"), all = T)
colnames(ptms_pro_pairs_sup)[ncol(ptms_pro_pairs_sup)] <- "source.depod"
ptms_pro_pairs_sup <- ptms_pro_pairs_sup[!is.na(ptms_pro_pairs_sup$GENE) & !is.na(ptms_pro_pairs_sup$SUB_GENE),]

write.table(x = ptms_pro_pairs_sup,
            file = paste0(makeOutDir(), "omnipath_phosphositeplus_depod_enzyme_substrate_protein_level_union.txt"),
            row.names = F, col.names = T, quote = F, sep = "\t")


# compile size-level ------------------------------------------------------
## input Omnipath PTM sites
op <- read_delim("~/Box Sync/pan3can_shared_data/Phospho_databases/OmniPath/ptms.txt","\t", escape_double = FALSE, trim_ws = TRUE)
## annotate omnipath site level to gene name
op <- merge(op, omnipath_all_gene_uniprot_mapped, by.x = c("enzyme"), by.y = c("From"), all.x = T)
colnames(op)[ncol(op)] <- "GENE"
op <- merge(op, omnipath_all_gene_uniprot_mapped, by.x = c("substrate"), by.y = c("From"), all.x = T)
colnames(op)[ncol(op)] <- "SUB_GENE"
op$SUB_MOD_RSD <- paste0(op$residue_type, ptms$residue_offset)
op_site_pairs <- cbind(data.frame(KINASE = op$GENE, KIN_ACC_ID = op$enzyme, GENE = op$GENE, KIN_ORGANISM = "human",
                                        SUBSTRATE = op$SUB_GENE, SUB_GENE_ID = "NULL", SUB_ACC_ID = op$substrate, SUB_GENE = op$SUB_GENE, SUB_ORGANISM = "human",
                                        SUB_MOD_RSD = op$SUB_MOD_RSD, SITE_GRP_ID = "NULL", SITE_...7_AA = "NULL", 
                                        networkin_score = Inf, Source = "OmniPath"))

write.table(x = op,
            file = paste0(makeOutDir(), "omnipath_genename_annotated.csv"),
            row.names = F, col.names = T, quote = F, sep = ",")


## input NetKIN predicted enzyme-substrate pairs
networkin <- read_delim("~/Box Sync/pan3can_shared_data/Phospho_databases/NetworKIN/networkin_human_predictions_3.1.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
networkin <- networkin[networkin$networkin_score >= 2,]
## get kinase ensembl id
networkin$kinase_ensembl <- str_split_fixed(string = networkin$string_path, pattern = ",", 2)[,1]
## annotate kinase name
ensembl = useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl",host="www.ensembl.org")
# filters = listFilters(ensembl)
mapTab = getBM(attributes = c("hgnc_symbol", "ensembl_peptide_id"), filters = "ensembl_peptide_id", values = unique(networkin$kinase_ensembl), mart = ensembl, uniqueRows=FALSE)
networkin <- merge(networkin, mapTab, by.x = c("kinase_ensembl"), by.y = c("ensembl_peptide_id"), all.x = T)
networkin <- networkin[!is.na(networkin$hgnc_symbol),]
## correct substrate name
substrate_names <- as.vector(networkin$substrate_name)
substrate_names2correct <- substrate_names[grepl(pattern = "ENSP00", x = substrate_names)]
mapTab = getBM(attributes = c("hgnc_symbol", "ensembl_peptide_id"), filters = "ensembl_peptide_id", values = substrate_names2correct, mart = ensembl, uniqueRows=FALSE)
tmp <- as.vector(mapTab$hgnc_symbol)
names(tmp) <- as.vector(mapTab$ensembl_peptide_id)
substrate_names[substrate_names %in% substrate_names2correct] <- tmp[substrate_names2correct]

### annotate residue type
residue_type <- vector("character", length = nrow(networkin))
for (residue in tolower(unique(op$residue_type))) {
  residue_type[grepl(x = as.vector(networkin$sequence), pattern = residue, ignore.case = F)] <- toupper(residue)
}
nk_site_pairs <- cbind(data.frame(KINASE = networkin$id, KIN_ACC_ID = "NULL", GENE = networkin$hgnc_symbol, KIN_ORGANISM = "human",
                                  SUBSTRATE = "NULL", SUB_GENE_ID = "NULL", SUB_ACC_ID = "NULL", SUB_GENE = substrate_names, SUB_ORGANISM = "human",
                                  SUB_MOD_RSD = paste0(residue_type, networkin$position), SITE_GRP_ID = "NULL", SITE_...7_AA = "NULL", 
                                  networkin_score = networkin$networkin_score, Source = "NetKIN"))


## input DEPOD
# ps_orig <- load_ks_table(protein = "phosphotase")
# mapTab = getBM(attributes = c("hgnc_symbol", "uniprotswissprot"), filters = "uniprotswissprot", values = unique(c(as.vector(ps_orig$Phosphatase_UniProtAC_human), as.vector(ps_orig$Substrate_UniProtAC_ref))), mart = ensembl, uniqueRows=FALSE)
# ps_orig <- merge(ps_orig, mapTab, by.x = c("Phosphatase_UniProtAC_human"), by.y = c("uniprotswissprot"), all.x = T)
# colnames(ps_orig)[ncol(ps_orig)] <- "GENE"
# ps_orig <- merge(ps_orig, mapTab, by.x = c("Substrate_UniProtAC_ref"), by.y = c("uniprotswissprot"), all.x = T)
# colnames(ps_orig)[ncol(ps_orig)] <- "SUB_GENE"
# ps_orig[is.na(ps_orig$SUB_GENE),]
# 
# ps_orig_pro_pairs$modification <- "dephosphorylation"
# ps_orig_pro_pairs$source <- "DEPOD"
# ps_orig_pro_pairs <- unique(ps_orig[, c("Phosphatase_UniProtAC_human", "Substrate_UniProtAC_ref")])
# ps_site_pairs <- cbind(data.frame(KINASE = networkin$id, KIN_ACC_ID = "NULL", GENE = networkin$hgnc_symbol, KIN_ORGANISM = "human",
#                                   SUBSTRATE = "NULL", SUB_GENE_ID = "NULL", SUB_ACC_ID = "NULL", SUB_GENE = substrate_names, SUB_ORGANISM = "human",
#                                   SUB_MOD_RSD = paste0(residue_type, networkin$position), SITE_GRP_ID = "NULL", SITE_...7_AA = "NULL", 
#                                   networkin_score = networkin$networkin_score, Source = "DEPOD"))
# 
# 
# ## input phosphositeplus
# psp_orig <- load_ks_table(protein = "kinase")


## write out table
ptms_site_pairs_sup <- rbind(op_site_pairs, nk_site_pairs)
write.table(x = ptms_site_pairs_sup,
            file = paste0(makeOutDir(), "omnipath_networkin_enzyme_substrate_site_level_union.csv"),
            row.names = F, col.names = T, quote = F, sep = ",")

# write.table(x = ptms_site_pairs_sup,
#             file = paste0(makeOutDir(), "omnipath_phosphositeplus_depod_enzyme_substrate_site_level_union.csv"),
#             row.names = F, col.names = T, quote = F, sep = ",")

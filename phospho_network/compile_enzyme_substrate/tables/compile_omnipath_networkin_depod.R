# Yige Wu @ WashU 2018 Apr
# compile protein-level and site-level ks pairs from Omnipath + networkin + DEPOD

# souce ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/cptac2p_analysis_shared.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
source("https://bioconductor.org/biocLite.R")
# biocLite("devtools")
# BiocInstaller::biocLite('grimbough/biomaRt')
library(biomaRt)

# compile size-level for OmniPath ------------------------------------------------------
## input Omnipath PTM sites
op <- read.delim(file = "./cptac2p/analysis_results/preprocess_files/tables/convert_omnipath_site2cptac/op_cptac_converted.txt", sep = ",", header = T)
op <- data.frame(op)
op_site_pairs <- data.frame(KINASE = op$GENE, KIN_ACC_ID = op$enzyme, GENE = op$GENE, KIN_ORGANISM = "human",
                            SUBSTRATE = op$SUB_GENE, SUB_GENE_ID = "NULL", SUB_ACC_ID = op$substrate, SUB_GENE = op$SUB_GENE, SUB_ORGANISM = "human",
                            SUB_MOD_RSD = paste0(op$residue_type, op$residue_offset), SITE_GRP_ID = "NULL", SITE_...7_AA = "NULL", 
                            networkin_score = Inf, Source = op$sources)

# compile size-level for networkin ------------------------------------------------------
## input NetKIN predicted enzyme-substrate pairs
networkin <- read_csv("./cptac2p/analysis_results/preprocess_files/tables/convert_networkin_site2cptac/networkin_cptac_converted.txt")
## get kinase read.delim id
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
                                  SUB_MOD_RSD = paste0(residue_type, networkin$position), SITE_GRP_ID = "NULL", SITE_...7_AA = networkin$sequence, 
                                  networkin_score = networkin$networkin_score, Source = "NetKIN"))


# compile size-level ------------------------------------------------------
## write out table
ptms_site_pairs_sup <- rbind(op_site_pairs, nk_site_pairs)
# write.table(x = ptms_site_pairs_sup,
#             file = paste0(makeOutDir(resultD = resultD), "omnipath_networkin_enzyme_substrate_site_level_union.csv"),
#             row.names = F, col.names = T, quote = F, sep = ",")

## try to organize pairs into the same row
# ptms_site_pairs_sup <- fread(input = "./cptac2p/analysis_results/phospho_network/compile_enzyme_substrate/tables/compile_omnipath_networkin/omnipath_networkin_enzyme_substrate_site_level_union.csv", data.table = F)
ptms_site_pairs_sup$pair <- paste0(ptms_site_pairs_sup$GENE, ":", ptms_site_pairs_sup$SUB_GENE, ":", ptms_site_pairs_sup$SUB_MOD_RSD)
ptms_site_pairs_sup <- ptms_site_pairs_sup[order(ptms_site_pairs_sup$pair),]
ptms_site_pairs_sup_nodup <- ptms_site_pairs_sup[!(ptms_site_pairs_sup$pair %in% ptms_site_pairs_sup$pair[duplicated(ptms_site_pairs_sup$pair)]),]
ptms_site_pairs_sup_dup <- ptms_site_pairs_sup[(ptms_site_pairs_sup$pair %in% ptms_site_pairs_sup$pair[duplicated(ptms_site_pairs_sup$pair)]),]
ptms_site_pairs_sup_dup_other <- ptms_site_pairs_sup_dup[ptms_site_pairs_sup_dup$Source != "NetKIN",]
ptms_site_pairs_sup_dup_other <- ptms_site_pairs_sup_dup_other[!(duplicated(paste0(ptms_site_pairs_sup_dup_other$pair, ":", ptms_site_pairs_sup_dup_other$Source))),]

ptms_site_pairs_sup_dup_other_new <- NULL
for (pair in unique(ptms_site_pairs_sup_dup_other$pair)) {
  tab_tmp <- ptms_site_pairs_sup_dup_other[ptms_site_pairs_sup_dup_other$pair == pair,]
  tab_tmp2 <- tab_tmp[1,]
  sources <- unlist(str_split(string = tab_tmp$Source, pattern = ";"))
  sources_cat <- paste0(sources, collapse = ";")
  tab_tmp2$Source <- sources_cat
  ptms_site_pairs_sup_dup_other_new <- rbind(ptms_site_pairs_sup_dup_other_new, tab_tmp2)
}

# Sometimes there are multiple scores for the same pair of kinase and substrate phosphosite (one reason is that one kinase can have multiple functional domain which has different preference for the same substrate, which result in different scoring), we chose the larger score.
ptms_site_pairs_sup_dup_networkin <- ptms_site_pairs_sup_dup[ptms_site_pairs_sup_dup$Source == "NetKIN",]
ptms_site_pairs_sup_dup_networkin <- ptms_site_pairs_sup_dup_networkin[order(ptms_site_pairs_sup_dup_networkin$pair, ptms_site_pairs_sup_dup_networkin$networkin_score, decreasing = T),]
ptms_site_pairs_sup_dup_networkin <- ptms_site_pairs_sup_dup_networkin[!duplicated(ptms_site_pairs_sup_dup_networkin$pair),]
which(is.na(ptms_site_pairs_sup_dup_networkin))

ptms_site_pairs_sup_nodup_new <- ptms_site_pairs_sup_dup_other_new
ptms_site_pairs_sup_nodup_new$networkin_score <- NULL
ptms_site_pairs_sup_nodup_new <- merge(ptms_site_pairs_sup_nodup_new, ptms_site_pairs_sup_dup_networkin[, c("GENE", "SUB_GENE", "SUB_MOD_RSD", "pair", "networkin_score")], by = c("GENE", "SUB_GENE", "SUB_MOD_RSD", "pair"), all = T)
which(is.na(ptms_site_pairs_sup_nodup_new$SUB_GENE))

ptms_site_pairs_sup_nodup_new <- rbind(ptms_site_pairs_sup_nodup_new, ptms_site_pairs_sup_nodup)
ptms_site_pairs_sup_nodup_new <- ptms_site_pairs_sup_nodup_new[order(ptms_site_pairs_sup_nodup_new$pair),]
ptms_site_pairs_sup_nodup_new$Source[is.na(ptms_site_pairs_sup_nodup_new$Source)] <- "NetKIN"
which(is.na(ptms_site_pairs_sup_nodup_new$SUB_GENE))

tmp <- as.vector(ptms_site_pairs_sup_nodup_new$networkin_score)
tmp[is.na(tmp)] <- Inf
ptms_site_pairs_sup_nodup_new$networkin_score <- tmp
tmp <- as.vector(ptms_site_pairs_sup_nodup_new$Source)
tmp2 <- tmp[!is.infinite(ptms_site_pairs_sup_nodup_new$networkin_score) & ptms_site_pairs_sup_nodup_new$Source != "NetKIN"]
tmp[!is.infinite(ptms_site_pairs_sup_nodup_new$networkin_score) & ptms_site_pairs_sup_nodup_new$Source != "NetKIN"] <- paste0(tmp2, ";", "NetKIN")
ptms_site_pairs_sup_nodup_new$Source <- tmp
which(is.na(ptms_site_pairs_sup_nodup_new$SUB_GENE))

write.table(x = ptms_site_pairs_sup_nodup_new,
            file = paste0(makeOutDir(resultD = resultD), "omnipath_networkin_enzyme_substrate_site_level_union_source_cleaned.csv"),
            row.names = F, col.names = T, quote = F, sep = ",")


# incorporate depod ------------------------------------------------------------------
# ptms_site_pairs_sup_nodup_new <- fread(input = paste0("./cptac2p/analysis_results/phospho_network/compile_enzyme_substrate/tables/compile_omnipath_networkin_depod/omnipath_networkin_enzyme_substrate_site_level_union_source_cleaned.csv"), data.table = F)
depod_site_pairs2add_cptac <- fread(input = paste0("./cptac2p/analysis_results/preprocess_files/tables/convert_depod_site2cptac/depod_site_pairs2add_cptac.txt"), data.table = F)
depod_site_pairs_in_cptac <- fread(input = paste0("./cptac2p/analysis_results/preprocess_files/tables/convert_depod_site2cptac/depod_site_pairs_in_cptac.txt"), data.table = F)
depod_pro_pairs <- fread(input = paste0("./cptac2p/analysis_results/preprocess_files/tables/convert_depod_site2cptac/depod_pro_pairs.txt"), data.table = F)

sources <- as.vector(ptms_site_pairs_sup_nodup_new$Source)
sources2add_depod <- sources[ptms_site_pairs_sup_nodup_new$pair %in% depod_site_pairs_in_cptac$pair]
sources[ptms_site_pairs_sup_nodup_new$pair %in% depod_site_pairs_in_cptac$pair] <- paste0(sources2add_depod, ";DEPOD")
ptms_site_pairs_sup_nodup_new$Source <- sources
depod_site_pairs2add_cptac <- unique(depod_site_pairs2add_cptac)
ptms_depod <- data.frame(GENE = depod_site_pairs2add_cptac$GENE,
                         SUB_GENE = depod_site_pairs2add_cptac$SUB_GENE,
                         SUB_MOD_RSD = depod_site_pairs2add_cptac$SUB_MOD_RSD,
                         pair = depod_site_pairs2add_cptac$pair,
                         KINASE = depod_site_pairs2add_cptac$Phosphatase_Gene,
                         KIN_ACC_ID = depod_site_pairs2add_cptac$Phosphatase_UniProtAC_human,
                         KIN_ORGANISM = "human",
                         SUBSTRATE = depod_site_pairs2add_cptac$Substrate_Gene,
                         SUB_GENE_ID = "NULL",
                         SUB_ACC_ID = depod_site_pairs2add_cptac$Substrate_UniProtAC_ref,
                         SUB_ORGANISM = "human",
                         SITE_GRP_ID = "NULL",
                         SITE_...7_AA = "NULL",
                         Source = "DEPOD",
                         networkin_score = Inf)

ptms_site_pairs_sup_nodup_wdepod <- rbind(ptms_depod, ptms_site_pairs_sup_nodup_new)
which(is.na(ptms_site_pairs_sup_nodup_new$SUB_GENE))
which(is.na(ptms_site_pairs_sup_nodup_wdepod$SUB_GENE))

write.table(x = ptms_site_pairs_sup_nodup_wdepod,
            file = paste0(makeOutDir(resultD = resultD), "omnipath_networkin_enzyme_substrate_site_level_union_source_cleaned_addDEPOD.csv"),
            row.names = F, col.names = T, quote = F, sep = ",")

ptms_depod_pro <- data.frame(GENE = depod_pro_pairs$GENE,
                         SUB_GENE = depod_pro_pairs$SUB_GENE,
                         SUB_MOD_RSD = "NULL",
                         pair = "NULL",
                         KINASE = depod_pro_pairs$Phosphatase_Gene,
                         KIN_ACC_ID = depod_pro_pairs$Phosphatase_UniProtAC_human,
                         KIN_ORGANISM = "human",
                         SUBSTRATE = depod_pro_pairs$Substrate_Gene,
                         SUB_GENE_ID = "NULL",
                         SUB_ACC_ID = depod_pro_pairs$Substrate_UniProtAC_ref,
                         SUB_ORGANISM = "human",
                         SITE_GRP_ID = "NULL",
                         SITE_...7_AA = "NULL",
                         Source = "DEPOD",
                         networkin_score = Inf)
ptms_site_pairs_sup_nodup_wdepod2 <- rbind(ptms_site_pairs_sup_nodup_wdepod, ptms_depod_pro)
op <- read_csv(paste0(ppnD, "compile_enzyme_substrate/compile_omnipath/omnipath_genename_annotated.csv"))
op$enzyme_type <- NA
op$enzyme_type[op$modification == "dephosphorylation"] <- "phosphatase"
phosphatase1 <- unique(op$GENE[op$enzyme_type == "phosphatase" & !is.na(op$enzyme_type)])
phosphatase2 <- unique(as.vector(ptms_site_pairs_sup$GENE)[grepl(x = ptms_site_pairs_sup$Source, pattern = "DEPOD")])
phosphatases <- unique(c(phosphatase1, phosphatase2))
ptms_site_pairs_sup_nodup_wdepod2$enzyme_type <- ifelse(ptms_site_pairs_sup_nodup_wdepod2$GENE %in% phosphatases, "phosphatase", "kinase")
which(is.na(ptms_site_pairs_sup_nodup_wdepod2$SUB_GENE))

write.table(x = ptms_site_pairs_sup_nodup_wdepod2,
            file = paste0(makeOutDir(resultD = resultD), "omnipath_networkin_enzyme_substrate_site_level_union_source_cleaned_addDEPOD_extended.csv"),
            row.names = F, col.names = T, quote = F, sep = ",")


# Yige Wu @ WashU 2018 Jul
# convert site-level/protein level phosphatase-substrate pairs from DEPOD by matching first uniprot ID then genomic coordinations with CPTAC data

# source -------------------------------------------------------------------
library(data.table)
library(readr)
source('/Users/yigewu/Box Sync/cptac2p_analysis/preprocess_files/preprocess_files.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
## for convert uniprot ID to ensembl peptide ID
source("https://bioconductor.org/biocLite.R")
# biocLite("devtools")
# BiocInstaller::biocLite('grimbough/biomaRt')
library(biomaRt)

## for convert ensembl peptide ID + residue to genomic coord
library(httr)
library(jsonlite)
library(xml2)

# inputs ------------------------------------------------------------------
op <- fread(input = "cptac2p/analysis_results/preprocess_files/tables/convert_omnipath2genomic_coord/op_w.phosphosite.availability.txt", data.table = F)
op <- data.frame(op)

ptms_site_pairs_sup <- read_csv(paste0(ppnD, "compile_enzyme_substrate/tables/compile_omnipath_networkin/omnipath_networkin_enzyme_substrate_site_level_union.csv"))
## input CPTAC gene name with uniprot ID mapping
cptac_uniprots <- fread(input = "./Ding_Lab/Projects_Current/CPTAC/cptac_shared/UniProt.human.20171228.RISnrNF.259contams.fasta.crossRefs.RefSeq.20171003_Human_ucsc_hg38_cpdb_mito_259contamsnr.fasta.categories.tsv", data.table = F)

## input uniprot ID to gene name mapping for depod
pp_uniprots <-  fread(input = paste0(pan3can_shared_dataD, "Phospho_databases/DEPOD/DEPOD_201612_human_phosphatase-protein_substrate_to_Kuan-lin.Phosphatase_UniProtAC_human.txt"), data.table = F)
sub_uniprots_mapped <-  fread(input = paste0(pan3can_shared_dataD, "Phospho_databases/DEPOD/DEPOD_201612_human_phosphatase-protein_substrate_to_Kuan-lin.Substrate_UniProtAC_ref.mapped.txt"), data.table = F)
sub_uniprots_manual <-  fread(input = paste0(pan3can_shared_dataD, "Phospho_databases/DEPOD/DEPOD_201612_human_phosphatase-protein_substrate_to_Kuan-lin.Substrate_UniProtAC_ref.mapped.manual.txt"), data.table = F)
sub_uniprots <- rbind(sub_uniprots_mapped, sub_uniprots_manual)

# input phosphosite genomic coordinates (hg38) -----------------------------------
pho_head_geno_pos <- fread(input = paste0("./cptac2p/analysis_results/preprocess_files/tables/convert_cptac_pho2genomic_coord/pho_head_cancers_single_site_wpid.txt"), data.table = F)

# input 201612 version of depod data --------------------------------------
depod <- fread(input = paste0(pan3can_shared_dataD, "Phospho_databases/DEPOD/DEPOD_201612_human_phosphatase-protein_substrate_to_Kuan-lin.csv"), data.table = F)

# convert uniprot ID to ensembl peptide ID ------------------------------------------------------------------
ensembl = useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl",host="www.ensembl.org")
filters = listFilters(ensembl)

# change depod to one gene-gene pair per line -----------------------------


# annotate uniprot ID with gene name --------------------------------------
## annotate with CPTAC uniprot IDs
depod <- merge(depod, cptac_uniprots[, c("accession_number", "geneSymbol")], by.x = c("Substrate_UniProtAC_ref"), by.y = c("accession_number"), all.x = T)
colnames(depod)[ncol(depod)] <- "geneSymbol_by_cptac_uniprot"

## annotate with uniprot gene name mapping
depod <- merge(depod, pp_uniprots, by.x = c("Phosphatase_UniProtAC_human"), by.y = c("From"), all.x = T)
which(is.na(depod$To))
colnames(depod)[ncol(depod)] <- "GENE"

## kuan's substrate mapping have mouse proteins
sub_geneTab <- 
depod <- merge(depod, sub_uniprots, by.x = c("Substrate_UniProtAC_ref"), by.y = c("From"), all.x = T)
which(is.na(depod$To))
colnames(depod)[ncol(depod)] <- "SUB_GENE.w.nonhuman"

depod_sub_gene_mapped <- depod[!is.na(depod$SUB_GENE),]
depod_sub_gene_na <- depod[is.na(depod$SUB_GENE),]

## manual check for gene names of unmapped ones
depod_sub_gene_na$Substrate_Gene
### RCML not in CPTAC uniprot list
### for tubulin I took all the genes in the tubulin family
tubulin_genes <- unique(cptac_uniprots$geneSymbol[str_sub(as.vector(cptac_uniprots$geneSymbol), start = 1, end = 3) == "TUB"])
num_tubulin_genes <- length(tubulin_genes)
depod_sub_gene_tubulin_orig <- depod_sub_gene_na[depod_sub_gene_na$Substrate_Gene == "Tubulin",]
depod_sub_gene_tubulin <- depod_sub_gene_tubulin_orig[rep(1:nrow(depod_sub_gene_tubulin_orig), vector(mode = "numeric", length = nrow(depod_sub_gene_tubulin_orig)) + num_tubulin_genes),]
depod_sub_gene_tubulin$SUB_GENE.w.nonhuman <- rep(tubulin_genes, nrow(depod_sub_gene_tubulin_orig))
## for mannose 6-phosphate glycoproteins, there's Mannose-6-Phosphate Receptor, but disregard for uncertainty
## for BCR-ABL, split to BCR and ABL1 separately
depod_sub_gene_bcrabl_orig <- depod_sub_gene_na[depod_sub_gene_na$Substrate_Gene == "p210 BCR-ABL",]
depod_sub_gene_bcrabl <- depod_sub_gene_bcrabl_orig[rep(1:nrow(depod_sub_gene_bcrabl_orig), 2),]
depod_sub_gene_bcrabl$SUB_GENE.w.nonhuman <- rep(c("BCR", "ABL1"), vector(mode = "numeric", length = 2) + nrow(depod_sub_gene_bcrabl_orig))

depod <- rbind(depod_sub_gene_tubulin, depod_sub_gene_bcrabl, depod_sub_gene_mapped)

## I kepped all the non-human substrates
depod$SUB_GENE <- toupper(depod$SUB_GENE.w.nonhuman)
depod$is_human <- (as.vector(depod$SUB_GENE) == as.vector(depod$SUB_GENE.w.nonhuman))
depod <- depod[depod$is_human,]
depod <- data.frame(depod)

## annotate whether there's any phosphosite measurement of the substrate
depod$Substrate_in_cptac_pho <- (depod$SUB_GENE %in% pho_head_geno_pos$SUBSTRATE)
table(depod$Substrate_in_cptac_pho)
## ~70% in cptac

## check how many pairs are in DEPOD but not in existing omnipath + networkin
depod$pair_pro <- paste0(depod$GENE, ":", depod$SUB_GENE)
depod_cptac1 <- depod[depod$Substrate_in_cptac_pho,]

depod_pro_pairs <- depod[depod$DephosphoSite == "N/A" | depod$DephosphoSite == "",]
depod_site_pairs <- depod[!(depod$DephosphoSite == "N/A" | depod$DephosphoSite == ""),]

ptms_site_pairs_sup$pair_pro <- paste0(ptms_site_pairs_sup$GENE, ":", ptms_site_pairs_sup$SUB_GENE)
depod$pair_pro[!(depod$pair_pro %in% ptms_site_pairs_sup$pair_pro)]
## only ~10% of the pairs are in the existing table

# ## after observation of the format of DephosphoSite, edit the original data frame
# depod_site_pairs_regular <- depod_site_pairs[!grepl(x = as.vector(depod_site_pairs$DephosphoSite), pattern = "in"),]
# depod_site_pairs_irregular <- depod_site_pairs[grepl(x = as.vector(depod_site_pairs$DephosphoSite), pattern = "in"),]
# depod_site_pairs_irregular_keep <- depod_site_pairs_irregular[depod_site_pairs_irregular$DephosphoSite != "Ser-75, Ser-99 in Q92934; Thr-117 in Q61337",]
# depod_site_pairs_irregular2edit <- depod_site_pairs_irregular[depod_site_pairs_irregular$DephosphoSite == "Ser-75, Ser-99 in Q92934; Thr-117 in Q61337",]
# depod_site_pairs_irregular_edited1 <- depod_site_pairs_irregular2edit; depod_site_pairs_irregular_edited1$DephosphoSite <- "Thr-117 in Q61337"
# depod_site_pairs_irregular_edited2 <- depod_site_pairs_irregular2edit; depod_site_pairs_irregular_edited2$DephosphoSite <- "Ser-75, Ser-99 in Q92934"; depod_site_pairs_irregular_edited2$Substrate_UniProtAC_ref <- "Q92934"
# depod_site_pairs <- rbind(depod_site_pairs_irregular_edited1, depod_site_pairs_irregular_edited2, depod_site_pairs_irregular_keep, depod_site_pairs_regular)

# split site level table and reshape --------------------------------------
depod_site_pairs_new <- NULL
for (i in c(1:nrow(depod_site_pairs))) {
  tab_tmp <- depod_site_pairs[i,]
  sites_string <- as.character(tab_tmp$DephosphoSite)
  sites <- unlist(str_split(sites_string, pattern = ', '))
  sites <- unlist(str_split(sites, pattern = ','))
  sites <- unlist(str_split(sites, pattern = '\\[|\\]'))
  sites <- unlist(str_split(sites, pattern = 'and '))
  sites <- unlist(str_split(sites, pattern = ';'))
  site_last <- sites[length(sites)]
  sites[length(sites)] <- unlist(str_split(string = site_last, pattern = " "))[1]
  site_first <- sites[1]
  site_first <- unlist(str_split(string = site_first, pattern = " "))
  site_first <- site_first[length(site_first)]
  sites[1] <- site_first
  sites <- unlist(str_split(sites, pattern = ' '))
  sites <- sites[sites != ""]
  sites <- unique(sites)
  sites <- str_split_fixed(string = sites, pattern = "-", 2)

  tab_tmp2 <- tab_tmp[vector(mode = "numeric", length = nrow(sites))+1,]
  tab_tmp2$amino_acid <- sites[,1]
  tab_tmp2$offset <- sites[,2]
  depod_site_pairs_new <- rbind(depod_site_pairs_new, tab_tmp2)
}
table(depod_site_pairs_new$amino_acid)

residues <- vector(mode = "character", length = nrow(depod_site_pairs_new))
residues[depod_site_pairs_new$amino_acid == "Ser"] <- "S"
residues[depod_site_pairs_new$amino_acid == "Thr"] <- "T"
residues[depod_site_pairs_new$amino_acid == "Tyr"] <- "Y"
residues[depod_site_pairs_new$amino_acid== "His"] <- "H"
depod_site_pairs_new$residue <- residues
depod_site_pairs_new$SUB_MOD_RSD <- paste0(depod_site_pairs_new$residue, depod_site_pairs_new$offset)
write.table(x = depod_site_pairs_new,
            file = paste0(makeOutDir(resultD = resultD), "depod_site_pairs_new.txt"),
            row.names = F, col.names = T, quote = F, sep = "\t")

# determine which pairs to add ot the original table ----------------------
# depod_pro_pairs2add <- depod_pro_pairs[!(depod_pro_pairs$pair_pro %in% ptms_site_pairs_sup$pair_pro),]
depod_pro_pairs2add <- depod_pro_pairs
## add around 500 pairs, 400 new pairs

depod_site_pairs_new$pair <- paste0(depod_site_pairs_new$pair_pro, ":", depod_site_pairs_new$SUB_MOD_RSD)
ptms_site_pairs_sup$pair <- paste0(ptms_site_pairs_sup$pair_pro, ":", ptms_site_pairs_sup$SUB_MOD_RSD)
depod_site_pairs2add <- depod_site_pairs_new[depod_site_pairs_new$Substrate_in_cptac_pho & !(depod_site_pairs_new$pair %in% ptms_site_pairs_sup$pair),]

## annoate the site
pho_head_geno_pos$site <- paste0(pho_head_geno_pos$SUBSTRATE, ":", pho_head_geno_pos$SUB_MOD_RSD)
depod_site_pairs2add$site <- paste0(depod_site_pairs2add$SUB_GENE, ":", depod_site_pairs2add$SUB_MOD_RSD)
depod_site_pairs2add$site_original_in_cptac_pho <- (depod_site_pairs2add$site %in% pho_head_geno_pos$site)
depod2add_sites <- unique(depod_site_pairs2add[, c("SUB_GENE", "offset", "Substrate_UniProtAC_ref", "residue")])

# check whether the uniprot ID of the substrate is the same as cpt --------
## input omnipath+networkin table with genomic coordinations
op_substrate_siteTab_wpid <- fread(input = "cptac2p/analysis_results/preprocess_files/tables/convert_omnipath2genomic_coord/op_substrate_siteTab_wpid.txt", data.table = F)
op_substrate_pids <- unique(op_substrate_siteTab_wpid[, c("SUB_GENE", "matched_ens_pid")])

depod2add_sites_wpid <- merge(depod2add_sites, op_substrate_pids, all.x = T)
depod2add_sites_wpid_mapped <- depod2add_sites_wpid[!is.na(depod2add_sites_wpid$matched_ens_pid),]
depod_unmapped_uniprots <- unique(depod2add_sites_wpid$Substrate_UniProtAC_ref[is.na(depod2add_sites_wpid$matched_ens_pid)])

## convert substrate uniprot IDs to ensembl peptide ID
depod_unmapped_pids = getBM(attributes = c("ensembl_peptide_id", "uniprotswissprot", "hgnc_symbol"),
               filters = "uniprotswissprot",
               values = depod_unmapped_uniprots, mart = ensembl, uniqueRows=T)

depod2add_sites_wpid_unmapped <- depod2add_sites_wpid[is.na(depod2add_sites_wpid$matched_ens_pid),]
depod2add_sites_wpid_unmapped$matched_ens_pid <- NULL
depod2add_sites_wpid_unmapped <- merge(depod2add_sites_wpid_unmapped, depod_unmapped_pids[, c("uniprotswissprot", "ensembl_peptide_id")], by.x = c("Substrate_UniProtAC_ref"), by.y = c("uniprotswissprot"), all.x = T)
## some NAs are due residual non-human uniprot IDs
depod2add_sites_wpid_unmapped <- depod2add_sites_wpid_unmapped[!is.na(depod2add_sites_wpid_unmapped$ensembl_peptide_id),]
colnames(depod2add_sites_wpid_unmapped)[ncol(depod2add_sites_wpid_unmapped)] <- "matched_ens_pid"
depod2add_sites_wpid <- rbind(depod2add_sites_wpid_mapped, depod2add_sites_wpid_unmapped[, colnames(depod2add_sites_wpid_unmapped)])
colnames(depod2add_sites_wpid)[ncol(depod2add_sites_wpid)] <- "ensembl_peptide_id"

## get the peptide sequence and length for each ensembl peptide IDs, sort by length
op_substrate_uniprots_mapTab <- fread(input = "cptac2p/analysis_results/preprocess_files/tables/convert_omnipath2genomic_coord/op_substrate_uniprots_mapTab.txt", data.table = F)
op_substrate_pid_seq <- unique(op_substrate_uniprots_mapTab[, c("ensembl_peptide_id", "peptide_seq", "peptide_len")])

depod2add_sites_wpid <- merge(depod2add_sites_wpid, op_substrate_pid_seq, by = c("ensembl_peptide_id"), all.x = T)
depod2add_sites_wpid_wseq <- depod2add_sites_wpid[!is.na(depod2add_sites_wpid$peptide_seq),]
depod2add_sites_wpid_nseq <- depod2add_sites_wpid[is.na(depod2add_sites_wpid$peptide_seq),]

peptide_seq <- sapply(as.vector(depod2add_sites_wpid_nseq$ensembl_peptide_id), FUN = function(pid) {
  r <- GET(paste("http://rest.ensembl.org", paste0("/sequence/id/", pid, "?"), sep = ""), content_type("application/json"))
  stop_for_status(r)
  seq <- fromJSON(toJSON(content(r)))$seq
  print(pid)
  return(seq)
})
depod2add_sites_wpid_nseq$peptide_seq <- peptide_seq
depod2add_sites_wpid_nseq$peptide_len <- sapply(as.vector(depod2add_sites_wpid_nseq$peptide_seq), FUN = function(s) {
  len <- nchar(s)
  return(len)
})
depod2add_sites_wpid_wseq <- rbind(depod2add_sites_wpid_wseq, depod2add_sites_wpid_nseq)
depod2add_sites_wpid_wseq <- depod2add_sites_wpid_wseq[order(depod2add_sites_wpid_wseq$SUB_GENE, depod2add_sites_wpid_wseq$peptide_len, decreasing = T),]
write.table(x = depod2add_sites_wpid_wseq,
            file = paste0(makeOutDir(resultD = resultD), "depod2add_sites_wpid_wseq.txt"),
            row.names = F, col.names = T, quote = F, sep = ",")

## do another round of matching of the suitable peptide id
matched_ens_pid <- vector(mode = "character", length = nrow(depod2add_sites_wpid_wseq))
for (i in 1:nrow(depod2add_sites_wpid_wseq)) {
  gene <- depod2add_sites_wpid_wseq$SUB_GENE[i]
  rsd <- depod2add_sites_wpid_wseq$residue[i]
  offset <- depod2add_sites_wpid_wseq$offset[i]
  ens_pids_rows <- which(as.vector(depod2add_sites_wpid_wseq$SUB_GENE) == gene)
  for (j in ens_pids_rows) {
    rsd_j <- substr(x = as.vector(depod2add_sites_wpid_wseq$peptide_seq)[j], start = offset, stop = offset)
    if (!is.na(rsd_j) & rsd_j == rsd) {
      matched_ens_pid[i] <- as.vector(depod2add_sites_wpid_wseq$ensembl_peptide_id)[j]
      break
    }
  }
}
depod2add_sites_wpid_wseq$matched_ens_pid <- matched_ens_pid
depod2add_sites_wpid_wseq_rsd_mapped <- depod2add_sites_wpid_wseq
depod2add_sites_wpid_wseq_rsd_mapped <- depod2add_sites_wpid_wseq_rsd_mapped[depod2add_sites_wpid_wseq_rsd_mapped$matched_ens_pid != "",]
depod2add_sites_wpid_wseq_rsd_mapped$ensembl_peptide_id <- NULL
depod2add_sites_wpid_wseq_rsd_mapped$matched_ens_pid_new <- NULL
depod2add_sites_wpid_wseq_rsd_mapped$peptide_seq <- NULL
depod2add_sites_wpid_wseq_rsd_mapped$peptide_len <- NULL

coords <- lapply(1:nrow(depod2add_sites_wpid_wseq_rsd_mapped), FUN = function(n, pids, offsets) {
  r <- GET(paste("http://rest.ensembl.org", paste0("/map/translation/", pids[n], "/", offsets[n], "..", offsets[n], "?"), sep = ""), content_type("application/json"))
  stop_for_status(r)
  tmp <- fromJSON(toJSON(content(r)))$mappings
  print(pids[n])
  return(tmp)
}, pids = as.vector(depod2add_sites_wpid_wseq_rsd_mapped$matched_ens_pid), offsets = as.vector(depod2add_sites_wpid_wseq_rsd_mapped$offset))

## TODO: resolve incorrect coordinates, e.g. i = 19
df <- matrix(NA, nrow=length(coords), ncol = ncol(coords[[1]]))
for (i in 1:length(coords)) {
  if (nrow(coords[[i]]) == 1) {
    df[i,] <- matrix(unlist(coords[[i]]), nrow = 1, byrow = T)
  }
}
df <- data.frame(df,stringsAsFactors=FALSE)
colnames(df) <- colnames(coords[[1]])
depod2add_sites_wpid_wseq_rsd_mapped <- cbind(depod2add_sites_wpid_wseq_rsd_mapped, df)

# edit the residue offset according to the genomic-matched cptac phosphosites ----------------
residue_offset_cptac <- vector(mode = "numeric", length = nrow(depod2add_sites_wpid_wseq_rsd_mapped)) + NA
for (i in 1:nrow(depod2add_sites_wpid_wseq_rsd_mapped)) {
  matched_row <- which(pho_head_geno_pos$seq_region_name == depod2add_sites_wpid_wseq_rsd_mapped$seq_region_name[i] & pho_head_geno_pos$start == depod2add_sites_wpid_wseq_rsd_mapped$start[i])
  if (length(matched_row) > 0) {
    residue_offset_cptac[i] <- pho_head_geno_pos$offset[matched_row[1]]
  }
}
depod2add_sites_wpid_wseq_rsd_mapped$residue_offset_cptac <- residue_offset_cptac
write.table(x = depod2add_sites_wpid_wseq_rsd_mapped,
            file = paste0(makeOutDir(resultD = resultD), "depod2add_sites_wpid_wseq_rsd_mapped.txt"),
            row.names = F, col.names = T, quote = F, sep = "\t")


depod_site_pairs2add <- merge(depod_site_pairs2add, depod2add_sites_wpid_wseq_rsd_mapped[, c("SUB_GENE", "residue", "offset", "residue_offset_cptac")], by = c("SUB_GENE", "residue", "offset"))
depod_site_pairs2add <- data.frame(depod_site_pairs2add)
depod_site_pairs2add_cptac <- depod_site_pairs2add[!is.na(depod_site_pairs2add$residue_offset_cptac),]


# write -------------------------------------------------------------------
write.table(x = depod_pro_pairs,
            file = paste0(makeOutDir(resultD = resultD), "depod_pro_pairs.txt"),
            row.names = F, col.names = T, quote = F, sep = "\t")

## merge
depod_site_pairs2add_cptac$offset <- depod_site_pairs2add_cptac$residue_offset_cptac
depod_site_pairs2add_cptac$residue_offset_cptac <- NULL
depod_site_pairs2add_cptac$SUB_MOD_RSD <- paste0(depod_site_pairs2add_cptac$residue, depod_site_pairs2add_cptac$offset)
write.table(x = depod_site_pairs2add_cptac,
            file = paste0(makeOutDir(resultD = resultD), "depod_site_pairs2add_cptac.txt"),
            row.names = F, col.names = T, quote = F, sep = "\t")

depod_site_pairs_in_cptac <- depod_site_pairs_new[depod_site_pairs_new$Substrate_in_cptac_pho & (depod_site_pairs_new$pair %in% ptms_site_pairs_sup$pair),]
write.table(x = depod_site_pairs_in_cptac,
            file = paste0(makeOutDir(resultD = resultD), "depod_site_pairs_in_cptac.txt"),
            row.names = F, col.names = T, quote = F, sep = "\t")

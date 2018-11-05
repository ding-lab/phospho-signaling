# Yige Wu @ WashU 2018 Jan
## shared parameters and functions for phospho_network

# library -----------------------------------------------------------------
library(stringr)
library(reshape)
library(grid)
require(plyr)
library(data.table)
library(readr)
library(rstudioapi)


# source ------------------------------------------------------------------
# folders <- strsplit(x = rstudioapi::getSourceEditorContext()$path, split = "\\/")[[1]]
# folder_num <- which(folders == "cptac2p_analysis") -1
# baseD <- paste0(strsplit(paste(folders[1:folder_num], collapse = "/"), split = "\\.")[[1]][1], "/")
# setwd(baseD)

folders <- strsplit(x = rstudioapi::getSourceEditorContext()$path, split = "\\/")[[1]]
baseD <- paste0(paste0(folders[1:which(folders == "Box Sync")], collapse = "/"), "/")
setwd(baseD)
source("./cptac2p_analysis/cptac2p_analysis_shared.R")

# static parameters -------------------------------------------------------
ppnD <- paste0(resultD, "phospho_network/")

## significance level for FDR from regression model
fdr_pp <- 0.1
fdr_pk <- 0.05

reg_sig <- c(fdr_pp, fdr_pk)
names(reg_sig) <- c("phosphatase", "kinase")

# functions ---------------------------------------------------------------

makeOutDir = function(resultD) {
  folders <- strsplit(x = rstudioapi::getSourceEditorContext()$path, split = "\\/")[[1]]
  folder_num <- which(folders == "cptac2p_analysis") + 1
  resultDnow <- paste(strsplit(paste(folders[folder_num:length(folders)], collapse = "/"), split = "\\.")[[1]][1], sep = "/")
  resultDnow <- paste0(resultD, resultDnow, "/")
  dir.create(resultDnow)
  resultDnow_son <- resultDnow
  dirs2make <- NULL
  while (!dir.exists(resultDnow_son)) {
    tmp <- strsplit(resultDnow_son, split = "\\/")[[1]]
    resultDnow_parent <-paste(tmp[-length(tmp)], collapse = "/")
    dir.create(resultDnow_parent)
    dir.create(resultDnow_son)
    dir.create(resultDnow)
    if (!dir.exists(resultDnow_son)) {
      dirs2make[length(dirs2make) + 1] <- resultDnow_son
    }
    resultDnow_son <- resultDnow_parent
  }
  
  if (length(dirs2make) > 0){
    for (i in 1:length(dirs2make)) {
      dir.create(dirs2make[i])
    }
  } 
  return(resultDnow)
}

sampID2partID = function(sampleID_vector, sample_map) {
  sampleIDs = sampleID_vector
  partIDs <- NULL
  for (sampID in sampleIDs) {
    partID <- as.character(unique(sample_map$Participant.ID[grepl(x = sample_map$Specimen.Label, pattern = sampID) & !is.na(sample_map$Specimen.Label) & !is.na(sample_map$Participant.ID)]))
    if (!is.null(partID)) {
      partIDs <- c(partIDs, partID)
    } else {
      partIDs <- c(partIDs, NA)
    }
  }
  return(partIDs)
}

sampID2clinicalMSI = function(sampleID_vector, sample_map, subtype_map) {
  sampleIDs = sampleID_vector
  partIDs <- sampID2partID(sampleIDs, sample_map)
  subtypes <- NULL
  for (partID in partIDs) {
    subtype <- unique(as.character(subtype_map$`Microsatellite Instability (Abnormal @ >33% loci tested)`[subtype_map$`Participant ID` == partID]))
    subtype <- subtype[!is.na(subtype)]
    if (subtype == 'Yes') {
      subtypes <- c(subtypes, "MSI")
    } else if (subtype == "No") {
      subtypes <- c(subtypes, "MSS")
    } else {
      subtypes <- c(subtypes, 'other')
    }
  }
  return(subtypes)
}

sampID2MSI = function(sampleID_vector, sample_map, subtype_map) {
  sampleIDs = sampleID_vector
  partIDs <- sampID2partID(sampleIDs, sample_map)
  subtypes <- NULL
  for (partID in partIDs) {
    subtype <- unique(as.character(subtype_map$MSI_type[subtype_map$Sample == partID]))
    subtype <- subtype[!is.na(subtype)]
    if (length(subtype) > 0) {
      subtypes <- c(subtypes, subtype)
    } else {
      subtypes <- c(subtypes, "no_genomic_data")
    }
  }
  return(subtypes)
}

sampID2pam50fromclinTab = function(sampleID_vector, subtype_map) {
  sampleIDs = sampleID_vector
  subtypes <- NULL
  for (sampID in sampleIDs) {
    subtype <- unique(as.character(subtype_map$`PAM50 Call`[grepl(x = subtype_map$`Specimen Label`, pattern = sampID)]))
    subtype <- subtype[!is.na(subtype)]
    if (length(subtype) > 0) {
      subtypes <- c(subtypes, subtype)
    } else {
      subtypes <- c(subtypes, NA)
    }
  }
  return(subtypes)
}

sampID2pam50 = function(sampleID_vector, pam50_map, sample_map) {
  sampleIDs = sampleID_vector
  partIDs <- sampID2partID(sampleIDs, sample_map)
  subtypes <- NULL
  for (partID in partIDs) {
    subtype <- unique(as.character(pam50_map$Call[pam50_map$V1 == partID]))
    subtype <- subtype[!is.na(subtype)]
    if (length(subtype) > 0) {
      subtypes <- c(subtypes, subtype)
    } else {
      subtypes <- c(subtypes, "no_RNA-seq")
    }
  }
  return(subtypes)
}

#function for normalize the variables for regression
range01 <- function(x, ...){(x - min(x, ...)) / (max(x, ...) - min(x, ...))}

markVadRetro = function(regression, sig_thres, enzyme_type) {
  ## Usage: input the regression result table and the significance level
  ## mark table with retrospective validated results (FDR < significance level and coefficient > 0)
  retro_table <- fread(input = "~/Box Sync/pan3can_shared_data/analysis_results/regression/table/kinase_substrate_regression_trans_edited.txt",
                       data.table = F)
  retro_table <- markSigKS(retro_table, sig_thres = sig_thres, enzyme_type = enzyme_type)
  retro_table_vad <- retro_table[retro_table$fdr_sig & retro_table$coef_sig,]
  table <- markSigKS(regression, sig_thres = sig_thres, enzyme_type = enzyme_type)
  table[, "vad_retro_pairwise"] <- (paste(table$KINASE, table$SUBSTRATE, sep = ":") %in% paste(retro_table_vad$KINASE, retro_table_vad$SUBSTRATE, sep = ":"))
  table[, "vad_retro_sitewise"] <- (table$pair %in% retro_table_vad$pair)
  return(table)
}

markSigKS = function(regression, sig_thres, enzyme_type) {
  ## Usage: input the regression result table and the significance level
  ## mark table with significant cis and trans pairs
  table <- regression
  table_cis <- table[table$SELF == "cis",]
  table_cis$fdr_sig <- (table_cis$FDR_pro_kin < sig_thres)
  table_trans <- table[table$SELF == "trans",]
  table_trans$fdr_sig <- (table_trans$FDR_pho_kin < sig_thres)
  if (enzyme_type == 'kinase') {
    table_cis$coef_sig <- (table_cis$coef_pro_kin > 0)
    table_trans$coef_sig <- (table_trans$coef_pho_kin > 0)
  } else if (enzyme_type == "phosphatase") {
    table_cis$coef_sig <- (table_cis$coef_pro_kin < 0)
    table_trans$coef_sig <- (table_trans$coef_pho_kin < 0)
  }
  table <- as.data.frame(rbind(table_cis, table_trans))
  return(table)
}

markSigCan = function(regression, sig_thres, enzyme_type) {
  ## Usage: input the regression result table and the significance level
  ## mark table with whether showing significant in other cancer types
  table <- regression
  table <- markSigKS(table, sig_thres = sig_thres, enzyme_type = enzyme_type)
  table$ks_pair <- paste(table$KINASE, table$SUBSTRATE, sep = ':')
  
  sig_pairs <- vector('list', 3)
  sig_ks_pairs <- vector('list', 3)
  for (cancer in unique(table$Cancer)) {
    table_can <- table[table$Cancer == cancer & table$fdr_sig,]
    sig_pairs[[cancer]] <- as.vector(table_can$pair)
    sig_ks_pairs[[cancer]] <- as.vector(table_can$ks_pair)
  }
  for (cancer in unique(table$Cancer)) {
    table[, paste0('sig_', cancer)] <- (table$pair %in% sig_pairs[[cancer]])
    table[, paste0('ks_sig_', cancer)] <- (table$ks_pair %in% sig_ks_pairs[[cancer]])
  }
  
  table$shared3can <- (table$ks_sig_BRCA & table$ks_sig_OV & table$ks_sig_CO)
  table$uniq_BRCA <- (table$ks_sig_BRCA & !table$ks_sig_OV & !table$ks_sig_CO)
  table$uniq_OV <- (!table$ks_sig_BRCA & table$ks_sig_OV & !table$ks_sig_CO)
  table$uniq_CO <- (!table$ks_sig_BRCA & !table$ks_sig_OV & table$ks_sig_CO)
  return(table)
}


markSigSiteCan = function(regression, sig_thres, enzyme_type) {
  ## Usage: input the regression result table and the significance level
  ## mark table with whether showing significant in other cancer types
  if ('Cancer' %in% colnames(regression)) {
    table <- regression
    table$pair <- paste0(table$GENE, ":", table$SUB_GENE, ":", table$SUB_MOD_RSD)
    table <- markSigKS(table, sig_thres = sig_thres, enzyme_type = enzyme_type)
    sig_pairs <- vector('list', 3); examined_pairs <- vector('list', 3)
    for (cancer in unique(table$Cancer)) {
      table_can <- table[table$Cancer == cancer,]
      examined_pairs[[cancer]] <- as.vector(table_can$pair[!is.na(table_can$fdr_sig)])
      table_can_sig <- table[table$Cancer == cancer & table$fdr_sig & table$coef_sig,]
      sig_pairs[[cancer]] <- as.vector(table_can_sig$pair)
    }
    for (cancer in unique(table$Cancer)) {
      table[, paste0('sig_', cancer)] <- NA
      examined <- table$pair %in% examined_pairs[[cancer]]
      table[examined, paste0('sig_', cancer)] <- (table$pair[examined] %in% sig_pairs[[cancer]])
    }
    
    ## ks pairs need to be examined in all 3 cancer before marking them as shared or unique
    examined3can <- (table$pair %in% examined_pairs[["BRCA"]]) & (table$pair %in% examined_pairs[["OV"]]) & (table$pair %in% examined_pairs[["CO"]])
    table$shared3can <-  NA; table$shared3can[examined3can] <- (table$sig_BRCA[examined3can] & table$sig_OV[examined3can] & table$sig_CO[examined3can])
    table$uniq_BRCA <-  NA; table$uniq_BRCA[examined3can] <- (table$sig_BRCA[examined3can] & !table$sig_OV[examined3can] & !table$sig_CO[examined3can])
    table$uniq_OV <- NA; table$uniq_OV[examined3can] <- (!table$sig_BRCA[examined3can] & table$sig_OV[examined3can] & !table$sig_CO[examined3can])
    table$uniq_CO <-NA; table$uniq_CO[examined3can]  <- (!table$sig_BRCA[examined3can]  & !table$sig_OV[examined3can]  & table$sig_CO[examined3can] )
    return(table)
  } else {
    print('no column named Cancer!')
  }
}

remove_outliers <- function(x, out_thres, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- out_thres * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}

formatPhosphosite = function(phosphosite_vector, gene_vector) {
  data <- phosphosite_vector
  gene <- gene_vector
  # pho_rsd_split <- data.frame(SUB_MOD_RSD = sapply(as.vector(data), FUN = function(p) strsplit(x = p, split = ":")[[1]][1]))
  pho_rsd_split <- as.data.frame(str_split_fixed(data, ":", 2))
  colnames(pho_rsd_split) <- c("transcript", "SUB_MOD_RSD")
  pho_rsd_split$SUB_MOD_RSD <- toupper(pho_rsd_split$SUB_MOD_RSD)
  pho_rsd_split$SUBSTRATE <- gene
  return(pho_rsd_split)
}

load_ks_table <- function(protein) {
  ## Usage: k_s_table <- load_ks_table({kinase/phosphatase})
  if ( protein == "kinase" ) {
    ### read in the kinase/substrate table/ phosphorylation data ###
    k_s_table = read.delim(paste(pan3can_shared_dataD,"Phospho_databases/PhosphositePlus/data/Kinase_Substrate_Dataset_human_final_hugoified.txt",sep=""))
    # k_s_table = read_table2(paste0(pan3can_shared_dataD,"Phospho_databases/PhosphositePlus/Aug_01_2018/Kinase_Substrate_Dataset"),  skip = 3)
  }
  
  if ( protein == "phosphatase" ) {
    ### read in the phosphatase/substrate table/ phosphorylation data ### 
    k_s_table <- read.csv(paste(pan3can_shared_dataD,"Phospho_databases/DEPOD/DEPOD_201612_human_phosphatase-protein_substrate_to_Kuan-lin.csv",sep = ""))
    colnames(k_s_table) <- c("Phosphatase_UniProtAC_human","GENE","Substrate_UniProtAC_ref","SUB_GENE","Substrate_Type","DephosphoSite","BioassayType", "PubMed_ID_rev")
  }
  return(k_s_table)
}


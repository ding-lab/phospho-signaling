# Yige Wu @ WashU Jan 2018
# functions and variables useful for preprocessing CPTAC CDAP data
## TODO: check for multiple samples (tumor or normal) for the same participant


# load packages -----------------------------------------------------------
library(limma)
library(readxl)
library(readr)
# library(dplyr)
library(stringr)
library(data.table)
library(WGCNA)

# inputs ------------------------------------------------------------------
setwd("/Users/yigewu/Box Sync/cptac2p_analysis/preprocess_files/")
baseD = "/Users/yigewu/Box\ Sync/"
inputD <- "/Users/yigewu/Box Sync/cptac2p/cptac_shared/"
resultD <- paste0(baseD, "analysis_results/preprocess_files/")
base_ov_pnnl <- paste(inputD, "4_CPTAC2_Ovarian_Prospective_Collection_PNNL", sep = "")
base_co_pnnl <- paste(inputD, "2_CPTAC2_Colon_Prospective_Collection_PNNL", sep = "")
base_brca_bi <- paste(inputD, "5_CPTAC2_Breast_Prospective_Collection_BI", sep = "")
base_ov_jhu <- paste(inputD, "3_CPTAC2_Ovarian_Prospective_Collection_JHU", sep = "")

prefix <- NULL
prefix["BRCA"] <- "BRCA_BI"
prefix["OV"] <- "OV_PNNL"
prefix["CO"] <- "CO_PNNL"
prefix["UCEC"] <- "UCEC_PNNL"
prefix["CCRCC"] <- "CCRCC_JHU"

raw <- matrix(NA, nrow = 5, ncol = 2, dimnames = list(row = c("BRCA","OV","CO", "UCEC", "CCRCC"),
                                                      col = c("protein", "phosphoprotein")))
raw["BRCA","protein"] <- paste(base_brca_bi,"/SummaryReports/CPTAC2_Breast_Prospective_Collection_BI_Proteome.tmt10.tsv", sep = "")
raw["BRCA", "phosphoprotein"] <- paste(base_brca_bi,"/SummaryReports/CPTAC2_Breast_Prospective_Collection_BI_Phosphoproteome.phosphosite.tmt10.tsv", sep = "")
raw["OV","protein"] <- paste(base_ov_pnnl,"/SummaryReports/CPTAC2_Ovarian_Prospective_Collection_PNNL_Proteome.tmt10.tsv", sep = "")
raw["OV", "phosphoprotein"] <- paste(base_ov_pnnl,"/SummaryReports/CPTAC2_Ovarian_Prospective_Collection_PNNL_Phosphoproteome.phosphosite.tmt10.tsv", sep = "")
raw["CO","protein"] <- paste(base_co_pnnl,"/SummaryReports/CPTAC2_Colon_Prospective_Collection_PNNL_Proteome.tmt10.tsv", sep = "")
raw["CO", "phosphoprotein"] <- paste(base_co_pnnl,"/SummaryReports/CPTAC2_Colon_Prospective_Collection_PNNL_Phosphoproteome.phosphosite.tmt10.tsv", sep = "")

datatypes <- c("PRO", "PHO", "collapsed_PHO")


# Functions ---------------------------------------------------------------
sortColumn = function(expression) {
  data <- expression
  idcols <- colnames(data)[colnames(data) %in% c("Gene", "Phosphosite")]
  samplecols <- sort(colnames(data)[!(colnames(data) %in% c("Gene", "Phosphosite"))])
  data <- data[,c(idcols, samplecols)]
  return(data)
}

unshared_pro = function(Pro.m){
  Pro <- Pro.m
  Pro = Pro[!(Pro$Gene %in% c('Mean', 'Median', 'StdDev')),]
  row.names(Pro) = make.names(Pro$Gene, unique=T)
  Pro = Pro[,grepl("Unshared", colnames(Pro))]
  colnames(Pro) = sub(".Unshared.Log.Ratio","",colnames(Pro))
  colnames(Pro) = sub(" Unshared Log Ratio","",colnames(Pro))
  colnames(Pro) = sub("^X","",colnames(Pro))
  print('-Extracting unshared log ratio for proteome data done!')
  Pro <- sortColumn(Pro)
  return(Pro)
}

aveRep = function(expression, clinical.m){
  print('-Averaging biological replicates')
  data <- expression
  clin <- clinical.m
  
  # average duplicates, while one parent specimen label per participant
  label <- str_split_fixed(colnames(data), "\\.", 2)
  dup_labels <- label[duplicated(label[,1]),1]
  if (length(dup_labels) > 0) {
    print('--Biological duplicate samples:')
    print(colnames(data)[label[,1] %in% dup_labels])
    dup_labels <- unique(dup_labels)
    dup_ave <- NULL
    for (i in dup_labels) {
      dup <- colnames(data)[label[,1] == i]
      dupave <- rowMeans(data[,dup], na.rm = T)
      dup_ave <- cbind(dup_ave, dupave)
    }
    colnames(dup_ave) <- dup_labels
    data_nodup <- data[,!(label[,1] %in% dup_labels)]
    data <- cbind(dup_ave, data_nodup)
  } else {
    print('no biological duplicate samples')
  }
  data <- sortColumn(data)
  return(data)
}

k1SamplePerPatient = function(expression, clinical.m){
  ## Usage: expression data have to be tumor-only or normal only samples
  print('-keep one tumor sample for multiple tumor samples of one patient or multiple one normal samples for one patient')
  data <- expression
  clin <- clinical.m
  
  label <- str_split_fixed(colnames(data), "\\.", 2)
  dup_labels <- label[duplicated(label[,1]),1]
  if (length(dup_labels) > 0) {
    print('--averaging samples:')
    print(colnames(data)[label[,1] %in% dup_labels])
    dup_labels <- unique(dup_labels)
    dup_ave <- NULL
    for (i in dup_labels) {
      dup <- colnames(data)[label[,1] == i]
      dupave <- rowMeans(data[,dup], na.rm = T)
      dup_ave <- cbind(dup_ave, dupave)
    }
    colnames(dup_ave) <- dup_labels
    data_nodup <- data[,!(label[,1] %in% dup_labels)]
    data <- cbind(dup_ave, data_nodup)
  } else {
    print('no multiple tumor/normal samples for the same patient')
  }
  data <- sortColumn(data)
  return(data)
}


get_tumor_raw = function(expression, clinical.m){
  ## Usage: for extracting tumor samples for proteomics data using unprocessed clinical info file
  print('-Extracting tumor samples from expression data')
  data <- expression
  clin <- clinical.m
  
  # get tumor only parent specimen label in clinical data
  tumor_p <- clin$Parent.Specimen.Label[grepl("Tumor", clin$Parent.Specimen.Title)]
  
  # get parent specimen label in expression data
  p <-   str_split_fixed(colnames(data), "_", 2)[,1]
  
  # keep tumor only specimens
  data.t <- data[,(p %in% tumor_p)]
  
  # check multiple samples for one participant
  participant_ids <- NULL
  for (i in 1:ncol(data.t)) {
    specimen_label <- colnames(data.t)[i]
    participant_id <- clin$Participant.ID[clin$Specimen.Label == specimen_label]
    participant_ids <- c(participant_ids, as.character(participant_id))
  }
  participant_ids_dup <- participant_ids[duplicated(participant_ids)]
  if (length(participant_ids_dup) > 0) {
    print(paste0('--multiple samples for one participant: ', participant_ids_dup))
  } else {
    print('--no multiple samples for one participant')
  }
  
  data.t <- sortColumn(data.t)
  return(data.t)
}

get_tumor = function(expression, clinical.m){
  ## Usage: for extracting tumor samples for proteomics data using processed clinical info file (with tumor_normal column)
  print('-Extracting tumor samples from expression data')
  data <- expression
  clin <- clinical.m
  
  # get tumor only parent specimen label in clinical data
  tumor_p <- clin$Parent.Specimen.Label[grepl("Tumor", clin$tumor_normal)]
  
  # get parent specimen label in expression data
  p <-   str_split_fixed(colnames(data), "_", 2)[,1]
  
  # keep tumor only specimens
  data.t <- data[,colnames(data)[p %in% tumor_p]]
  
  # check multiple samples for one participant
  participant_ids <- NULL
  for (i in 1:ncol(data.t)) {
    specimen_label <- colnames(data.t)[i]
    participant_id <- clin$Participant.ID[clin$Specimen.Label == specimen_label]
    participant_id <- unique(participant_id[!is.na(participant_id)])
    participant_ids <- c(participant_ids, as.character(participant_id))
  }
  participant_ids_dup <- participant_ids[duplicated(participant_ids)]
  participant_ids_dup <- participant_ids_dup[!is.na(participant_ids_dup)]
  if (length(participant_ids_dup) > 0) {
    print(paste0('--multiple samples for one participant: ', participant_ids_dup))
    print('---keep one tumor sample with the most non-NA values for multiple tumor samples of one patient')
    for (partID in participant_ids_dup) {
      print(partID)
      print(colnames(data.t))
      print(participant_ids)
      sampIDs <- colnames(data.t)[participant_ids == partID]
      print(sampIDs)
      nonNAnum <- colSums(!is.na(data.t[, sampIDs]))
      sampID <- sampIDs[which.max(nonNAnum)]
      print(paste0('----keep ', sampID, " for patient ", partID))
      ## take out unselected samples
      data.t <- data.t[, colnames(data.t)[!(colnames(data.t) %in% sampIDs[sampIDs %in% c(sampID)])]]
    }
  } else {
    print('--no multiple samples for one participant')
  }
  
  data.t <- sortColumn(data.t)
  return(data.t)
}

get_normal_raw = function(expression, clinical.m){
  print('-Extracting normal samples from expression data')
  data <- expression
  clin <- clinical.m
  
  # get normal only parent specimen label in clinical data
  normal_p <- clin$Parent.Specimen.Label[grepl("Normal", clin$Parent.Specimen.Title)]
  
  # get parent specimen label in expression data
  p <-   str_split_fixed(colnames(data), "_", 2)[,1]
  
  # keep normal only specimens and no replicates
  data.n <- data[,(p %in% normal_p)]
  # data.n <- data.frame(data.n)
  
  # check multiple samples for one participant
  participant_ids <- NULL
  for (i in 1:ncol(data.n)) {
    specimen_label <- colnames(data.n)[i]
    participant_id <- clin$Participant.ID[clin$Specimen.Label == specimen_label]
    participant_ids <- c(participant_ids, as.character(participant_id))
  }
  participant_ids_dup <- participant_ids[duplicated(participant_ids)]
  if (length(participant_ids_dup) > 0) {
    print(paste0('--multiple samples for one participant: ', participant_ids_dup))
  } else {
    print('--no multiple samples for one participant')
  }
  data.n <- sortColumn(data.n)
  return(data.n)
}

get_normal = function(expression, clinical.m){
  print('-Extracting normal samples from expression data')
  data <- expression
  clin <- clinical.m
  
  # get normal only parent specimen label in clinical data
  normal_p <- clin$Parent.Specimen.Label[grepl("Normal", clin$tumor_normal)]
  
  # get parent specimen label in expression data
  p <-   str_split_fixed(colnames(data), "_", 2)[,1]
  
  # keep normal only specimens and no replicates
  data.n <- data[,colnames(data)[p %in% normal_p]]
  # data.n <- data.frame(data.n)
  
  # check multiple samples for one participant
  participant_ids <- NULL
  for (i in 1:ncol(data.n)) {
    specimen_label <- colnames(data.n)[i]
    participant_id <- clin$Participant.ID[clin$Specimen.Label == specimen_label]
    participant_ids <- c(participant_ids, as.character(participant_id))
  }
  participant_ids_dup <- participant_ids[duplicated(participant_ids)]
  participant_ids_dup <- participant_ids_dup[!is.na(participant_ids_dup)]
  if (length(participant_ids_dup) > 0) {
    print(paste0('--multiple samples for one participant: ', participant_ids_dup))
    print('---keep one normal sample with the most non-NA values for multiple normal samples of one patient')
    for (partID in participant_ids_dup) {
      sampIDs <- colnames(data.n)[participant_ids == partID]
      nonNAnum <- colSums(!is.na(data.n[, sampIDs]))
      sampID <- sampIDs[which.max(nonNAnum)]
      print(paste0('----keep ', sampID, " for patient ", partID))
      ## take out unselected samples
      data.t <- data.n[, colnames(data.n)[!(colnames(data.n) %in% sampIDs[sampIDs %in% c(sampID)])]]
    }
  } else {
    print('--no multiple samples for one participant')
  }
  data.n <- sortColumn(data.n)
  return(data.n)
}


print_format_pro <- function(expression) {
  data <- expression
  data2p <- data.frame(Gene = rownames(expression))
  data2p <- cbind(data2p, data)
  return(data2p)
}

unshared_pho = function(Pho.m, colPhospho){
  data <- Pho.m
  col <- colPhospho
  data = data[, colnames(data)[grepl(pattern = "Ratio", x = colnames(data))]]
  data <- cbind(col, data)
  colnames(data) = sub(" Log Ratio", "", colnames(data))
  print('-Extracting unshared log ratio for phosphoproteome data done!')
  data <- sortColumn(data)
  return(data)
}

print_format_pho <- function(expression, colPhospho) {
  data <- expression
  col <- colPhospho
  data2p <- cbind(col, data)
  data2p <- sortColumn(data2p)
  return(data2p)
}

normalize_by_sample = function(m){
  m = as.matrix(m)
  m[!is.finite(m)] = NA
  for (i in 1:ncol(m)){
    mean_tmp <-  mean(m[,i], na.rm=T)
    sd_tmp <- sd(m[,i], na.rm=T)
    if ( !is.nan(mean_tmp) & !is.na(sd_tmp) ) {
      m[,i] = (m[,i] - mean_tmp)/sd_tmp
    } else {
      m[,i] = NA
    }
  }
  return(m)
}

# add to clinical data whether a certain data type exist for each sample
# expreType = "protein"/"phosphoprotein"
anno_clinical <- function(expression, clinical.m, expreType) {
  print('-Annotating clinical info table:')
  clin <- clinical.m
  
  # add in the parent specimen label
  s <-   unique(str_split_fixed(colnames(expression)[!(colnames(expression) %in% c("Gene", "Phosphosite"))], "\\.", 2)[,1])

  if (!(expreType %in% colnames(clin))) {
    clin[,expreType] <- FALSE
  }
  for (i in s) {
    tmp <- strsplit(i, "_", 3)[[1]]
    specimen_label <- paste0(tmp[!(tmp %in% c("1", "2"))], collapse = "_")
    p <- tmp[1]
    specimen_title <- unique(na.omit(clin[clin$Parent.Specimen.Label==p, "Parent.Specimen.Title"]))
    if (!(specimen_label %in% clin$Specimen.Label) & !is.na(specimen_label)) {
      print(paste0("--Add ", specimen_label))
      clin[nrow(clin) + 1,] <- NA
      clin[nrow(clin),"Specimen.Label"] <- specimen_label
      tmp <- as.vector(clin[, expreType])
      tmp[clin$Specimen.Label==specimen_label]  <- TRUE
      clin[, expreType]  <- tmp
      
      if ( length(p) > 0 ) {
        tmp <- as.vector(clin[, "Parent.Specimen.Label"])
        tmp[clin$Specimen.Label==specimen_label]  <- p
        clin[, "Parent.Specimen.Label"]  <- tmp
        
        if ( p %in% clin$Parent.Specimen.Label ) {
          partID <- unique(na.omit(clin[clin$Parent.Specimen.Label==p, "Participant.ID"]))
          if (length(partID) == 0) {
            print(paste0('--no Participant.ID: ', p))
          } else {
            tmp <- as.vector(clin[, "Participant.ID"])
            tmp[clin$Specimen.Label==specimen_label]  <- partID
            clin[, "Participant.ID"]  <- tmp
          }
          if (length(specimen_title) == 0) {
            print(paste0('--no Parent.Specimen.Title: ', p))
          } else {
            tmp <- as.vector(clin[, "Parent.Specimen.Title"])
            tmp[clin$Specimen.Label==specimen_label]  <- specimen_title
            clin[, "Parent.Specimen.Title"]  <- tmp
          }
        }
      } else {
        print(paste0('--no Parent.Specimen.Label: ', p))
      }
      
    } else if (!is.na(specimen_label)) {
      print(paste0("--Add ", specimen_label))
      tmp <- as.vector(clin[, expreType])
      tmp[clin$Specimen.Label==specimen_label]  <- TRUE
      clin[, expreType]  <- tmp
      if (length(specimen_title) == 0) {
        print(paste0('--no Parent.Specimen.Title: ', p))
      }
    }
  }
  clin$tumor_normal <- NA
  clin$tumor_normal[grepl('Normal', as.vector(clinical$Parent.Specimen.Title),ignore.case = T)] <- 'Normal'
  clin$tumor_normal[grepl('Tumor', as.vector(clinical$Parent.Specimen.Title),ignore.case = T)] <- 'Tumor'
  clin <- unique(clin)
  return(clin)
}

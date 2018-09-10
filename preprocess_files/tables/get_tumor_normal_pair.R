# Yige Wu @ WashU 2018 Feb
# get the paired tumor-normal info
## TODO: add unmatch normal samples


# input -------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/preprocess_files/preprocess_files.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')

clinical <- read.delim(paste0(inputD, "Specimen_Data_20161005_Yige_20180307.txt"), sep = "\t")
clinical <- data.frame(clinical)


# functions ---------------------------------------------------------------
sampID2partID <- function(sampIDs, clinical) {
  clinical <- clinical[!is.na(clinical$tumor_normal),]
  partIDs <- NULL
  for (sampID in sampIDs) {
    partID <- unique(as.character(clinical$Participant.ID[clinical$Specimen.Label == sampID]))
    if (length(partID) != 0 ) {
      partIDs <- c(partIDs, partID)
    } else {
      partIDs <- c(partIDs, NA)
      print(paste0(sampID, "not in clinical info!"))
    }
  }
  return(partIDs)
}

tumorSamp2normal <- function(tumorSampIDs, clinical) {
  clinical <- clinical[!is.na(clinical$tumor_normal)& !is.na(clinical$protein) & clinical$protein,]
  tumorPartIDs <-sampID2partID(tumorSampIDs, clinical)
  normalSampIDs <- NULL
  for (partID in tumorPartIDs) {
    normalSampID <- unique(as.character(clinical$Specimen.Label[clinical$Participant.ID== partID & clinical$tumor_normal == "Normal"]))
    if (length(normalSampID) != 0) {
      normalSampIDs <- c(normalSampIDs, normalSampID)
    } else {
      normalSampIDs <- c(normalSampIDs, NA)
      # print(paste0(partID, "has no normal sample!"))
    }
  }
  ids <- data.frame(partID = tumorPartIDs, tumorSampIDs = tumorSampIDs, normalSampIDs = normalSampIDs)
  rownames(ids) <- ids$tumorSampIDs
  return(ids)
}

tumor_normal_match <- list()
for (cancer in c("BRCA","OV","CO","UCEC")) {
# for (cancer in c("UCEC")) {
  Pro.n_aveRep <- fread(paste(inputD, cancer, "/", prefix[cancer], "_PRO_formatted_normalized_replicate_averaged.txt",sep=""),
                        data.table = F)
  Pro.n_tumor <- get_tumor(Pro.n_aveRep, clinical.m = clinical)
  
  tumor_normal_match_tmp <- tumorSamp2normal(as.vector(colnames(Pro.n_tumor)), clinical = clinical)
  tumor_normal_match[[cancer]] <- tumor_normal_match_tmp
  write.table(x = tumor_normal_match_tmp, file = paste0(makeOutDir(), cancer, "_tumor_normal_matched_specimenIDs.txt"),
              row.names = F, quote = F)
}

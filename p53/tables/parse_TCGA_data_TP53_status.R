# Yige Wu @ WashU 2018 Aug
## For formating CCRCC data in data freeze to be in concord with CDAP proteomics data and genomics data


# source ------------------------------------------------------------------
wd <- getwd()
if (wd != "/Users/yigewu/Box Sync") {
  setwd("/Users/yigewu/Box Sync")
}
source("./cptac2p_analysis/phospho_network/phospho_network_shared.R")
library(readr)

# set variables -----------------------------------------------------------
pipeline_type <- "TCGA"
cancers2process <- c("BRCA", "OV", "CO", "UCEC", "CCRCC", "LIHC")
sample_type <- "tumor"
gene_altered <- "TP53"

for (cancer in cancers2process) {

# input sample list -------------------------------------------------------


# input maf ---------------------------------------------------------------


# generate mutation matrix ------------------------------------------------
## get the patient ids with deletion
  
  ## get the patient ids with truncations
  
  ## make classification (simple)
  
  ## make classification complex
  
  
  # process RPPA --------------------------------------------------------------
  ## input 
  files_tmp <- list.files(path = paste0(cBioPortalD, cancer, "/"), recursive = T)
  files_tmp
  files_tmp <- files_tmp[grepl(x = files_tmp, pattern = "data_rppa", ignore.case = T) & grepl(x = files_tmp, pattern = "tcga", ignore.case = T)]
  files_tmp
  if (length(files_tmp) > 1) {
    stop("multiple RNA files!")
  }
  if (length(files_tmp) == 0) {
    files_tmp <- paste0("./Ding_Lab/Projects_Current/TP53_shared_data/resources/TCGA/RPPA/", cancer, ".rppa.txt")
    if (cancer == "CO") {
      files_tmp <- paste0("./Ding_Lab/Projects_Current/TP53_shared_data/resources/TCGA/RPPA/", "COADREAD", ".rppa.txt")
    }
    if(!file.exists(files_tmp)) {
      stop("no RPPA data!")
    }
    data_rppa <- read_delim(files_tmp, 
                            "\t", escape_double = FALSE, trim_ws = TRUE)
  } else {
    data_rppa <- read_delim(paste0(cBioPortalD,  cancer, "/", files_tmp), 
                            "\t", escape_double = FALSE, trim_ws = TRUE)
  }

  partIDs <- colnames(data_rppa)[!(colnames(data_rppa) %in% c("Composite.Element.REF"))]
  
  ## divide protein and phospho
  rppa_entries <- as.vector(data_rppa$Composite.Element.REF)
  rppa_entries %>%
    length()
  
  ### get protein
  rppa_entries_pro <- rppa_entries[!grepl(x = rppa_entries, pattern = "_p")]
  rppa_entries_pro
  rppa_entries_pro %>% length()
  rppa_entries_pro <- data.frame(Composite.Element.REF = rppa_entries_pro, Gene = str_split_fixed(string = rppa_entries_pro, pattern = "\\|", 2)[,1])
  rppa_entries_pro %>%
    head
  rppa_data_pro <- cbind(rppa_entries_pro, data.frame(data_rppa[data_rppa$Composite.Element.REF %in% rppa_entries_pro$Composite.Element.REF, partIDs]))
  rppa_data_pro %>% head
  rppa_data_pro$Composite.Element.REF <- NULL
  
  colnames(rppa_data_pro) <- c("Gene", partIDs)
  rppa_data_pro %>% nrow()
  expresson_type <- "PRO"
  write.table(x = rppa_data_pro, file = paste0(makeOutDir(resultD = resultD), cancer, "_", expresson_type, "_", sample_type, "_", pipeline_type, "_", "RPPA", "_", "partID", ".txt"), quote = F, row.names = F, sep = "\t")
  
  ## get phospho data
  rppa_entries_pho <- rppa_entries[grepl(x = rppa_entries, pattern = "_p")]
  rppa_entries_pho
  rppa_entries_pho %>% length()
  rppa_entries_pho <- data.frame(Composite.Element.REF = rppa_entries_pho, 
                                 Gene = str_split_fixed(string = rppa_entries_pho, pattern = "\\|", 2)[,1],
                                 Phosphosite = str_split_fixed(string = str_split_fixed(string = rppa_entries_pho, pattern = "_p", 2)[,2], pattern = "-", 2)[,1])
  rppa_entries_pho %>%
    head
  rppa_data_pho <- cbind(rppa_entries_pho, data.frame(data_rppa[data_rppa$Composite.Element.REF %in% rppa_entries_pho$Composite.Element.REF, partIDs]))
  rppa_data_pho %>% head
  rppa_data_pho$Composite.Element.REF <- NULL
  colnames(rppa_data_pho) <- c("Gene", "Phosphosite", partIDs)
  expresson_type <- "PHO"
  write.table(x = rppa_data_pho, file = paste0(makeOutDir(resultD = resultD), cancer, "_", expresson_type, "_", sample_type, "_", pipeline_type, "_", "RPPA", "_", "partID", ".txt"), quote = F, row.names = F, sep = "\t")
  
}

# Yige Wu @ WashU 2018 Jan
## shared parameters and functions for cptac2p_analysis

#  source -----------------------------------------------------------------
library(data.table)
library(readxl)

# heavily-used data paths ---------------------------------------------------------
## base directory at Box Sync
folders <- strsplit(x = rstudioapi::getSourceEditorContext()$path, split = "\\/")[[1]]
baseD <- paste0(paste0(folders[1:which(folders == "Box Sync")], collapse = "/"), "/")
setwd(baseD)
# baseD = "/Users/yigewu/Box\ Sync/"
# setwd(baseD)

## CPTAC2 prospective and CPTAC proteomics data
cptac_sharedD <- "./Ding_Lab/Projects_Current/CPTAC/cptac_shared/"

## CPTAC2 prospective genomics data
cptac2pD <- "./Ding_Lab/Projects_Current/CPTAC/CPTAC_Prospective_Samples/"
cptac2p_genomicD <- "./Ding_Lab/Projects_Current/CPTAC/CPTAC_Prospective_Samples/"

## CPTAC2 retrospective proteomics/genomics data
pan3can_shared_dataD <- paste0(cptac2pD, "pan3can_shared_data/")

## processed data directory
resultD <- "./cptac2p/analysis_results/"
preprocess_filesD <- paste0(resultD, "preprocess_files/")

# set variables -----------------------------------------------------------
cancers_sort <- c("BRCA", "OV", "CO")
prefix <- NULL
prefix["BRCA"] <- "BRCA_BI"
prefix["OV"] <- "OV_PNNL"
prefix["CO"] <- "CO_PNNL"
prefix["UCEC"] <- "UCEC_PNNL"
prefix["CCRCC"] <- "CCRCC_JHU"
prefix["LUAD"] <- "LUAD_BI"

datatypes <- c("PRO", "PHO", "collapsed_PHO")

maf_files <- list.files(path = paste0(cptac2p_genomicD, "Somatic/"))
#list.files(path = paste0(cptac2p_genomicD, "FPKM/08092018/"), recursive = T)
fpkm_files <- paste0(cptac2p_genomicD, "FPKM/08092018/", 
                     c("fpkm_prospective_breastcancer_110617_swap_corrected.csv",
                       "fpkm_ov_allgenes.011718.csv",
                       "fpkm_crc_allgenes.v1.0.10232017.csv"))
names(fpkm_files) <- cancers_sort



# gene lists --------------------------------------------------------------
## SMGs
SMGs <- list()

SMGs[["BRCA"]] <- c("PIK3CA", "TP53", "MAP3K1", "MAP2K4", "GATA3", "MLL3", "CDH1", "PTEN", "PIK3R1", "AKT1", "RUNX1", "CBFB", "TBX3", "NCOR1", "CTCF", "FOXA1", "SF3B1", "CDKN1B", "RB1", "AFF2", "NF1", "PTPN22", "PTPRD")

SMGs[["OV"]] <- c("TP53", "NF1", "BRCA1", "BRCA2", "RB1", "CDK12")

SMGs[["CO"]] <- c("APC", "TP53", "KRAS", "PIK3CA", "FBXW7", "SMAD4", "TCF7L2", "NRAS", "BRAF")

## significant SCNA
CNAs <- list()
cancer_synonyms <- c("BRCA", "OV", "CRC", "All_cancer")
names(cancer_synonyms) <- c(cancers_sort, "PANCAN")
for (cancer in c(cancers_sort, "PANCAN")) {
  CNAs[[cancer]] <- list()
  for (cna_type in c("amplification", "deletion")) {
    if (cancer == "PANCAN") {
      cna_tab <- readxl::read_xlsx(path = "./Ding_Lab/Projects_Current/Gene Lists/copy_number_alteration/ng.2760-S5.xlsx", 
                                   sheet = paste0(substr(x = cna_type, start = 1, stop = 3), "_genes.", cancer_synonyms[cancer]))
    } else {
      cna_tab <- readxl::read_xlsx(path = "./Ding_Lab/Projects_Current/Gene Lists/copy_number_alteration/ng.2760-S5.xlsx", 
                                   sheet = paste0(substr(x = cna_type, start = 1, stop = 3), "_genes.annotated.", cancer_synonyms[cancer], ".txt"))
    }
    cna_tab.m <- unlist(cna_tab[4:nrow(cna_tab), 2:ncol(cna_tab)], use.names = F)
    cna_tab.m <- cna_tab.m[!is.na(cna_tab.m) & cna_tab.m != "NA"]
    cna_tab.m <- str_split_fixed(string = cna_tab.m, pattern = "\\-", 2)[,1]
    CNAs[[cancer]][[cna_type]] <- cna_tab.m
  }
}


# oncogenes and TSGs ------------------------------------------------------
driver_genes <- read_excel("./Ding_Lab/Projects_Current/TCGA_data/gene_lists/mmc1.xlsx", 
                           sheet = "Table S1", skip = 3)
driver_genes <- data.frame(driver_genes)
oncogenes <- driver_genes$Gene[grepl(x = driver_genes$Tumor.suppressor.or.oncogene.prediction..by.20.20.., pattern = "oncogene")]
tsgs <- driver_genes$Gene[grepl(x = driver_genes$Tumor.suppressor.or.oncogene.prediction..by.20.20.., pattern = "tsg")]
oncogenes <- c(oncogenes, "CTNND1", "PIK3R1")
tsgs <- tsgs[!(tsgs %in% c("CDH1", "CTNND1", "PIK3R1"))]
oncogenes <- unique(oncogenes)
tsgs <- unique(tsgs)


# functions ---------------------------------------------------------------
loadMSIMap <- function() {
  msi_score <- fread(input = "./Ding_Lab/Projects_Current/MSI/MSI_shared_data/MSI_score/CPTAC.MSI.score.tsv", data.table = F)
  msi_score$MSI_type <- ""
  msi_score$MSI_type[msi_score$Score >= 3.5] <- "MSI-H"
  msi_score$MSI_type[msi_score$Score < 3.5 & msi_score$Score >= 1.0] <- "MSI-L"
  msi_score$MSI_type[msi_score$Score < 1.0] <- "MSS"
  return(msi_score)
}

loadPAM50Map <- function() {
  file <- fread(paste0(cptac_sharedD, "BRCA/fpkm_pam50_brca_080718_pam50scores_swap_corrected.txt"), data.table = F)
  return(file)
}

# loadpam50fromclinTab <- function() {
#   file <- read_excel(paste0(cptac_sharedD, "5_CPTAC2_Breast_Prospective_Collection_BI/proteome-data-v1.01011/data/20171106_CPTAC2_ProspectiveBreastCancer_Sample_Annotations_Table_v79.xlsx"))
#   return(file)
# }

loadSampMap <- function() {
  file <- fread(input = paste0(cptac_sharedD, "Specimen_Data_20161005_Yige_20180815.txt"), data.table = F, sep = "\t")
  return(file)
}

loadMaf <- function(cancer, maf_files) {
  maf <- fread(input = paste0("./Ding_Lab/Projects_Current/CPTAC/CPTAC_Prospective_Samples/Somatic/", maf_files[grepl(x = maf_files, pattern = cancer)]), data.table = F)
  return(maf)
}

loadProteinNormalizedTumor <- function(cancer) {
  prefix <- NULL
  prefix["BRCA"] <- "BRCA_BI"
  prefix["OV"] <- "OV_PNNL"
  prefix["CO"] <- "CO_PNNL"
  prefix["UCEC"] <- "UCEC_PNNL"
  prefix["CCRCC"] <- "CCRCC_JHU"
  pro <- fread(input = paste0("./Ding_Lab/Projects_Current/CPTAC/cptac_shared/", cancer, "/", prefix[cancer], "_PRO_formatted_normalized_noControl.txt"), data.table = F)
  return(pro)
}

loadPhosphoproteinNormalizedTumor <- function(cancer) {
  prefix <- NULL
  prefix["BRCA"] <- "BRCA_BI"
  prefix["OV"] <- "OV_PNNL"
  prefix["CO"] <- "CO_PNNL"
  prefix["UCEC"] <- "UCEC_PNNL"
  prefix["CCRCC"] <- "CCRCC_JHU"
  pro <- fread(input = paste0("./Ding_Lab/Projects_Current/CPTAC/cptac_shared/", cancer, "/", prefix[cancer], "_collapsed_PHO_formatted_normalized_replicate_averaged_Tumor.txt"), data.table = F)
  return(pro)
}

loadPhosphositeNormalizedTumor <- function(cancer) {
  prefix <- NULL
  prefix["BRCA"] <- "BRCA_BI"
  prefix["OV"] <- "OV_PNNL"
  prefix["CO"] <- "CO_PNNL"
  prefix["UCEC"] <- "UCEC_PNNL"
  prefix["CCRCC"] <- "CCRCC_JHU"
  pro <- fread(input = paste0("./Ding_Lab/Projects_Current/CPTAC/cptac_shared/", cancer, "/", prefix[cancer], "_PHO_formatted_normalized_noControl.txt"), data.table = F)
  return(pro)
}

loadGeneList = function(gene_type, cancer, is.soft.limit) {
  if (gene_type == "RTK") {
    RTK_file = read.table(paste(pan3can_shared_dataD,"reference_files/RTKs_list.txt",sep = ""))
    genelist <- as.vector(t(RTK_file))
  }
  driver_table <- read_excel("./Ding_Lab/Projects_Current/TCGA_data/gene_lists/mmc1.xlsx", 
                                sheet = "Table S1", skip = 3)
  
  driver_table <- data.frame(driver_table)
  if (cancer != "CO") {
    driver_cancer = driver_table[driver_table$Cancer == cancer, ]
    ## oncogenes for OV: RAB25
  } else if (cancer == "CO") {
    driver_cancer = driver_table[driver_table$Cancer == "COADREAD", ]
  }
  if (gene_type == "driver") {
    genelist <- as.vector(driver_cancer$Gene)
  }
  if (gene_type == "tsg" || gene_type == "oncogene") {
    if (is.soft.limit == "soft") {
      genelist <- as.vector(driver_cancer$Gene)[grepl(pattern = gene_type, x = as.vector(driver_cancer$Tumor.suppressor.or.oncogene.prediction..by.20.20..))]
    } else {
      genelist <- as.vector(driver_cancer$Gene)[driver_cancer$Tumor.suppressor.or.oncogene.prediction..by.20.20.. == gene_type]
    }
  }
  
  return(genelist)
}

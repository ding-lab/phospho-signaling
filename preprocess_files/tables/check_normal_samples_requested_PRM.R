# Yige Wu @ WashU Jan 2018
# annotated the samples to request for PRM with whether the same patient have tumor samples with global proteomics data

# source -----------------------------------------------------------
library(data.table)
library(readr)
library(readxl)
source('/Users/yigewu/Box Sync/cptac2p_analysis/preprocess_files/preprocess_files.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')


# inputs ------------------------------------------------------------------
## input files logging which samples have proteomics data
clinical <- read.delim(paste0(inputD, "Specimen_Data_20161005_Yige_20180307.txt"), sep = "\t")
clinical <- data.frame(clinical)
clinical_pro <- clinical[clinical$protein,]
## input normal samples Sherri picked for PRM
normal_tab <- read_excel("~/Box Sync/cptac2p/sherri/BRCA_ BCR inventory_ 03 22 18.xlsx")
normal_chosen <- read_excel("~/Box Sync/cptac2p/sherri/BRCA_ BCR inventory_ 03 22 18_sherri_picked.xlsx")

# annotate whether any tumor sample from the same patient has proteomics data----------------------------------------------------------------
normal_tab$global.tumor <- (normal_tab$`Subject ID` %in% clinical_pro$Participant.ID[clinical_pro$tumor_normal == "Tumor"])
normal_tab$global.normal <- (normal_tab$`Subject ID` %in% clinical_pro$Participant.ID[clinical_pro$tumor_normal == "Normal"])
normal_chosen$global.tumor <- (normal_chosen$`Subject ID` %in% clinical_pro$Participant.ID[clinical_pro$tumor_normal == "Tumor"])
normal_chosen$global.normal <- (normal_chosen$`Subject ID` %in% clinical_pro$Participant.ID[clinical_pro$tumor_normal == "Normal"])

write.table(normal_tab, row.names = F, quote=F, sep = ',', 
            file=paste(makeOutDir(), "BRCA_BCR_inventory_03_22_18_global_availability_annotated.csv",sep=""))
write.table(normal_chosen, row.names = F, quote=F, sep = ',', 
            file=paste(makeOutDir(), "BRCA_BCR_inventory_03_22_18_sherri_picked_global_availability_annotated.csv",sep=""))

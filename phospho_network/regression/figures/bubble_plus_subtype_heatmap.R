# Yige Wu @ WashU 2018 Jan
# plotting bubble chart and subtype phosphorylation/protein expression heatmap side by side for BRCA dataset

# choose kinase or phosphotase, significance level, outlier threshold and least sample number-------------------------
# protein <- "kinase"
#protein <- "phosphotase"
top <- 50 # choose top n rows for fdr for model1

# library -----------------------------------------------------------------
library(grid)
library(dplyr)
library(gridExtra)
library(gtable)
library(readxl)

# input -------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
source(paste0(baseD, "cptac2p_analysis/pan3can_aes.R")) # aes for general purposes; it should be one directory out of the working directory
resultDnow <- makeOutDir()
RTKs <- loadGeneList("RTK")
var_list <- c('pro_kin', 'pho_kin')
names(var_list) <- c('cis', 'trans')
subtypes_list <- c('pam50', 'MSI_type', 'tumor_normal_WBasal'); names(subtypes_list) <- c("BRCA", "CO", "OV")

# bubble chart, top XX significant result ---------------------
#for (cancer in c("BRCA")) {
#for (cancer in c("CO", "BRCA", "OV")) {
for (cancer in c("CO")) {
  pho_subtype_mean <- fread(input = paste0(inputD, cancer, "/", prefix[cancer], "_PHO_formatted_normalized_",subtypes_list[cancer],"_mean.txt"), data.table = F)
  pro_subtype_mean <- fread(input = paste0(inputD, cancer, "/", prefix[cancer], "_PRO_formatted_normalized_",subtypes_list[cancer],"_mean.txt"), data.table = F)
  phog_subtype_mean <- fread(input = paste0(inputD, cancer, "/", prefix[cancer], "_PHO_formatted_normalized_collapsed_",subtypes_list[cancer],"_mean.txt"), data.table = F)
  
  for (protein in c('kinase', 'phosphotase')) {
    tn = paste0(resultD,"regression/tables/",protein,"_substrate_regression_cptac2p_3can.txt")
    table_3can <- read.delim(file = tn)
    table_3can <- markSigKS(table_3can, sig_thres = sig, protein_type = protein)
    
    for (self in c("cis","trans")) {
      t0 <- table_3can[table_3can$SELF == self & table_3can$fdr_sig & table_3can$coef_sig & table_3can$Cancer==cancer,]
      var <- var_list[self]
      fdr_var <- paste("FDR_",var,sep = "");  coef_var <- paste("coef_",var,sep = "")
      
      if (self == "cis") {
        table_kinase <- pro_subtype_mean
      } else {
        table_kinase <- phog_subtype_mean
      }
      table_substrate <- pho_subtype_mean
      subtype_sort <- colnames(table_substrate)[-1]
      subtype_sort <- c(subtype_sort[!grepl("Normal", subtype_sort, ignore.case = T)],
                        subtype_sort[grepl("Normal", subtype_sort, ignore.case = T)])
      if (nrow(t0) > 0){
        ## draw the top significant FDR pairs
        t0 <- t0[order(t0[,fdr_var]), c("KINASE","SUBSTRATE","SUB_MOD_RSD",fdr_var, coef_var,"Cancer","pair")]
        colnames(t0) <- c("KINASE","SUBSTRATE","SUB_MOD_RSD","fdr","coef","Cancer","pair")
        t1 <- t0 %>%
          group_by(KINASE) %>%
          top_n(n = 3, wt = fdr )
        len <- min(top,nrow(t1))
        table_bubble <- t1[1:len,]
        table_bubble$sig <- (table_bubble$fdr <= sig)
        if (nrow(table_bubble) > 0){
          bubble_heatmap(cancer, protein, self, 'most_sig', 
                         table_bubble, table_substrate, table_kinase, resultDnow, subtype_sort)
        } else {
          print('top not enough data')
        }
        
        ## draw pairs with given kinase list
        table_bubble <- t0[t0$KINASE %in% RTKs, ]
        table_bubble$sig <- (table_bubble$fdr <= sig)
        if (nrow(table_bubble) > 0){
          bubble_heatmap(cancer, protein, self, 'RTK', 
                         table_bubble, table_substrate, table_kinase, resultDnow, subtype_sort)
        } else {
          print('RTK not enough data')
        }
      }
    }
  }
}



# show kinase or substrate == RTK across all cancer typesß------------------------------------------------------
GeneList_bubble_heatmap(resultD = resultD, inputD = inputD, resultDnow = resultDnow, sig_thres = sig, prefix = prefix, 
                        genelist = RTKs, 
                        id = "RTK-")

# show kinase or substrate == SMGs across all cancer typesß------------------------------------------------------
GeneList_bubble_heatmap(resultD = resultD, inputD = inputD, resultDnow = resultDnow, sig_thres = sig, prefix = prefix, 
                        genelist =c("PIK3CA", "MAP3K1", "GATA3", "TP53", "CDH1", "MAP2K4"), 
                        id = "LumA-SMG-")

# show kinase or substrate == breast germline predisposition variant genes across all cancer typesß------------------------------------------------------
GeneList_bubble_heatmap(resultD = resultD, inputD = inputD, resultDnow = resultDnow, sig_thres = sig, prefix = prefix, 
                        genelist =c("ATM", "BRCA1", "BRCA2", "BRIP1", "CHEK2", "NBN", "PTEN", "RAD51C","TP53"), 
                        id = "BRCA-germline-")

# show kinase or substrate == breast CNV genes across all cancer typesß------------------------------------------------------
GeneList_bubble_heatmap(resultD = resultD, inputD = inputD, resultDnow = resultDnow, sig_thres = sig, prefix = prefix, 
                        genelist = c("PIK3CA", "EGFR", "FOXA1", "ERBB2", "MLL3", "PTEN", "RB1", "MAP2K4"), 
                        id = "BRCA-CNV-")

# provide kinase list, show top cis regulating kinases (retrospective) across all cancer typesß------------------------------------------------------
genelist <- c("WNK1", "ERBB2", "RPS6KA4", "MARK2", "NEK9", 
              "RIPK2","CDC42BPA", "GTF2F1","MARK3","PAK1",
              "PAK4", "PRKCD")
genelist_id <- "retro_top_cis_kinases"



# see the shared and unique ks pairs in all 3 cancers ---------------------
GeneList_bubble_heatmap(resultD = resultD, inputD = inputD, resultDnow = resultDnow, sig_thres = sig, prefix = prefix, 
                        genelist = c("CAMK2G", "GTF2F1", "PRKCD", "PRKDC", "CASK", "NEK9", "NME2", "PAK1", "PKN1" ,"PRKAA1"), 
                        id = "uniq_CO-cis")

GeneList_bubble_heatmap(resultD = resultD, inputD = inputD, resultDnow = resultDnow, sig_thres = sig, prefix = prefix, 
                        genelist = c("EEF2K", "PAK4", "WNK1", "MARK3", "NEK9", "RIPK2", "SRC"), 
                        id = "uniq_OV-cis")

GeneList_bubble_heatmap(resultD = resultD, inputD = inputD, resultDnow = resultDnow, sig_thres = sig, prefix = prefix, 
                        genelist = c("BRAF", "RPS6KA4", "MAP2K4", "MARK3", "PAK1", "ABL1", "ATM", "CAMK2G"), 
                        id = "uniq_BRCA-cis")

GeneList_bubble_heatmap(resultD = resultD, inputD = inputD, resultDnow = resultDnow, sig_thres = sig, prefix = prefix, 
                        genelist = c("EEF2K", "EGFR", "PAK4", "PRKCD", "RIPK2", "AKT1", "CDC42BPA", "CDK9", "ERBB2", "PAK2"), 
                        id = "shared3can-cis-")

# see TP53 related ones ---------------------------------------------------

GeneList_bubble_heatmap(resultD = resultD, inputD = inputD, resultDnow = resultDnow, sig_thres = 1, prefix = prefix, 
                        genelist = c("TP53"), 
                        id = "TP53-")


# see shared trans ones ---------------------------------------------------
protein <- "kinase"
tn = paste0(resultD, "signal_hub/tables/overlap_diffexp_regression/", protein,"_substrate_regression_overlap_cptac2p_3can_diffPhosphositeAnnotated.txt")
table_3can <- fread(tn, data.table = F)
table_3can <- markSigSiteCan(regression = table_3can, sig_thres = sig, protein_type = protein)
for (pairset in c('shared3can', 'uniq_BRCA', 'uniq_OV', 'uniq_CO')) {
  kspairList_bubble_heatmap(resultD = resultD, inputD = inputD, resultDnow = resultDnow, sig_thres = sig, prefix = prefix, 
                          pairlist = as.vector(table_3can$pair[table_3can[, pairset]]), 
                          id = paste0(pairset, "-trans-", protein, "-substrate-"))
}

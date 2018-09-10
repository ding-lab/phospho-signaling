# Yige Wu @ WashU 2018 Feb
# differential expression of protein and phosphoprotein (missing values imputed) using SAM

# inputs ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/cptac2p_analysis_shared.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/preprocess_files/get_tumor_normal_pair.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')

library(sSeq)
resultDnow <- makeOutDir()
clinical <- read.delim(paste0(inputD, "Specimen_Data_20161005_Yige_20180307.txt"), sep = "\t")

# function ----------------------------------------------------------------
run2pairedsamrbypair <- function(exp_t, exp_n, expression, identifier, clinical, df_tn_match, testType) {
  ## df_tn_match is a data frame matching tumors and normals
  ## construct vector differentiating tumor and normal
  y = c(1, -1)
  x = as.matrix(cbind(exp_t, exp_n))

    ## keep only the tumors and normals
  data=list(x=x, y=y,
            geneid = as.vector(expression[, identifier]), 
            genenames = as.vector(expression[, "Gene"]),
            logged2=TRUE)
  samfit<-samr(data,  resp.type="Two class paired", nperms=100, testStatistic = "wilcoxon")
  return(samfit)
}

run2pairedsamr <- function(expression, identifier, clinical, df_tn_match, sig_thres, outputID, outputDir, testType) {
  ## df_tn_match is a data frame matching tumors and normals
  ## construct vector differentiating tumor and normal
  tumorSampIDs <- as.vector(df_tn_match$tumorSampIDs[!is.na(df_tn_match$normalSampIDs)])
  normalSampIDs <- as.vector(df_tn_match[tumorSampIDs, "normalSampIDs"])
  
  y = c(1:length(tumorSampIDs), -(1:length(tumorSampIDs)))
  x = as.matrix(expression[,c(tumorSampIDs, normalSampIDs)])
  
  ## keep only the tumors and normals
  data=list(x=x, y=y,
            geneid = as.vector(expression[, identifier]), 
            genenames = as.vector(expression[, "Gene"]),
            logged2=TRUE)
  samfit<-samr(data,  resp.type="Two class paired", nperms=100, testStatistic = testType)
  return(samfit)
}

datatype2id <- c("Gene", "Phosphosite", "Gene")
names(datatype2id) <- datatypes
cancers <- c("BRCA", "OV", "CO", "UCEC")
# paired tumor vs normal ---------------------------------------------------------------
for (testType in c("wilcoxon")) {
  diffexp <- vector('list', 3)
  sam.obj <- vector('list', 3)
  
  for (datatype in datatypes[1]) {
    #for (datatype in c("PHO")) {
    diffexp[[datatype]]  <- vector('list', length = length(cancers))
    sam.obj[[datatype]]  <- vector('list', length = length(cancers))
    
    for (cancer in c("OV")) {
    # for (cancer in cancers) {
      diffexp[[datatype]][[cancer]] <- NULL
      sam.obj[[datatype]][[cancer]] <- NULL
      
      ## input protein or phosphoprotein data
      exp_data <- fread(paste0(inputD, cancer, "/", prefix[cancer], "_", datatype, "_formatted_normalized_replicate_averaged_imputed.txt"),
                        data.table = F)
      tn_match <- tumor_normal_match[[cancer]]
      diffexp[[datatype]][[cancer]][["genes"]] <- exp_data$Gene
      diffexp[[datatype]][[cancer]][["foldchange"]] <- NULL
      diffexp[[datatype]][[cancer]][["pvalue"]] <- NULL
      diffexp[[datatype]][[cancer]][["fdr"]] <- NULL
      samobj <- run2pairedsamr(expression = exp_data, identifier = datatype2id[datatype], 
                                  clinical = clinical, df_tn_match = tumor_normal_match[[cancer]], 
                                  sig_thres = sig, outputID = paste0(cancer, "_", datatype), outputDir = resultDnow, testType = testType)
      
 
    }
  }
}

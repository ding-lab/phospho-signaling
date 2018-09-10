# Yige Wu @ WashU 2018 Feb
# differential expression of protein and phosphoprotein (missing values imputed) using SAM

# inputs ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/preprocess_files/get_tumor_normal_pair.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
library(samr)
resultDnow <- makeOutDir()
clinical <- read.delim(paste0(inputD, "Specimen_Data_20161005_Yige_20180307.txt"), sep = "\t")

## input MSI scores
msi_score <- read_delim("~/Box Sync/MSI_CPTAC/Data/MSIsensor_Score_qgao/CPTAC.MSI.score.tsv",
                        "\t", escape_double = FALSE, trim_ws = TRUE)
# function ----------------------------------------------------------------
run2unpairedSAM <- function(expr1, expr2, geneids, genenames, sig_thres, outputID, outputDir, powerfig_id) {
  ## construct vector differentiating two groups
  y = c(rep(2, ncol(expr1)), rep(1, ncol(expr2)))
  x = cbind(expr1, expr2)
  
  ## keep only the tumors and normals
  samfit<-SAM(x, y, resp.type="Two class unpaired", nperms = 100,
              geneid = geneids, genenames = genenames,
              fdr.output = sig_thres,
              logged2 = T, testStatistic = "wilcoxon")
  ## plot power results
  samr.assess.samplesize.obj <- samr.assess.samplesize(samfit$samr.obj, samfit$samr.obj, 
                                                       dif = 1, samplesize.factors=c(1))
  fn = paste0(outputDir, outputID, "_diffexp_",powerfig_id,"2unpaired_wilcoxon_FDR_", sig_thres, '_diffnum_', as.character(samfit$siggenes.table$ngenes.up + samfit$siggenes.table$ngenes.lo),'.pdf')
  pdf(fn, height = 4, width = 5)
  print(samr.assess.samplesize.plot(samr.assess.samplesize.obj))
  dev.off()
  return(samfit)
}

run3unpairedSAM <- function(expr1, expr2, expr3, geneids, genenames, sig_thres, outputID, outputDir, powerfig_id) {
  ## construct vector differentiating two groups
  y = c(rep(1, ncol(expr3)), rep(2, ncol(expr2), rep(3, ncol(expr1))))
  x = cbind(expr3, expr2, expr1)
  
  ## keep only the tumors and normals
  samfit<-SAM(x, y, resp.type="Multiclass", nperms = 100,
              geneid = geneids, genenames = genenames,
              fdr.output = sig_thres,
              logged2 = T, testStatistic = "wilcoxon")
  return(samfit)
}

run2pairedSAM <- function(expression, identifier, clinical, df_tn_match, sig_thres, outputID, outputDir, testType) {
  ## df_tn_match is a data frame matching tumors and normals
  ## construct vector differentiating tumor and normal
  tumorSampIDs <- as.vector(df_tn_match$tumorSampIDs[!is.na(df_tn_match$normalSampIDs)])
  normalSampIDs <- as.vector(df_tn_match[tumorSampIDs, "normalSampIDs"])
  
  y = c(1:length(tumorSampIDs), -(1:length(tumorSampIDs)))
  x = as.matrix(expression[,c(tumorSampIDs, normalSampIDs)])

  ## keep only the tumors and normals
  samfit<-SAM(x,y,resp.type="Two class paired", nperms = 100,
              geneid = expression[, identifier], genenames = expression[, "Gene"],
              fdr.output = sig_thres,
              logged2 = T, testStatistic = testType)
  
  ## plot power results
  samr.assess.samplesize.obj <- samr.assess.samplesize(samfit$samr.obj, samfit$samr.obj, 
                                                       dif = 1, samplesize.factors=c(1))
  fn = paste0(outputDir, outputID, "_diffexp_FDR_", sig_thres, '_diffnum_', as.character(samfit$siggenes.table$ngenes.up + samfit$siggenes.table$ngenes.lo),"_", testType,'.pdf')
  pdf(fn, height = 4, width = 5)
  print(samr.assess.samplesize.plot(samr.assess.samplesize.obj))
  dev.off()
  return(samfit)
}

filterSAM <- function(samfit) {
  ## filter by fold change
  siggenes.table <- samfit$siggenes.table
  genes.up.filtered <- as.data.frame(siggenes.table$genes.up)
  if (nrow(genes.up.filtered) > 0){
    genes.up.filtered <- genes.up.filtered[as.numeric(as.vector(genes.up.filtered$`Fold Change`)) > 2,]
    genes.up.filtered$diffexp_type <- "up"
  }
  genes.lo.filtered <- as.data.frame(siggenes.table$genes.lo)
  
  if (nrow(genes.lo.filtered) > 0) {
    genes.lo.filtered <- genes.lo.filtered[as.numeric(as.vector(genes.lo.filtered$`Fold Change`)) < 1/2,]
    genes.lo.filtered$diffexp_type <- "down"
  }
  genes.diffexp <- rbind(genes.up.filtered, genes.lo.filtered)
  return(genes.diffexp)
}

datatype2id <- c("Gene", "Phosphosite", "Gene")
names(datatype2id) <- datatypes

# paired tumor vs normal ---------------------------------------------------------------
for (testType in c("standard", "wilcoxon")) {
  for (sig in c(0.1, 0.05, 0.2)) {
  # for (sig in c(1)) {
    diffexp <- vector('list', 3)
    for (datatype in datatypes) {
      #for (datatype in c("PHO")) {
      diffexp[[datatype]]  <- NULL
      #for (cancer in c("OV")) {
      for (cancer in c("BRCA", "OV", "CO", "UCEC")) {
        diffexp[[datatype]][[cancer]] <- NULL
        ## input protein or phosphoprotein data
        exp_data <- fread(paste0(inputD, cancer, "/", prefix[cancer], "_", datatype, "_formatted_normalized_replicate_averaged_imputed.txt"),
                          data.table = F)
        ## run SAM and filter results by fold change
        samobj <- run2pairedSAM(expression = exp_data, identifier = datatype2id[datatype], 
                                clinical = clinical, df_tn_match = tumor_normal_match[[cancer]], 
                                sig_thres = sig, outputID = paste0(cancer, "_", datatype), outputDir = resultDnow, testType = testType)
        
        ## store the result
        diffexp_tmp <- filterSAM(samfit = samobj)
        colnames(diffexp_tmp)[1:2] <- c("Gene", "Phosphosite");
        diffexp_tmp$Cancer <- cancer
        write.table(diffexp_tmp, row.names = F, quote=F, sep = '\t', 
                    file=paste0(makeOutDir(), cancer, "_", datatype, "_diffexp_paired_FDR", sig, "_", testType, ".txt"))
        
        diffexp_3can <- diffexp[[datatype]]
        diffexp_3can <- rbind(diffexp_3can, diffexp_tmp)
        diffexp[[datatype]] <- diffexp_3can
      }
      # write out ---------------------------------------------------------------
      write.table(diffexp[[datatype]], row.names = F, quote=F, sep = '\t',
      file=paste0(resultDnow, datatype, "_diffexp_4can_paired_FDR", sig, "_", testType, ".txt"))
    }
  }
}


# MSI-H tumors vs other tumors --------------------------------------------
for (sig in c(0.05, 0.1, 0.2)) {
  cancers <- c("BRCA", "CO")
  diffexp <- vector('list', length = length(cancers))
  for (datatype in datatypes) {
    diffexp[[datatype]]  <- NULL
    for (cancer in cancers) {
      diffexp[[datatype]][[cancer]] <- NULL
      ## input protein or phosphoprotein data
      exp_data <- fread(paste0(inputD, cancer, "/", prefix[cancer], "_", datatype, "_formatted_normalized_replicate_averaged_imputed.txt"),
                        data.table = F)
      ## get MSI-H tumors and the rest
      exp_t <- get_tumor(expression = exp_data, clinical.m = clinical)
      samples_t <- colnames(exp_t)
      patients_t <- sampID2partID(sampleID_vector = samples_t, sample_map = clinical)
      exp_t_msih <- exp_t[, patients_t %in% as.vector(msi_score$Sample[msi_score$Score>=3.5])]
      exp_t_other <- exp_t[, !(patients_t %in% as.vector(msi_score$Sample[msi_score$Score>=3.5]))]
      
      ## run SAM and filter results by fold change
      samobj <- run2unpairedSAM(expr1 = exp_t_msih, expr2 = exp_t_other, 
                                geneids = exp_data[,datatype2id[datatype]], genenames = as.vector(exp_data$Gene),
                                sig_thres = sig, outputID = paste0(cancer, "_", datatype), outputDir = resultDnow, powerfig_id = "tumor_MSIH_vs_tumor_rest")
      
      ## store the result
      diffexp_tmp <- filterSAM(samfit = samobj)
      colnames(diffexp_tmp)[1:2] <- c("Gene", "Phosphosite");
      diffexp_tmp$Cancer <- cancer
      diffexp_3can <- diffexp[[datatype]]
      diffexp_3can <- rbind(diffexp_3can, diffexp_tmp)
      diffexp[[datatype]] <- diffexp_3can
    }
    # write out ---------------------------------------------------------------
    write.table(diffexp[[datatype]], row.names = F, quote=F, sep = '\t', 
                file=paste0(resultDnow, datatype, "_diffexp_tumor_MSIH_vs_tumor_rest_2unpaired_FDR", sig, ".txt"))
  }
}


# MSI-H tumors vs other tumors and normals --------------------------------
for (sig in c(0.05, 0.1, 0.2)) {
  cancers <- c("BRCA", "CO")
  diffexp <- vector('list', length = length(cancers))
  for (datatype in datatypes) {
    diffexp[[datatype]]  <- NULL
    for (cancer in cancers) {
      diffexp[[datatype]][[cancer]] <- NULL
      ## input protein or phosphoprotein data
      exp_data <- fread(paste0(inputD, cancer, "/", prefix[cancer], "_", datatype, "_formatted_normalized_replicate_averaged_imputed.txt"),
                        data.table = F)
      ## get MSI-H tumors and the rest
      exp_t <- get_tumor(expression = exp_data, clinical.m = clinical)
      samples_t <- colnames(exp_t)
      patients_t <- sampID2partID(sampleID_vector = samples_t, sample_map = clinical)
      exp_t_msih <- exp_t[, patients_t %in% as.vector(msi_score$Sample[msi_score$Score>=3.5])]
      exp_t_other <- exp_t[, !(patients_t %in% as.vector(msi_score$Sample[msi_score$Score>=3.5]))]
      exp_n <- get_normal(expression = exp_data, clinical.m = clinical)
      
      ## run SAM and filter results by fold change
      samobj <- run3unpairedSAM(expr1 = exp_t_msih, expr2 = exp_t_other, expr3 = exp_n,
                                geneids = exp_data[,datatype2id[datatype]], genenames = as.vector(exp_data$Gene),
                                sig_thres = sig, outputID = paste0(cancer, "_", datatype), outputDir = resultDnow, powerfig_id = "tumor_MSIH_vs_tumor_rest&normal")
      
      ## store the result
      diffexp_tmp <- filterSAM(samfit = samobj)
      colnames(diffexp_tmp)[1:2] <- c("Gene", "Phosphosite");
      diffexp_tmp$Cancer <- cancer
      diffexp_3can <- diffexp[[datatype]]
      diffexp_3can <- rbind(diffexp_3can, diffexp_tmp)
      diffexp[[datatype]] <- diffexp_3can
    }
    # write out ---------------------------------------------------------------
    write.table(diffexp[[datatype]], row.names = F, quote=F, sep = '\t', 
                file=paste0(resultDnow, datatype, "_diffexp_tumor_MSIH_vs_tumor_rest&normal_3unpaired_FDR", sig, ".txt"))
  }
}




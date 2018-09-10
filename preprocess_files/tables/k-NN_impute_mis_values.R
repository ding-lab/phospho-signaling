# Yige Wu @ WashU 2018 Feb
# Missing value imputation for protein and phosphosite data
# For algorithms that could not handle missing values in the data—e.g., k-means clustering, marker selections using SAM85 and GSEA18— missing values were imputed using k-nearest neighbor (k-NN) imputation


# load libraries ----------------------------------------------------------
library(pamr)
source('/Users/yigewu/Box Sync/cptac2p_analysis/preprocess_files/preprocess_files.R')


# loop around 3 cancers ---------------------------------------------------
# for (cancer in c("BRCA", "OV", "CO", "UCEC")) {
for (cancer in c("UCEC")) {
  for (datatype in datatypes) {
    ## input protein and phosphoprotein data
    exp_data <- fread(paste0(inputD, cancer, "/", prefix[cancer], "_", datatype, "_formatted_normalized_replicate_averaged.txt"),
                      data.table = F)
    samples <- colnames(exp_data)[!(colnames(exp_data) %in% c("Gene", "Phosphosite"))]
    headers <- colnames(exp_data)[(colnames(exp_data) %in% c("Gene", "Phosphosite"))]
    ## discard protein and phosphosites with >50% missing data
    mis50 <- rowSums(is.na(exp_data[, samples]))
    mis50_filtered <- exp_data[mis50 <= 0.5*(length(samples)),]
    ## inpute missing values using k-nearest neighbor imputation
    ### remove proteins and phosphosites with more than 50% missing data
    exp_imputed <- pamr.knnimpute(data = list(x = as.matrix(mis50_filtered[, samples]), y = samples,
                                  k = 10, rowmax = 0.5, colmax = 0.8, maxp = 1500))
    exp_imputed2p <- data.frame(mis50_filtered[,headers]);
    colnames(exp_imputed2p) <- headers
    exp_imputed2p <- cbind(exp_imputed2p, exp_imputed[['x']])
    ## write out
    write.table(exp_imputed2p, row.names = F, quote=F, sep = '\t', file=paste(inputD, cancer, "/", prefix[cancer], "_", datatype, "_formatted_normalized_replicate_averaged_imputed.txt",sep=""))
  }
}



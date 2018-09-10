# Yige Wu @ WashU March 2018
# analyze cohort level proteomics data and convert to different matrices in a sample-gene format

# source ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/preprocess_files/preprocess_files.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')


# inputs ------------------------------------------------------------------
unfactorize = function(df){
  for(i in which(sapply(df, class) == "factor")) df[[i]] = as.numeric(as.character(df[[i]]))
  return(df)
}

ecdf_fun = function(x,perc) ecdf(x)(perc)

expression_effect = function(m){ 
  cat("##### EXPRESSION ANALYSIS #####\n")
  minNum = 5
  m = as.matrix(m)
  num = nrow(m)
  m2 = as.matrix(m[rowSums(!is.na(m)) >= minNum, ])
  num_NA= nrow(m2)
  cat(paste("Original number of markers:", num, "; NA filtered:", num_NA, "\n", sep=" "))
  
  # initiate tables
  outlier = matrix(data = NA,nrow=dim(m2)[1],ncol=dim(m2)[2])
  row.names(outlier) = row.names(m2)
  colnames(outlier) = colnames(m2)
  exp_score = outlier
  exp_quantile = outlier
  
  # gene-wise expression score and quantile score
  for (i in 1:nrow(m2)){
    #IQR = quantile(m2[i,], probs=0.75, na.rm=T) - quantile(m2[i,], probs=0.25, na.rm=T) 
    exp_score[i,] = m2[i,]#(m2[i,] - quantile(m2[i,], probs=0.50, na.rm=T))/IQR
    exp_quantile[i,] = ecdf_fun(m2[i,],m2[i,])
  }
  
  return(list("exp_score"=exp_score, "exp_quantile"=exp_quantile))
}

# process per cancer in cptac2 prospective ------------------------------------------------------
exp_quantile_tables = vector("list")
exp_score_tables = vector("list")
for (cancer in c("CO", "OV", "BRCA")) {
  ## input proteomics tumor only data
  pro.t <- fread(input = paste0(inputD, cancer, "/", prefix[cancer], "_PRO_formatted_normalized_noControl.txt"), 
                 data.table = F)
  ## derive expression quantile
  exp_table <- pro.t[, !(colnames(pro.t) %in% c("Gene"))]
  rownames(exp_table) <- pro.t$Gene
  exp_results = expression_effect(exp_table)
  
  exp_quantile_table_m = melt(exp_results$exp_quantile)
  colnames(exp_quantile_table_m) <- c("Gene", "Specimen.ID", "qt")
  exp_quantile_table_m$cancer = cancer
  
  ## derive expression score
  exp_score_table_m = melt(exp_results$exp_score)
  colnames(exp_score_table_m) <- c("Gene", "Specimen.ID", "score")
  exp_score_table_m$cancer = cancer
  
  ## store in RDS
  exp_quantile_tables[[cancer]] = exp_quantile_table_m
  exp_score_tables[[cancer]] = exp_score_table_m
}

resultDnow <- makeOutDir()
saveRDS(object = exp_quantile_tables, file = paste0(resultDnow, "cptac2p_protein_quantile.RDS"))
rm(exp_quantile_tables)
saveRDS(object = exp_score_tables, file = paste0(resultDnow, "cptac2p_protein_score.RDS"))
rm(exp_score_tables)

# process per cancer in CPTAC3  ------------------------------------------------------
exp_quantile_tables = vector("list")
exp_score_tables = vector("list")
for (cancer in c("UCEC")) {
  ## input proteomics tumor only data
  pro.t <- fread(input = paste0(inputD, cancer, "/", prefix[cancer], "_PRO_formatted_normalized_noControl.txt"), 
                 data.table = F)
  ## derive expression quantile
  exp_table <- pro.t[, !(colnames(pro.t) %in% c("Gene"))]
  rownames(exp_table) <- pro.t$Gene
  exp_results = expression_effect(exp_table)
  
  exp_quantile_table_m = melt(exp_results$exp_quantile)
  colnames(exp_quantile_table_m) <- c("Gene", "Specimen.ID", "qt")
  exp_quantile_table_m$cancer = cancer
  
  ## derive expression score
  exp_score_table_m = melt(exp_results$exp_score)
  colnames(exp_score_table_m) <- c("Gene", "Specimen.ID", "score")
  exp_score_table_m$cancer = cancer
  
  ## store in RDS
  exp_quantile_tables[[cancer]] = exp_quantile_table_m
  exp_score_tables[[cancer]] = exp_score_table_m
}

resultDnow <- makeOutDir()
saveRDS(object = exp_quantile_tables, file = paste0(resultDnow, "cptac3_protein_quantile.RDS"))
rm(exp_quantile_tables)
saveRDS(object = exp_score_tables, file = paste0(resultDnow, "cptac3_protein_score.RDS"))
rm(exp_score_tables)

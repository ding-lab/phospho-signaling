# Yige Wu @ WashU March 2018
# analyze cohort level phosphoproteomics data and convert to different matrices in a sample-gene format

# source ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/preprocess_files/preprocess_files.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/expression_matrices/expression_shared.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')


# inputs ------------------------------------------------------------------

# process per cancer in cptac2 prospective ------------------------------------------------------
exp_quantile_tables = vector("list")
exp_score_tables = vector("list")
for (cancer in c("CO", "OV", "BRCA")) {
  ## input proteomics tumor only data
  pro.t <- fread(input = paste0(inputD, cancer, "/", prefix[cancer], "_collapsed_PHO_formatted_normalized_replicate_averaged_Tumor.txt"), 
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
saveRDS(object = exp_quantile_tables, file = paste0(resultDnow, "cptac2p_phosphoprotein_quantile.RDS"))
rm(exp_quantile_tables)
saveRDS(object = exp_score_tables, file = paste0(resultDnow, "cptac2p_phosphoprotein_score.RDS"))
rm(exp_score_tables)

# process per cancer in CPTAC3  ------------------------------------------------------
exp_quantile_tables = vector("list")
exp_score_tables = vector("list")
for (cancer in c("UCEC")) {
  ## input proteomics tumor only data
  pro.t <- fread(input = paste0(inputD, cancer, "/", prefix[cancer], "_collapsed_PHO_formatted_normalized_replicate_averaged_Tumor.txt"), 
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
saveRDS(object = exp_quantile_tables, file = paste0(resultDnow, "cptac3_phosphoprotein_quantile.RDS"))
rm(exp_quantile_tables)
saveRDS(object = exp_score_tables, file = paste0(resultDnow, "cptac3_phosphoprotein_score.RDS"))
rm(exp_score_tables)

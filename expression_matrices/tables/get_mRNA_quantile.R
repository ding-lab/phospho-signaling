# Yige Wu @ WashU March 2018
# analyze cohort level proteomics data and convert to different matrices in a sample-gene format

# source ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/preprocess_files/preprocess_files.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/expression_matrices/expression_shared.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')

# inputs ------------------------------------------------------------------
fpkm_sampID_map <- fread(input = paste0(cptac_sharedD, "DCC/manifest-cptac2-final.txt"), data.table = F)
fpkm_sampID_map$rna_id <- str_split_fixed(string = fpkm_sampID_map$file_name, pattern = "\\.", 4)[,1]


# process per cancer in cptac2 prospective ------------------------------------------------------
exp_quantile_tables = vector("list")
exp_score_tables = vector("list")
for (cancer in c("BRCA")) {
  ## input proteomics tumor only data
  fpkm.t <- fread(input = fpkm_files[cancer], data.table = F)
  rna_ids <- colnames(fpkm.t[,-1])
  
  ## input sample map and transform into specimen IDs
  if (any(grepl("-R", rna_ids)) | any(grepl("_R", rna_ids)))  {
    Participant.IDs <- NULL
    for (rna_id in rna_ids) {
      tmp <- unique(fpkm_sampID_map$cases.submitter_id[fpkm_sampID_map$rna_id == rna_id])
      tmp <- tmp[1]
      Participant.IDs <- c(Participant.IDs, tmp)    
    }
    id_map <- data.frame(rna_id = rna_ids, Participant.ID = Participant.IDs)
  } else {
    id_map <- data.frame(rna_id = rna_ids, Participant.ID = rna_ids)
  }

  ## derive expression quantile
  exp_table <- fpkm.t[, !(grepl(pattern = "gene", x = colnames(fpkm.t), ignore.case = T))]
  rownames(exp_table) <- fpkm.t$gene
  exp_results = expression_effect(exp_table)
  
  exp_quantile_table_m = melt(exp_results$exp_quantile)
  colnames(exp_quantile_table_m) <- c("Gene", "rna_id", "qt")
  exp_quantile_table_m <- merge(exp_quantile_table_m, id_map, all.x = T) 
  exp_quantile_table_m$rna_id <- NULL
  exp_quantile_table_m$Specimen.ID <- NULL
  exp_quantile_tables[[cancer]] = exp_quantile_table_m
  rm(exp_quantile_table_m)
  
  ## derive expression score
  exp_score_table_m = melt(exp_results$exp_score)
  colnames(exp_score_table_m) <- c("Gene", "rna_id", "score")
  exp_score_table_m <- merge(exp_score_table_m, id_map, all.x = T) 
  exp_score_table_m$rna_id <- NULL
  exp_score_table_m$Specimen.ID <- NULL
  exp_score_tables[[cancer]] = exp_score_table_m
  rm(exp_score_table_m)
  rm(exp_results)
}

resultDnow <- makeOutDir(resultD = resultD)
saveRDS(object = exp_quantile_tables, file = paste0(makeOutDir(resultD = resultD), "cptac2p_mRNA_quantile.RDS"))
rm(exp_quantile_tables)
saveRDS(object = exp_score_tables, file = paste0(resultDnow, "cptac2p_mRNA_score.RDS"))
rm(exp_score_tables)

# process per cancer in CPTAC3  ------------------------------------------------------
# exp_quantile_tables = vector("list")
# exp_score_tables = vector("list")
# for (cancer in c("UCEC", "CCRC")) {
#   ## input proteomics tumor only data
#   fpkm.t <- fread(input = paste0("MSI_CPTAC/Data/fpkm/", cancer, ".fpkm"), data.table = F)
#   fpkm.t <- fpkm.t[!duplicated(fpkm.t$gene),]
#   ## derive expression quantile
#   exp_table <- fpkm.t[, !(grepl(pattern = "gene", x = colnames(fpkm.t), ignore.case = T))]
#   rownames(exp_table) <- fpkm.t$gene
#   exp_results = expression_effect(exp_table)
#   
#   exp_quantile_table_m = melt(exp_results$exp_quantile)
#   colnames(exp_quantile_table_m) <- c("Gene", "rna_id", "qt")
#   exp_quantile_table_m$Participant.ID <- str_split_fixed(string = exp_quantile_table_m$rna_id, pattern = "__", 2)[,2]
#   exp_quantile_table_m$cancer = cancer
#   
#   ## derive expression score
#   exp_score_table_m = melt(exp_results$exp_score)
#   colnames(exp_score_table_m) <- c("Gene", "rna_id", "score")
#   exp_score_table_m$Participant.ID <- str_split_fixed(string = exp_score_table_m$rna_id, pattern = "__", 2)[,2]
#   exp_score_table_m$cancer = cancer
#   
#   ## store in RDS
#   exp_quantile_tables[[cancer]] = exp_quantile_table_m
#   exp_score_tables[[cancer]] = exp_score_table_m
# }
# 
# resultDnow <- makeOutDir(resultD = resultD)
# saveRDS(object = exp_quantile_tables, file = paste0(resultDnow, "cptac3_mRNA_quantile.RDS"))
# rm(exp_quantile_tables)
# saveRDS(object = exp_score_tables, file = paste0(resultDnow, "cptac3_mRNA_score.RDS"))
# rm(exp_score_tables)

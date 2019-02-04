# Yige Wu @ WashU 2019 Feb
## gather the regression results of all the iterations for  from the same cohort downsized to the same numbers (same batch)
## see which iteration give the largest overlap with known kinase-substrate relationships

# source ------------------------------------------------------------------
wd <- getwd()
if (wd != "/Users/yigewu/Box Sync") {
  setwd("/Users/yigewu/Box Sync")
}
source('./cptac2p_analysis/phospho_network/phospho_network_shared.R')
library(dplyr)

# set variables -----------------------------------------------------------
reg_nonNA <- 20
reg_sig <- c(0.05, 0.05); names(reg_sig) <- c("kinase", "phosphatase")
names_stats <- c("size", "cancer", "num_regulated_direct", "num_regulated",
                 "is.present_batch", "file_prefix", 
                  "file_path")
num_iteration2stop_overlap <- 5
# inputs ------------------------------------------------------------------

# get the per downsize statistics -----------------------------------------
## get the table of regression results from downsized cohorts
downsize_fd <- "./cptac2p/analysis_results/phospho_network/regression/tables/regression_downsize/"
downsize_fns <- list.files(path = downsize_fd, 
                           pattern = "regression", recursive = T)
downsize_fns <- downsize_fns[!grepl(x = downsize_fns, pattern = "backups")]
downsize_fns

## group file names into each cohort
## group file names into each downsize number

downsize_fn_prefixs <- unique(str_split_fixed(string = downsize_fns, pattern = "_iteration", 2)[,1])
downsize_fn_prefixs %>% head()

summary_mat_allBatches <- NULL
for (downsize_fn_prefix in downsize_fn_prefixs) {
  downsize_fns2process <- downsize_fns[grepl(x = downsize_fns, pattern = downsize_fn_prefix)]
  downsize_fps2process <- paste0(downsize_fd, downsize_fns2process)
  downsize_fps2process <- downsize_fps2process[order(downsize_fps2process)]

  summary_mat <- matrix(data = NA, nrow = length(downsize_fps2process), ncol = length(names_stats),
                        dimnames = list(c(downsize_fps2process), names_stats))
  summary_mat
  
  for (downsize_fp in downsize_fps2process) {
    ## input table
    tab_tmp <- fread(input = downsize_fp, data.table = F)
    ## annotate with enzyme type
    tab_tmp <- annotate_enzyme_type(regression = tab_tmp, kinases = kinases, phosphatases = phosphatases)
    
    ## change nonNA number and mark significance
    tab_tmp <- change_regression_nonNA(regression = tab_tmp, reg_nonNA = reg_nonNA, reg_sig = reg_sig)
    
    ## annotate with direct kinase-substrate relations
    tab_tmp <- annotate_ks_source(regression = tab_tmp)
    tab_tmp %>% head()
    
    regulated_pairs <- as.vector(unique(tab_tmp$pair[tab_tmp$regulated]))
    regulated_pairs %>% head()
    regulated_direct_pairs <- as.vector(unique(tab_tmp$pair[tab_tmp$regulated & tab_tmp$is.direct]))
    regulated_direct_pairs %>% head()
    ## record stats
    summary_mat[downsize_fp, "size"] <- as.numeric(str_split(string = downsize_fn_prefix, pattern = "downsize|\\/")[[1]][2])
    summary_mat[downsize_fp, "cancer"] <- str_split(string = downsize_fn_prefix, pattern = "downsize|\\/")[[1]][3]
    summary_mat[downsize_fp, "num_regulated"] <- length(regulated_pairs)
    summary_mat[downsize_fp, "num_regulated_direct"] <- length(regulated_direct_pairs)
    summary_mat[downsize_fp, "file_path"] <- downsize_fp
    summary_mat[downsize_fp, "file_prefix"] <- str_split(string = downsize_fn_prefix, pattern = "\\/|_downsize")[[1]][3]
  }
  

  summary_mat <- data.frame(summary_mat, row.names = NULL)
  summary_mat$num_regulated_direct <- as.numeric(as.vector(summary_mat$num_regulated_direct))
  summary_mat$is.present_batch <- ifelse(summary_mat$num_regulated_direct == max(summary_mat$num_regulated_direct), T, F)
  summary_mat %>% head()

  ## get one iteration that is closest to the overlap and filter by overlapped pairs
  tab2w_input_fp <- as.character(unique(summary_mat$file_path[summary_mat$is.present_batch]))[1]
  tab2w_tmp <- fread(input = tab2w_input_fp, data.table = F)

  ## gather summary stats
  summary_mat_allBatches <- rbind(summary_mat_allBatches, summary_mat)
}
summary_mat_allBatches <- summary_mat_allBatches[!is.na(summary_mat_allBatches$cancer),]

# process original regression sup table-----------------------------------------------------------------
## input original regression sup table
original_fp <- paste0(ppnD, "regression/tables/", "change_regression_nonNA/", 
                      "regression", "_", "cptac2p_cptac3", "_", "tumor", "_reg_nonNA", reg_nonNA, ".txt")
original_tab <- fread(input = original_fp, data.table = F)

summary_mat <- matrix(data = NA, nrow = length(unique(original_tab$Cancer)), ncol = length(colnames(summary_mat_allBatches)),
                      dimnames = list(c(unique(original_tab$Cancer)), colnames(summary_mat_allBatches)))
summary_mat
for (cancer in unique(original_tab$Cancer)) {
  tab_tmp <- original_tab[original_tab$Cancer == cancer,]
  tab_tmp <- change_regression_nonNA(regression = tab_tmp, reg_nonNA = reg_nonNA, reg_sig = reg_sig)
  file_prefix_tmp <- unique(as.character(summary_mat_allBatches$file_prefix[summary_mat_allBatches$cancer == cancer]))
  original_new_fp <- paste0(makeOutDir(resultD = resultD), 
                            file_prefix_tmp, ".txt" )
  if (!file.exists(original_new_fp)) {
    write.table(x = tab_tmp, file = original_new_fp, quote = F, sep = "\t", row.names = F)
    
  }
  regulated_pairs <- unique(tab_tmp$pair[tab_tmp$regulated])
  tab_tmp <- annotate_ks_source(regression = tab_tmp)
  regulated_direct_pairs <- unique(tab_tmp$pair[tab_tmp$regulated & tab_tmp$is.direct])
  
  summary_mat[cancer, "size"] <- max(tab_tmp$Size)
  summary_mat[cancer, "cancer"] <- cancer
  summary_mat[cancer, "num_regulated"] <- length(regulated_pairs)
  summary_mat[cancer, "num_regulated_direct"] <- length(regulated_direct_pairs)
  summary_mat[cancer, "file_path"] <- original_new_fp
  summary_mat[cancer, "file_prefix"] <- file_prefix_tmp
  summary_mat[cancer, "is.present_batch"] <- T
  
}
summary_mat_allBatches <- summary_mat_allBatches[!is.na(summary_mat_allBatches$is.present_batch),]
## gather summary stats
summary_mat_allBatches <- rbind(summary_mat_allBatches, data.frame(summary_mat))
summary_mat_allBatches <- unique(summary_mat_allBatches)
## clean up
summary_mat_allBatches <- summary_mat_allBatches[!(grepl(x = summary_mat_allBatches$file_path, pattern = "OV") & grepl(x = summary_mat_allBatches$file_path, pattern = "downsize83")),]
summary_mat_allBatches$size <- as.numeric(as.vector(summary_mat_allBatches$size))
## write out summary stats

write.table(x = summary_mat_allBatches, file = paste0(makeOutDir(resultD = resultD), "summary_mat_allBatches.txt"), quote = F, sep = "\t", row.names = F)

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
names_stats <- c("size", "cancer", "num_regulated_overlap", "num_regulated", "num_iterations4overlap",
                 "is.present_batch", "file_prefix", 
                 "num_regulated_gain2overlap", 
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
  downsize_fps2process <- downsize_fps2process[1:min(length(downsize_fps2process), num_iteration2stop_overlap)]
  
  summary_mat <- matrix(data = NA, nrow = length(downsize_fps2process), ncol = length(names_stats),
                        dimnames = list(c(downsize_fps2process), names_stats))
  summary_mat
  downsize_sig_pairs_overlap <- NULL
  downsize_sig_pairs <- list()
  for (downsize_fp in downsize_fps2process) {
    ## input table
    tab_tmp <- fread(input = downsize_fp, data.table = F)
    ## annotate with enzyme type
    tab_tmp <- annotate_enzyme_type(regression = tab_tmp, kinases = kinases, phosphatases = phosphatases)
    
    ## change nonNA number and mark significance
    tab_tmp <- change_regression_nonNA(regression = tab_tmp, reg_nonNA = reg_nonNA, reg_sig = reg_sig)
    
    ## record significant pairs
    sig_pairs <- as.vector(unique(tab_tmp$pair[tab_tmp$regulated]))
    downsize_sig_pairs[[downsize_fp]] <- sig_pairs
    
    ## get the overlap of all iterations
    if (!is.null(downsize_sig_pairs_overlap)) {
      downsize_sig_pairs_overlap <- intersect(downsize_sig_pairs_overlap, sig_pairs)
    } else {
      downsize_sig_pairs_overlap <- sig_pairs
    }
    print(length(intersect(downsize_sig_pairs_overlap, sig_pairs)))
  }
  downsize_sig_pairs_overlap %>% head()
  
  ## say we have decided 5 iterations are goog enough to get the overlap
  for (downsize_fp in names(downsize_sig_pairs)) {
    ## go through all iteration belong to the same batch and get the difference to the overlap
    sig_pairs <- downsize_sig_pairs[[downsize_fp]]
    sig_pairs_gain <- setdiff(sig_pairs, downsize_sig_pairs_overlap)

    summary_mat[downsize_fp, "size"] <- as.numeric(str_split(string = downsize_fn_prefix, pattern = "downsize|\\/")[[1]][2])
    summary_mat[downsize_fp, "cancer"] <- str_split(string = downsize_fn_prefix, pattern = "downsize|\\/")[[1]][3]
    summary_mat[downsize_fp, "num_regulated"] <- length(sig_pairs)
    summary_mat[downsize_fp, "num_regulated_overlap"] <- length(downsize_sig_pairs_overlap)
    summary_mat[downsize_fp, "num_regulated_gain2overlap"] <- length(sig_pairs_gain)
    summary_mat[downsize_fp, "num_iterations4overlap"] <- length(names(downsize_sig_pairs))
    summary_mat[downsize_fp, "file_path"] <- downsize_fp
    summary_mat[downsize_fp, "file_prefix"] <- str_split(string = downsize_fn_prefix, pattern = "\\/")[[1]][3]
    
  }
  summary_mat <- data.frame(summary_mat, row.names = NULL)
  summary_mat$num_regulated_gain2overlap <- as.numeric(as.vector(summary_mat$num_regulated_gain2overlap))
  summary_mat$is.present_batch <- ifelse(summary_mat$num_regulated_gain2overlap == min(summary_mat$num_regulated_gain2overlap), T, F)
  
  summary_mat %>% head()

  ## get one iteration that is closest to the overlap and filter by overlapped pairs
  tab2w_input_fp <- as.character(unique(summary_mat$file_path[summary_mat$is.present_batch]))[1]
  tab2w_tmp <- fread(input = tab2w_input_fp, data.table = F)
  tab2w_tmp$is.iterations_overlap <- (tab2w_tmp$pair %in% downsize_sig_pairs_overlap)

  ## write table
  write.table(x = tab2w_tmp, file = paste0(makeOutDir(resultD = resultD), 
                                           str_split(string = downsize_fn_prefix, pattern = "\\/")[[1]][3],
                                           ".txt"), quote = F, row.names = F, col.names = F)
  
  ## gather summary stats
  summary_mat_allBatches <- rbind(summary_mat_allBatches, summary_mat)
}

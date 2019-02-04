# Yige Wu @ WashU 2019 Feb

# source ------------------------------------------------------------------
baseD = "/Users/yigewu/Box\ Sync/"
wd <- getwd()
if (wd != baseD) {
  setwd(baseD)
}
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_plotting.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
library(ggrepel)
library(dplyr)
# set variables -----------------------------------------------------------
# num_nonNA2test <- seq(from = 5, to = 120, by = 5)
# num_nonNA2test <- seq(from = 5, to = 120, by = 5)
reg_nonNA <- 20
downsize_list <- c(83)
cancers2process <- c("BRCA", "OV", "CO", "UCEC", "CCRCC")

# inputs ------------------------------------------------------------------
## input summary stats that tell you which regression result to use
summary_mat_allBatches <- fread("./cptac2p/analysis_results/phospho_network/regression/tables/gather_regression_downsize_tables_and_stats/summary_mat_allBatches.txt", data.table = F)

# Plot how we choose the samples during down-sizing -----------------------
for (size in unique(summary_mat_allBatches$size)) {
  ## create the data frame for plotting
  tab2p <- summary_mat_allBatches
  ## filtering
  tab2p <- tab2p[tab2p$size == size,]
  
  tab2p$y <- as.numeric(as.vector(tab2p$num_regulated_direct))
  tab2p$x <- get_order_by_id(value_vector = tab2p$y, id_vector = tab2p$cancer)
  tab2p$is.top <- get_top_by_id(value_vector = tab2p$y, id_vector = tab2p$cancer, num_top = 1)
  ## order the lines 
  
  tab2p %>% head
  
  ## plotting
  p <- ggplot()
  p <- p + geom_point(data = tab2p, mapping = aes(x = x, y = y, group = cancer))
  p <- p + geom_line(data = tab2p, mapping = aes(x = x, y = y, group = cancer, color = cancer))
  p <- p + geom_text_repel(data = tab2p[tab2p$is.top,], mapping = aes(x = x, y = y, label = cancer), linetype = 2)
  p <- p + theme_bw()
  p <- p + scale_color_manual(values = color_cancers2)
  p <- p + ylab("number of correlated directed kinase-substrate pairs\n(supported by experimental evidence)") + xlab(paste0("iteration of random sampling\n(", size, " samples)"))
  p <- p + scale_x_continuous(breaks = 1:max(tab2p$x))
  p <- p + guides(color = F)
  p
  ggsave(filename = paste0(makeOutDir(resultD = resultD), "downsize", size, "_correlation_overlap_direct_ks_pairs.pdf"),
         width = 5, height = 4)
}

# Plot how the # regulated change after down-sizing without limited to commonly detected pairs -----------------------
## get the pairs detected in all the indicated cancer types above (>= 20)
### make sure the sample size is the same
summary_mat_allBatches <- data.frame(summary_mat_allBatches)
for (size in c(83)) {
  pairs_detected_overlap <- NULL
  for (cancer in cancers2process) {
    ## input the regression result table the represent the downsampled cohort
    
    fp_tmp <- summary_mat_allBatches$file_path[(summary_mat_allBatches$size < (size + 2) & summary_mat_allBatches$size > (size - 2)) & summary_mat_allBatches$cancer == cancer & summary_mat_allBatches$is.present_batch ]
    fp_tmp <- unique(as.character(fp_tmp))[1]
    if (is.na(fp_tmp) | !file.exists(fp_tmp)) {
      next()
    }
    tab_tmp <- fread(input = fp_tmp, data.table = F)
    tab_tmp <- change_regression_nonNA(regression = tab_tmp, reg_nonNA = reg_nonNA, reg_sig = reg_sig)
    
    detected_pairs <- unique(tab_tmp$pair)
    if (!is.null(pairs_detected_overlap)) {
      pairs_detected_overlap <- unique(intersect(pairs_detected_overlap, detected_pairs))
    } else {
      pairs_detected_overlap <- detected_pairs
    }
  }
}
length(pairs_detected_overlap)
detected_overlap_summary_mat <- NULL
for (size in c(83, 96, 105, 123)) {
  ## only take these pairs and summarise how many within them are regualted
  names_stats <- c("size", "cancer", "num_regulated", "num_pairs_detected_overlap")
  ## the number of regulated pairs for each downsampled size
  summary_mat_tmp <- matrix(data = NA, nrow = length(cancers2process), ncol = length(names_stats),
                            dimnames = list(c(cancers2process), names_stats))
  summary_mat_tmp
  for (cancer in cancers2process) {
    ## input the regression result table the represent the downsampled cohort
    fp_tmp <- summary_mat_allBatches$file_path[(summary_mat_allBatches$size < (size + 2) & summary_mat_allBatches$size > (size - 2)) & summary_mat_allBatches$cancer == cancer & summary_mat_allBatches$is.present_batch ]
    fp_tmp <- unique(as.character(fp_tmp))[1]
    if (is.na(fp_tmp) | !file.exists(fp_tmp)) {
      next()
    }
    tab_tmp <- fread(input = fp_tmp, data.table = F)
    
    ## get the regulated pairs
    tab_tmp <- change_regression_nonNA(regression = tab_tmp, reg_nonNA = reg_nonNA, reg_sig = reg_sig)
    tab_tmp <- tab_tmp[tab_tmp$pair %in% pairs_detected_overlap,]
    head(tab_tmp)
    fp_tmp <- summary_mat_allBatches$file_prefix[summary_mat_allBatches$size == size & summary_mat_allBatches$cancer == cancer & summary_mat_allBatches$is.present_batch ]
    fp_tmp <- unique(as.character(fp_tmp))[1]
    
    file_path_tmp <- paste0(makeOutDir(resultD = resultD), "regression_", cancer, "_size", size, "_detected_in_", paste0(cancers2process, collapse = "_"),".txt")
    if (!file.exists(file_path_tmp)) {
      write.table(x = tab_tmp, file = file_path_tmp, quote = F, sep = "\t", row.names = F)
    }
    
    regulated_pairs <- unique(tab_tmp$pair[tab_tmp$regulated])
    summary_mat_tmp[cancer, "size"] <- size
    summary_mat_tmp[cancer, "cancer"] <- cancer
    summary_mat_tmp[cancer, "num_regulated"] <- length(regulated_pairs)
    summary_mat_tmp[cancer, "num_pairs_detected_overlap"] <- length(pairs_detected_overlap)
    
  }
  detected_overlap_summary_mat <- rbind(detected_overlap_summary_mat, data.frame(summary_mat_tmp, row.names = NULL))
}
detected_overlap_summary_mat <- detected_overlap_summary_mat[!is.na(detected_overlap_summary_mat$size),]
detected_overlap_summary_mat

## create the data frame for plotting
tab2p <- detected_overlap_summary_mat
## filtering

tab2p$y <- as.numeric(as.vector(tab2p$num_regulated))
tab2p$x <- as.numeric(as.vector(tab2p$size))
## order the lines 

tab2p %>% head

## plotting
p <- ggplot()
p <- p + geom_point(data = tab2p, mapping = aes(x = x, y = y, group = cancer))
p <- p + geom_line(data = tab2p, mapping = aes(x = x, y = y, group = cancer, color = cancer))
p <- p + geom_text_repel(data = tab2p[tab2p$size == 83,], mapping = aes(x = x, y = y, label = cancer, color = cancer), force = 2)
p <- p + theme_bw()
p <- p + scale_color_manual(values = color_cancers2)
p <- p + ylab("number of correlated kinase-substrate pairs") + xlab(paste0("number of samples\n(chosen by largest overlap between correlated and direct kinase target among 10 iterations)"))
# p <- p + scale_x_continuous(breaks = 1:max(tab2p$x))
p <- p + guides(color = F)
p
ggsave(filename = paste0(makeOutDir(resultD = resultD), "num_regulated_vs_downsample_size_limit_to_detected_in_all.pdf"),
       width = 5, height = 4)


# Plot how the # regulated change after down-sizing without limited to commonly detected pairs -----------------------
## get the pairs detected in all the indicated cancer types above (>= 20)
### make sure the sample size is the same
summary_mat_allBatches <- data.frame(summary_mat_allBatches)
detected_overlap_summary_mat <- NULL
for (size in c(83, 96, 105, 123)) {
  pairs_detected_overlap <- NULL
  for (cancer in cancers2process) {
    ## input the regression result table the represent the downsampled cohort

    fp_tmp <- summary_mat_allBatches$file_path[(summary_mat_allBatches$size < (size + 2) & summary_mat_allBatches$size > (size - 2)) & summary_mat_allBatches$cancer == cancer & summary_mat_allBatches$is.present_batch ]
    fp_tmp <- unique(as.character(fp_tmp))[1]
    if (is.na(fp_tmp) | !file.exists(fp_tmp)) {
      next()
    }
    tab_tmp <- fread(input = fp_tmp, data.table = F)
    tab_tmp <- change_regression_nonNA(regression = tab_tmp, reg_nonNA = reg_nonNA, reg_sig = reg_sig)
    
    detected_pairs <- unique(tab_tmp$pair)
    if (!is.null(pairs_detected_overlap)) {
      pairs_detected_overlap <- intersect(pairs_detected_overlap, detected_pairs)
    } else {
      pairs_detected_overlap <- detected_pairs
    }
  }

  ## only take these pairs and summarise how many within them are regualted
  names_stats <- c("size", "cancer", "num_regulated", "num_pairs_detected_overlap")
  ## the number of regulated pairs for each downsampled size
  summary_mat_tmp <- matrix(data = NA, nrow = length(cancers2process), ncol = length(names_stats),
                                     dimnames = list(c(cancers2process), names_stats))
  summary_mat_tmp
  for (cancer in cancers2process) {
    ## input the regression result table the represent the downsampled cohort
    fp_tmp <- summary_mat_allBatches$file_path[(summary_mat_allBatches$size < (size + 2) & summary_mat_allBatches$size > (size - 2)) & summary_mat_allBatches$cancer == cancer & summary_mat_allBatches$is.present_batch ]
    fp_tmp <- unique(as.character(fp_tmp))[1]
    if (is.na(fp_tmp) | !file.exists(fp_tmp)) {
      next()
    }
    tab_tmp <- fread(input = fp_tmp, data.table = F)

    ## get the regulated pairs
    tab_tmp <- change_regression_nonNA(regression = tab_tmp, reg_nonNA = reg_nonNA, reg_sig = reg_sig)
    tab_tmp <- tab_tmp[tab_tmp$pair %in% pairs_detected_overlap,]
    head(tab_tmp)
    fp_tmp <- summary_mat_allBatches$file_prefix[summary_mat_allBatches$size == size & summary_mat_allBatches$cancer == cancer & summary_mat_allBatches$is.present_batch ]
    fp_tmp <- unique(as.character(fp_tmp))[1]

    file_path_tmp <- paste0(makeOutDir(resultD = resultD), "regression_", cancer, "_size", size, "_detected_in_", paste0(cancers2process, collapse = "_"),".txt")
    if (!file.exists(file_path_tmp)) {
      # write.table(x = tab_tmp, file = file_path_tmp, quote = F, sep = "\t", row.names = F)
    }

    regulated_pairs <- unique(tab_tmp$pair[tab_tmp$regulated])
    summary_mat_tmp[cancer, "size"] <- size
    summary_mat_tmp[cancer, "cancer"] <- cancer
    summary_mat_tmp[cancer, "num_regulated"] <- length(regulated_pairs)
    summary_mat_tmp[cancer, "num_pairs_detected_overlap"] <- length(pairs_detected_overlap)

  }
  detected_overlap_summary_mat <- rbind(detected_overlap_summary_mat, data.frame(summary_mat_tmp, row.names = NULL))
}
detected_overlap_summary_mat <- detected_overlap_summary_mat[!is.na(detected_overlap_summary_mat$size),]
detected_overlap_summary_mat

## create the data frame for plotting
tab2p <- detected_overlap_summary_mat
## filtering

tab2p$y <- as.numeric(as.vector(tab2p$num_regulated))
tab2p$x <- as.numeric(as.vector(tab2p$size))
## order the lines

tab2p %>% head

## plotting
p <- ggplot()
p <- p + geom_point(data = tab2p, mapping = aes(x = x, y = y, group = cancer))
p <- p + geom_line(data = tab2p, mapping = aes(x = x, y = y, group = cancer, color = cancer))
p <- p + geom_text_repel(data = tab2p[tab2p$size == 83,], mapping = aes(x = x, y = y, label = cancer, color = cancer), force = 2)
p <- p + theme_bw()
p <- p + scale_color_manual(values = color_cancers2)
p <- p + ylab("number of correlated kinase-substrate pairs") + xlab(paste0("number of samples\n(chosen by largest overlap between\ncorrelated and direct kinase target)"))
# p <- p + scale_x_continuous(breaks = 1:max(tab2p$x))
p <- p + guides(color = F)
p
ggsave(filename = paste0(makeOutDir(resultD = resultD), "num_regulated_vs_downsample_size_not_limit_to_detected_in_all.pdf"),
       width = 5, height = 4)

# Yige Wu @ WashU 2018 Jan
# draw a grid showing the distribution of correlated kinase-substrate pairs are distributed across oncogenic pathways

# source ------------------------------------------------------------------
baseD = "/Users/yigewu/Box\ Sync/"
wd <- getwd()
if (wd != baseD) {
  setwd(baseD)
}
source("./cptac2p_analysis/phospho_network/phospho_network_shared.R")
source("./cptac2p_analysis/phospho_network/phospho_network_plotting.R")
library(pheatmap)

# set variables -----------------------------------------------------------
reg_nonNA <- 20
fdr_thres <- c(0.05, 0.1); names(fdr_thres) <- c("kinase", "phosphatase")
SELF = "trans"
num_top <- 10
enzyme_type <- "kinase"
# enzyme_type <- "phosphatase"

for (enzyme_type in c("kinase")) {
  summary_tab_tmp <- NULL
  for (cancer in cancers2process) {
    ## input the regression table for only pairs detectable in all cancers
    file_path_tmp <- paste0("./cptac2p/analysis_results/phospho_network/regression/figures/plot_downsize_affect_regression/", "regression_", cancer, "_size", size, "_detected_in_", paste0(cancers2process, collapse = "_"),".txt")
    tab_tmp <- fread(input = file_path_tmp, data.table = F)
    tab_tmp <- tab_tmp[order(tab_tmp$FDR_pho_kin, tab_tmp$pair, decreasing = F),]
    tab_tmp$GENE <- tab_tmp$KINASE
    tab_tmp$SUB_GENE <- tab_tmp$SUBSTRATE
    tab_tmp %>% head()
    tab_tmp <- tab_tmp[!duplicated(tab_tmp$pair),]
    rownames(tab_tmp) <- tab_tmp$pair
    ## filtering
    tab_tmp <- tab_tmp[tab_tmp$enzyme_type == enzyme_type & tab_tmp$SELF == SELF,]
    
    ## mark whether each cancer type is significantly correlated for each pair
    if (is.null(summary_tab_tmp)) {
      summary_tab_tmp <- tab_tmp[,c("pair", "SELF", "GENE", "SUB_GENE", "SUB_MOD_RSD")]
      rownames(summary_tab_tmp) <- summary_tab_tmp$pair
    }
    summary_tab_tmp[, paste0("regulated_", cancer)] <- as.numeric(tab_tmp[as.vector(summary_tab_tmp$pair), "regulated"])
  }
  regression <- NULL
  for (cancer in cancers2process) {
    ## input the regression table for only pairs detectable in all cancers
    file_path_tmp <- paste0("./cptac2p/analysis_results/phospho_network/regression/figures/plot_downsize_affect_regression/", "regression_", cancer, "_size", size, "_detected_in_", paste0(cancers2process, collapse = "_"),".txt")
    tab_tmp <- fread(input = file_path_tmp, data.table = F)
    tab_tmp <- tab_tmp[order(tab_tmp$FDR_pho_kin, tab_tmp$pair, decreasing = F),]
    tab_tmp$GENE <- tab_tmp$KINASE
    tab_tmp$SUB_GENE <- tab_tmp$SUBSTRATE
    tab_tmp %>% head()
    tab_tmp <- tab_tmp[!duplicated(tab_tmp$pair),]
    rownames(tab_tmp) <- tab_tmp$pair
    
    tab_tmp$regulated_uniq <- ifelse(tab_tmp$pair %in% summary_tab_tmp$pair[summary_tab_tmp[,paste0("regulated_", cancer)] == 1], T, F)
    
    if (is.null(regression)) {
      regression <- tab_tmp
      ordered_cols <- colnames(regression)
    } else {
      regression <- rbind(regression, tab_tmp[,ordered_cols])
      
    }
  }
  regression %>% tail
  
  tab_pairs <- data.frame(table(regression[, c("GENE", "Cancer", "regulated_uniq")]))
  tab_pairs %>%
    tail()

  tab_pairs.wratio <- merge(tab_pairs[tab_pairs$regulated_uniq == "TRUE",],
                            tab_pairs[tab_pairs$regulated_uniq == "FALSE", c("GENE", "Cancer", "Freq")],
                            by = c("GENE", "Cancer"), suffixes = c(".regulated_uniqT", ".regulated_uniqF"), all.y = T)
  tab_pairs.wratio %>%
    tail()
  tab_pairs.wratio$ratio.regulated_uniqT <- tab_pairs.wratio$Freq.regulated_uniqT/(tab_pairs.wratio$Freq.regulated_uniqT +  tab_pairs.wratio$Freq.regulated_uniqF)
  tab_pairs.wratio %>%
    head()

  # plotting the top ratio kinases alone----------------------------------------------------------------
  ## compile data frame for plotting
  tab2p <- tab_pairs.wratio
  tab2p$y <- tab2p$Cancer
  tab2p$x <- tab2p$GENE
  tab2p$fill <- tab2p$ratio.regulated_uniqT
  
  ## filtering
  tab2p <- tab2p[!is.na(tab2p$x) & !is.na(tab2p$y),]
  tab2p %>%
    head()
  boxplot(tab2p$fill)
  # if (enzyme_type == "kinase") {
  #   top_ratios <- get_top_by_id(value_vector = tab2p$ratio.regulated_uniqT[tab2p$Freq.regulated_uniqT >= 5], id_vector = tab2p$Cancer[tab2p$Freq.regulated_uniqT >= 5], num_top = num_top)
  # }
  ### get the top 3 of each cancer type
  tab2p$is.top <- get_top_by_id(value_vector = tab2p$fill, id_vector = tab2p$Cancer, num_top = num_top)
  
  
  tab2p <- tab2p[!is.na(tab2p$fill) & tab2p$fill > 0,]
  
  ## ordering
  tab2p$y <- order_cancer(tab2p$y)
  
  df_value <- dcast(data = tab2p, y ~ x, value.var = "fill")
  mat_value <- as.matrix(df_value[,-1])
  rownames(mat_value) <- df_value$y
  head(mat_value)
  mat_value[is.na(mat_value)] <- 0
  
  df_label <- dcast(data = tab2p, y ~ x, value.var = "Freq.regulated_uniqT")
  mat_label <- as.matrix(df_label[,-1])
  rownames(mat_label) <- df_label$y
  if (enzyme_type == "kinase") {
    mat_value[mat_label < 5] <- 0
    mat_value <- mat_value[, colSums(mat_label, na.rm = T) >= 5]
    mat_label <- mat_label[, colSums(mat_label, na.rm = T) >= 5]
  }
  
  mat_label[is.na(mat_label)] <- ""
  
  
  fn = paste(makeOutDir(resultD = resultD), enzyme_type, '_with_top', num_top, '_ratio.pdf',sep = "")
  
  if (enzyme_type == "kinase") {
    my_heatmap <- pheatmap(mat_value, 
                           color = c("white", col_paletteR(100)),
                           na_col = "white", 
                           display_numbers = mat_label,
                           cluster_rows=F,
                           show_colnames = T)
    save_pheatmap_pdf(x = my_heatmap, 
                      filename = fn, 
                      width = 9, height = 3)
  } else {
    my_heatmap <- pheatmap(mat_value, 
                           color = c("white", col_paletteB(100)),
                           na_col = "white", 
                           display_numbers = mat_label,
                           cluster_rows=F,
                           show_colnames = T)
    save_pheatmap_pdf(x = my_heatmap, 
                      filename = fn, 
                      width = 6, height = 3)
  }
  
  
  
}

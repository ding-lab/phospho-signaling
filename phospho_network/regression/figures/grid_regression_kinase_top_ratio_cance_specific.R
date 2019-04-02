# Yige Wu @ WashU 2018 Jan
# draw a grid showing the distribution of correlated kinase-substrate pairs are distributed across oncogenic pathways

# source ------------------------------------------------------------------
baseD = "/Users/yigewu/Box\ Sync/"
setwd(baseD)

source("./cptac2p_analysis/phospho_network/phospho_network_shared.R")
source("./cptac2p_analysis/phospho_network/phospho_network_plotting.R")
library(pheatmap)

# set variables -----------------------------------------------------------
reg_nonNA <- 20
# fdr_thres <- c(0.05, 0.05); names(fdr_thres) <- c("kinase", "phosphatase")
fdr_thres <- c(0.05, 0.05); names(fdr_thres) <- c("kinase", "phosphatase")

SELF = "trans"
num_top <- 25
enzyme_type <- "kinase"
# enzyme_type <- "phosphatase"
size <- 83
cancers2process <- c("BRCA", "OV", "CO", "UCEC", "CCRCC")

file_path_tmp <- paste0("./cptac2p/analysis_results/phospho_network/regression/tables/generate_regression_regulated_uniq_marked/",  "regression_size", size, "_FDR", fdr_thres[1], "_detected_in_", paste0(cancers2process, collapse = "_"),".txt")
regression <- fread(input = file_path_tmp, data.table = F)
regression <- adjust_regression_by_nonNA(regression = regression, reg_nonNA = 20, reg_sig = fdr_thres)

regression$pair_pro <- paste0(regression$GENE, ":", regression$SUB_GENE)
regression <- regression %>%
  filter(pair_pro %in% omnipath_tab$pair_pro[!(omnipath_tab$Source %in% c("MIMP", "NetKIN", "PhosphoNetworks"))])

for (enzyme_type in c("kinase")) {
  source("./cptac2p_analysis/phospho_network/regression/tables/generate_regression_regulated_uniq_marked.R")
  # file_path_tmp <- paste0("./cptac2p/analysis_results/phospho_network/regression/figures/grid_regression_kinase_top_ratio_cance_specific/", "regression_", enzyme_type, "_substrate_size", size, "_detected_in_", paste0(cancers2process, collapse = "_"),".txt")
  # regression <- fread(input = file_path_tmp, data.table = F)
  
  for (col2evaluate in c("GENE", "SUB_GENE")) {
    tab_pairs <- data.frame(table(regression[, c(col2evaluate, "Cancer", "regulated_uniq")]))
    tab_pairs %>%
      tail()
    
    tab_pairs.wratio <- merge(tab_pairs[tab_pairs$regulated_uniq == "TRUE",],
                              tab_pairs[tab_pairs$regulated_uniq == "FALSE", c(col2evaluate, "Cancer", "Freq")],
                              by = c(col2evaluate, "Cancer"), suffixes = c(".regulated_uniqT", ".regulated_uniqF"), all.y = T)
    tab_pairs.wratio %>%
      tail()
    tab_pairs.wratio$ratio.regulated_uniqT <- tab_pairs.wratio$Freq.regulated_uniqT/(tab_pairs.wratio$Freq.regulated_uniqT +  tab_pairs.wratio$Freq.regulated_uniqF)
    tab_pairs.wratio %>%
      head()
    
    for (cancer_tmp in cancers2process) {
      tab4GSEA <- tab_pairs.wratio %>%
        filter(Cancer == cancer_tmp) %>%
        filter(ratio.regulated_uniqT > 0)
      tab4GSEA <- tab4GSEA[, c(col2evaluate, "ratio.regulated_uniqT")]
      write.table(x = tab4GSEA, file = paste0(makeOutDir(resultD = resultD), cancer_tmp, "_", col2evaluate, "_regulated_uniq_ratio2detected.txt"), quote = F, row.names = F, col.names = F)
    }
    
    # plotting the top ratio kinases alone----------------------------------------------------------------
    ## compile data frame for plotting
    tab2p <- tab_pairs.wratio
    tab2p$y <- tab2p$Cancer
    tab2p$x <- tab2p[,col2evaluate]
    tab2p$fill <- tab2p$ratio.regulated_uniqT
    
    ## filtering
    tab2p <- tab2p[!is.na(tab2p$x) & !is.na(tab2p$y),]
    tab2p %>%
      head()
    boxplot(tab2p$fill)
    tab2p$is.top <- F
    
    if (enzyme_type == "kinase") {
      tab2p <- tab2p[tab2p$GENE != "AKT1S1",]
      tab2p$is.top[tab2p$Freq.regulated_uniqF >= 5] <- get_top_by_id(value_vector = tab2p$ratio.regulated_uniqT[tab2p$Freq.regulated_uniqF >= 5], id_vector = tab2p$Cancer[tab2p$Freq.regulated_uniqF >= 5], num_top = num_top)
    } else {
      tab2p$is.top <- get_top_by_id(value_vector = tab2p$fill, id_vector = tab2p$Cancer, num_top = num_top)
      
    }
    stop("")
    ### get the top 3 of each cancer type
    
    ## filtering
    tab2p <- tab2p[tab2p[, col2evaluate] %in% tab2p[tab2p$is.top, col2evaluate],]
    
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
    
    mat_label[is.na(mat_label)] <- ""
    
    fn = paste(makeOutDir(resultD = resultD), enzyme_type, '_with_top', num_top, "_", col2evaluate, '_regulated_uniq_ratio.pdf',sep = "")
    
    if (enzyme_type == "kinase") {
      my_heatmap <- pheatmap(mat_value, 
                             color = c("white", col_paletteR(100)),
                             na_col = "white", 
                             display_numbers = mat_label,
                             cluster_rows=F,
                             show_colnames = T)
      save_pheatmap_pdf(x = my_heatmap, 
                        filename = fn, 
                        width = 6, height = 3)
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
}
}
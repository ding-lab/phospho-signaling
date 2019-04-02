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
# fdr_thres <- c(0.05, 0.1); names(fdr_thres) <- c("kinase", "phosphatase")
fdr_thres <- c(0.05, 0.05); names(fdr_thres) <- c("kinase", "phosphatase")
# stop("phosphatase threshold!")
SELF = "trans"
num_top <- 10
num_top <- 15
# num_top <- 5
qt_top <- 0.5
# enzyme_type <- "kinase"
enzyme_type <- "phosphatase"


# input regression table --------------------------------------------------
regression_sup <- fread(input = paste0(ppnD, "regression/tables/annotate_regression_with_mut_impact/", 
                                       "regression_cptac2p_cptac3_tumor_reg_nonNA20_mut_impact_cancer_specificity_annotated.txt"),
                        data.table = F)
regression_sup <- regression_sup %>%
  filter(pair_pro %in% omnipath_tab$pair_pro[!(omnipath_tab$Source %in% c("NetKIN", "PhosphoNetworks", "MIMP"))] | pair_pro %in% psp_tab$pair_pro)

regression_sup <- adjust_regression_by_nonNA(regression = regression_sup, reg_nonNA = 20, reg_sig = reg_sig)
regression_sup <- annotate_ks_source(regression = regression_sup)

drug_genes <- fread(input = "./Ding_Lab/Projects_Current/CPTAC/pan3can_shared_data/reference_files/gene_drug_list/Premed_raw_databases/drugBank/drug_list.txt", data.table = F, col.names = "gene")

for (SELF_tmp in c("trans")) {
  # first plot site-level comparison with known ones ------------------------
  for (enzyme_type_tmp in c("kinase", "phosphatase")) {
    regression <- regression_sup %>%
      filter(enzyme_type == enzyme_type_tmp) %>%
      filter(SELF == SELF_tmp)
    
    # get table of how many pairs are annotated to different pathways ---------
    tab_pairs <- data.frame(table(regression[, c("GENE", "Cancer", "regulated")]))
    tab_pairs %>%
      head()
    
    tab_pairs.wratio <- merge(tab_pairs[tab_pairs$regulated == "TRUE",], 
                              tab_pairs[tab_pairs$regulated == "FALSE", c("GENE", "Cancer", "Freq")],
                              by = c("GENE", "Cancer"), suffixes = c(".regulatedT", ".regulatedF"), all.y = T)
    tab_pairs.wratio %>%
      head()
    tab_pairs.wratio$ratio.regulatedT <- tab_pairs.wratio$Freq.regulatedT/(tab_pairs.wratio$Freq.regulatedT +  tab_pairs.wratio$Freq.regulatedF)
    tab_pairs.wratio %>%
      head()
    
    # Plot range of the number of substrate phoshposites per kinase -----------------------------------
    tab2p <- tab_pairs.wratio
    tab2p <- tab2p[(tab2p$Freq.regulatedF + tab2p$Freq.regulatedT) > 0,]
    # tab2p <- data.frame(table(tab2p[, c("Freq.regulatedT", "Cancer")]))
    tab2p %>%
      head()
    tab2p$x <- tab2p$Freq.regulatedT
    tab2p$Cancer <- order_cancer(tab2p$Cancer)
    
    ## filtering
    tab2p <- tab2p[tab2p$x < 100 & tab2p$x > 0,]
    
    p <- ggplot(data = tab2p)
    p <- p + geom_histogram(mapping = aes(x = x, fill = Cancer, group = Cancer), binwidth = 5)
    p <- p + scale_fill_manual(values = color_cancers2)
    p <- p + facet_grid(Cancer~.)
    p <- p + theme_bw()
    p
    fn = paste(makeOutDir(resultD = resultD), SELF_tmp, "_",  enzyme_type_tmp, '_barplot_num_substrates.pdf',sep = "")
    ggsave(filename = fn, width = 6, height = 4)
    
    
    # plotting tyrosine vs ST -------------------------------------------------
    tab2p <- regression
    tab2p$GENE_STY <- ifelse(tab2p$GENE %in% TK_list, "Y", "ST")
    tab2p$RSD <- ifelse(substr(tab2p$SUB_MOD_RSD, start = 1, stop = 1) == "Y", "Y", "ST")
    tab2p <- data.frame(table(tab2p[tab2p$regulated, c("Cancer", "GENE_STY", "RSD")]))
    tab2p %>%
      head()
    tab2p$num_pairs_per_cancer <- sapply(tab2p$Cancer, FUN = function(cancer, tab2p) sum(tab2p$Freq[tab2p$Cancer == cancer]), tab2p)
    tab2p$group <- paste0(tab2p$GENE_STY, "_kinase->", tab2p$RSD)
    tab2p$y <- tab2p$Freq/tab2p$num_pairs_per_cancer
    tab2p$x <- order_cancer(tab2p$Cancer)
    
    p <- ggplot(data = tab2p)
    p <- p + geom_bar(mapping = aes(x = x, y = y, fill = group), stat = "identity")
    p <- p + theme_bw()
    p <- p + coord_flip()
    p <- p + theme(axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position = "right")
    p
    fn = paste(makeOutDir(resultD = resultD), SELF_tmp, "_",  enzyme_type_tmp, '_barplot_ratio_STY.pdf',sep = "")
    ggsave(filename = fn, width = 6, height = 4)
    
    # plotting the top ratio kinases----------------------------------------------------------------
    ## compile data frame for plotting
    summary(tab_pairs.wratio$ratio.regulatedT)
    
    tab2p <- tab_pairs.wratio
    tab2p$y <- tab2p$Cancer
    tab2p$x <- tab2p$GENE
    tab2p$fill <- tab2p$ratio.regulatedT
    
    ## filtering
    tab2p <- tab2p[!is.na(tab2p$x) & !is.na(tab2p$y),]
    tab2p %>%
      head()
    
    tab2p$is.top <- F
    tab2p$is.top[tab2p$Freq.regulatedT >= 5] <- get_top_by_id(value_vector = tab2p$ratio.regulatedT[tab2p$Freq.regulatedT >= 5], id_vector = tab2p$Cancer[tab2p$Freq.regulatedT >= 5], num_top = num_top)
    # top_ratios <- get_top_by_id(value_vector = tab2p$ratio.regulatedT, id_vector = tab2p$Cancer, num_top = num_top)
    ### get the top 3 of each cancer type
    tab2p <- tab2p[tab2p$GENE %in% tab2p$GENE[tab2p$is.top],]
    
    tab2p <- tab2p[!is.na(tab2p$ratio.regulatedT) & tab2p$ratio.regulatedT > 0,]
    
    ## ordering
    tab2p$y <- order_cancer(tab2p$y)
    
    df_value <- dcast(data = tab2p, y ~ x, value.var = "fill")
    mat_value <- as.matrix(df_value[,-1])
    rownames(mat_value) <- df_value$y
    head(mat_value)
    mat_value[is.na(mat_value)] <- 0
    
    df_label <- dcast(data = tab2p, y ~ x, value.var = "Freq.regulatedT")
    mat_label <- as.matrix(df_label[,-1])
    rownames(mat_label) <- df_label$y
    mat_value[mat_label < 5] <- 0
    mat_value <- mat_value[, colSums(mat_label, na.rm = T) >= 5]
    mat_label <- mat_label[, colSums(mat_label, na.rm = T) >= 5]
    
    mat_label[is.na(mat_label)] <- ""
    
    
    fn = paste(makeOutDir(resultD = resultD), SELF_tmp, "_",  enzyme_type_tmp, '_with_top', num_top, '_ratio.pdf',sep = "")
    my_heatmap <- pheatmap(mat_value, 
                           color = c("white", col_paletteR(100)),
                           na_col = "white", 
                           display_numbers = mat_label,
                           cluster_rows=F,
                           show_colnames = T)
    save_pheatmap_pdf(x = my_heatmap, 
                      filename = fn, 
                      width = 10, height = 3)
    
    
    # plotting the kinases with ratio of regulated substrate within top 25% quantile----------------------------------------------------------------
    ## compile data frame for plotting
    tab2p <- tab_pairs.wratio
    tab2p$y <- tab2p$Cancer
    tab2p$x <- tab2p$GENE
    tab2p$fill <- tab2p$ratio.regulatedT
    
    ## filtering
    tab2p <- tab2p[!is.na(tab2p$x) & !is.na(tab2p$y),]
    tab2p %>%
      head()
    
    tab2p$is.top <- F
    tab2p$is.top[tab2p$Freq.regulatedT >= 5] <- get_top_QT_by_id(value_vector = tab2p$ratio.regulatedT[tab2p$Freq.regulatedT >= 5], id_vector = tab2p$Cancer[tab2p$Freq.regulatedT >= 5], qt = qt_top)
    
    # top_ratios <- get_top_by_id(value_vector = tab2p$ratio.regulatedT, id_vector = tab2p$Cancer, num_top = num_top)
    ### get the top 3 of each cancer type
    tab2p <- tab2p[tab2p$GENE %in% tab2p$GENE[tab2p$is.top],]
    
    tab2p <- tab2p[!is.na(tab2p$ratio.regulatedT) & tab2p$ratio.regulatedT > 0,]
    
    ## ordering
    tab2p$y <- order_cancer(tab2p$y)
    
    df_value <- dcast(data = tab2p, y ~ x, value.var = "fill")
    mat_value <- as.matrix(df_value[,-1])
    rownames(mat_value) <- df_value$y
    head(mat_value)
    mat_value[is.na(mat_value)] <- 0
    
    mat_colnames <- colnames(mat_value)
    mat_colnames[mat_colnames %in% drug_genes$gene] <- paste0(mat_colnames[mat_colnames %in% drug_genes$gene], "*")
    colnames(mat_value) <- mat_colnames
    
    df_label <- dcast(data = tab2p, y ~ x, value.var = "Freq.regulatedT")
    mat_label <- as.matrix(df_label[,-1])
    rownames(mat_label) <- df_label$y
    mat_value[mat_label < 5] <- 0
    mat_value <- mat_value[, colSums(mat_label, na.rm = T) >= 5]
    mat_label <- mat_label[, colSums(mat_label, na.rm = T) >= 5]
    
    mat_label[is.na(mat_label)] <- ""
    
    
    fn = paste(makeOutDir(resultD = resultD), SELF_tmp, "_",  enzyme_type_tmp, '_with_topQT', qt_top, '_ratio.pdf',sep = "")
    my_heatmap <- pheatmap(mat_value, 
                           color = c("white", col_paletteR(100)),
                           na_col = "white", 
                           display_numbers = mat_label,
                           cluster_rows=F,
                           show_colnames = T)
    if (enzyme_type == "kinase") {
      save_pheatmap_pdf(x = my_heatmap, 
                        filename = fn, 
                        width = 10, height = 3)
    }
    if (enzyme_type == "kinase") {
      save_pheatmap_pdf(x = my_heatmap, 
                        filename = fn, 
                        width = 3, height = 2.5)
    }

    
    
  }
}

#
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
fdr_thres <- c(0.05, 0.2); names(fdr_thres) <- c("kinase", "phosphatase")
enzyme_type <- "kinase"
# enzyme_type <- "phosphatase"

SELF = "trans"
num_top <- 10

# inputs ------------------------------------------------------------------
regression <- fread(input = paste0(ppnD, "regression/tables/", "change_regression_nonNA/", 
                                   "regression", "_", "cptac2p_cptac3", "_", "tumor", "_reg_nonNA", reg_nonNA, ".txt"), data.table = F)
## filtering
regression <- regression[regression$enzyme_type == enzyme_type & regression$SELF == SELF,]
regression <- markSigKS(regression = regression, sig_thres = fdr_thres[enzyme_type], enzyme_type = enzyme_type)

# get table of how many pairs are annotated to different pathways ---------
tab_pairs <- data.frame(table(regression[, c("SUB_GENE", "Cancer", "regulated")]))
tab_pairs %>%
  head()

tab_pairs.wratio <- merge(tab_pairs[tab_pairs$regulated == "TRUE",], 
                          tab_pairs[tab_pairs$regulated == "FALSE", c("SUB_GENE", "Cancer", "Freq")],
                          by = c("SUB_GENE", "Cancer"), suffixes = c(".regulatedT", ".regulatedF"), all.y = T)
tab_pairs.wratio %>%
  head()
tab_pairs.wratio$ratio.regulatedT <- tab_pairs.wratio$Freq.regulatedT/(tab_pairs.wratio$Freq.regulatedT +  tab_pairs.wratio$Freq.regulatedF)
tab_pairs.wratio %>%
  tail()


# compile data frame for plotting -----------------------------------------
tab2p <- tab_pairs.wratio
tab2p$y <- tab2p$Cancer
tab2p$x <- tab2p$SUB_GENE
tab2p$fill <- tab2p$ratio.regulatedT

## filtering
tab2p <- tab2p[!is.na(tab2p$x) & !is.na(tab2p$y),]
tab2p %>%
  head()
boxplot(tab2p$fill)
if (enzyme_type == "kinase") {
  tab2p_for_top <- tab2p[tab2p$Freq.regulatedT >= 5,]
  tab2p_for_top$top <- get_top_by_id(value_vector = tab2p_for_top$ratio.regulatedT, id_vector = tab2p_for_top$Cancer, num_top = num_top)
  tab2p_for_top <- tab2p_for_top[tab2p_for_top$top,]
} else {
  tab2p$top <- get_top_by_id(value_vector = tab2p$ratio.regulatedT, id_vector = tab2p$Cancer, num_top = num_top)
  tab2p_for_top <- tab2p[tab2p$top,]
}
### get the top 3 of each cancer type
tab2p <- tab2p[tab2p$SUB_GENE %in% tab2p_for_top$SUB_GENE,]
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
if (enzyme_type == "kinase") {
  mat_value[mat_label < 5] <- 0
  mat_value <- mat_value[, colSums(mat_label, na.rm = T) >= 5]
  mat_label <- mat_label[, colSums(mat_label, na.rm = T) >= 5]
}

mat_label[is.na(mat_label)] <- ""


fn = paste(makeOutDir(resultD = resultD), "substrate_of_" , enzyme_type, '_with_top', num_top, '_ratio.pdf',sep = "")

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
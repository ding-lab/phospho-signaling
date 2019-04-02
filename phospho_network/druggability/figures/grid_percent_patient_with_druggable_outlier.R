# Yige Wu @ WashU 2018 Jan
# draw a grid showing the distribution of correlated kinase-substrate pairs are distributed across oncogenic pathways

# source ------------------------------------------------------------------
baseD = "/Users/yigewu/Box\ Sync/"
setwd(baseD)
source('./cptac2p_analysis/preprocess_files/preprocess_files_shared.R')
source("./cptac2p_analysis/phospho_network/phospho_network_shared.R")
source("./cptac2p_analysis/phospho_network/phospho_network_plotting.R")
library(pheatmap)

# set variables -----------------------------------------------------------
reg_nonNA <- 20
num_top <- 10
outlier_sd <- 1.5
cancers2process <- c("BRCA", "OV", "CO", "UCEC", "CCRCC")
enzyme_type <- "kinase"
reg_sig <- c(0.05, 0.05); names(reg_sig) <- c("kinase", "phosphatase")

# input druggable gene list -----------------------------------------------
# drug_genes <- fread(input = "./Ding_Lab/Projects_Current/CPTAC/pan3can_shared_data/reference_files/gene_drug_list/Premed_raw_databases/drugBank/all_target_ids_all.txt.human.tsv_hugoified.tsv", data.table = F)
# drug_genes %>%
#   head
# drug_genes <- as.vector(drug_genes$V3)

drug_genes <- read.delim(file = "./Ding_Lab/Projects_Current/CPTAC/pan3can_shared_data/reference_files/gene_drug_list/Premed_raw_databases/drugBank/drug_list.txt", header = F, col.names = "gene")
drug_genes <- as.vector(drug_genes$gene)

# input regression --------------------------------------------------------
regression <- fread(input = paste0(ppnD, "regression/tables/annotate_regression_with_mut_impact/", 
                                   "regression_cptac2p_cptac3_tumor_reg_nonNA20_mut_impact_cancer_specificity_annotated.txt"),
                    data.table = F)
regression <- regression %>%
  filter(pair_pro %in% omnipath_tab$pair_pro[!(omnipath_tab$Source %in% c("NetKIN", "PhosphoNetworks", "MIMP"))] | pair_pro %in% psp_tab$pair_pro)

regression <- adjust_regression_by_nonNA(regression = regression, reg_nonNA = 20, reg_sig = reg_sig)
regression <- annotate_ks_source(regression = regression)


# calculate % cases -------------------------------------------------------
esscore_tab_outlier_drug <- NULL
num_partIDs <- NULL
for (cancer in cancers2process) {
  esscore_tab <- fread(input = paste0("./cptac2p/analysis_results/phospho_network/es_score/table/calculate_es_pair_score/kinase/esscore_tab_", cancer, "_kinase_reg_nonNA20.txt"), data.table = F)
  partIDs <- colnames(esscore_tab); partIDs <- partIDs[!(partIDs %in% c("GENE", "SUB_GENE", "SUB_MOD_RSD", "SELF", "pair", "site_id"))]
  esscore_tab <- esscore_tab[esscore_tab$GENE %in% drug_genes, ]
  esscore_tab <- annotate_ks_source(regression = esscore_tab)
  esscore_tab <- esscore_tab %>%
    mutate(pair_cancer = paste0(pair, ":", cancer)) %>%
    # filter(is.direct == T) %>%
    filter(pair_cancer %in% regression$pair_cancer[regression$regulated])
  
  esscore_tab <- esscore_tab[!duplicated(esscore_tab$pair),]
  pairs2generate_head <- esscore_tab[, c("GENE", "SUB_GENE", "SUB_MOD_RSD", "SELF", "pair")]
  pairs2generate_head %>%
    head()
  esscore_tab_scaled <- esscore_tab[, partIDs]
  esscore_tab_scaled <- scale_by_row(esscore_tab_scaled)
  esscore_tab_outlier <- cbind(pairs2generate_head, (esscore_tab_scaled > outlier_sd))
  colnames(esscore_tab_outlier) <- c(colnames(pairs2generate_head), partIDs)
  
  esscore_tab_outlier.m <- melt.data.frame(esscore_tab_outlier, 
                                           id.vars = c("GENE", "SUB_GENE", "SUB_MOD_RSD", "SELF", "pair"))
  esscore_tab_outlier.m %>%
    head()
  esscore_tab_outlier.m <- esscore_tab_outlier.m %>%
    filter(!is.na(value) & value == T)
  
  esscore_tab_outlier.m$cancer <- cancer
  esscore_tab_outlier.m %>%
    head()
  esscore_tab_outlier_drug <- rbind(esscore_tab_outlier_drug, esscore_tab_outlier.m)
  num_partIDs[cancer] <- length(partIDs)
}
write.table(x = esscore_tab_outlier_drug, file = paste0(makeOutDir(resultD = resultD), "esscore_tab_outlier_drug_genes.txt"), quote = F, sep = "\t", row.names = F)




# for all druggable genes cis/trans tgt showing ratio ---------
tab_pairs <- esscore_tab_outlier_drug %>%
  select(GENE, cancer, variable) %>%
  unique() %>%
  select(GENE, cancer) %>%
  table %>%
  data.frame()
tab_pairs %>%
  head()
tmp <- data.frame(num_partIDs = num_partIDs, cancer = names(num_partIDs))
tab_pairs.wratio <- merge(tab_pairs, tmp, by = c("cancer"), all.x = T)
tab_pairs.wratio %>%
  head()
tab_pairs.wratio$ratio.outlier <- tab_pairs.wratio$Freq/(tab_pairs.wratio$num_partIDs)
tab_pairs.wratio %>%
  head()

## compile data frame for plotting
tab2p <- tab_pairs.wratio
tab2p$y <- tab2p$cancer
tab2p$x <- tab2p$GENE
tab2p$fill <- as.numeric(as.vector(tab2p$ratio.outlier))

## filtering
tab2p <- tab2p[!is.na(tab2p$x) & !is.na(tab2p$y),]
tab2p %>%
  head()
boxplot(tab2p$fill)

tab2p <- tab2p[!is.na(tab2p$ratio.outlier) & tab2p$ratio.outlier > 0,]

## ordering
tab2p$y <- order_cancer(tab2p$y)
tab2p$percentage <- paste0(format(100*tab2p$fill, digits = 1))
df_value <- dcast(data = tab2p, y ~ x, value.var = "fill")
mat_value <- as.matrix(df_value[,-1])
rownames(mat_value) <- df_value$y
head(mat_value)
mat_value[is.na(mat_value)] <- 0

df_label <- dcast(data = tab2p, y ~ x, value.var = "percentage")
mat_label <- as.matrix(df_label[,-1])
rownames(mat_label) <- df_label$y
# if (enzyme_type == "kinase") {
#   mat_value[mat_label < 5] <- 0
#   # mat_value <- mat_value[, colSums(mat_label, na.rm = T) >= 5]
#   # mat_label <- mat_label[, colSums(mat_label, na.rm = T) >= 5]
# }

mat_label[is.na(mat_label)] <- ""
col_order4 <- c("AKT1", "MTOR", "MAPK1", "MAP2K1")
# mat_value <- mat_value[, c(col_order4, colnames(mat_value)[!(colnames(mat_value) %in% col_order4)])]
# mat_label <- mat_label[, c(col_order4, colnames(mat_value)[!(colnames(mat_value) %in% col_order4)])]
fn = paste(makeOutDir(resultD = resultD), "druggable_", enzyme_type, '_outlier_ratio_in_patients_show_percent.pdf',sep = "")

my_heatmap <- pheatmap(mat_value, 
                       color = c("white", col_paletteR(100)),
                       na_col = "white", 
                       display_numbers = mat_label,
                       cluster_rows=F,
                       cluster_cols = F,
                       show_colnames = T)
save_pheatmap_pdf(x = my_heatmap, 
                  filename = fn, 
                  width = 8, height = 2)









# for all druggable genes cis/trans tgt ---------
tab_pairs <- esscore_tab_outlier_drug %>%
  select(GENE, cancer, variable) %>%
  unique() %>%
  select(GENE, cancer) %>%
  table %>%
  data.frame()
tab_pairs %>%
  head()
tmp <- data.frame(num_partIDs = num_partIDs, cancer = names(num_partIDs))
tab_pairs.wratio <- merge(tab_pairs, tmp, by = c("cancer"), all.x = T)
tab_pairs.wratio %>%
  head()
tab_pairs.wratio$ratio.outlier <- tab_pairs.wratio$Freq/(tab_pairs.wratio$num_partIDs)
tab_pairs.wratio %>%
  head()

## compile data frame for plotting
tab2p <- tab_pairs.wratio
tab2p$y <- tab2p$cancer
tab2p$x <- tab2p$GENE
tab2p$fill <- as.numeric(as.vector(tab2p$ratio.outlier))

## filtering
tab2p <- tab2p[!is.na(tab2p$x) & !is.na(tab2p$y),]
tab2p %>%
  head()
boxplot(tab2p$fill)

tab2p <- tab2p[!is.na(tab2p$ratio.outlier) & tab2p$ratio.outlier > 0,]

## ordering
tab2p$y <- order_cancer(tab2p$y)

df_value <- dcast(data = tab2p, y ~ x, value.var = "fill")
mat_value <- as.matrix(df_value[,-1])
rownames(mat_value) <- df_value$y
head(mat_value)
mat_value[is.na(mat_value)] <- 0

df_label <- dcast(data = tab2p, y ~ x, value.var = "Freq")
mat_label <- as.matrix(df_label[,-1])
rownames(mat_label) <- df_label$y
# if (enzyme_type == "kinase") {
#   mat_value[mat_label < 5] <- 0
#   # mat_value <- mat_value[, colSums(mat_label, na.rm = T) >= 5]
#   # mat_label <- mat_label[, colSums(mat_label, na.rm = T) >= 5]
# }

mat_label[is.na(mat_label)] <- ""
col_order4 <- c("AKT1", "MTOR", "MAPK1", "MAP2K1")
# mat_value <- mat_value[, c(col_order4, colnames(mat_value)[!(colnames(mat_value) %in% col_order4)])]
# mat_label <- mat_label[, c(col_order4, colnames(mat_value)[!(colnames(mat_value) %in% col_order4)])]
fn = paste(makeOutDir(resultD = resultD), "druggable_", enzyme_type, '_outlier_ratio_in_patients.pdf',sep = "")

my_heatmap <- pheatmap(mat_value, 
                       color = c("white", col_paletteR(100)),
                       na_col = "white", 
                       display_numbers = mat_label,
                       cluster_rows=F,
                       cluster_cols = F,
                       show_colnames = T)
save_pheatmap_pdf(x = my_heatmap, 
                  filename = fn, 
                  width = 6, height = 2)






# for all druggable genes cis/trans separate ---------
tab_pairs <- esscore_tab_outlier_drug %>%
  select(GENE, cancer, SELF, variable) %>%
  unique() %>%
  select(GENE, cancer, SELF) %>%
  table %>%
  data.frame()
tab_pairs %>%
  head()
tmp <- data.frame(num_partIDs = num_partIDs, cancer = names(num_partIDs))
tab_pairs.wratio <- merge(tab_pairs, tmp, by = c("cancer"), all.x = T)
tab_pairs.wratio %>%
  head()
tab_pairs.wratio$ratio.outlier <- tab_pairs.wratio$Freq/(tab_pairs.wratio$num_partIDs)
tab_pairs.wratio %>%
  head()

for (SELF in c("cis", "trans")) {
  ## compile data frame for plotting
  tab2p <- tab_pairs.wratio
  tab2p <- tab2p[tab2p$SELF == SELF,]
  tab2p$y <- tab2p$cancer
  tab2p$x <- tab2p$GENE
  tab2p$fill <- as.numeric(as.vector(tab2p$ratio.outlier))
  
  ## filtering
  tab2p <- tab2p[!is.na(tab2p$x) & !is.na(tab2p$y),]
  tab2p %>%
    head()
  boxplot(tab2p$fill)
  
  tab2p <- tab2p[!is.na(tab2p$ratio.outlier) & tab2p$ratio.outlier > 0,]
  
  ## ordering
  tab2p$y <- order_cancer(tab2p$y)
  
  df_value <- dcast(data = tab2p, y ~ x, value.var = "fill")
  mat_value <- as.matrix(df_value[,-1])
  rownames(mat_value) <- df_value$y
  head(mat_value)
  mat_value[is.na(mat_value)] <- 0
  
  df_label <- dcast(data = tab2p, y ~ x, value.var = "Freq")
  mat_label <- as.matrix(df_label[,-1])
  rownames(mat_label) <- df_label$y
  # if (enzyme_type == "kinase") {
  #   mat_value[mat_label < 5] <- 0
  #   # mat_value <- mat_value[, colSums(mat_label, na.rm = T) >= 5]
  #   # mat_label <- mat_label[, colSums(mat_label, na.rm = T) >= 5]
  # }
  
  mat_label[is.na(mat_label)] <- ""
  col_order4 <- c("AKT1", "MTOR", "MAPK1", "MAP2K1")
  # mat_value <- mat_value[, c(col_order4, colnames(mat_value)[!(colnames(mat_value) %in% col_order4)])]
  # mat_label <- mat_label[, c(col_order4, colnames(mat_value)[!(colnames(mat_value) %in% col_order4)])]
  fn = paste(makeOutDir(resultD = resultD), "druggable_", enzyme_type, '_', SELF, '_outlier_ratio_in_patients.pdf',sep = "")
  
  my_heatmap <- pheatmap(mat_value, 
                         color = c("white", col_paletteR(100)),
                         na_col = "white", 
                         display_numbers = mat_label,
                         cluster_rows=F,
                         cluster_cols = F,
                         show_colnames = T)
  save_pheatmap_pdf(x = my_heatmap, 
                    filename = fn, 
                    width = ifelse(SELF == "cis", 4, 6), height = 2)
}


# for only druggable oncogenes ---------
tab_pairs <- esscore_tab_outlier_drug %>%
  filter(GENE %in% oncogenes) %>%
  select(GENE, cancer, SELF, variable) %>%
  unique() %>%
  select(GENE, cancer, SELF) %>%
  table %>%
  data.frame()
tab_pairs %>%
  head()
tmp <- data.frame(num_partIDs = num_partIDs, cancer = names(num_partIDs))
tab_pairs.wratio <- merge(tab_pairs, tmp, by = c("cancer"), all.x = T)
tab_pairs.wratio %>%
  head()
tab_pairs.wratio$ratio.outlier <- tab_pairs.wratio$Freq/(tab_pairs.wratio$num_partIDs)
tab_pairs.wratio %>%
  head()


for (SELF in c("cis", "trans")) {
  ## compile data frame for plotting
  tab2p <- tab_pairs.wratio
  tab2p <- tab2p[tab2p$SELF == SELF,]
  tab2p$y <- tab2p$cancer
  tab2p$x <- tab2p$GENE
  tab2p$fill <- as.numeric(as.vector(tab2p$ratio.outlier))
  
  ## filtering
  tab2p <- tab2p[!is.na(tab2p$x) & !is.na(tab2p$y),]
  tab2p %>%
    head()
  boxplot(tab2p$fill)
  
  tab2p <- tab2p[!is.na(tab2p$ratio.outlier) & tab2p$ratio.outlier > 0,]
  
  ## ordering
  tab2p$y <- order_cancer(tab2p$y)
  
  df_value <- dcast(data = tab2p, y ~ x, value.var = "fill")
  mat_value <- as.matrix(df_value[,-1])
  rownames(mat_value) <- df_value$y
  head(mat_value)
  mat_value[is.na(mat_value)] <- 0
  
  df_label <- dcast(data = tab2p, y ~ x, value.var = "Freq")
  mat_label <- as.matrix(df_label[,-1])
  rownames(mat_label) <- df_label$y
  # if (enzyme_type == "kinase") {
  #   mat_value[mat_label < 5] <- 0
  #   # mat_value <- mat_value[, colSums(mat_label, na.rm = T) >= 5]
  #   # mat_label <- mat_label[, colSums(mat_label, na.rm = T) >= 5]
  # }
  
  mat_label[is.na(mat_label)] <- ""
  col_order4 <- c("AKT1", "MTOR", "MAPK1", "MAP2K1")
  mat_value <- mat_value[, c(col_order4, colnames(mat_value)[!(colnames(mat_value) %in% col_order4)])]
  mat_label <- mat_label[, c(col_order4, colnames(mat_value)[!(colnames(mat_value) %in% col_order4)])]
  fn = paste(makeOutDir(resultD = resultD), "druggable_oncogene_", enzyme_type, '_', SELF, '_outlier_ratio_in_patients.pdf',sep = "")
  
  my_heatmap <- pheatmap(mat_value, 
                         color = c("white", col_paletteR(100)),
                         na_col = "white", 
                         display_numbers = mat_label,
                         cluster_rows=F,
                         cluster_cols = F,
                         show_colnames = T)
  save_pheatmap_pdf(x = my_heatmap, 
                    filename = fn, 
                    width = ifelse(SELF == "cis", 3, 2.5), height = 2)
}

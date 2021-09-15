# Yige Wu @ WashU 2019 Dec
## make table for inputting into Protein Paint

# source ------------------------------------------------------------------
source("~/Box/Ding_Lab/Projects_Current/PanCan_Phospho-signaling/phospho-signaling_analysis/phospho_network/phospho_network_shared.R")
packages = c(
  "plyr",
  "rDGIdb"
)

for (pkg_name_tmp in packages) {
  if (!(pkg_name_tmp %in% installed.packages()[,1])) {
    BiocManager::install(pkgs = pkg_name_tmp, update = F)
  }
  library(package = pkg_name_tmp, character.only = T)
}

# set run id --------------------------------------------------------------
version_tmp <- 2
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)


# set  output directory ---------------------------------------------------
dir_out <- makeOutDir()
dir_out_tmp <- paste0(dir_out, run_id, "/")
dir.create(dir_out_tmp)

# input the overall KS result table ---------------------------------------
ks_result_tab <- fread(input = "./Ding_Lab/Projects_Current/PanCan_Phospho-signaling/Resources/GBM/regression_cis_trans_regulation_kinase_phosphatase_GBM_tumor_nonNA20.20191104.v2.1_pathway_annotated.txt", data.table = F)

# input phosphosite-level outlier table -------------------------------
phosphosite_outlier_tab <- fread(input = "./Ding_Lab/Projects_Current/PanCan_Phospho-signaling/Resources/GBM/phosphosite_outlier_table.20191006.v1.csv", data.table = F)
phosphosite_outlier_tab <- data.frame(phosphosite_outlier_tab)
colnames_header <- c("V1")
colnames_case <- colnames(phosphosite_outlier_tab)[!(colnames(phosphosite_outlier_tab) %in% colnames_header)]
header_df <- phosphosite_outlier_tab %>%
  dplyr::rename(gene_phosphosite_dirty = V1) %>%
  select(gene_phosphosite_dirty) %>%
  mutate(gene = str_split_fixed(string = gene_phosphosite_dirty, pattern = "-", n = 2)[,1]) %>%
  mutate(phosphosite_dirty = str_split_fixed(string = gene_phosphosite_dirty, pattern = "-", n = 2)[,2]) %>%
  mutate(phosphosite = gsub(x = phosphosite_dirty, pattern = "[sty]", replacement = "", ignore.case = F)) %>%
  mutate(gene_phosphosite = paste0(gene, "_", phosphosite))

phosphosite_outlier_df <- cbind(header_df, phosphosite_outlier_tab[,colnames_case])

# make combined table -----------------------------------------------------
## filter for significant associated kinase-substrate pairs
ks_result_path_tmp <- ks_result_tab %>%
  filter(ENZ_pathway != "Unknown") %>%
  filter(SUB_pathway != "Unknown") %>%
  filter(FDR_pho_kin < 0.05) %>%
  filter((coef_pho_kin > 0 & enzyme_type == "kinase") | (coef_pho_kin < 0 & enzyme_type == "phosphatase"))
## get all the phosphosites involved in the above ks pairs
phosphosites_tmp <- unique(c(as.vector(ks_result_path_tmp$ENZ_phosphosite), as.vector(ks_result_path_tmp$SUB_phosphosite)))
phosphosites_tmp
## filter the outlier status table by these phosphosite to make it smaller
phosphosite_outlier_df_tmp <- phosphosite_outlier_df[phosphosite_outlier_df$gene_phosphosite %in% phosphosites_tmp,]

ks_result_path_tmp$ES_phosphosite.num_outlier <- sapply(X = 1:nrow(ks_result_path_tmp), 
                                                        FUN = function(i, phosphosite_pair_tab, phospho_outlier_tab, cases) {
                                                          p_enz <- phosphosite_pair_tab[i, "ENZ_phosphosite"]
                                                          p_sub <- phosphosite_pair_tab[i, "SUB_phosphosite"]
                                                          p_enz_outlier <- phospho_outlier_tab[phospho_outlier_tab$gene_phosphosite == p_enz, cases]
                                                          p_sub_outlier <- phospho_outlier_tab[phospho_outlier_tab$gene_phosphosite == p_sub, cases]
                                                          
                                                          p_enz_sub_outlier <- (p_enz_outlier & p_sub_outlier)
                                                          num_enz_sub_outlier <- length(p_enz_sub_outlier[!is.na(p_enz_sub_outlier) & p_enz_sub_outlier == T])
                                                          return(num_enz_sub_outlier)
                                                        }, 
                                                        phosphosite_pair_tab = ks_result_path_tmp, 
                                                        phospho_outlier_tab = phosphosite_outlier_df_tmp, 
                                                        cases = colnames_case)

ks_result_path_tmp$ES_phosphosite.frac_outlier <- ks_result_path_tmp$ES_phosphosite.num_outlier/ks_result_path_tmp$Size

pairs_outlier <- ks_result_path_tmp$pair[ks_result_path_tmp$ES_phosphosite.num_outlier > 0]
genes_outlier <- unique(c(as.vector(ks_result_path_tmp$GENE[ks_result_path_tmp$pair %in% pairs_outlier]), as.vector(ks_result_path_tmp$SUB_GENE[ks_result_path_tmp$pair %in% pairs_outlier])))
pairs_outlier_first_degree <- ks_result_path_tmp$pair[(ks_result_path_tmp$GENE %in% genes_outlier | ks_result_path_tmp$SUB_GENE %in% genes_outlier) & ks_result_path_tmp$ES_phosphosite.num_outlier == 0]


# write the network file for only pairs showing outlier -------------------
## write network file
network_tab <- ks_result_path_tmp %>%
  filter(pair %in% pairs_outlier) %>%
  group_by(GENE, SUB_GENE, SUB_phosphosite) %>%
  dplyr::summarize(coef_pho_kin_ave = mean(coef_pho_kin), ES_phosphosite.frac_outlier.max = max(ES_phosphosite.frac_outlier)) %>%
  mutate(coef_pho_kin_ave.abs = abs(coef_pho_kin_ave))

write.table(x = network_tab, 
            file = paste0(dir_out_tmp, "Cytoscape.Network_File.GBM_KS.", "Phospho_Outlier_Pair_Only.", run_id, ".txt"), 
            sep = "\t", quote = F, row.names = F)

## write node file
genes_tmp <- unique(c(as.vector(network_tab$GENE), as.vector(network_tab$SUB_GENE)))
### annotate druggability in DGIdb
dgi_result <- queryDGIdb(genes_tmp)
dgi_result_by_gene <- byGene(dgi_result)
druggable_genes_tmp <- dgi_result_by_gene$Gene[dgi_result_by_gene$DistinctDrugCount > 0]

### annotate druggability in DGIdb
dgi_result_by_gene <- byGene(dgi_result)
druggable_genes_tmp <- dgi_result_by_gene$Gene[dgi_result_by_gene$DistinctDrugCount > 0]

node_tab <- data.frame(gene_symbol = genes_tmp, 
                       gene_druggable = (genes_tmp %in% druggable_genes_tmp),
                       gene_ks_type = ifelse(genes_tmp %in% ks_result_path_tmp$GENE[ks_result_path_tmp$enzyme_type == "kinase"], "kinase",
                                             ifelse(genes_tmp %in% ks_result_path_tmp$GENE[ks_result_path_tmp$enzyme_type == "phosphatase"], "phosphatase", "substrate")))
write.table(x = node_tab, 
            file = paste0(dir_out_tmp, "Cytoscape.Node_Annotation_File.GBM_KS.", "Phospho_Outlier_Pair_Only.", run_id, ".txt"), 
            sep = "\t", quote = F, row.names = F)

write.table(x = dgi_result_by_gene, 
            file = paste0(dir_out_tmp, "DGIdb_Result_By_Gene.", "Phospho_Outlier_Pair_Only.", run_id, ".txt"), 
            sep = "\t", quote = F, row.names = F)

dgi_result_summary <- resultSummary(dgi_result)
write.table(x = dgi_result_summary, 
            file = paste0(dir_out_tmp, "DGIdb_Result_Summary.", "Phospho_Outlier_Pair_Only.", run_id, ".txt"), 
            sep = "\t", quote = F, row.names = F)

dgi_detailed_result <- detailedResults(dgi_result)
write.table(x = dgi_detailed_result, 
            file = paste0(dir_out_tmp, "DGIdb_Detailed_Results.", "Phospho_Outlier_Pair_Only.", run_id, ".txt"), 
            sep = "\t", quote = F, row.names = F)

# plot the patients showing KS phospho outlier--------------------------------------------------------------------
## filter for top 10 KS pairs showing the higher fraction of patients showing outlier KS
ks_result_tmp <- ks_result_path_tmp %>%
  filter(pair %in% pairs_outlier) %>%
  top_n(n = 10, wt = ES_phosphosite.num_outlier)
## for KS pairs showing outlier, get all involving phosphosite
phosphosites_in_pairs_outlier <- unique(c(as.vector(ks_result_tmp$ENZ_phosphosite), as.vector(ks_result_tmp$SUB_phosphosite)))
## filter the phosphosite outlier status result for these phosphosites
plot_df <- phosphosite_outlier_df[phosphosite_outlier_df$gene_phosphosite %in% phosphosites_in_pairs_outlier,]
## group the substrate phosphosites by kinase phosphosite
plot_df$kinase_phosphosite_group <- mapvalues(x = plot_df$gene_phosphosite, from = ks_result_tmp$SUB_phosphosite, to = ks_result_tmp$ENZ_phosphosite)
## order the data frame by associated kinase phosphosite
plot_df <- plot_df %>%
  arrange(kinase_phosphosite_group)
plot_mat <- as.matrix(plot_df[,colnames_case])
rownames(plot_mat) <- plot_df$gene_phosphosite
## derive the numbers for dividing the heatmap rows
gaps_row_vec <- cumsum(as.vector(table(plot_df$kinase_phosphosite_group)))
## plot heatmap
library(pheatmap)
my_heatmap <- pheatmap(plot_mat, 
                       color = c("white", "black"),
                       gaps_row = gaps_row_vec,
                       scale = "none",
                       na_col = "grey",
                       show_colnames = T,
                       cluster_rows=F, 
                       cluster_cols=T)
save_pheatmap_png(x = my_heatmap, 
                  filename = fn, 
                  width = 1200, height = 350, res = 150)


# filter by substrate pathway and write table ---------------------------------------------
sub_path_tmp <- "RTK RAS"
sub_path_tmp <- "WNT"

ks_result_tab.w.outlier <- NULL
ks_outlier_sample_list <- list()
for (sub_path_tmp in unique(ks_result_tab$SUB_pathway)) {
  dir_out_tmp <- paste0(dir_out, sub_path_tmp, "/")
  dir.create(dir_out_tmp)
  
  ks_result_path_tmp <- ks_result_tab %>%
    filter(SUB_pathway == sub_path_tmp) %>%
    filter(ENZ_pathway != "Unknown") %>%
    filter(FDR_pho_kin < 0.05) %>%
    filter((coef_pho_kin > 0 & enzyme_type == "kinase") | (coef_pho_kin < 0 & enzyme_type == "phosphatase"))
  
  phosphosites_tmp <- unique(c(as.vector(ks_result_path_tmp$ENZ_phosphosite), as.vector(ks_result_path_tmp$SUB_phosphosite)))
  phosphosites_tmp
  
  phosphosite_outlier_df_tmp <- phosphosite_outlier_df[phosphosite_outlier_df$gene_phosphosite %in% phosphosites_tmp,]
  
  ks_result_path_tmp$ES_phosphosite.num_outlier <- sapply(X = 1:nrow(ks_result_path_tmp), 
                                                           FUN = function(i, phosphosite_pair_tab, phospho_outlier_tab, cases) {
                                                             p_enz <- phosphosite_pair_tab[i, "ENZ_phosphosite"]
                                                             p_sub <- phosphosite_pair_tab[i, "SUB_phosphosite"]
                                                             p_enz_outlier <- phospho_outlier_tab[phospho_outlier_tab$gene_phosphosite == p_enz, cases]
                                                             p_sub_outlier <- phospho_outlier_tab[phospho_outlier_tab$gene_phosphosite == p_sub, cases]
                                                             
                                                             p_enz_sub_outlier <- (p_enz_outlier & p_sub_outlier)
                                                             num_enz_sub_outlier <- length(p_enz_sub_outlier[!is.na(p_enz_sub_outlier) & p_enz_sub_outlier == T])
                                                             return(num_enz_sub_outlier)
                                                             }, 
                                                           phosphosite_pair_tab = ks_result_path_tmp, 
                                                           phospho_outlier_tab = phosphosite_outlier_df_tmp, 
                                                           cases = colnames_case)
  
  ks_result_path_tmp$ES_phosphosite.frac_outlier <- ks_result_path_tmp$ES_phosphosite.num_outlier/ks_result_path_tmp$Size
  
  genes_tmp <- unique(c(as.vector(ks_result_path_tmp$GENE), as.vector(ks_result_path_tmp$SUB_GENE)))
  dgi_result <- queryDGIdb(genes_tmp)
  
  dgi_result_by_gene <- byGene(dgi_result)
  druggable_genes_tmp <- dgi_result_by_gene$Gene[dgi_result_by_gene$DistinctDrugCount > 0]
  
  ## write network file
  network_tab <- ks_result_path_tmp %>%
    group_by(GENE, SUB_GENE, SUB_phosphosite) %>%
    dplyr::summarize(coef_pho_kin_ave = mean(coef_pho_kin), ES_phosphosite.frac_outlier.max = max(ES_phosphosite.frac_outlier)) %>%
    mutate(coef_pho_kin_ave.abs = abs(coef_pho_kin_ave))
  
  write.table(x = network_tab, 
              file = paste0(dir_out_tmp, "Cytoscape.Network_File.GBM_KS.", "SUB_pathway_", sub_path_tmp, ".", run_id, ".txt"), 
              sep = "\t", quote = F, row.names = F)
  
  ## write node file
  genes_tmp <- unique(c(as.vector(ks_result_path_tmp$GENE), as.vector(ks_result_path_tmp$SUB_GENE)))
  node_tab <- data.frame(gene_symbol = genes_tmp, 
                         gene_druggable = (genes_tmp %in% druggable_genes_tmp),
                         gene_ks_type = ifelse(genes_tmp %in% ks_result_path_tmp$GENE[ks_result_path_tmp$enzyme_type == "kinase"], "kinase",
                                               ifelse(genes_tmp %in% ks_result_path_tmp$GENE[ks_result_path_tmp$enzyme_type == "phosphatase"], "phosphatase", "substrate")))
  write.table(x = node_tab, 
        file = paste0(dir_out_tmp, "Cytoscape.Node_Annotation_File.GBM_KS.", "SUB_pathway_", sub_path_tmp, ".", run_id, ".txt"), 
        sep = "\t", quote = F, row.names = F)
  
  ## bind to new KS result table
  ks_result_path_tmp$ENZ.druggable <- (ks_result_path_tmp$GENE %in% druggable_genes_tmp)
  ks_result_path_tmp$SUB.druggable <- (ks_result_path_tmp$SUB_GENE %in% druggable_genes_tmp)
  ks_result_tab.w.outlier <- rbind(ks_result_path_tmp, ks_result_tab.w.outlier)
  
  ## record which samples have outlier phosphorlation in KS pair
  ks_outlier_samples_list_tmp <- lapply(X = 1:nrow(ks_result_path_tmp), 
                                        FUN = function(i, phosphosite_pair_tab, phospho_outlier_tab, cases) {
                                          p_enz <- phosphosite_pair_tab[i, "ENZ_phosphosite"]
                                          p_sub <- phosphosite_pair_tab[i, "SUB_phosphosite"]
                                          p_enz_outlier <- phospho_outlier_tab[phospho_outlier_tab$gene_phosphosite == p_enz, cases]
                                          p_sub_outlier <- phospho_outlier_tab[phospho_outlier_tab$gene_phosphosite == p_sub, cases]
                                          
                                          p_enz_sub_outlier <- (p_enz_outlier & p_sub_outlier)
                                          cases_enz_sub_outlier <- colnames_case[!is.na(p_enz_sub_outlier) & p_enz_sub_outlier == T]
                                          return(cases_enz_sub_outlier)
                                        }, 
                                        phosphosite_pair_tab = ks_result_path_tmp, 
                                        phospho_outlier_tab = phosphosite_outlier_df_tmp, 
                                        cases = colnames_case)
  names(ks_outlier_samples_list_tmp) <- ks_result_path_tmp$pair
  ks_outlier_sample_list <- c(ks_outlier_sample_list, ks_outlier_samples_list_tmp)
}

# calculate how many pairs are druggable ----------------------------------
ks_result_tab.w.outlier %>%
  filter(ES_phosphosite.num_outlier > 0) %>%
  filter(ENZ.druggable | SUB.druggable) %>%
  select(GENE, SUB_phosphosite) %>%
  unique() %>%
  nrow()

ks_result_tab.w.outlier %>%
  filter(ES_phosphosite.num_outlier > 0) %>%
  filter(ENZ.druggable | SUB.druggable) %>%
  select(GENE, SUB_GENE) %>%
  unique() %>%
  nrow()

ks_result_tab.w.outlier %>%
  filter(ES_phosphosite.num_outlier > 0) %>%
  filter(ENZ.druggable | SUB.druggable) %>%
  select(GENE) %>%
  unique() %>%
  nrow()

ks_result_tab.w.outlier %>%
  filter(ES_phosphosite.num_outlier > 0) %>%
  filter(ENZ.druggable | SUB.druggable) %>%
  select(SUB_GENE) %>%
  unique() %>%
  nrow()

druggable_pairs <- ks_result_tab.w.outlier$pair[(ks_result_tab.w.outlier$ENZ.druggable | ks_result_tab.w.outlier$SUB.druggable) & ks_result_tab.w.outlier$ES_phosphosite.num_outlier > 0]
length(unique(unlist(ks_outlier_sample_list[druggable_pairs])))
length(colnames_case)

# calculate the number of samples with outlier event in each druggable gene --------
druggable_enzymes <- unique(ks_result_tab.w.outlier$GENE[ks_result_tab.w.outlier$ENZ.druggable & ks_result_tab.w.outlier$ES_phosphosite.num_outlier > 0])
druggable_substrates <- unique(ks_result_tab.w.outlier$SUB_GENE[ks_result_tab.w.outlier$SUB.druggable & ks_result_tab.w.outlier$ES_phosphosite.num_outlier > 0])
druggable_substrates

druggable_genes <- c(setdiff(druggable_enzymes, druggable_substrates),
                     intersect(druggable_substrates, druggable_enzymes),
                     setdiff(druggable_substrates, druggable_enzymes))
numbers_case_outlier <- sapply(X = druggable_genes, FUN = function(g, case_list, pair_table) {
  pairs_tmp <- pair_table$pair[pair_table$GENE == g | pair_table$SUB_GENE == g]
  cases_tmp <- unique(unlist(case_list[pairs_tmp]))
  return(length(cases_tmp))
}, case_list = ks_outlier_sample_list, pair_table = ks_result_tab.w.outlier)

tab2p <- data.frame(gene_symbol = druggable_genes, number_case_outlier = numbers_case_outlier)
tab2p <- tab2p %>%
  top_n(wt = number_case_outlier, n = 20) %>%
  arrange(number_case_outlier)

tab2p$gene_symbol <- factor(tab2p$gene_symbol, levels = as.vector(tab2p$gene_symbol))

p <- ggplot()
p <- p + geom_bar(data = tab2p, mapping = aes(x = gene_symbol, y = number_case_outlier), stat = "identity")
p <- p + coord_flip()
p <- p + theme_bw()
pdf(file = paste0(dir_out, "Barplot.Num_Case_Outlier_KS_Phospho.ByGene.pdf"), width = 4, height = 6)
print(p)
dev.off()







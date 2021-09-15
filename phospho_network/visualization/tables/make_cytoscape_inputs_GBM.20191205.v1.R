# Yige Wu @ WashU 2019 Dec
## make table for inputting into Protein Paint

# source ------------------------------------------------------------------
source("~/Box/Ding_Lab/Projects_Current/PanCan_Phospho-signaling/phospho-signaling_analysis/phospho_network/phospho_network_shared.R")
library(plyr)

# set run id --------------------------------------------------------------
version_tmp <- 2
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)


# set  output directory ---------------------------------------------------
dir_out <- makeOutDir()


# input the overall KS result table ---------------------------------------
ks_result_tab <- fread(input = "./Ding_Lab/Projects_Current/PanCan_Phospho-signaling/Resources/GBM/regression_cis_trans_regulation_kinase_phosphatase_GBM_tumor_nonNA20.20191104.v2.1_pathway_annotated.txt", data.table = F)


# input phospho gene-level outlier fraction -------------------------------
pho_gene_fraction_tab <- fread(input = "./Ding_Lab/Projects_Current/PanCan_Phospho-signaling/Resources/GBM/fraction_table.csv", data.table = F)
pho_gene_fraction_cohort_tab <- 
pho_gene_fraction_cohort_tab <- data.frame(gene_symbol = pho_gene_fraction_tab$V1, 
                                           pho_outlier_frac_median = sapply(1:nrow(pho_gene_fraction_tab), FUN = function(i, mat) median(mat[i,], na.rm = T), mat = as.matrix(pho_gene_fraction_tab[,-1])),
                                           pho_outlier_frac_mean = sapply(1:nrow(pho_gene_fraction_tab), FUN = function(i, mat) mean(mat[i,], na.rm = T), mat = as.matrix(pho_gene_fraction_tab[,-1])))

quantile_level <- 0.5
iqr_factor <- 1.5

# input phospho data ------------------------------------------------------
pho_file <- fread(input = "./Ding_Lab/Projects_Current/CPTAC3-GBM/Data_Freeze/cptac3_gbm_data_freeze_v3.0.20191121/phosphoproteome/mssm/phosphoproteome_mssm_per_gene_clean.v3.0.20191121.tsv.gz", data.table = F)
pho_file %>% colnames()
colnames_header <- c("gene", "refseq_id", "peptide", "site")
colnames_case <- colnames(pho_file)[!(colnames(pho_file) %in% colnames_header)]
header_df <- pho_file %>%
  select(colnames_header) %>%
  mutate(phosphosite_dirty = str_split_fixed(string = site, pattern = "-", n = 2)[,2]) %>%
  mutate(phosphosite = gsub(x = phosphosite_dirty, pattern = "[sty]", replacement = "", ignore.case = F)) %>%
  mutate(gene_phosphosite = paste0(gene, "_", phosphosite))
rownames(pho_file) <- header_df$gene_phosphosite


# filter by substrate pathway and write table ---------------------------------------------
ks_result_tab %>% colnames()
ks_result_tab %>%
  select(SUB_pathway) %>%
  unique()

sub_path_tmp <- "RTK RAS"
sub_path_tmp <- "WNT"

for (sub_path_tmp in unique(ks_result_tab$SUB_pathway)) {
  dir_out_tmp <- paste0(dir_out, sub_path_tmp, "/")
  dir.create(dir_out_tmp)
  
  ks_result_path_tmp <- ks_result_tab %>%
    filter(SUB_pathway == sub_path_tmp) %>%
    filter(ENZ_pathway != "Unknown") %>%
    # filter(ENZ_phosphosite.is_regulatory == T) %>%
    # filter(SUB_phosphosite.is_regulatory == T) %>%
    filter(FDR_pho_kin < 0.05) %>%
    filter((coef_pho_kin > 0 & enzyme_type == "kinase") | (coef_pho_kin < 0 & enzyme_type == "phosphatase"))
  
  phosphosites_tmp <- unique(c(as.vector(ks_result_path_tmp$ENZ_phosphosite), as.vector(ks_result_path_tmp$SUB_phosphosite)))
  phosphosites_tmp
  
  pho_mat_tmp <- pho_file[header_df$gene_phosphosite %in% phosphosites_tmp, colnames_case]
  header_tmp <- header_df[header_df$gene_phosphosite %in% phosphosites_tmp,]
  
  pho_quantile <- sapply(1:nrow(pho_mat_tmp), FUN = function(i, mat) quantile(x = mat[i,], probs = quantile_level, na.rm = T), mat = as.matrix(pho_mat_tmp))
  pho_iqr <- sapply(1:nrow(pho_mat_tmp), FUN = function(i, mat) IQR(x = mat[i,], na.rm = T), mat = as.matrix(pho_mat_tmp))
  pho_outlier_thres <- pho_quantile + iqr_factor*pho_iqr
  
  pho_outlier_mat_tmp <- (pho_mat_tmp >= pho_outlier_thres)

  pho_outlier_df_tmp <- cbind(header_tmp, pho_outlier_mat_tmp)
  pho_outlier_df_tmp <- pho_outlier_df_tmp[!duplicated(header_tmp$gene_phosphosite),]
  
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
                                                           phospho_outlier_tab = pho_outlier_df_tmp, 
                                                           cases = colnames_case)
  
  ks_result_path_tmp$ES_phosphosite.frac_outlier <- ks_result_path_tmp$ES_phosphosite.num_outlier/ks_result_path_tmp$Size
  
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
  node_tab <- data.frame(gene_symbol = genes_tmp, gene_ks_type = ifelse(genes_tmp %in% ks_result_path_tmp$GENE[ks_result_path_tmp$enzyme_type == "kinase"], "kinase",
                                                                        ifelse(genes_tmp %in% ks_result_path_tmp$GENE[ks_result_path_tmp$enzyme_type == "phosphatase"], "phosphatase", "substrate")))
  write.table(x = node_tab, 
        file = paste0(dir_out_tmp, "Cytoscape.Node_Annotation_File.GBM_KS.", "SUB_pathway_", sub_path_tmp, ".", run_id, ".txt"), 
        sep = "\t", quote = F, row.names = F)
}


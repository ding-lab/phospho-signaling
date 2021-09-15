# Yige Wu @ WashU 2018 Aug
## For formating GBM manuscript data freeze v2.0

# set working directory ---------------------------------------------------
baseD = "~/Box/"
setwd(baseD)
source("./Ding_Lab/Projects_Current/PanCan_Phospho-signaling/phospho-signaling_analysis/phospho_network/phospho_network_shared.R")

# set variables -----------------------------------------------------------
cancer <- "GBM"
pipeline_type <- "PGDAC"
pro_normalization_method <- "d4"
pho_normalization_method <- "d6"

# input bobo’s PTMcosmos to add peptide ID --------------------------------
uniprot_ensembl_map <- fread("./cptac2p/resources/PTMcosmos/Data dump/2019-01-11/all_uniprot_ensembl_id_mappings.tsv", data.table = F)
ptm_sites_map <- fread("./cptac2p/resources/PTMcosmos/Data dump/2019-01-11/all_ptm_sites.all_txs.tsv", data.table = F)


# rewrite proteomics/phosphoproteomics data split by tumor and normal ----------------------------------
for (sample_type in c("tumor")) {
  pro_data <- fread(input = "Ding_Lab/Projects_Current/CPTAC3-GBM/Data_Freeze/cptac3_gbm_data_freeze_v2.1.20190927/proteome/mssm/proteome_mssm_per_gene_clean.v2.1.20190927.tsv.gz", data.table = F)
  head(pro_data)
  col_names_origin <- colnames(pro_data)[!(colnames(pro_data) %in% c("gene", "refseq_id"))]
  col_names_keep <- col_names_origin
  data2w <- data.frame(Gene = pro_data$gene)
  data2w <- cbind(data2w, pro_data[, col_names_keep])
  expresson_type <- "PRO"
  write.table(x = data2w, file = paste0(makeOutDir(), cancer, "_", expresson_type, "_", sample_type, "_", pipeline_type, "_", pro_normalization_method, "_", "partID", ".txt"), quote = F, row.names = F, sep = "\t")
  
  pho_data <- fread(input = "Ding_Lab/Projects_Current/CPTAC3-GBM/Data_Freeze/cptac3_gbm_data_freeze_v2.1.20190927/phosphoproteome/mssm/phosphoproteome_mssm_per_gene_clean.v2.1.20190927.tsv.gz", data.table = F)
  head(pho_data)
  col_names_origin <- colnames(pho_data)[!(colnames(pho_data) %in% c("gene", "refseq_id", "peptide", "site"))]
  col_names_keep <- col_names_origin
  data2w <- data.frame(Gene = pho_data$gene, Phosphosite = str_split_fixed(string = pho_data$site, pattern = "-", n = 2)[,2])
  data2w %>% head()
  data2w <- cbind(data2w, pho_data[, col_names_keep])
  expresson_type <- "PHO"
  write.table(x = data2w, file = paste0(makeOutDir(), cancer, "_", expresson_type, "_", sample_type, "_", pipeline_type, "_", pho_normalization_method, "_", "partID", ".txt"), quote = F, row.names = F, sep = "\t")
}


# format copy number status -----------------------------------------------
cna_per_gene_tab <- fread(input = "Ding_Lab/Projects_Current/CPTAC3-GBM/Data_Freeze/cptac3_gbm_data_freeze_v2.0.20190905/somatic_cnv/wgs_somatic_cnv_per_gene.v2.0.20190905.tsv.gz", data.table = F)
col_names_origin <- colnames(cna_per_gene_tab)[!(colnames(cna_per_gene_tab) %in% c("symbol", "gene_id", "gene_id_version", "original_symbol"))]
col_names_keep <- col_names_origin
col_names_keep
cna_per_gene_tab <- cna_per_gene_tab[, c("symbol", col_names_keep)]
colnames(cna_per_gene_tab) <- c("gene", col_names_keep)
head(cna_per_gene_tab)
write.table(x = cna_per_gene_tab, 
            file = paste0(makeOutDir(), "somatic_CNA", ".", cancer, ".partID.txt"), quote = F, row.names = F, sep = "\t")


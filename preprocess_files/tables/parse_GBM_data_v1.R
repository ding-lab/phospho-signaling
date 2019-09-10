# Yige Wu @ WashU 2018 Aug
## For formating GBM WG data v1.0.20190802


# source ------------------------------------------------------------------
code_top_dir <- "~/Box/Ding_Lab/Projects_Current/PanCan_Phospho-signaling/phospho-signaling_analysis/"
path2phospho_network_shared <- paste0(code_top_dir, "phospho_network/phospho_network_shared.R")
source(path2phospho_network_shared)

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
  pro_data <- fread(input = "Ding_Lab/Projects_Current/CPTAC/PGDAC/GBM/v1.0.20190802/proteome_pnnl_per_gene_d4.v1.0.20190802.tsv", data.table = F)
  head(pro_data)
  data2w <- pro_data
  colnames(data2w)[1] <- "Gene"
  expresson_type <- "PRO"
  write.table(x = data2w, file = paste0(makeOutDir(), cancer, "_", expresson_type, "_", sample_type, "_", pipeline_type, "_", pro_normalization_method, "_", "partID", ".txt"), quote = F, row.names = F, sep = "\t")
  
  pho_data <- fread(input = "Ding_Lab/Projects_Current/CPTAC/PGDAC/GBM/v1.0.20190802/phosphoproteome_pnnl_d6.v1.0.20190802.tsv", data.table = F)
  head(pho_data)
  col_names_origin <- colnames(pho_data)[!(colnames(pho_data) %in% c("gene", "peptide", "site"))]
  col_names_keep <- col_names_origin
  data2w <- data.frame(Gene = pho_data$gene, Phosphosite = str_split_fixed(string = pho_data$site, pattern = "-", n = 2)[,2])
  data2w %>% head()
  data2w <- cbind(data2w, pho_data[, col_names_keep])
  expresson_type <- "PHO"
  write.table(x = data2w, file = paste0(makeOutDir(), cancer, "_", expresson_type, "_", sample_type, "_", pipeline_type, "_", pho_normalization_method, "_", "partID", ".txt"), quote = F, row.names = F, sep = "\t")
}

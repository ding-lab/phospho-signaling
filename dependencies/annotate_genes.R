# Yige Wu @ WashU 2018 Jan
## annotate columns in table

# set workding directory ------------------------------------------------
baseD = "/Users/yigewu/Box\ Sync/"
setwd(baseD)

# annotate druggable genes ------------------------------------------------
drug_genes <- fread(input = "./Ding_Lab/Projects_Current/CPTAC/pan3can_shared_data/reference_files/gene_drug_list/Premed_raw_databases/drugBank/drug_list.txt", data.table = F, col.names = "gene")

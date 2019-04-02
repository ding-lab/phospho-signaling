# Yige Wu @ WashU 2019 Feb
## 

# source ------------------------------------------------------------------
setwd(dir = "~/Box Sync/")
source("./cptac2p_analysis/phospho_network/phospho_network_shared.R")

# compile kinase/phosphatase-substrate table ---------------------------
## input enzyme-substrate table
omnipath_tab <- load_omnipath()
### take out pairs that are purely computationally predicted
omnipath_tab <- omnipath_tab[!(omnipath_tab$Source %in% c("NetKIN", "PhosphoNetworks", "MIMP")) ,]
### take unique combinations of enzymes and substrates
enzyme_sub_tab <- unique(omnipath_tab[, c("GENE", "SUB_GENE")])
### input psp
psp_tab <- load_psp()
enzyme_sub_tab <- unique(rbind(enzyme_sub_tab, unique(psp_tab[, c("GENE", "SUB_GENE")])))
enzyme_sub_tab_rev <- data.frame(GENE = enzyme_sub_tab$SUB_GENE, SUB_GENE = enzyme_sub_tab$GENE)
enzyme_sub_pair_tab <- unique(rbind(enzyme_sub_tab, enzyme_sub_tab_rev))


# compile complex table ---------------------------------------------------
## input complex table
complex_pair_tab <- fread(input = "./cptac2p/analysis_results/phospho_network/compile_enzyme_substrate/tables/parse_corum_signor_reactome/sup_complex_pair_uniq.txt", data.table = F)
### reformat complex table
complex_pair_tab <- data.frame(GENE = complex_pair_tab$geneA, SUB_GENE = complex_pair_tab$geneB)

## merge the final table to test
pair_tab_trans <- unique(rbind(enzyme_sub_tab, enzyme_sub_tab_rev, complex_pair_tab))
pair_tab_trans %>%
  head()
all_genes_involved <- unique(c(as.vector(pair_tab_trans$GENE), as.vector(pair_tab_trans$SUB_GENE)))
pair_tab_cis <- data.frame(GENE = all_genes_involved, SUB_GENE = all_genes_involved)
pair_tab <- unique(rbind(pair_tab_trans, pair_tab_cis))

## clean up
bad_strings <- c("DNA damage", "")
pair_tab <- pair_tab[!(pair_tab$GENE %in% bad_strings) & !(pair_tab$SUB_GENE %in% bad_strings),]


# input TF relations ------------------------------------------------------
tf_tab <- fread(input = "./Ding_Lab/Projects_Current/TP53_shared_data/resources/PPI/TF_interactions.txt", data.table = F)
tf_tab %>% head()
tf_pair_tab <- data.frame(GENE = tf_tab$source_genesymbol, SUB_GENE = tf_tab$target_genesymbol)

## double check the important TP53 downstreams are there
list_downstream <- fread(input = "./Ding_Lab/Projects_Current/TP53_shared_data/resources/PPI/list_downstream.txt", data.table = F, col.names = "Gene", header = F)
list_downstream %>% head()
list_downstream <- as.vector(list_downstream$Gene)
list_downstream[!(list_downstream %in% tf_pair_tab$SUB_GENE[tf_pair_tab$GENE == "TP53"])]
## ~10 are off, add back
tf_pair_tab <- unique(rbind(data.frame(GENE = "TP53", SUB_GENE = list_downstream),
                            tf_pair_tab))

tf_pair_tab_rev <- data.frame(GENE = tf_pair_tab$SUB_GENE, SUB_GENE = tf_pair_tab$GENE)

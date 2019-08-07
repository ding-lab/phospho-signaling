# Yige Wu @ WashU 2019 Apr
## test mutation impact on protein/phosphorylation within kinase-substrate pairs or protein complex pairs

# source ------------------------------------------------------------------
baseD = "/Users/yigewu/Box\ Sync/"
wd <- getwd()
if (wd != baseD) {
  setwd(baseD)
}
source("./cptac2p_analysis/phospho_network/phospho_network_shared.R")
source("./cptac2p_analysis/p53/TP53_shared.R")


# input enzyme-substrate table --------------------------------------------
ptms_site_pairs_sup <- load_omnipath()
ptms_site_pairs_sup <- data.frame(ptms_site_pairs_sup)
### take out pairs that are purely computationally predicted
ptms_site_pairs_sup <- ptms_site_pairs_sup[!(ptms_site_pairs_sup$Source %in% c("NetKIN", "PhosphoNetworks", "MIMP")) ,]
### take unique combinations of enzymes and substrates
enzyme_sub_tab <- unique(ptms_site_pairs_sup[, c("GENE", "SUB_GENE")])
enzyme_sub_tab_rev <- data.frame(GENE = enzyme_sub_tab$SUB_GENE, SUB_GENE = enzyme_sub_tab$GENE)


# input complex -----------------------------------------------------------
complex_pair_tab <- fread(input = "./cptac2p/analysis_results/phospho_network/compile_enzyme_substrate/tables/parse_corum_signor_reactome/sup_complex_pair_uniq.txt", data.table = F)
### reformat complex table
complex_pair_tab <- data.frame(GENE = complex_pair_tab$geneA, SUB_GENE = complex_pair_tab$geneB)


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

## input corum
corum_tab <- fread(input = "./cptac2p/analysis_results/phospho_network/compile_enzyme_substrate/tables/compile_corum_complex/corum_complex_pair_tab.txt", data.table = F)
## https://wustl.app.box.com/folder/64826080258
corum_tab$pair_pro <- paste0(corum_tab$geneA, ":",  corum_tab$geneB)


# Merge all protein pairs -------------------------------------------------
pair_tab_trans <- unique(rbind(enzyme_sub_tab, enzyme_sub_tab_rev, 
                               complex_pair_tab,
                               tf_pair_tab, tf_pair_tab_rev))
pair_tab_trans %>%
  head()
all_genes_involved <- unique(c(as.vector(pair_tab_trans$GENE), as.vector(pair_tab_trans$SUB_GENE)))
pair_tab_cis <- data.frame(GENE = all_genes_involved, SUB_GENE = all_genes_involved)
pair_tab <- unique(rbind(pair_tab_trans, pair_tab_cis))

## clean up
bad_strings <- c("DNA damage", "")
pair_tab <- pair_tab[!(pair_tab$GENE %in% bad_strings) & !(pair_tab$SUB_GENE %in% bad_strings),]

## filtering
pair_tab <- pair_tab[pair_tab$GENE == "TP53",]

write.table(x = pair_tab, file = paste0(makeOutDir(resultD = resultD), "TP53_interactor_pair_tab.txt"), quote = F, sep = "\t", row.names = F)
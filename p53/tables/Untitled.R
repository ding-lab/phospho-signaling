# Yige Wu @ WashU 2019 Jan
## 

# inputs pair table------------------------------------------------------------------
## input enzyme-substrate table
ptms_site_pairs_sup <- read_csv(paste0(ppnD, "compile_enzyme_substrate/tables/compile_omnipath_networkin_depod_signor_manual/Omnipath_NetworKin_DEPOD_SignorNotSiteMapped_manualAdded.csv"))
ptms_site_pairs_sup <- data.frame(ptms_site_pairs_sup)
### take out pairs that are purely computationally predicted
ptms_site_pairs_sup <- ptms_site_pairs_sup[!(ptms_site_pairs_sup$Source %in% c("NetKIN", "PhosphoNetworks", "MIMP")) ,]
### take unique combinations of enzymes and substrates
enzyme_sub_tab <- unique(ptms_site_pairs_sup[, c("GENE", "SUB_GENE")])
enzyme_sub_tab_rev <- data.frame(GENE = enzyme_sub_tab$SUB_GENE, SUB_GENE = enzyme_sub_tab$GENE)

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

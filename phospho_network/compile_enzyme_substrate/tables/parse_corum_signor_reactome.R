# Yige Wu @ WashU 2018 Nov
# parse complex information from corum, signor and reactome paired gene table



# source ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')

# inputs ------------------------------------------------------------------
## input signor pair table
signor_complex_pair_tab <- fread(input = "./cptac2p/analysis_results/phospho_network/compile_enzyme_substrate/tables/compile_signor_complex/signor_complex_pair_tab.txt", data.table = F)
signor_complex_pair_tab$pair_pro <- paste0(signor_complex_pair_tab$geneA, ":", signor_complex_pair_tab$geneB)
which(duplicated(signor_complex_pair_tab$pair_pro))
signor_complex_pair_uniq <- signor_complex_pair_tab[!duplicated(signor_complex_pair_tab$pair_pro), c("geneA", "geneB", "pair_pro")]

## input corum pair table
corum_complex_pair_tab <- fread(input = "./cptac2p/analysis_results/phospho_network/compile_enzyme_substrate/tables/compile_corum_complex/corum_complex_pair_tab.txt", data.table = F)
corum_complex_pair_tab$pair_pro <- paste0(corum_complex_pair_tab$geneA, ":", corum_complex_pair_tab$geneB)
corum_complex_pair_uniq <- corum_complex_pair_tab[!duplicated(corum_complex_pair_tab$pair_pro), c("geneA", "geneB", "pair_pro")]

# Business  ---------------------------------------------------------------
## rbind
sup_complex_pair_uniq <- unique(rbind(corum_complex_pair_uniq, signor_complex_pair_uniq))
sup_complex_pair2add <- data.frame(geneA = c("TP53", "TP53", "TP53", "TP53", "TP53", "TP53", "TP53", "TP53", "TP53", "TP53", "TP53"), geneB = c("TP53BP1", "CREBBP", "KAT5", "SIRT1", "BAX", "PUMA", "NOXA", "CDKN1A", "SLC7A11", "GLS2", "TIGAR"))
sup_complex_pair2addrev <- data.frame(geneA = sup_complex_pair2add$geneB, geneB = sup_complex_pair2add$geneA)
sup_complex_pair2add <- rbind(sup_complex_pair2add, sup_complex_pair2addrev)
sup_complex_pair2add$pair_pro <- paste0(sup_complex_pair2add$geneA, ":", sup_complex_pair2add$geneB)
sup_complex_pair_uniq <- unique(rbind(sup_complex_pair_uniq, sup_complex_pair2add))
nrow(sup_complex_pair_uniq)
write.table(x = sup_complex_pair_uniq, file = paste0(makeOutDir(resultD = resultD), "sup_complex_pair_uniq.txt"), quote = F, sep = "\t", row.names = F)

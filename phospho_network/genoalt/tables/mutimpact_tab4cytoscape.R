# Yige Wu @ WashU 2018 July
## output table for cytoscape visualization

# source ------------------------------------------------------------------
source("./cptac2p_analysis/phospho_network/phospho_network_plotting.R")
source("./cptac2p_analysis/phospho_network/phospho_network_shared.R")

# set variables -----------------------------------------------------------
sig_thres <- 0.1

# inputs ------------------------------------------------------------------
mutimpact <- fread(input = paste0(ppnD, "genoalt/tables/merge_mutation_impact/mutation_impact.txt"), data.table = F)
mutimpact$log10_pvalue <- -log10(mutimpact$p_value)
mutimpact$SUB_MOD_RSD <- toupper(mutimpact$Phosphosite)
mutimpact$pair <- paste0(mutimpact$Mutated_Gene, ":", mutimpact$Substrate_Gene, ":", mutimpact$SUB_MOD_RSD)

driver_genes <- read.delim(file = "./TCGA_data/reference_files/Consensus.Genelist.full.txt")

ptms_site_pairs_sup <- read_csv(paste0(ppnD, "compile_enzyme_substrate/tables/compile_omnipath_networkin/omnipath_networkin_enzyme_substrate_site_level_union.csv"))

ptms_site_pairs_sup <- read_csv(paste0(ppnD, "compile_enzyme_substrate/tables/compile_omnipath_networkin/omnipath_networkin_enzyme_substrate_site_level_union_source_cleaned.csv"))
mutimpact <- merge(mutimpact, ptms_site_pairs_sup[, c("pair", "Source", "networkin_score")], all.x = T)
mutimpact_sig <- mutimpact[mutimpact$p_value < sig_thres,]
mutimpact_sig_pho <- mutimpact_sig[mutimpact_sig$Phosphosite != "Protein",]
mutimpact_sig_pho_source <- mutimpact_sig_pho[!is.na(mutimpact_sig_pho$Source),]
mutimpact_sig_pho_source_predict <- mutimpact_sig_pho_source[mutimpact_sig_pho_source$Source %in% c("NetKIN", "MIMP"),]
mutimpact_sig_pho_source_experiment <- mutimpact_sig_pho_source[!(mutimpact_sig_pho_source$Source %in% c("NetKIN", "MIMP")),]

# write network file for BRCA ---------------------------------------------
for (cancer in c("BRCA", "CO")) {
  network <- mutimpact_sig_pho_source_experiment[mutimpact_sig_pho_source_experiment$Cancer == cancer,]
  write.table(x = network, file = paste0(makeOutDir(resultD = resultD), cancer,  "mutimpact_sig_pho_source_experiment.txt"), 
              quote = F, row.names = F, col.names = T, sep = "\t")
  ## write node file 
  node_tab <- data.frame(shared_name = unique(c(as.vector(network$Mutated_Gene), as.vector(network$Substrate_Gene))))
  node_tab$is_enzyme <- (node_tab$shared_name %in% mutimpact$Mutated_Gene)
  node_tab$enzyme_driver_type <- ""
  node_tab$enzyme_driver_type[node_tab$shared_name %in% driver_genes$Gene[grepl(pattern = "oncogene", x = driver_genes$Tumor.suppressor.or.oncogene.prediction..by.20.20..)]] <- "oncogene"
  node_tab$enzyme_driver_type[node_tab$shared_name %in% driver_genes$Gene[grepl(pattern = "tsg", x = driver_genes$Tumor.suppressor.or.oncogene.prediction..by.20.20..)]] <- "tsg"
  node_tab$enzyme_is_druggable <- (node_tab$shared_name %in% mutimpact$Mutated_Gene[mutimpact$enzyme_is_druggable])
  node_tab$cancer <- cancer
  write.table(x = node_tab, file = paste0(makeOutDir(resultD = resultD), cancer,  "mutimpact_sig_pho_source_experiment_node_tab.txt"), 
              quote = F, row.names = F, col.names = T, sep = "\t")
}
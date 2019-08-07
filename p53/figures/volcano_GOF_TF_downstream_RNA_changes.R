# Yige Wu @ WashU 2019 Apr
## plot the expression changes in the transcriptional targets of TFs interacting with mutant p53 that enable GOF


# source ------------------------------------------------------------------
baseD = "/Users/yigewu/Box\ Sync/"
setwd(baseD)
source("./cptac2p_analysis/phospho_network/phospho_network_shared.R")
source("./cptac2p_analysis/phospho_network/phospho_network_plotting.R")
source("./cptac2p_analysis/p53/TP53_shared.R")


# input result table for comparing GOF missense vs other missense ---------
mut_impact_tab <- fread(input = "./cptac2p/analysis_results/p53/tables/test_TP53_missense_impact_TF_downstream/TP53_missense_mut_impact_tab.txt", data.table = F)

# input proteins known to interact with p53 -------------------------------
TP53_pair_tab <- fread("./cptac2p/analysis_results/p53/tables/compile_TP53_interactor_pair_tab/TP53_interactor_pair_tab.txt", data.table = F)


# input TF-gene table -----------------------------------------------------
tf_tab <- fread(input = "./Ding_Lab/Projects_Current/TP53_shared_data/resources/PPI/TF_interactions.txt", data.table = F)

# construct table for scatterplot -----------------------------------------
significant_affected_targets <- mut_impact_tab %>%
  filter(fdr < 0.05) %>%
  filter(affected_exp_type == "RNA") %>%
  select(SUB_GENE) %>%
  unique()

significant_affected_targets <- significant_affected_targets$SUB_GENE

## filter for the transcriptional targets of GOF-mediating TFs
TFs2p <- tf_tab %>%
  filter(target_genesymbol %in% significant_affected_targets) %>%
  filter(source_genesymbol %in% TP53_GOF_interacting_TFs) %>%
  select(source_genesymbol) %>%
  unique() %>%
  unlist(use.names = F)
 
TFs2merge <- tf_tab %>%
  filter(source_genesymbol %in% TFs2p) %>%
  select(source_genesymbol, target_genesymbol)

## add TFs into result table
mut_impact_tab.w.TF <- merge(mut_impact_tab, TFs2merge, by.x = c("SUB_GENE"), by.y = c("target_genesymbol"))

# ## filter out affected genes known to directly interact with p53
# mut_impact_tab.w.TF.notTP53interactor <- mut_impact_tab.w.TF %>%
#   filter(!(SUB_GENE %in% TP53_pair_tab$SUB_GENE))

# Plot per cancer ---------------------------------------------------------
for (cancer_tmp in c("BRCA", "OV", "CO", "UCEC", "LIHC")) {
  ## get the genes to highlight: significant by both gene and protein level changes
  genes_RNA <- mut_impact_tab %>%
    filter(cancer == cancer_tmp) %>%
    filter(affected_exp_type == "RNA") %>%
    filter(SELF == "trans") %>%
    filter(fdr < 0.05) %>%
    select(SUB_GENE)
  
  genes_PRO <- mut_impact_tab %>%
    filter(cancer == cancer_tmp) %>%
    filter(affected_exp_type == "PRO") %>%
    filter(SELF == "trans") %>%
    filter(fdr < 0.05) %>%
    select(SUB_GENE)
  
  genes2highlight <- intersect(genes_RNA$SUB_GENE, genes_PRO$SUB_GENE)
  
  for (expression_type_tmp in c("RNA", "PRO")) {
    tab2p <- mut_impact_tab %>%
      filter(cancer == cancer_tmp) %>%
      filter(affected_exp_type == expression_type_tmp) %>%
      filter(SELF == "trans") %>%
      mutate(point_alpha = ifelse(SUB_GENE %in% genes2highlight, 0.8, 0.6)) %>%
      mutate(log10_FDR = -log10(fdr))
    
    cap <- min(max(abs(tab2p$meddiff), na.rm = T), 3)
    tab2p$meddiff_capped <- tab2p$meddiff
    tab2p$meddiff_capped[tab2p$meddiff > cap] <- cap
    tab2p$meddiff_capped[tab2p$meddiff < (-cap)] <- (-cap)
    
    p <- ggplot()
    p <- p + geom_vline(xintercept = 0, linetype = 2, color = "grey50")
    p <- p + geom_point(data = tab2p %>%
                        filter(fdr >= 0.05), 
                        mapping = aes(x = meddiff, y = log10_FDR), color = "grey50", alpha = 0.4, position = pos)
    p <- p + geom_point(data = tab2p %>%
                          filter(fdr < 0.05) %>%
                          filter(meddiff > 0), 
                        mapping = aes(x = meddiff, y = log10_FDR), color = "red", position = pos)
    p <- p + geom_text_repel(data = tab2p %>% 
                               filter(fdr < 0.05) %>%
                               filter(meddiff > 0), 
                             mapping = aes(x = meddiff, y = log10_FDR, label = SUB_GENE, alpha = point_alpha), color = "red", position = pos)
    p <- p + geom_point(data = tab2p %>%
                          filter(fdr < 0.05) %>%
                          filter(meddiff < 0), 
                        mapping = aes(x = meddiff, y = log10_FDR), color = "blue", alpha = 0.8, position = pos)
    p <- p + geom_text_repel(data = tab2p %>% 
                               filter(fdr < 0.05) %>%
                               filter(meddiff < 0), 
                             mapping = aes(x = meddiff, y = log10_FDR, label = SUB_GENE, alpha = point_alpha), color = "blue", position = pos)
    p = p + theme_nogrid()
    p <- p + ylab(paste0("-log10(FDR)")) + xlab(paste0(expression_type_tmp, " level difference\nMedian(GOF Missense) - Median(other Missense)"))
    p <- p + guides(alpha = F)
    p
    fn <- paste0(makeOutDir(resultD = resultD), cancer_tmp, "_", expression_type_tmp, "_", "GOF_missense_vs_other_missense.pdf")
    ggsave(filename = fn, width = 8, height = 6)
    
  }
}


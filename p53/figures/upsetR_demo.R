# Yige Wu @ WashU 2019 Feb
# look at the intersections of affected expression types (RNA, PRO, PHO) across genes


# source -------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
source(paste0(baseD, "cptac2p_analysis/pan3can_aes.R")) # aes for general purposes
library(UpSetR)


# Input mutational impact table -------------------------------------------
mut_impac_tab <- fread(input = "./cptac2p/analysis_results/p53/tables/test_mut_impact_proteome_TP53/TP53_mut_impact_proteome_RNA_cptac2p_cptac3_tab.txt", data.table = F)
mut_impac_tab %>%
  head

# adjust fdr --------------------------------------------------------------
mut_impac_tab$affected_exp_type <- mut_impac_tab$SUB_MOD_RSD
mut_impac_tab$affected_exp_type[!(mut_impac_tab$SUB_MOD_RSD %in% c("RNA", "PRO"))] <- "PHO"

mut_impac_tab$fdr <- FDR_by_id_columns(p_vector = mut_impac_tab$p, 
                                       id_columns = c("variant_class", "SELF", "cancer", "affected_exp_type",
                                                      "SUB_GENE.is_downstream", "SUB_GENE.is_kinase", "SUB_GENE.is_phosphatase", "SUB_GENE.is_complex_corum"), 
                                       df = mut_impac_tab)

# filter FDR --------------------------------------------------------------
mut_impac_tab <- mut_impac_tab[mut_impac_tab$fdr < 0.05,]
mut_impac_tab %>%
  nrow()
# sort data frame for plotting --------------------------------------------
variant_class2filter <- c("not_silent")
fdr2filter <- 0.05
tab2p <- mut_impac_tab
tab2p <- tab2p %>%
  filter(variant_class %in% variant_class2filter) %>%
  filter(fdr < fdr2filter) %>%
  mutate(SUB_GENE_cancer = paste0(SUB_GENE, "_", cancer)) %>%
  select(SUB_GENE_cancer, affected_exp_type) %>%
  unique() %>%
  table
tab2p <- as.matrix(tab2p)

gene_cancers <- str_split_fixed(string = rownames(tab2p), pattern = "_", n = 2)
tab2p <- data.frame(Gene = gene_cancers[,1], cancer = gene_cancers[,2], 
                    RNA = tab2p[, "RNA"], PRO = tab2p[, "PRO"], PHO = tab2p[, "PHO"])

tab2p %>% nrow()
tab2p %>% head()

# plotting ----------------------------------------------------------------
for (cancer_tmp in unique(tab2p$cancer)) {
  tab2p_tmp <- tab2p %>%
    filter(cancer == cancer_tmp)
  if (nrow(unique(tab2p_tmp[, c("RNA", "PRO", "PHO")])) > 1) {
    fn <- paste0(makeOutDir(resultD = resultD), cancer_tmp, "_TP53_", fdr2filter, "FDR", fdr2filter, ".pdf")
    pdf(fn, height = 3, width = 4, onefile = F)
    upset(tab2p_tmp, 
          sets=c("RNA", "PRO", "PHO"), 
          keep.order=T, 
          order.by = "freq")
    dev.off()
  }
}


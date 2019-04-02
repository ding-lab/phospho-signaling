# Yige Wu @ WashU Mar 2019
## alluvial plot showing kinase connecting to pathways

# source ------------------------------------------------------------------
baseD = "/Users/yigewu/Box\ Sync/"
setwd(baseD)
source("./cptac2p_analysis/phospho_network/phospho_network_shared.R")
source("./cptac2p_analysis/phospho_network/phospho_network_plotting.R")
source('/Users/yigewu/Box Sync/cptac2p_analysis/dependencies/tables/load_TCGA_pathways.R')
library(alluvial)

# set variables -----------------------------------------------------------


# input regression --------------------------------------------------------
regression <- fread(input = paste0(ppnD, "regression/tables/annotate_regression_with_mut_impact/", 
                                   "regression_cptac2p_cptac3_tumor_reg_nonNA20_mut_impact_cancer_specificity_annotated.txt"),
                    data.table = F)
reg_sig <- c(0.05, 0.05); names(reg_sig) <- c("kinase", "phosphatase")
regression %>% nrow()
regression <- regression %>%
  filter(pair_pro %in% omnipath_tab$pair_pro[!(omnipath_tab$Source %in% c("NetKIN", "PhosphoNetworks", "MIMP"))] | pair_pro %in% psp_tab$pair_pro)
regression %>% nrow()

regression <- adjust_regression_by_nonNA(regression = regression, reg_nonNA = 20, reg_sig = reg_sig)
regression <- annotate_ks_source(regression = regression)
table(regression$Cancer)

regression %>%
  select(SUB_GENE) %>%
  unique() %>%
  nrow()

regression %>%
  select(GENE) %>%
  unique() %>%
  nrow()

regression %>%
  filter(regulated == T) %>%
  select(SUB_GENE) %>%
  unique() %>%
  nrow()

regression %>%
  filter(regulated == T) %>%
  select(GENE) %>%
  unique() %>%
  nrow()

# annotate substratres ----------------------------------------------------
sub_genes2pathways <- map2TCGApathwaways(gene_list = unique(regression$SUB_GENE), pathway_list = tcga_pathways_pluskegg_and_pathway)

sub_genes2pathways <- data.frame(SUB_GENE = rep(x = names(sub_genes2pathways), sapply(X = sub_genes2pathways, FUN = function(x) length(x))), 
                                 SUB_GENE.path = unlist(sub_genes2pathways, use.names = F))


# construct table for alluvial plot ---------------------------------------


# plot selected kinases by cancer ---------------------------------------------------
for (cancer_tmp in unique(regression$Cancer)) {
  tab2p <- regression %>%
    filter(regulated == T) %>%
    filter(Cancer == cancer_tmp) %>%
    filter(enzyme_type == "kinase") %>%
    select(GENE, SUB_GENE) %>%
    unique
  tab2p <- merge(tab2p, sub_genes2pathways, by = c("SUB_GENE"), all.x = T)  
  tmp <- as.vector(tab2p$SUB_GENE.path)
  tmp[is.na(tmp)] <- "other"
  tab2p$SUB_GENE.pat <- tmp
  
  enzyme_sort_tab <- data.frame(table(tab2p[, c("GENE")]))
  enzyme_sort_tab <- enzyme_sort_tab[order(enzyme_sort_tab$Freq, decreasing = T),]
  
  tab2p_v2 <- data.frame(table(tab2p[, c("GENE", "SUB_GENE.path")]))
  tab2p_v2$color <- tcga_path_colors[as.vector(tab2p_v2$SUB_GENE.path)]
  tab2p_v2$GENE <- factor(tab2p_v2$GENE, levels = enzyme_sort_tab$Var1)
  
  tab2p_v2 <- tab2p_v2 %>%
    filter(GENE %in% enzyme_sort_tab$Var1[1:10])
  
  fn <- paste0(makeOutDir(resultD = resultD), cancer_tmp, "_selected.pdf")
  pdf(file = fn, width = 5, height = 8, useDingbats = F)
  alluvial(tab2p_v2[, 1:2], freq=tab2p_v2$Freq,
           col = tab2p_v2$color,
           border = NA,
           hide = tab2p_v2$Freq == 0,
           blocks= "bookends",
           axis_labels = c("Kinase", "Substrate Pathway"),
           cex = 0.7
  )
  dev.off()
  
}


# plot all kinases ---------------------------------------------------
tab2p_v2 <- data.frame(table(tab2p[, c("GENE", "SUB_GENE.path")]))
tab2p_v2$color <- tcga_path_colors[as.vector(tab2p_v2$SUB_GENE.path)]
tab2p_v2$GENE <- factor(tab2p_v2$GENE, levels = enzyme_sort_tab$Var1)

# tab2p_v2 <- tab2p_v2 %>%
#   filter(GENE %in% enzyme_sort_tab$Var1[1:10])

alluvial(tab2p_v2[, 1:2], freq=tab2p_v2$Freq,
         col = tab2p_v2$color,
         border = NA,
         hide = tab2p_v2$Freq == 0,
         blocks= "bookends",
         cex = 0.7
)

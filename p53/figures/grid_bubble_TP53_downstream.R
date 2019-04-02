# Yige Wu @ WashU Mar 2019
## to show the downstream effects of TP53 mutations

# source ------------------------------------------------------------------
setwd(dir = "Box Sync/")
source("./cptac2p_analysis/phospho_network/phospho_network_plotting.R")
source("./cptac2p_analysis/phospho_network/phospho_network_shared.R")
source("./cptac2p_analysis/p53/TP53_shared.R")


# input mutational associations -------------------------------------------
mut_impact_tab <- fread(input = "./cptac2p/analysis_results/p53/tables/test_mut_impact_proteome_TP53/TP53_mut_impact_proteome_RNA_cptac2p_cptac3_tab.txt", data.table = F)

# filter for downstrea ----------------------------------------------------
mut_impact_tab_downstream <- mut_impact_tab %>%
  filter(SUB_GENE.is_downstream == T)

mut_impact_tab %>%
  filter(cancer == "CCRCC" & variant_class == "truncation")

mut_impact_tab %>%
  filter(SUB_GENE == "ESR1")

## make sure every cancer type is presented in each facet
mut_impact_tab_downstream2add <- mut_impact_tab_downstream %>%
  filter(cancer == "CCRCC" & variant_class == "missense") %>%
  mutate(variant_class = "truncation") %>%
  mutate(meddiff = NA) %>%
  mutate(fdr = NA)

mut_impact_tab_downstream <- rbind(mut_impact_tab_downstream, mut_impact_tab_downstream2add)

# set variables -----------------------------------------------------------
gene_altered <- "TP53"
fdr_thres2p <- 0.05
sig_colors <- c("black", "white")
names(sig_colors) <- c(paste0('FDR<', fdr_thres2p), paste0('FDR>', fdr_thres2p))
expression_cap <- 1.5
expression_types2p <- c("RNA", "PRO")

# input TF relation table ----------------------------------------------------------
# TF_tab <- fread(input = "./Ding_Lab/Projects_Current/TP53_shared_data/resources/PPI/TF_interactions.txt", data.table = F)
TF_tab <- fread(input = "./Ding_Lab/Projects_Current/TP53_shared_data/resources/PPI/TF_interactions_TP53_manual.txt", data.table = F)
TF_tab <- TF_tab %>%
  mutate(TF_effect = Effect)

TF_tab %>%
  filter(source_genesymbol == "TP53" & target_genesymbol == "ESR1")

TF_tab %>%
  filter(source_genesymbol == gene_altered) %>%
  select(TF_effect) %>%
  table()

for (TF_effect_tmp in unique(TF_tab$TF_effect)) {
  tab2p <- mut_impact_tab_downstream %>%
    filter(SUB_GENE %in% TF_tab$target_genesymbol[TF_tab$source_genesymbol == gene_altered & TF_tab$TF_effect == TF_effect_tmp]) %>%
    filter(affected_exp_type %in% expression_types2p)
  
  ## only plot genes with significant results
  affected_genes2p <- unique(tab2p$SUB_GENE[tab2p$fdr_sig & order(tab2p$meddiff)])
  tab2p <- tab2p %>%
    filter(SUB_GENE %in% affected_genes2p)

  tab2p$meddiff_capped <- add_cap(x = tab2p$meddiff, cap = expression_cap)
  tab2p$log10_FDR <- -log10(x = tab2p$fdr)
  tab2p$cancer <- order_cancer_rev(x = tab2p$cancer)
  tab2p$affected_exp_type <- factor(x = tab2p$affected_exp_type, levels = expression_types2p)
  tab2p$variant_class <- order_variant_class(x = tab2p$variant_class)
  tab2p$SUB_GENE <- factor(x = tab2p$SUB_GENE, levels = affected_genes2p)
  
  p <- ggplot()
  p <- p + geom_point(data = tab2p, mapping = aes(x = SUB_GENE, y = cancer, fill = meddiff_capped, size = log10_FDR,
                                               colour = ifelse(fdr_sig == TRUE, paste0('FDR<', fdr_thres2p), paste0('FDR>', fdr_thres2p))), 
                      alpha = 0.6, shape = 22)
  p <- p + facet_grid(affected_exp_type + variant_class ~ ., drop=F, space = "free",scales = "free")
  p <- p + scale_color_manual(values = sig_colors)
  p = p + scale_fill_gradientn(name= "Mut-WT\nExpression\nFold Change", na.value=NA, colours=RdBu1024, limit=c(-expression_cap,expression_cap))
  p <- p + xlab("Affected Gene") + ylab("Cancer Type")
  p <- p + guides(colour = guide_legend(title = "FDR"))
  p <- p + theme_bw()
  p <- p + theme(panel.spacing.x = unit(0, "lines"), panel.spacing.y = unit(0, "lines"))
  p <- p + theme(axis.text.y = element_text(size = 10, face = "bold"))
  p <- p + theme(axis.text.x = element_text(size = 10, face = "bold", angle = 90, hjust = 0.5, vjust = 0.5))
  p <- p + theme(panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "grey80"))
  p <- p + theme(strip.text.x = element_text(angle = 90))
  p
  fn <- paste0(makeOutDir(resultD = resultD), gene_altered, "_" , TF_effect_tmp, "_targets_expression_changes.pdf")
  if (TF_effect_tmp == "Unknown") {
    ggsave(filename = fn,
           width = 10, height = 6)
  }
  if (TF_effect_tmp == "Activation" | TF_effect_tmp == "Inhibition") {
    ggsave(filename = fn,
           width = 10, height = 6)
  }  

}


# Yige Wu @ WashU 2018 Mar
## scatter plot for comparing correlations of prospective BRCA PRM and TMT protein-level abundance with mRNA abundance

# souce -------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/preprocess_files/preprocess_files.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_plotting.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
library(ggrepel)

# inputs ------------------------------------------------------------------
## cancer type
cancer <- "BRCA"

coef_names <- c("product moment correlation coefficient", "rho")
names(coef_names) <- c("pearson", "spearman")



# scatterplot -------------------------------------------------------------
for (method in c("spearman")) {
  sig <- 0.1
  cor_prm <- fread(input = paste0("./cptac2p/cptac_shared/analysis_results/validation/tables/cor_prm&tmt_vs_mRNA/", cancer, "_", method, "_correlation_gene_level_prm_mRNA.txt"))
  cor_prm <- data.frame(cor_prm)
  cor_tmt <- fread(input = paste0("./cptac2p/cptac_shared/analysis_results/validation/tables/cor_prm&tmt_vs_mRNA/", cancer, "_", method, "_correlation_gene_level_tmt_mRNA_limittoprmgenes.txt"))
  cor_tmt <- data.frame(cor_tmt)
  int_cols <- intersect(colnames(cor_prm), colnames(cor_tmt))
  cor_df <- merge(cor_tmt[,int_cols], cor_prm[,int_cols], 
                  by = c("Gene", "ProteinNames", "GroupName", "is_kinase", "panel"), suffixes = c(".tmt",".prm"), all.y = T)
  cor_df <- unique(cor_df)
  for (panel in unique(cor_df$panel)[1]) {
    cor_df_panel <- cor_df[cor_df$panel == panel,]
    p = ggplot(cor_df_panel,aes(x=rho.prm, y=rho.tmt, color = GroupName))
    p = p + geom_point(alpha=0.5, stroke = 0)
    p = p + geom_text_repel(aes(rho.prm, rho.tmt, label= as.character(Gene), color = GroupName, segment.color = GroupName),
                            force = 1, 
                            segment.size = 0.5, segment.alpha = 0.2, 
                            size=1.5,alpha=0.8)
    p = p + geom_abline(slope = 1, intercept = 0, linetype = 2, alpha = 0.8)
    p = p + geom_vline(xintercept = 0, color = 'grey', alpha = 0.5)
    p = p + geom_hline(yintercept = 0, color = 'grey', alpha = 0.5)
    p <- p + facet_grid(.~panel, scales = "free_x", space = "free_x")
    p = p + labs(x = paste0(method, "'s ", coef_names[method], "(PRM~mRNA)"), 
                 y=paste0(method, "'s ", coef_names[method], "(global~mRNA)"), 
                 color = "kinase family")
    p = p + theme_nogrid()
    p = p + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
    p = p + theme(title = element_text(size = 8),
                  legend.key.size = unit(1, "cm"), legend.text = element_text(size = 5), 
                  legend.background = element_rect(fill = NA), legend.box.background = element_rect(fill = NA),
                  legend.key = element_rect(fill = NA), legend.key.height = unit(0.5, "cm"),
                  legend.title = element_text(size = 5), legend.justification = c(1, 0), legend.position = c(1, 0))
    p
    resultDnow <- makeOutDir()
    fn = paste0(resultDnow, cancer, "_", panel, "_", method, "_correlations_gene_level_prm_or_tmt~mRNA_facet_kinase_metabolic", "_fdr",sig,".pdf")
    ggsave(file=fn, height=6, width=6, useDingbats=FALSE)
  }
}


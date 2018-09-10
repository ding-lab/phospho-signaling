# Yige Wu @ WashU 2018 Mar
## volcano plot for correlation of prospective BRCA PRM/TMT protein-level abundance with mRNA abundance
## match parent specimen labels
## TODO: plot per-gene/per-sample density curves before correlation for both MS and PRM data

# souce -------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/preprocess_files/preprocess_files.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_plotting.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
library(ggrepel)

# inputs ------------------------------------------------------------------
## input PRM protein-level values
prm_parent_labeled <- fread(input = "./cptac2p/cptac_shared/analysis_results/validation/tables/cor_prm_ms_match_specimen/BRCA_prm_trplicatemean_median_parent_labeled.txt", data.table = F)
## input TMT protein-level values
ms_parent_labeled <- fread(input = "./cptac2p/cptac_shared/analysis_results/validation/tables/cor_prm_ms_match_specimen/BRCA_tmt_match_prm_parent_labeled.txt", data.table = F)

## cancer type
cancer <- "BRCA"

# prm~mRNA ----------------------------------------------------------------
coef_names <- c("product moment correlation coefficient", "rho")
names(coef_names) <- c("pearson", "spearman")
for (method in c("spearman", "pearson")) {
  cor_df <- fread(input = paste0("./cptac2p/cptac_shared/analysis_results/validation/tables/cor_prm&tmt_vs_mRNA/", cancer, "_", method, "_correlation_gene_level_prm_mRNA.txt"))
  cor_df <- data.frame(cor_df)
  cor_df <- cor_df[!is.na(cor_df$panel),]
  for (panel in unique(cor_df$panel)) {
    cor_df_panel <- cor_df[cor_df$panel == panel,]
    for (sig in c(0.05, 0.1)) {
      sink(file = paste0(resultDnow, cancer, "_", panel, "_", method, "_correlation_gene_level_prm_mRNA", "_fdr",sig,".txt"))
      cat(paste0(length(which(!is.na(cor_df_panel$rho))), " proteins in ", panel, " were tested for correlation btw protein and mRNA abundance\n"))
      cat(paste0(length(which(cor_df_panel$rho > 0)), " proteins in ", panel, " has positive correlation btw protein and mRNA abundance\n"))
      cat(paste0(length(which(cor_df_panel$rho > 0 & cor_df_panel$fdr <= sig)), " proteins in ", panel, " has significant positive correlation btw protein and mRNA abundance\n"))
      sink()
      closeAllConnections()
      p = ggplot(cor_df_panel,aes(x=rho, y=log_fdr, color = GroupName))
      p = p + geom_point(alpha=0.5, stroke = 0)
      p = p + geom_text_repel(aes(rho, log_fdr, label= ifelse(log_fdr > -log10(sig), as.character(Gene), NA), color = GroupName, segment.color = GroupName),
                              force = 1, 
                              segment.size = 0.5, segment.alpha = 0.2, 
                              size=1.5,alpha=0.8)
      p = p + geom_vline(xintercept = 0, color = 'grey')
      p = p + labs(x = paste0(method, "'s ", coef_names[method]), y="-log10(FDR)", color = paste0(panel, " group"))
      p = p + ggtitle(label = paste0(panel, "(", length(which(!is.na(cor_df_panel$rho))), " proteins)" ), 
                      subtitle = paste0(method, " correlation of PRM protein abundance and mRNA abundance(N=", length(prm_mrna_partIDs_int), ")"))
      p = p + theme_nogrid()
      p = p + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
      p = p + theme(title = element_text(size = 8),
                    legend.key.size = unit(1, "cm"), legend.text = element_text(size = 5), 
                    legend.background = element_rect(fill = NA), legend.box.background = element_rect(fill = NA),
                    legend.key = element_rect(fill = NA), legend.key.height = unit(0.5, "cm"),
                    legend.title = element_text(size = 5), legend.justification = c(0, 1), legend.position = c(0, 1))
      p
      resultDnow <- makeOutDir()
      fn = paste0(resultDnow, cancer, "_", panel, "_", method, "_correlation_gene_level_prm_mRNA", "_fdr",sig,".pdf")
      ggsave(file=fn, height=6, width=6, useDingbats=FALSE)
    }
  }
}
# write.table(unique(cor_df$Gene[!(cor_df$is_kinase)]), file = paste0(makeOutDir(), "proteins_in_prm_notkinase.txt"), row.names = F, quote = F, col.names = F)

for (method in c("spearman", "pearson")) {
  cor_df <- fread(input = paste0("./cptac2p/cptac_shared/analysis_results/validation/tables/cor_prm&tmt_vs_mRNA/", cancer, "_", method, "_correlation_gene_level_prm_mRNA.txt"))
  cor_df$panel <- ifelse(cor_df$is_kinase, "Kinome_panel", "Metabolome_panel")
  for (sig in c(0.05, 0.1)) {
    p = ggplot(cor_df,aes(x=rho, y=log_fdr, color = GroupName))
    p = p + geom_point(alpha=0.5, stroke = 0)
    p = p + geom_text_repel(aes(rho, log_fdr, label= ifelse(log_fdr > -log10(sig), as.character(Gene), NA), color = GroupName, segment.color = GroupName),
                            force = 1, 
                            segment.size = 0.5, segment.alpha = 0.2, 
                            size=1.5,alpha=0.8)
    p = p + geom_vline(xintercept = 0, color = 'grey')
    p <- p + facet_grid(.~panel, scales = "free_x", space = "free_x")
    p = p + labs(x = paste0(method, "'s ", coef_names[method]), y="-log10(FDR)", color = "enzyme family")
    p = p + ggtitle(label = paste0(method, " correlation of PRM protein abundance and mRNA abundance(gene-level)"), subtitle = paste0("matching patient IDs(N=", length(prm_mrna_partIDs_int), ")"))
    p = p + theme_nogrid()
    p = p + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
    p = p + theme(title = element_text(size = 8))
    p
    resultDnow <- makeOutDir()
    fn = paste0(resultDnow, cancer, "_", method, "_correlation_gene_level_prm_mRNA_facet_kinase_metabolic", "_fdr",sig,".pdf")
    ggsave(file=fn, height=5, width=10, useDingbats=FALSE)
  }
}


# tmt~mRNA only prm genes----------------------------------------------------------------
#for (method in c("spearman")) {
for (method in c("spearman", "pearson")) {
  cor_df <- fread(input = paste0("./cptac2p/cptac_shared/analysis_results/validation/tables/cor_prm&tmt_vs_mRNA/", cancer, "_", method, "_correlation_gene_level_tmt_mRNA_limittoprmgenes.txt"))
  cor_df <- data.frame(cor_df)
  for (panel in unique(cor_df$panel)) {
    cor_df_panel <- cor_df[cor_df$panel == panel,]
    for (sig in c(0.05, 0.1)) {
      sink(paste0(resultDnow, cancer, "_", panel, "_", method, "_correlation_gene_level_tmt_mRNA_limittoprmgenes", "_fdr",sig,".txt"))
      cat(paste0(length(which(!is.na(cor_df_panel$rho))), " proteins in ", panel, " were tested for correlation btw protein and mRNA abundance\n"))
      cat(paste0(length(which(cor_df_panel$rho > 0 & !is.na(cor_df_panel$rho))), " proteins in ", panel, " has positive correlation btw protein and mRNA abundance\n"))
      cat(paste0(length(which(cor_df_panel$rho > 0 & cor_df_panel$fdr <= sig & !is.na(cor_df_panel$rho))), " proteins in ", panel, " has significant positive correlation btw protein and mRNA abundance\n"))
      sink()
      closeAllConnections()
      p = ggplot(cor_df_panel,aes(x=rho, y=log_fdr, color = GroupName))
      p = p + geom_point(alpha=0.5, stroke = 0)
      p = p + geom_text_repel(aes(rho, log_fdr, label= ifelse(log_fdr > -log10(sig), as.character(Gene), NA), color = GroupName, segment.color = GroupName),
                              force = 1, 
                              segment.size = 0.5, segment.alpha = 0.2, 
                              size=1.5,alpha=0.8)
      p = p + geom_vline(xintercept = 0, color = 'grey')
      p = p + labs(x = paste0(method, "'s ", coef_names[method]), y="-log10(FDR)", color = paste0(panel, " group"))
      p = p + ggtitle(label = paste0(panel, "(", length(which(!is.na(cor_df_panel$rho))), " proteins)" ), 
                      subtitle = paste0(method, " correlation of TMT10 protein abundance and mRNA abundance(N=", length(prm_mrna_partIDs_int), ")"))
      p = p + theme_nogrid()
      p = p + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
      p = p + theme(title = element_text(size = 8),
                    legend.key.size = unit(1, "cm"), legend.text = element_text(size = 5), 
                    legend.background = element_rect(fill = NA), legend.box.background = element_rect(fill = NA),
                    legend.key = element_rect(fill = NA), legend.key.height = unit(0.5, "cm"),
                    legend.title = element_text(size = 5), legend.justification = c(0, 1), legend.position = c(0, 1))
      p
      resultDnow <- makeOutDir()
      fn = paste0(resultDnow, cancer, "_", panel, "_", method, "_correlation_gene_level_tmt_mRNA_limittoprmgenes", "_fdr",sig,".pdf")
      ggsave(file=fn, height=6, width=6, useDingbats=FALSE)
    }
  }
}



for (method in c("spearman", "pearson")) {
  cor_df <- fread(input = paste0("./cptac2p/cptac_shared/analysis_results/validation/tables/cor_prm&tmt_vs_mRNA/", cancer, "_", method, "_correlation_gene_level_tmt_mRNA_limittoprmgenes.txt"))
  cor_df$panel <- ifelse(cor_df$is_kinase, "Kinome_panel", "Metabolome_panel")
  for (sig in c(0.05, 0.1)) {
    cor_df <- cor_df[!(cor_df$fdr > sig & cor_df$rho > quantile(x = cor_df$rho, probs = 0.9, na.rm = T)),]
    
    p = ggplot(cor_df,aes(x=rho, y=log_fdr, color = GroupName))
    p = p + geom_point(alpha=0.5, stroke = 0)
    p = p + geom_text_repel(aes(rho, log_fdr, label= ifelse(log_fdr > -log10(sig), as.character(Gene), NA), color = GroupName, segment.color = GroupName),
                            force = 1, 
                            segment.size = 0.5, segment.alpha = 0.2, 
                            size=1.5,alpha=0.8)
    p = p + geom_vline(xintercept = 0, color = 'grey')
    p <- p + facet_grid(.~panel, scales = "free_x", space = "free_x")
    p = p + labs(x = paste0(method, "'s ", coef_names[method]), y="-log10(FDR)", color = "enzyme family")
    p = p + ggtitle(label = paste0(method, " correlation of TMT10 protein abundance and mRNA abundance(gene-level)"), subtitle = paste0("only examine proteins detected by PRM, matching patient IDs(N=", length(tmt_mrna_partIDs_int), ")"))
    p = p + theme_nogrid()
    p = p + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
    p = p + theme(title = element_text(size = 8))
    p
    resultDnow <- makeOutDir()
    fn = paste0(resultDnow, cancer, "_", method, "_correlation_gene_level_tmt_mRNA_limittoprmgenes_facet_kinase_metabolic", "_fdr",sig,".pdf")
    ggsave(file=fn, height=5, width=10, useDingbats=FALSE)
  }
}






# Yige Wu @ WashU 2018 Mar
## correlate prospective BRCA PRM peptide-level abundance with MS/MS protein abundance from CDAP
## match parent specimen labels
## TODO: plot per-gene/per-sample density curves before correlation for both MS and PRM data

# souce -------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/preprocess_files/preprocess_files.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_plotting.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')


# inputs ------------------------------------------------------------------
cancer <- "BRCA"
prm_mean_median <- fread(input = paste0("./cptac2p/cptac_shared/analysis_results/preprocess_files/format_RPM_peptide/", "prm_kinase_median_per_triplicate_specimen_labeled.txt"),
                         data.table = F)
Pro.n_tumor2p <- fread(input = paste0(inputD, cancer, "/", prefix[cancer], "_PRO_formatted_normalized_noControl.txt"),
                       data.table = F)
clinical <- read.delim(paste0(inputD, "Specimen_Data_20161005_Yige_20180127.txt"), sep = "\t")
clinical <- data.frame(clinical)

# format ------------------------------------------------------------------
## transform MS protein columns to patient IDs
ms_sampIDs <- colnames(Pro.n_tumor2p)
ms_parent_sampIDs <- str_split_fixed(string = ms_sampIDs, pattern = "_", 2)[,1]
prm_sampIDs <- colnames(prm_mean_median)
prm_parent_sampIDs <- str_split_fixed(string = prm_sampIDs, pattern = "_", 2)[,1]

## get the intersect of patient IDs
parent_sampIDs_int <- intersect(ms_parent_sampIDs, prm_parent_sampIDs)
genes_int <- intersect(Pro.n_tumor2p$Gene, prm_mean_median$Gene_Name)
prm_parent_labeled <- prm_mean_median; colnames(prm_parent_labeled) <- prm_parent_sampIDs
ms_parent_labeled <- Pro.n_tumor2p; colnames(ms_parent_labeled) <- ms_parent_sampIDs
# write.table(prm_parent_labeled, row.names = F, quote=F, sep = '\t', file=paste(resultDnow, cancer, "_prm_trplicatemean_median_parent_labeled.txt",sep=""))
# write.table(ms_parent_labeled, row.names = F, quote=F, sep = '\t', file=paste(resultDnow, cancer, "_tmt_match_prm_parent_labeled.txt",sep=""))

library(ggrepel)
coef_names <- c("product moment correlation coefficient", "rho")
names(coef_names) <- c("pearson", "spearman")
for (method in c("pearson", "spearman")) {
  cor_stats <- lapply(1:length(genes_int), FUN = function(n, m1, m2, gs, ps, mt) cor.test(x = as.numeric(as.vector(m1[m1$Gene == gs[n], ps])), 
                                                                                          y = as.numeric(as.vector(m2[m2$Gene == gs[n], ps])), 
                                                                                          method = mt),
         m1 = prm_parent_labeled, m2 = ms_parent_labeled, gs = genes_int, ps = parent_sampIDs_int, mt = method)
  cor_df <- data.frame(Gene = genes_int,
                       pvalue = sapply(1:length(genes_int), FUN = function(n, l) l[[n]]$p.value, l = cor_stats),
                       rho = sapply(1:length(genes_int), FUN = function(n, l) l[[n]]$estimate[[1]], l = cor_stats))
  cor_df$fdr <- p.adjust(cor_df$pvalue,method = "fdr")
  cor_df$log_fdr <- -log10(cor_df$fdr)
  write.table(cor_df, row.names = F, quote=F, sep = '\t', file=paste(resultDnow, cancer, "_", method,"_correlation_prm_ms_match_parent_specimen.txt",sep=""))
  
  for (sig in c(0.05, 0.1)) {
    p = ggplot(cor_df,aes(x=rho, y=log_fdr))
    p = p + geom_point(alpha=0.5, stroke = 0 , color = "grey")
    p = p + geom_text_repel(aes(rho, log_fdr, label= ifelse(log_fdr > -log10(sig), as.character(Gene), NA)),
                            force = 1, segment.color = '#cccccc', segment.size = 0.5, segment.alpha = 0.2, 
                            size=1.5,alpha=0.8)
    # p = p + geom_text(aes(label= ifelse(-log10(fdr) > -log10(0.05), as.character(Gene), NA), color = "red"),size=1.5,alpha=0.8)
    p = p + geom_vline(xintercept = 0, color = 'grey')
    p = p + labs(x = paste0(method, "'s ", coef_names[method]), y="-log10(FDR)")
    p = p + ggtitle(label = paste0(method, " correlation of PRM and TMT detected abundance"), subtitle = paste0("matching parent specimen labels(N=", length(parent_sampIDs_int), ")"))
    p = p + theme_nogrid()
    p = p + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
    p = p + theme(title = element_text(size = 8))
    p = p + theme(legend.position = "none")
    p
    resultDnow <- makeOutDir()
    fn = paste0(resultDnow, method,'_correlation_prm_ms_match_parent_specimen_for_kinases_in_', cancer, "_fdr",sig,".pdf")
    ggsave(file=fn, height=6, width=6, useDingbats=FALSE)
  }
}

## correlation
## direction
## across peptide per protein
## within triplicate
## across triplicate
## PRM captured but TMT missed


# Yige Wu @ WashU 2018 Mar
## correlate prospective BRCA PRM peptide-level abundance with MS/MS protein abundance from CDAP
## TODO: plot per-gene/per-sample density curves before correlation for both MS and PRM data

# souce -------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/preprocess_files/preprocess_files.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_plotting.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')


# inputs ------------------------------------------------------------------
cancer <- "BRCA"
prm_mean_median <- fread(input = paste0("./cptac2p/cptac_shared/analysis_results/preprocess_files/format_RPM_peptide/", "prm_kinase_median_per_triplicate.txt"),
                         data.table = F)
Pro.n_tumor2p <- fread(input = paste0(inputD, cancer, "/", prefix[cancer], "_PRO_formatted_normalized_noControl.txt"),
                       data.table = F)
clinical <- read.delim(paste0(inputD, "Specimen_Data_20161005_Yige_20180127.txt"), sep = "\t")
clinical <- data.frame(clinical)

# format ------------------------------------------------------------------
## transform MS protein columns to patient IDs
sampIDs <- colnames(Pro.n_tumor2p[,-1])
partIDs <- sampID2partID(sampleID_vector = sampIDs, sample_map = clinical)

colnames(Pro.n_tumor2p)[-1] <- partIDs

## get the intersect of patient IDs
partIDs_int <- intersect(partIDs, colnames(prm_mean_median))
genes_int <- intersect(Pro.n_tumor2p$Gene, prm_mean_median$Gene_Name)
prm_int <- prm_mean_median[, partIDs_int]
ms_int <- Pro.n_tumor2p[, partIDs_int]

for (method in c("pearson", "spearman")) {
  cor_stats <- lapply(1:length(genes_int), FUN = function(n, m1, m2, gs, ps) cor.test(x = as.numeric(as.vector(m1[m1$Gene_Name == gs[n], ps])), y = as.numeric(as.vector(m2[m2$Gene == gs[n], ps])), method = "spearman"),
         m1 = prm_mean_median, m2 = Pro.n_tumor2p, gs = genes_int, ps = partIDs_int)
  cor_df <- data.frame(Gene = genes_int,
                       pvalue = sapply(1:length(genes_int), FUN = function(n, l) l[[n]]$p.value, l = cor_stats),
                       rho = sapply(1:length(genes_int), FUN = function(n, l) l[[n]]$estimate[[1]], l = cor_stats))
  cor_df$fdr <- p.adjust(cor_df$pvalue,method = "fdr")
  
  p = ggplot(cor_df,aes(x=rho, y=-log10(fdr)))
  p = p + geom_point(alpha=0.5 ,  stroke = 0 , color = "grey")
  p = p + geom_text(aes(label= ifelse(-log10(fdr) > 2, as.character(Gene), NA), color = "red"),size=1.5,alpha=0.5)
  p = p + theme_bw() #+ theme_nogrid()
  p = p + geom_vline(xintercept = 0, color = 'grey')
  p = p + theme(axis.title = element_text(size=10), axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
  p = p + labs(x = paste0(method, " correlation coefficient"), y="-log10(FDR)")
  p
  resultDnow <- makeOutDir()
  fn = paste0(resultDnow, method,'_correlation_prm_ms_kinases_', cancer, ".pdf")
  ggsave(file=fn, height=6, width=6, useDingbats=FALSE)
  
}




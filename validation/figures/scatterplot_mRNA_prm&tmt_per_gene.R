# Yige Wu @ WashU 2018 Mar
## scatter plot for comparing correlations of prospective BRCA PRM and TMT protein-level abundance with mRNA abundance

# souce -------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/preprocess_files/preprocess_files.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_plotting.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
library(ggrepel)

# inputs ------------------------------------------------------------------
cancer <- "BRCA"
## input PRM protein-level values
prm_PPI_labeled <- fread(input = "./cptac2p/cptac_shared/analysis_results/preprocess_files/format_RPM_peptide/prm_kinase_median_per_triplicate_PPI_labeled.txt", data.table = F)
## input sample mapping info
clinical <- read.delim(paste0(inputD, "Specimen_Data_20161005_Yige_20180127.txt"), sep = "\t")
clinical <- data.frame(clinical)

## input TMT protein-level values
tmt_specimen_labeled <- fread(input = paste0(inputD, cancer, "/", prefix[cancer], "_PRO_formatted_normalized_noControl.txt"),
                              data.table = F)
tmt_specimen_labeled.t <- get_tumor(expression = tmt_specimen_labeled, clinical.m = clinical)
colnames(tmt_specimen_labeled.t) <- sampID2partID(sampleID_vector = colnames(tmt_specimen_labeled.t), sample_map = clinical)
tmt_PPI_labeled.t <- data.frame(gene = tmt_specimen_labeled$Gene)
tmt_PPI_labeled.t <- cbind(tmt_PPI_labeled.t, tmt_specimen_labeled.t)
## input mRNA abundance
fpkm <- fread(input = paste0("./cptac2p/FPKM/BRCA/11062017/fpkm_prospective_breastcancer_110617.csv"), data.table = F)
## log transform mRNA abundance
fpkm.log2 <- data.frame(gene = fpkm$gene)
fpkm.log2 <- cbind(fpkm.log2, log2(as.matrix(fpkm[, !(colnames(fpkm) %in% c("gene"))])+0.01))
## input PRM peptide-level values
prm_raw.f <- fread(input = "./cptac2p/cptac_shared/analysis_results/preprocess_files/format_RPM_peptide/prm_genename_annotated_specimen_labeled.txt", data.table = F)

## input kinase family info
kinase_family_per_synonym <- fread(input = "./cptac2p/cptac_shared/analysis_results/validation/tables/cor_prm_ms_match_specimen/kinase_family_per_synonym,txt", data.table = F)
## get intersecting genes and patient IDs
prm_mrna_genes_int <- intersect(prm_PPI_labeled$Gene_Name, fpkm$gene)
prm_mrna_partIDs_int <- intersect(colnames(prm_PPI_labeled), colnames(fpkm))
tmt_mrna_genes_int <- intersect(tmt_PPI_labeled.t$gene, fpkm$gene)
tmt_mrna_partIDs_int <- intersect(colnames(tmt_PPI_labeled.t), colnames(fpkm))
prm_tmt_mrna_genes_int <- intersect(prm_mrna_genes_int, tmt_mrna_genes_int)


# faceted scatterplot remove x&y outliers -------------------------------------------------------
resultDnow <- makeOutDir()
subdir <- paste0(resultDnow, "rm_x&y_outlier/")
dir.create(subdir)
for (gene in prm_tmt_mrna_genes_int) {
  df_prm <- data.frame(protein = as.numeric(as.vector(prm_PPI_labeled[prm_PPI_labeled$Gene == gene, prm_mrna_partIDs_int])), 
                         fpkm.log2 = as.numeric(as.vector(fpkm.log2[fpkm.log2$gene == gene, prm_mrna_partIDs_int])),
                         technology = "PRM")
  df_prm_rm_outlier <- df_prm[df_prm$fpkm.log2 <= quantile(df_prm$fpkm.log2, probs = 0.9, na.rm = T) & df_prm$fpkm.log2 >= quantile(df_prm$fpkm.log2, probs = 0.1, na.rm = T) & df_prm$protein <= quantile(df_prm$protein, probs = 0.9, na.rm = T),]
  
  df_tmt <- data.frame(protein = as.numeric(as.vector(tmt_PPI_labeled.t[tmt_PPI_labeled.t$gene == gene, tmt_mrna_partIDs_int])), 
                        fpkm.log2 = as.numeric(as.vector(fpkm.log2[fpkm.log2$gene == gene, tmt_mrna_partIDs_int])),
                         technology = "global")
  df_tmt_rm_outlier <- df_tmt[df_tmt$fpkm.log2 <= quantile(df_tmt$fpkm.log2, probs = 0.9, na.rm = T) & df_tmt$fpkm.log2 >= quantile(df_tmt$fpkm.log2, probs = 0.1, na.rm = T) & df_tmt$protein <= quantile(df_tmt$protein, probs = 0.9, na.rm = T),]
  
  df_rm_outlier <- rbind(df_prm_rm_outlier, df_tmt_rm_outlier)
  df_rm_outlier <- df_rm_outlier[!is.na(df_rm_outlier$fpkm.log2),]
  p = ggplot(df_rm_outlier, aes(x=fpkm.log2, y=protein))
  p = p + geom_point(alpha=0.5, stroke = 0)
  p <- p + facet_grid(technology~., scales = "free", space = "fixed", drop = T)
  p = p + labs(x = paste0("mRNA abundance(log2(FPKM+0.01))"), 
               y=paste0("PRM/global protein abundance"))
  p = p + theme_nogrid()
  p = p + theme(axis.title = element_text(size=10), 
                axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), 
                axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
  p = p + theme(title = element_text(size = 8))
  p
  resultDnow <- makeOutDir()
  fn = paste0(subdir, cancer, "_", gene, "_prm_or_tmt~mRNA.pdf")
  ggsave(file=fn, height=6, width=4, useDingbats=FALSE)
}

# faceted scatterplot remove x outliers -------------------------------------------------------
resultDnow <- makeOutDir()
subdir <- paste0(resultDnow, "rm_x_outlier/")
dir.create(subdir)
for (gene in prm_tmt_mrna_genes_int) {
  df_prm <- data.frame(protein = as.numeric(as.vector(prm_PPI_labeled[prm_PPI_labeled$Gene == gene, prm_mrna_partIDs_int])), 
                       fpkm.log2 = as.numeric(as.vector(fpkm.log2[fpkm.log2$gene == gene, prm_mrna_partIDs_int])),
                       technology = "PRM")
  df_prm_rm_outlier <- df_prm[df_prm$fpkm.log2 <= quantile(df_prm$fpkm.log2, probs = 0.9, na.rm = T) & df_prm$fpkm.log2 >= quantile(df_prm$fpkm.log2, probs = 0.1, na.rm = T) ,]
  
  df_tmt <- data.frame(protein = as.numeric(as.vector(tmt_PPI_labeled.t[tmt_PPI_labeled.t$gene == gene, tmt_mrna_partIDs_int])), 
                       fpkm.log2 = as.numeric(as.vector(fpkm.log2[fpkm.log2$gene == gene, tmt_mrna_partIDs_int])),
                       technology = "global")
  df_tmt_rm_outlier <- df_tmt[df_tmt$fpkm.log2 <= quantile(df_tmt$fpkm.log2, probs = 0.9, na.rm = T) & df_tmt$fpkm.log2 >= quantile(df_tmt$fpkm.log2, probs = 0.1, na.rm = T) ,]
  df_rm_outlier <- rbind(df_prm_rm_outlier, df_tmt_rm_outlier)
  df_rm_outlier <- df_rm_outlier[!is.na(df_rm_outlier$fpkm.log2),]
  p = ggplot(df_rm_outlier, aes(x=fpkm.log2, y=protein))
  p = p + geom_point(alpha=0.5, stroke = 0)
  p <- p + facet_grid(technology~., scales = "free_y", space = "fixed")
  p = p + labs(x = paste0("mRNA abundance(log2(FPKM+0.01))"), 
               y=paste0("PRM/global protein abundance"))
  p = p + theme_nogrid()
  p = p + theme(axis.title = element_text(size=10), 
                axis.text.x = element_text(colour="black", size=6,angle=90, vjust=0.5), 
                axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
  p = p + theme(title = element_text(size = 8))
  p
  fn = paste0(subdir, cancer, "_", gene, "_prm_or_tmt~mRNA.pdf")
  ggsave(file=fn, height=6, width=4, useDingbats=FALSE)
}




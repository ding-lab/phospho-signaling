# Yige Wu @ WashU 2018 Mar
## correlate prospective BRCA PRM protein-level abundance with MS/MS protein abundance from CDAP


# souce -------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/preprocess_files/preprocess_files.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')


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

## input metabolome annotation info
metabolism_group <- read_excel(path = "./cptac2p/cptac_shared/analysis_results/validation/tables/annotate_prm_raw/proteins_in_prm_notkinase_DAVID_functional_annotation.xlsx")
metabolism_group2m <- metabolism_group[, c("ID","panel", "GroupName", "Gene Name")]
colnames(metabolism_group2m) <- c("Gene", "panel", "GroupName", "ProteinNames")


# prm~mRNA ----------------------------------------------------------------
## get intersecting genes and patient IDs
prm_mrna_genes_int <- intersect(prm_PPI_labeled$Gene_Name, fpkm$gene)
prm_mrna_partIDs_int <- intersect(colnames(prm_PPI_labeled), colnames(fpkm))

for (method in c("spearman", "pearson")) {
  cor_stats <- lapply(1:length(prm_mrna_genes_int), FUN = function(n, m1, m2, gs, ps, mt) cor.test(x = as.numeric(as.vector(m1[m1$Gene == gs[n], ps])), 
                                                                                          y = as.numeric(as.vector(m2[m2$gene == gs[n], ps])), 
                                                                                          method = mt),
                      m1 = prm_PPI_labeled, m2 = fpkm.log2, gs = prm_mrna_genes_int, ps = prm_mrna_partIDs_int, mt = method)
  cor_df <- data.frame(Gene = prm_mrna_genes_int,
                       pvalue = sapply(1:length(prm_mrna_genes_int), FUN = function(n, l) l[[n]]$p.value, l = cor_stats),
                       rho = sapply(1:length(prm_mrna_genes_int), FUN = function(n, l) l[[n]]$estimate[[1]], l = cor_stats),
                       num_samples = sapply(1:length(prm_mrna_genes_int), FUN = function(n, m1, m2, gs, ps) length(which(!is.na(m1[m1$Gene == gs[n], ps]) & !is.na(m2[m2$gene == gs[n], ps]))),
                                           m1 = prm_PPI_labeled, m2 = fpkm.log2, gs = prm_mrna_genes_int, ps = prm_mrna_partIDs_int))
  cor_df$fdr <- p.adjust(cor_df$pvalue,method = "fdr")
  cor_df$log_fdr <- -log10(cor_df$fdr)
  ## add in the long protein name to identify the kinases
  protein_names <- unique(prm_raw.f[, c("Protein.Name", "Gene_Name")])
  cor_df <- merge(cor_df, protein_names, by.x = c("Gene"), by.y = c("Gene_Name"), all.x = T)
  ## add in the group name for kinases
  cor_df <- merge(cor_df, kinase_family_per_synonym[, c("synonym", "GroupName", "ProteinNames")], by.x = c("Gene"), by.y = c("synonym"), all.x = T)
  cor_df$is_kinase <- (grepl(pattern = "kinase", x = as.vector(cor_df$ProteinNames)) | !is.na(cor_df$ProteinNames))
  cor_df.k <- cor_df[cor_df$is_kinase,]; cor_df.k$panel <- "Kinome_panel"
  ## add in the group name for metabolic enzymes
  cor_df.m <- cor_df[!cor_df$is_kinase,]; cor_df.m$GroupName <- NULL; cor_df.m$ProteinNames <- NULL
  cor_df.m <- merge(cor_df.m, metabolism_group2m, by = c("Gene"), all.x = T)
  ## merge
  cor_df <- rbind(cor_df.k, cor_df.m[, colnames(cor_df.k)])
  resultDnow <- makeOutDir()
  write.table(cor_df, row.names = F, quote=F, sep = '\t', file=paste(resultDnow, cancer, "_", method,"_correlation_gene_level_prm_mRNA.txt",sep=""))
}


# tmt~mRNA ----------------------------------------------------------------
## get intersecting genes and patient IDs
tmt_mrna_genes_int <- intersect(tmt_PPI_labeled.t$gene, fpkm$gene)
tmt_mrna_partIDs_int <- intersect(colnames(tmt_PPI_labeled.t), colnames(fpkm))

for (method in c("spearman", "pearson")) {
  cor_stats <- lapply(1:length(prm_mrna_genes_int), FUN = function(n, m1, m2, gs, ps, mt) try(cor.test(x = as.numeric(as.vector(m1[m1$gene == gs[n], ps])), 
                                                                                                       y = as.numeric(as.vector(m2[m2$gene == gs[n], ps])), 
                                                                                                       method = mt), silent = T),
                      m1 = tmt_PPI_labeled.t, m2 = fpkm.log2, gs = prm_mrna_genes_int, ps = tmt_mrna_partIDs_int, mt = method)
  cor_df <- data.frame(Gene = prm_mrna_genes_int,
                       pvalue = sapply(1:length(prm_mrna_genes_int), FUN = function(n, l) {
                         if ("p.value" %in% names(l[[n]])) {
                           return(l[[n]]$p.value)
                         } else {
                           return(NA)
                         }
                        }, l = cor_stats),
                       rho = sapply(1:length(prm_mrna_genes_int), FUN = function(n, l) {
                         if ("estimate" %in% names(l[[n]])) {
                           return(l[[n]]$estimate[[1]])
                         } else {
                           return(NA)
                         }
                       }, l = cor_stats),
                       num_samples = sapply(1:length(prm_mrna_genes_int), FUN = function(n, m1, m2, gs, ps) length(which(!is.na(m1[m1$Gene == gs[n], ps]) & !is.na(m2[m2$gene == gs[n], ps]))),
                                            m1 = tmt_PPI_labeled.t, m2 = fpkm.log2, gs = prm_mrna_genes_int, ps = tmt_mrna_partIDs_int))
  cor_df$fdr <- p.adjust(cor_df$pvalue,method = "fdr")
  cor_df$log_fdr <- -log10(cor_df$fdr)
  ## add in the group name for kinases
  cor_df <- merge(cor_df, kinase_family_per_synonym[, c("synonym", "GroupName", "ProteinNames")], by.x = c("Gene"), by.y = c("synonym"), all.x = T)
  cor_df$is_kinase <- (grepl(pattern = "kinase", x = as.vector(cor_df$ProteinNames)) | !is.na(cor_df$ProteinNames))
  cor_df.k <- cor_df[cor_df$is_kinase,]; cor_df.k$panel <- "Kinome_panel"
  ## add in the group name for metabolic enzymes
  cor_df.m <- cor_df[!cor_df$is_kinase,]; cor_df.m$GroupName <- NULL; cor_df.m$ProteinNames <- NULL
  cor_df.m <- merge(cor_df.m, metabolism_group2m, by = c("Gene"), all.x = T)
  ## merge
  cor_df <- rbind(cor_df.k, cor_df.m[, colnames(cor_df.k)])
  resultDnow <- makeOutDir()
  write.table(cor_df, row.names = F, quote=F, sep = '\t', file=paste(resultDnow, cancer, "_", method,"_correlation_gene_level_tmt_mRNA_limittoprmgenes.txt",sep=""))
}

# for (method in c("spearman")) {
#   cor_stats <- lapply(1:length(prm_mrna_genes_int), FUN = function(n, m1, m2, gs, ps, mt) try(cor.test(x = as.numeric(as.vector(m1[m1$gene == gs[n], ps])), 
#                                                                                                        y = as.numeric(as.vector(m2[m2$gene == gs[n], ps])), 
#                                                                                                        method = mt), silent = T),
#                       m1 = tmt_PPI_labeled.t, m2 = fpkm.log2, gs = tmt_mrna_genes_int, ps = tmt_mrna_partIDs_int, mt = method)
#   cor_df <- data.frame(Gene = prm_mrna_genes_int,
#                        pvalue = sapply(1:length(prm_mrna_genes_int), FUN = function(n, l) {
#                          if ("p.value" %in% names(l[[n]])) {
#                            return(l[[n]]$p.value)
#                          } else {
#                            return(NA)
#                          }
#                        }, l = cor_stats),
#                        rho = sapply(1:length(prm_mrna_genes_int), FUN = function(n, l) {
#                          if ("estimate" %in% names(l[[n]])) {
#                            return(l[[n]]$p.value)
#                          } else {
#                            return(NA)
#                          }
#                        }, l = cor_stats))
#   cor_df$fdr <- p.adjust(cor_df$pvalue,method = "fdr")
#   cor_df$log_fdr <- -log10(cor_df$fdr)
#   ## add in the group name for kinases
#   cor_df <- merge(cor_df, kinase_family_per_synonym[, c("synonym", "GroupName", "ProteinNames")], by.x = c("Gene"), by.y = c("synonym"), all.x = T)
#   cor_df$is_kinase <- grepl(pattern = "kinase", x = as.vector(cor_df$ProteinNames))
#   resultDnow <- makeOutDir()
#   write.table(cor_df, row.names = F, quote=F, sep = '\t', file=paste(resultDnow, cancer, "_", method,"_correlation_gene_level_tmt_mRNA.txt",sep=""))
# }
# 
# 



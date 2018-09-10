# Yige Wu @ WashU 2018 Mar
## correlate prospective BRCA PRM protein-level abundance with MS/MS protein abundance from CDAP
## match parent specimen labels
## TODO: plot per-gene/per-sample density curves before correlation for both MS and PRM data

# souce -------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/preprocess_files/preprocess_files.R')
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
resultDnow <- makeOutDir()
write.table(prm_parent_labeled, row.names = F, quote=F, sep = '\t', file=paste(resultDnow, cancer, "_prm_trplicatemean_median_parent_labeled.txt",sep=""))
write.table(ms_parent_labeled, row.names = F, quote=F, sep = '\t', file=paste(resultDnow, cancer, "_tmt_match_prm_parent_labeled.txt",sep=""))

## write table of correlation
## input PRM peptide-level values
prm_raw.f <- fread(input = "./cptac2p/cptac_shared/analysis_results/preprocess_files/format_RPM_peptide/prm_genename_annotated_specimen_labeled.txt", data.table = F)
## input kinase family info, use GroupName
kinase_family <- read_delim("~/Box Sync/pan3can_shared_data/Gene_family/knome_render_peerj-01-126-s001_wgene.txt","\t", escape_double = FALSE, trim_ws = TRUE)
kinase_family <- data.frame(kinase_family)
### reformat kinase family info to have one Synonyms per line
synonyms <- lapply(1:nrow(kinase_family), FUN = function(n, m) unlist(str_split(string = as.vector(m$Synonyms[n]), pattern = ",")), m = kinase_family)
num_synonyms <- sapply(synonyms, FUN = function(x) length(x))
kinase_family_per_synonym <- kinase_family[rep(1:nrow(kinase_family), times = num_synonyms),]
kinase_family_per_synonym$synonym <- unlist(synonyms)
# write.table(kinase_family_per_synonym, row.names = F, quote=F, sep = '\t', file=paste(makeOutDir(), "kinase_family_per_synonym,txt",sep=""))

## input metabolome annotation info
metabolism_group <- read_excel(path = "./cptac2p/cptac_shared/analysis_results/validation/tables/annotate_prm_raw/proteins_in_prm_notkinase_DAVID_functional_annotation.xlsx")
metabolism_group2m <- metabolism_group[, c("ID","panel", "GroupName", "Gene Name")]
colnames(metabolism_group2m) <- c("Gene", "panel", "GroupName", "ProteinNames")


for (method in c("pearson", "spearman")) {
  cor_stats <- lapply(1:length(genes_int), FUN = function(n, m1, m2, gs, ps, mt) cor.test(x = as.numeric(as.vector(m1[m1$Gene == gs[n], ps])), 
                                                                                          y = as.numeric(as.vector(m2[m2$Gene == gs[n], ps])), 
                                                                                          method = mt),
                      m1 = prm_parent_labeled, m2 = ms_parent_labeled, gs = genes_int, ps = parent_sampIDs_int, mt = method)
  cor_df <- data.frame(Gene = genes_int,
                       pvalue = sapply(1:length(genes_int), FUN = function(n, l) l[[n]]$p.value, l = cor_stats),
                       rho = sapply(1:length(genes_int), FUN = function(n, l) l[[n]]$estimate[[1]], l = cor_stats),
                       num_samples = sapply(1:length(genes_int), FUN = function(n, m1, m2, gs, ps) length(which(!is.na(m1[m1$Gene == gs[n], ps]) & !is.na(m2[m2$Gene == gs[n], ps]))),
                                            m1 = prm_parent_labeled, m2 = ms_parent_labeled, gs = genes_int, ps = parent_sampIDs_int))
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
  write.table(cor_df, row.names = F, quote=F, sep = '\t', file=paste(resultDnow, cancer, "_", method,"_correlation_prm_ms_match_parent_specimen.txt",sep=""))
}



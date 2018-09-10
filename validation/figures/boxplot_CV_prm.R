# Yige Wu @ WashU 2018 Mar
## showing the coefficient of variantion for all the proteins detected by PRM

# souce -------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/preprocess_files/preprocess_files.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_plotting.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')


# inputs ------------------------------------------------------------------
cancer <- "BRCA"
prm_raw.f <- fread(input = paste0("cptac2p/cptac_shared/analysis_results/preprocess_files/format_RPM_peptide/prm_CV_genename_annotated_specimen_labeled.txt"),
                         data.table = F)
# prm_raw.f.cv <- prm_raw.f[,c(colnames(prm_raw.f)[grepl(pattern = "CV", x = colnames(prm_raw.f))])]
prm_raw.f.cv <- prm_raw.f[,!(colnames(prm_raw.f) %in% c("Protein.Name","Peptide.Modified.Sequence", "X.", "Gene_Name"))]
cv_nonzeroorna <- rowSums(!is.na(prm_raw.f.cv) & (prm_raw.f.cv != 0 ))
cv_median <- sapply(1:nrow(prm_raw.f), FUN = function(n, m) median(as.numeric(as.vector(m[n,])), na.rm = T), m = prm_raw.f.cv)
prm_raw.f$num_nonzeroorna <- cv_nonzeroorna
prm_raw.f$cv_median <- cv_median
prm_raw.m <- melt(cbind(prm_raw.f[,c("Protein.Name","Peptide.Modified.Sequence","X.","Gene_Name", "num_nonzeroorna", "cv_median")], prm_raw.f.cv), 
                  id.vars = c("Protein.Name","Peptide.Modified.Sequence","X.","Gene_Name", "num_nonzeroorna", "cv_median"))

prm_raw.m$Peptide.Modified.Sequence <- factor(prm_raw.m$Peptide.Modified.Sequence, levels = as.vector(prm_raw.m$Peptide.Modified.Sequence)[order(cv_median)])
num_peptide <- data.frame(table(prm_raw.f[,c("Gene_Name")]))
prm_raw.f <- merge(prm_raw.f, num_peptide, by.x = c("Gene_Name"), by.y = c("Var1"), all.x = T)

# statistics for CV -------------------------------------------------------
length(unique(prm_raw.f$Gene_Name[prm_raw.f$cv_median > 0.2]))
nrow(prm_raw.f[prm_raw.f$cv_median > 0.2,])
length(unique(prm_raw.f$Gene_Name[prm_raw.f$cv_median > 0.2 & prm_raw.f$Freq == 1]))

nrow(prm_raw.f[prm_raw.f$cv_median > 0.2 & prm_raw.f$Freq == 1,])
length(unique(prm_raw.f$Gene_Name[prm_raw.f$cv_median < 0.2]))
nrow(prm_raw.f[prm_raw.f$cv_median < 0.2,])
length(unique(prm_raw.f$Gene_Name[prm_raw.f$cv_median < 0.2 & prm_raw.f$Freq == 1]))
nrow(prm_raw.f[prm_raw.f$cv_median < 0.2 & prm_raw.f$Freq == 1,])
length(unique(prm_raw.f$Gene_Name[prm_raw.f$cv_median > 0.3]))
nrow(prm_raw.f[prm_raw.f$cv_median > 0.3,])


# single peptide proteins --------------------------------------------------------------------
df <- prm_raw.m[prm_raw.m$Gene_Name %in% num_peptide$Var1[num_peptide$Freq==1] & prm_raw.m$num_nonzeroorna > 0,]
## annotate data issues
df$CV_issue <- ifelse(df$cv_median > 0.2, "median CV > 20%", "normal")
df$sample_size_issue <- ifelse(df$num_nonzeroorna < 5, "non-zero samples < 5", "normal")
df$value <- as.numeric(as.vector(df$value))

p <- ggplot(data = df)
p <- p + geom_hline(yintercept = 0.2, color = "grey")
p <- p + geom_boxplot(mapping = aes(x = Peptide.Modified.Sequence, y = value, color = Gene_Name, group = Gene_Name), 
                      width=0.2, alpha = 0.8, size = 0.5, 
                      outlier.size = 0.8, outlier.stroke = 0, outlier.alpha = 0.5)
p <- p + scale_x_discrete(breaks=as.vector(df$Peptide.Modified.Sequence),
                          labels=as.vector(df$Gene_Name))
p <- p + geom_text(mapping = aes(x = Peptide.Modified.Sequence, y = -0.1, label = num_nonzeroorna), size = 2)
p <- p + ylab("coefficient of variation (CV, the ratio of standard deviation to mean)")
p <- p + ggtitle(label = paste0("PRM measurement reproducibility for ", length(unique(df$Gene_Name)), " proteins with only one target peptide"))
p <- p + facet_grid(.~CV_issue, scales = "free_x", space = "free_x")
p <- p + theme_nogrid()
p <- p + theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 5, angle = 45), 
               legend.position = "none", plot.title = element_text(size=10, face = "bold"))
p
fn = paste0(makeOutDir(),'prm_single_peptide_per_protein.pdf')
ggsave(file=fn, height=5, width=20, useDingbats=FALSE)

# multiple peptide proteins --------------------------------------------------------------------
df <- prm_raw.m[prm_raw.m$Gene_Name %in% num_peptide$Var1[num_peptide$Freq>1] & prm_raw.m$num_nonzeroorna > 0,]
## annotate data issues
df$CV_issue <- ifelse(df$cv_median > 0.2, "median CV > 20%", "normal")
df$sample_size_issue <- ifelse(df$num_nonzeroorna < 5, "non-zero samples < 5", "normal")

df_multi <- df
df_multi$Peptide.Modified.Sequence <- reorder(df_multi$Peptide.Modified.Sequence, df_multi$Gene_Name)
df_multi$value <- as.numeric(as.vector(df_multi$value))
p <- ggplot(data = df_multi)
p <- p + geom_hline(yintercept = 0.2, color = "grey", linetype = 2)
p <- p + geom_boxplot(mapping = aes(x = Peptide.Modified.Sequence, y = value, color = Gene_Name), 
                      width=0.2, alpha = 0.8, size = 0.5, 
                      outlier.size = 0.8, outlier.stroke = 0, outlier.alpha = 0.5)
p <- p + scale_x_discrete(breaks=as.vector(df$Peptide.Modified.Sequence),
                          labels=as.vector(df$Gene_Name))
p <- p + geom_text(mapping = aes(x = Peptide.Modified.Sequence, y = -0.1, label = num_nonzeroorna), size = 1)
p <- p + ylab("coefficient of variation (CV, the ratio of standard deviation to mean)")
p <- p + ggtitle(label = paste0("PRM measurement reproducibility for peptides of ", length(unique(df_multi$Gene_Name))," protein with multiple peptides measured"))
p <- p + facet_grid(.~CV_issue, scales = "free_x", space = "free_x")
p <- p + theme_nogrid()
p <- p + theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 5, angle = 45), 
               legend.position = "none", plot.title = element_text(size=10, face = "bold"))
p
fn = paste0(makeOutDir(), 'prm_multi_peptide_per_protein','.pdf')
ggsave(file=fn, height=5, width=25, useDingbats=FALSE)

## per gene per plot
dir.create(paste0(makeOutDir(), "multi_peptide_proteins/"))
for (gene in unique(df$Gene_Name)) {
  df_gene <- df[df$Gene_Name == gene,]
  p <- ggplot(data = df_gene)
  p <- p + geom_hline(yintercept = 0.2, color = "grey")
  p <- p + geom_boxplot(mapping = aes(x = Peptide.Modified.Sequence, y = value, color = Gene_Name), 
                        width=0.2, alpha = 0.8, size = 0.5, 
                        outlier.size = 0.8, outlier.stroke = 0, outlier.alpha = 0.5)
  p <- p + scale_x_discrete(breaks=as.vector(df$Peptide.Modified.Sequence),
                            labels=as.vector(df$Gene_Name))
  p <- p + geom_text(mapping = aes(x = Peptide.Modified.Sequence, y = -0.1, label = num_nonzeroorna), size = 2)
  p <- p + geom_text(mapping = aes(x = Peptide.Modified.Sequence, y = max(df_gene$value, na.rm = T)*0.8, label = Peptide.Modified.Sequence), 
                     size = 2, angle = -90, color = "grey", alpha = 0.5, nudge_x = -0.1)
  p <- p + ylab("coefficient of variation (CV, the ratio of standard deviation to mean)")
  p <- p + ggtitle(label = paste0("PRM measurement reproducibility for peptides of ", gene, " protein"))
  p <- p + facet_grid(.~CV_issue+sample_size_issue, scales = "free_x", space = "free_x")
  p <- p + theme_nogrid()
  p <- p + theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 5, angle = 45), 
                 legend.position = "none", plot.title = element_text(size=10))
  p
  fn = paste0(makeOutDir(), "multi_peptide_proteins/", 'prm_multi_peptide_per_protein_', gene,'.pdf')
  ggsave(file=fn, height=5, width=5, useDingbats=FALSE)
}








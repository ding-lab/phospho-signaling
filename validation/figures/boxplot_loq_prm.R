# Yige Wu @ WashU 2018 Mar
## showing the limit of quantification for all the proteins detected by PRM

# souce -------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/preprocess_files/preprocess_files.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_plotting.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
library(ggrepel)
# inputs ------------------------------------------------------------------
cancer <- "BRCA"
prm_raw.f <- fread(input = paste0("cptac2p/cptac_shared/analysis_results/validation/tables/annotate_prm_raw/BRCA_prm_genename_annotated_specimen_lod_loq_labeled.txt"),
                   data.table = F)
tmp <- as.matrix(prm_raw.f[, colnames(prm_raw.f)[grepl(pattern = "Average", x = colnames(prm_raw.f))]])
tmp[tmp==0] <- NA
prm_raw.f$Average_min <- rowMins(tmp, na.rm = T)

# box plot LOQ ------------------------------------------------------------
prm_raw.f$GroupName[prm_raw.f$panel == "Kinome_panel" & is.na(prm_raw.f$GroupName)] <- "Other"

for (panel in unique(prm_raw.f$panel)) {
  df <- prm_raw.f[prm_raw.f$panel == panel,]
  df <- df[!is.na(df$CPTAC.LLOQ),]
  df$loq_remove_outlier <- df$CPTAC.LLOQ
  iqr_tmp <- IQR(df$CPTAC.LLOQ, na.rm = T)
  outliers <- data.frame(GroupName = unique(df$GroupName), 
                         med = sapply(unique(df$GroupName), FUN = function(g, m) quantile(m$CPTAC.LLOQ[m$GroupName == g], probs = 0.5, na.rm = T)[1], m = df),
                         upper = sapply(unique(df$GroupName), FUN = function(g, m) quantile(m$CPTAC.LLOQ[m$GroupName == g], probs = 0.75, na.rm = T)[1] + 1.5*(IQR(m$CPTAC.LLOQ[m$GroupName == g], na.rm = T)), m = df),
                         lower = sapply(unique(df$GroupName), FUN = function(g, m) quantile(m$CPTAC.LLOQ[m$GroupName == g], probs = 0.25, na.rm = T)[1] - 1.5*(IQR(m$CPTAC.LLOQ[m$GroupName == g], na.rm = T)), m = df))
  df <- merge(df, outliers, by = c("GroupName"), all.x = T)
  df$loq_remove_outlier[df$loq_remove_outlier <= df$upper & df$loq_remove_outlier >= df$lower] <- NA
  df$GroupName <- factor(df$GroupName, levels = as.vector(df$GroupName)[order(as.vector(df$med))])
  
  p <- ggplot(data = df)
  p <- p + geom_boxplot(mapping = aes(x = GroupName, y = CPTAC.LLOQ), 
                        width=0.2, alpha = 0.8, size = 0.5, 
                        outlier.size = 0.8, outlier.stroke = 0, outlier.alpha = 0.5)
  p = p + geom_text_repel(aes(GroupName, loq_remove_outlier, label= Gene_Name, color = GroupName, segment.color = GroupName),
                          force = 1, 
                          segment.size = 0.5, segment.alpha = 0.2, 
                          size=1.5, alpha=0.8)
  # p <- p + geom_text(mapping = aes(x = GroupName, y = loq_remove_outlier, label = Gene_Name, color = GroupName), size = 1, alpha = 0.5 , position = position_jitter()) # position = position_jitter()
  p <- p + theme_nogrid()
  p <- p + ylab("LOQ (fmol/ul)") + xlab("Gene Group")
  p <- p + theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 5, angle = 30, hjust = 1),  
                 legend.position = "none", plot.title = element_text(size=10, face = "bold"))
  p
  fn = paste0(makeOutDir(),'prm_all_LOQ_per_peptide', panel,'.pdf')
  ggsave(file=fn, height=5, width=5, useDingbats=FALSE)
  
  p_zoom <- p + coord_cartesian(ylim = c(0, quantile(x = df$CPTAC.LLOQ, probs = 0.98, na.rm = T)))
  p_zoom
  fn = paste0(makeOutDir(),'prm_all_LOQ_per_peptide', panel,'.qt98zoom.pdf')
  ggsave(file=fn, height=5, width=5, useDingbats=FALSE)
}

for (panel in unique(prm_raw.f$panel)) {
  df <- prm_raw.f[prm_raw.f$panel == panel & !is.na(prm_raw.f$CPTAC.LLOQ) & !is.na(prm_raw.f$Average_min),]
  df$Average_min_rm_outlier <- df$Average_min
  iqr_tmp <- IQR(df$Average_min, na.rm = T)
  upper <- quantile(df$Average_min, probs = 0.5, na.rm = T) + 1.5*iqr_tmp
  lower <- quantile(df$Average_min, probs = 0.5, na.rm = T) - 1.5*iqr_tmp
  df$Average_min_rm_outlier[df$Average_min_rm_outlier >= upper | df$Average_min_rm_outlier <= lower] <- NA
  
  p <- ggplot(data = df)
  p <- p + geom_point(mapping = aes(x = CPTAC.LLOQ, y = Average_min_rm_outlier, color = GroupName), size = 1, alpha = 0.5, stroke = 0)
  p = p + geom_text_repel(data = df[df$Average_min <= df$CPTAC.LLOQ,], mapping = aes(x = CPTAC.LLOQ, y = Average_min_rm_outlier, label= Gene_Name, color = GroupName, segment.color = GroupName),
                          force = 1,
                          segment.size = 0.5, segment.alpha = 0.2,
                          size=1.5, alpha=0.8)
  p = p + geom_abline(slope = 1, intercept = 0, linetype = 2, alpha = 0.8)
  p <- p + theme_nogrid()
  p <- p + ylim(c(0, 6)) + xlim(c(0, 6))
  p <- p + xlab("LOQ (limit of quantification, fmol/ul)") + ylab("lowest non-zero PRM-detected peptide abundance (fmol/ul)")
  p <- p + theme(axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),  axis.text.x = element_text(size = 5),
                 legend.key.size = unit(1, "cm"), legend.text = element_text(size = 5), 
                 legend.background = element_rect(fill = NA), legend.box.background = element_rect(fill = NA),
                 legend.key = element_rect(fill = NA),
                 legend.title = element_text(size = 5), legend.justification = c(1, 1), legend.position = c(1, 1),
                 plot.title = element_text(size=10, face = "bold"))
  p
  fn = paste0(makeOutDir(),'prm_all_LOQ_vs_lowest_per_peptide', panel,'.pdf')
  ggsave(file=fn, height=5, width=5)
  
}

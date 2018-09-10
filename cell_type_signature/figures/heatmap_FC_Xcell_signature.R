# Yige Wu @ WashU 2018 Apr
# take cell type specific mark genes from Xcell and test difference between prospective tumors and normals


# souce ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/cptac2p_analysis_shared.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/preprocess_files/get_tumor_normal_pair.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_plotting.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')

# inputs ------------------------------------------------------------------
pvalue_fisher <- fread(input = paste0(inputD, "analysis_results/cell_type_signature/tables/diffexp_Xcell_signature/pvalue_fisher.txt"), data.table = F)
xcell_wilcox <- readRDS(file = paste0(inputD, "analysis_results/cell_type_signature/tables/diffexp_Xcell_signature/xcell_wilcox.RDS"))

cancers <- c("BRCA", "OV", "CO")

# heatmap -----------------------------------------------------------------
for (cell_type in pvalue_fisher$cell_type) {
  all_genes <- names(xcell_wilcox[["OV"]][[cell_type]][["fdrs_ordered"]])
  
  pro <- NULL
  for (gene in all_genes[1:min(10, length(all_genes))]) {
    ## input protein abundance
    pro_g <- NULL
    for (cancer in cancers) {
      pro_t <- fread(paste0(inputD, cancer, "/", prefix[cancer], "_", "PRO", "_formatted_normalized_noControl.txt"),
                     data.table = F)
      pro_n <- fread(paste0(inputD, cancer, "/", prefix[cancer], "_", "PRO", "_formatted_normalized_Control.txt"),
                     data.table = F)
      
      if (gene %in% pro_t$Gene && gene %in% pro_n$Gene) {
        tn_match <- tumor_normal_match[[cancer]]
        samples_n <- as.vector(tn_match$normalSampIDs[!is.na(tn_match$normalSampIDs)])
        samples_t <- as.vector(tn_match$tumorSampIDs[!is.na(tn_match$normalSampIDs)])
        
        pro_t_g <- pro_t[pro_t$Gene == gene, c("Gene", samples_t)]
        pro_t_g.m <- melt(pro_t_g)
        colnames(pro_t_g.m)[ncol(pro_t_g.m)] <- "pro_t"
        
        
        pro_n_g <- pro_n[pro_n$Gene == gene, c("Gene", samples_n)]
        pro_n_g.m <- melt(pro_n_g)
        pro_t_g.m$pro_n <- pro_n_g.m$value
        pro_t_g.m$cancer <- cancer
        
        pro_t_g.m$fdr <- xcell_wilcox[[cancer]][[cell_type]][["fdrs_ordered"]][gene]
        pro_t_g.m$log2FC <- xcell_wilcox[[cancer]][[cell_type]][["log2fc_ordered"]][gene]
        pro_t_g.m$x <- min(pro_t_g.m$pro_n, na.rm = T)
        pro_t_g.m$y <- max(pro_t_g.m$pro_t, na.rm = T)
        
        pro_g <- rbind(pro_g, pro_t_g.m)
      }
    }
    pro <- rbind(pro, pro_g)
  }
  pro$Gene <- factor(pro$Gene, levels = all_genes[1:min(10, length(all_genes))])
  pro_sig <- pro[!is.na(pro$fdr) & pro$fdr < 0.05, c("x", "y", "fdr", "cancer", "Gene", "log2FC")]
  pro_sig <- unique(pro_sig)
  
  p <- ggplot(data = pro)
  p <- p + geom_point(mapping = aes(x = pro_n, y = pro_t, color = cancer), alpha = 0.8)
  p = p + geom_abline(slope = 1, intercept = 0, linetype = 2, alpha = 0.8)
  p <- p + facet_wrap(~Gene, nrow = 2, scales = "free")
  texts <- paste0("FDR=", signif(pro_sig$fdr, digits = 1), ",log2FC=", signif(pro_sig$log2FC, digits = 2))
  dat_text <- data.frame(label = texts, 
                         Gene  = as.vector(pro_sig$Gene),
                         cancer = as.vector(pro_sig$cancer),
                         x = pro_sig$x,
                         y = pro_sig$y)
  p <- p + geom_text(data    = dat_text, mapping = aes(x = 0.3*x, y = 0.9*y, label = label, color = cancer),
                     hjust   = -0.1, vjust   = -1, size = 4)
  p <- p + ggtitle(label = paste0(cell_type, " signature gene sets"),
                   subtitle = paste0("top significantly differentially expressed proteins"))
  p = p + theme_nogrid()
  p = p + labs(x = paste0("protein abundance of normal tissue(log2)"), 
               y=paste0("protein abundance of tumor tissue(log2)"))
  p
  resultDnow <- makeOutDir()
  fn = paste0(resultDnow, cell_type, "_signature.pdf")
  ggsave(file=fn, height=10, width=20, useDingbats=FALSE)
}


# Yige Wu @ WashU March 2018
# plot boxplot plot showing differential protein mRNA/protein/phosphoprotein per gene for MMR

# source ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/preprocess_files/preprocess_files.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_plotting.R')
library(ggpubr)

# input -------------------------------------------------------------------
## input MMR genes
mmr_gene <- read_delim("~/Box Sync/MSI_CPTAC/Data/mmr_gene.txt","\t", escape_double = FALSE, col_names = FALSE,trim_ws = TRUE)
mmr_gene <- mmr_gene$X1

## input MSI scores
msi_score <- read_delim("~/Box Sync/MSI_CPTAC/Data/MSIsensor_Score_qgao/CPTAC.MSI.score.tsv",
                        "\t", escape_double = FALSE, trim_ws = TRUE)
## input clinical info
clinical <- read.delim(paste0(inputD, "Specimen_Data_20161005_Yige_20180307.txt"), sep = "\t")
clinical <- data.frame(clinical)

# boxplots ------------------------------------------------------------
for (cancer in c("UCEC", "CO")) {
  msi_exps_list <- readRDS(file = paste0(inputD, "analysis_results/expression_matrices/tables/intergrate_msi_expression/mmr_gene_msi_exps_list.RDS"))
  
  for (msi_cutoff in c(0.8, 0.9, 1, 3.5)) {
    ## input mRNA abundance
    df_mRNA <- msi_exps_list[[cancer]][["mRNA"]][, c("Gene", "log2RPKM", "Score")]
    df_mRNA$MSI_status <- ifelse(df_mRNA$Score >= msi_cutoff, "MSI-H_tumors", "other_tumors")
    df_mRNA$expression_type <- "mRNA"
    colnames(df_mRNA) <- c("Gene", "score", "MSI_score", "MSI_status", "expression_type")
    ## input protein abundance
    df_protein <- msi_exps_list[[cancer]][["protein"]][, c("Gene", "score", "Score")]
    df_protein$MSI_status <- ifelse(df_protein$Score >= msi_cutoff, "MSI-H_tumors", "other_tumors")
    df_protein$expression_type <- "protein"
    colnames(df_protein) <- c("Gene", "score", "MSI_score", "MSI_status", "expression_type")
    ## input phosphoprotein abundance
    df_phosphoprotein <- msi_exps_list[[cancer]][["phosphoprotein"]][, c("Gene", "score", "Score")]
    df_phosphoprotein$MSI_status <- ifelse(df_phosphoprotein$Score >= msi_cutoff, "MSI-H_tumors", "other_tumors")
    df_phosphoprotein$expression_type <- "phosphoprotein"
    colnames(df_phosphoprotein) <- c("Gene", "score", "MSI_score", "MSI_status", "expression_type")
    df <- rbind(df_mRNA, df_protein, df_phosphoprotein)
    df$expression_type <- factor(df$expression_type, levels = c("mRNA", "protein", "phosphoprotein"))
    
    p <- ggboxplot(df, x = "Gene", y = "score",
                   color = "MSI_status", palette = "jco",
                   add = "jitter", xlab = "", ylab = "gene/protein abundance (log2 transformed)", font.label = list(size = 5, style = "plain"))
    p + stat_compare_means(aes(group = MSI_status), label = "p.signif")
    facet(p = p, facet.by = "expression_type", nrow = 3, scales = "free_y" )
    resultDnow <- makeOutDir()
    fn = paste0(resultDnow, cancer, '_MMR_mRNA&protein&phospho_abundance_between_MSI-H_cutoff', msi_cutoff, '_and_rest_in_', cancer , ".pdf")
    ggsave(file=fn, height=8, width=24, useDingbats=FALSE)
  }
}




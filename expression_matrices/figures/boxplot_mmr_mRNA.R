# Yige Wu @ WashU March 2018
# plot Violin plot showing differential mRNA expression per gene for MMR

# source ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/preprocess_files/preprocess_files.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_plotting.R')

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

## initiate list to record Mann-Whitney test results
mann_results <- vector("list")

# COAD --------------------------------------------------------------------
## input expression score(normalized log ratios) and take CO samples and intent gene set
cancer <- "CO"
exps <- readRDS(file = paste0(inputD, "analysis_results/expression_matrices/tables/get_mRNA_quantile/cptac2p_mRNA_score.RDS"))
exps <- exps[[cancer]]
exps <- exps[exps$Gene %in% mmr_gene,]

## transform sample IDs to patient IDs
exps$Participant.ID <- sampID2partID(sampleID_vector = exps$Specimen.ID, sample_map = clinical)
library(ggpubr)
## merge exps with MSI score
msi_exps <- merge(msi_score, exps, by.x = c("Sample"), by.y = c("Participant.ID"))
msi_exps$MSI_status <- ifelse(msi_exps$Score >= 3.5, paste0("MSI-H_tumors"), paste0("other_tumors"))
msi_exps$log2RPKM <- log2(msi_exps$score)
resultDnow <- makeOutDir()
dir.create(path = paste0(resultDnow, cancer))
## MLH1 not in colorectal proteomics
mann_results[[cancer]] <- vector("list")
for (gene in mmr_gene) {
  exps_gene <- msi_exps[msi_exps$Gene == gene & !is.na(msi_exps$log2RPKM),]
  exps_msih <- exps_gene$log2RPKM[exps_gene$MSI_status == "MSI-H_tumors"]
  exps_msil <- exps_gene$log2RPKM[exps_gene$MSI_status == "other_tumors"]
  
  if (length(exps_msih) > 4) {
    stat <- wilcox.test(x = exps_msih, y = exps_msil)
    mann_results[[cancer]][[gene]] <- list(p.value = stat$p.value)
    p <- ggplot(data = exps_gene)
    p <- p + geom_violin(mapping = aes(x = MSI_status, y = log2RPKM, fill = MSI_status), alpha = 0.5, color = NA)
    p <- p + scale_fill_manual(values = c("MSI-H_tumors" = set1[1], "other_tumors" = "grey"))
    p <- p + geom_boxplot(mapping = aes(x = MSI_status, y = log2RPKM), width=0.1, alpha = 0.8)
    p <- p + geom_point(mapping = aes(x = MSI_status, y = log2RPKM), alpha = 0.5, stroke = 0, color = "black", position = position_jitter())
    p <- p + ylab("mRNA abundance (log2RPKM)")
    p <- p + annotate("text", x = 1, y = 2 , label = paste0("Mann-Whitney, p = ", signif(stat$p.value, digits = 2)))
    p <- p + ggtitle(label = paste0(gene, " mRNA abundance between MSI-H tumor(N=", length(exps_msih),") and other tumors(N=", length(exps_msil), ")"))
    p <- p + theme_nogrid()
    p <- p + theme(axis.title.x = element_blank()) + theme(legend.position = "none") + theme(plot.title = element_text(size=8))
    p
    fn = paste0(resultDnow, cancer, "/", gene, '_mRNA_abundance_between_MSI-H_and_rest_in_', cancer, ".pdf")
    ggsave(file=fn, height=4, width=5, useDingbats=FALSE)
  } else {
    print(paste0(gene, " not enough data!"))
  }
}

# tmp <- unlist(mann_results[[cancer]])
# genes_sort <- as.vector(str_split_fixed(string = names(tmp[order(tmp)]), pattern = "\\.", 2)[,1])
# msi_exps$Gene <- factor(msi_exps$Gene, levels = genes_sort)
p <- ggboxplot(msi_exps, x = "Gene", y = "log2RPKM",
               color = "MSI_status", palette = "jco",
               add = "jitter", xlab = "", ylab = "mRNA abundance (log2RPKM)", font.label = list(size = 5, style = "plain") )
p + stat_compare_means(aes(group = MSI_status), label = "p.signif")
resultDnow <- makeOutDir()
fn = paste0(resultDnow, cancer, '_MMR_gene_mRNA_abundance_between_MSI-H_and_rest_in_', cancer, ".pdf")
ggsave(file=fn, height=5, width=24, useDingbats=FALSE)
msi_exps_list[[cancer]][["mRNA"]] <- msi_exps

# UCEC --------------------------------------------------------------------
cancer <- "UCEC"
exps <- readRDS(file = paste0(inputD, "analysis_results/expression_matrices/tables/get_mRNA_quantile/cptac3_mRNA_score.RDS"))
exps <- exps[[cancer]]
exps <- exps[exps$Gene %in% mmr_gene,]

## transform sample IDs to patient IDs
# exps$Participant.ID <- sampID2partID(sampleID_vector = exps$Specimen.ID, sample_map = clinical)
library(ggpubr)
## merge exps with MSI score
msi_exps <- merge(msi_score, exps, by.x = c("Sample"), by.y = c("Participant.ID"))
msi_exps$MSI_status <- ifelse(msi_exps$Score >= 3.5, paste0("MSI-H_tumors"), paste0("other_tumors"))
msi_exps$log2RPKM <- log2(msi_exps$score)
resultDnow <- makeOutDir()
dir.create(path = paste0(resultDnow, cancer))
mann_results[[cancer]] <- vector("list")
for (gene in mmr_gene) {
  exps_gene <- msi_exps[msi_exps$Gene == gene & !is.na(msi_exps$log2RPKM),]
  exps_msih <- exps_gene$log2RPKM[exps_gene$MSI_status == "MSI-H_tumors"]
  exps_msil <- exps_gene$log2RPKM[exps_gene$MSI_status == "other_tumors"]
  
  if (length(exps_msih) > 4) {
    stat <- wilcox.test(x = exps_msih, y = exps_msil)
    mann_results[[cancer]][[gene]] <- list(p.value = stat$p.value)
    p <- ggplot(data = exps_gene)
    p <- p + geom_violin(mapping = aes(x = MSI_status, y = log2RPKM, fill = MSI_status), alpha = 0.5, color = NA)
    p <- p + scale_fill_manual(values = c("MSI-H_tumors" = set1[1], "other_tumors" = "grey"))
    p <- p + geom_boxplot(mapping = aes(x = MSI_status, y = log2RPKM), width=0.1, alpha = 0.8)
    p <- p + geom_point(mapping = aes(x = MSI_status, y = log2RPKM), alpha = 0.5, stroke = 0, color = "black", position = position_jitter())
    p <- p + ylab("mRNA abundance (log2RPKM)")
    p <- p + annotate("text", x = 1, y = 2 , label = paste0("Mann-Whitney, p = ", signif(stat$p.value, digits = 2)))
    p <- p + ggtitle(label = paste0(gene, " mRNA abundance between MSI-H tumor(N=", length(exps_msih),") and other tumors(N=", length(exps_msil), ")"))
    p <- p + theme_nogrid()
    p <- p + theme(axis.title.x = element_blank()) + theme(legend.position = "none") + theme(plot.title = element_text(size=8))
    p
    
    fn = paste0(resultDnow, cancer, "/", gene, '_mRNA_abundance_between_MSI-H_and_rest_in_', cancer, ".pdf")
    ggsave(file=fn, height=4, width=5, useDingbats=FALSE)
  } else {
    print(paste0(gene, " not enough data!"))
  }
}

# tmp <- unlist(mann_results[[cancer]])
# genes_sort <- as.vector(str_split_fixed(string = names(tmp[order(tmp)]), pattern = "\\.", 2)[,1])
# msi_exps$Gene <- factor(msi_exps$Gene, levels = genes_sort)
msi_exps$log2RPKM_capped <- msi_exps$log2RPKM
msi_exps$log2RPKM_capped[msi_exps$log2RPKM < -5] <- -5
p <- ggboxplot(msi_exps, x = "Gene", y = "log2RPKM_capped",
               color = "MSI_status", palette = "jco",
               add = "jitter", xlab = "", ylab = "mRNA abundance (log2RPKM)", font.label = list(size = 5, style = "plain") )
p + stat_compare_means(aes(group = MSI_status), label = "p.signif")
resultDnow <- makeOutDir()
fn = paste0(resultDnow, cancer, '_MMR_gene_mRNA_abundance_between_MSI-H_and_rest_in_', cancer, ".pdf")
ggsave(file=fn, height=5, width=24, useDingbats=FALSE)

msi_exps_list[[cancer]][["mRNA"]] <- msi_exps


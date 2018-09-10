# Yige Wu @ WashU March 2018
# plot Violin plot showing differential protein expression per gene for MMR

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

## initiate list to record Mann-Whitney test results
mann_results <- vector("list")

# COAD & UCEC MMR genes--------------------------------------------------------------------
for (msi_thres in c(0.8, 3.5)) {
  for (cancer in c("CO", "UCEC")) {
    for (datatype in c("mRNA")) {
      exps <- readRDS(file = paste0(inputD, "analysis_results/expression_matrices/tables/intergrate_msi_expression/mmr_gene_msi_exps_list.RDS"))
      exps <- exps[[cancer]][["mRNA"]][, c("Sample", "Gene", "Score", "score")]; exps$expression_type <- "mRNA"
      colnames(exps) <- c("Sample", "Gene", "Score", "FPKM", "expression_type")
      exps$score <- log2(exps$FPKM + 0.001)
      exps <- exps[exps$Gene %in% mmr_gene,]
      exps$MSI_status <- ifelse(exps$Score >= msi_thres, paste0("MSI-H_tumors"), paste0("other_tumors"))
      mann_results[[cancer]] <- vector("list")
      dir.create(paste0(makeOutDir(), datatype))
      dir.create(paste0(makeOutDir(), datatype, "/", cancer))
      dir.create(paste0(makeOutDir(), datatype, "/", cancer, "/msi_thres", msi_thres))
      
      for (gene in mmr_gene) {
        exps_gene <- exps[exps$Gene == gene & !is.na(exps$score),]
        exps_msih <- exps_gene$score[exps_gene$MSI_status == "MSI-H_tumors"]
        exps_msil <- exps_gene$score[exps_gene$MSI_status == "other_tumors"]
        
        if (length(exps_msih) > 4) {
          stat <- wilcox.test(x = exps_msih, y = exps_msil)
          mann_results[[cancer]][[gene]] <- list(p.value = stat$p.value)
          p <- ggplot(data = exps_gene)
          p <- p + geom_violin(mapping = aes(x = MSI_status, y = score, fill = MSI_status), alpha = 0.5, color = NA)
          p <- p + scale_fill_manual(values = c("MSI-H_tumors" = set1[1], "other_tumors" = "grey"))
          p <- p + geom_boxplot(mapping = aes(x = MSI_status, y = score), width=0.1, alpha = 0.8)
          p <- p + geom_point(mapping = aes(x = MSI_status, y = score, color = Score), alpha = 0.5, stroke = 0, position = position_jitter()) # , color = "black"
          p <- p + ylab("mRNA abundance (log2(FPKM+0.001))")
          p <- p + annotate("text", x = 1, y = max(exps_gene$score) , label = paste0("Mann-Whitney, p = ", signif(stat$p.value, digits = 2)))
          p <- p + ggtitle(label = paste0(gene, " mRNA abundance between MSI-H tumor(N=", length(exps_msih),") and other tumors(N=", length(exps_msil), ")"))
          p <- p + theme_nogrid()
          p <- p + theme(axis.title.x = element_blank()) + theme(legend.position = "none") + theme(plot.title = element_text(size=8))
          p
          fn = paste0(makeOutDir(), datatype, "/", cancer, "/msi_thres", msi_thres, "/", gene, '_', datatype, '_abundance_between_MSI-H_and_rest_in_', cancer, ".pdf")
          ggsave(file=fn, height=4, width=5, useDingbats=FALSE)
        } else {
          print(paste0(gene, " not enough data!"))
        }
      }
      tmp <- unlist(mann_results[[cancer]])
      genes_sort <- as.vector(str_split_fixed(string = names(tmp[order(tmp)]), pattern = "\\.", 2)[,1])
      exps$Gene <- factor(exps$Gene, levels = genes_sort)
      p <- ggboxplot(exps, x = "Gene", y = "score",
                     color = "MSI_status", palette = "jco",
                     add = "jitter", xlab = "", ylab = "mRNA abundance (log2(FPKM+0.001))")
      p + stat_compare_means(aes(group = MSI_status), label = "p.signif")
      fn = paste0(makeOutDir(), datatype, "/", cancer, "/msi_thres", msi_thres, '_MMR_gene_', datatype, '_abundance_between_MSI-H_and_rest_in_', cancer, ".pdf")
      ggsave(file=fn, height=5, width=20, useDingbats=FALSE)
    }
  }
  
  for (cancer in c("CO", "UCEC")) {
    for (datatype in c( "protein", "phosphoprotein")) {
      exps <- readRDS(file = paste0(inputD, "analysis_results/expression_matrices/tables/intergrate_msi_expression/mmr_gene_msi_exps_list.RDS"))
      
      exps <- exps[[cancer]][[datatype]][, c("Sample", "Gene", "Score", "score")]; exps$expression_type <- datatype
      exps <- exps[exps$Gene %in% mmr_gene,]
      exps$MSI_status <- ifelse(exps$Score >= msi_thres, paste0("MSI-H_tumors"), paste0("other_tumors"))
      mann_results[[cancer]] <- vector("list")
      dir.create(paste0(makeOutDir(), datatype))
      dir.create(paste0(makeOutDir(), datatype, "/", cancer))
      dir.create(paste0(makeOutDir(), datatype, "/", cancer, "/msi_thres", msi_thres))
      
      for (gene in mmr_gene) {
        exps_gene <- exps[exps$Gene == gene & !is.na(exps$score),]
        exps_msih <- exps_gene$score[exps_gene$MSI_status == "MSI-H_tumors"]
        exps_msil <- exps_gene$score[exps_gene$MSI_status == "other_tumors"]
        
        if (length(exps_msih) > 0) {
          stat <- wilcox.test(x = exps_msih, y = exps_msil)
          mann_results[[cancer]][[gene]] <- list(p.value = stat$p.value)
          p <- ggplot(data = exps_gene)
          p <- p + geom_violin(mapping = aes(x = MSI_status, y = score, fill = MSI_status), alpha = 0.5, color = NA)
          p <- p + scale_fill_manual(values = c("MSI-H_tumors" = set1[1], "other_tumors" = "grey"))
          p <- p + geom_boxplot(mapping = aes(x = MSI_status, y = score), width=0.1, alpha = 0.8)
          p <- p + geom_point(mapping = aes(x = MSI_status, y = score), alpha = 0.5, stroke = 0, color = "black", position = position_jitter())
          p <- p + ylab(paste0(datatype, " abundance (log2 ratio)"))
          p <- p + annotate("text", x = 1, y = max(exps_gene$score) , label = paste0("Mann-Whitney, p = ", signif(stat$p.value, digits = 2)))
          p <- p + ggtitle(label = paste0(gene, " ", datatype, " abundance between MSI-H tumor(N=", length(exps_msih),") and other tumors(N=", length(exps_msil), ")"))
          p <- p + theme_nogrid()
          p <- p + theme(axis.title.x = element_blank()) + theme(legend.position = "none") + theme(plot.title = element_text(size=8))
          p
          fn = paste0(makeOutDir(), datatype, "/", cancer, "/msi_thres", msi_thres, "/", gene, '_', datatype, '_abundance_between_MSI-H_and_rest_in_', cancer, ".pdf")
          ggsave(file=fn, height=4, width=5, useDingbats=FALSE)
        } else {
          print(paste0(gene, " not enough data!"))
        }
      }
      tmp <- unlist(mann_results[[cancer]])
      genes_sort <- as.vector(str_split_fixed(string = names(tmp[order(tmp)]), pattern = "\\.", 2)[,1])
      exps$Gene <- factor(exps$Gene, levels = genes_sort)
      p <- ggboxplot(exps, x = "Gene", y = "score",
                     color = "MSI_status", palette = "jco",
                     add = "jitter", xlab = "", ylab = paste0(datatype, " abundance (log2 ratio)"))
      p + stat_compare_means(aes(group = MSI_status), label = "p.signif")
      fn = paste0(makeOutDir(), datatype, "/", cancer, "/msi_thres", msi_thres, '_MMR_gene_', datatype, '_abundance_between_MSI-H_and_rest_in_', cancer, ".pdf")
      ggsave(file=fn, height=5, width=20, useDingbats=FALSE)
    }
  }
}

# COAD & UCEC other genes -------------------------------------------------
genes_extra <- c("BLM", "AKT1", "AKT2", "AKT3")
for (msi_thres in c(0.8, 3.5)) {
  for (cancer in c("CO", "UCEC")) {
    for (datatype in c("mRNA")) {
      exps <- readRDS(file = paste0(inputD, "analysis_results/expression_matrices/tables/intergrate_msi_expression/mmr_gene_msi_exps_list.RDS"))
      exps <- exps[[cancer]][["mRNA"]][, c("Gene", "log2RPKM", "MSI_status")]; exps$expression_type <- "mRNA"
      colnames(exps) <- c("Gene", "score", "MSI_status", "expression_type")
      exps <- exps[exps$Gene %in% genes_extra,]
      mann_results[[cancer]] <- vector("list")
      dir.create(paste0(makeOutDir(), datatype))
      dir.create(paste0(makeOutDir(), datatype, "/", cancer))
      for (gene in genes_extra) {
        exps_gene <- exps[exps$Gene == gene & !is.na(exps$score),]
        exps_msih <- exps_gene$score[exps_gene$MSI_status == "MSI-H_tumors"]
        exps_msil <- exps_gene$score[exps_gene$MSI_status == "other_tumors"]
        
        if (length(exps_msih) > 4) {
          stat <- wilcox.test(x = exps_msih, y = exps_msil)
          mann_results[[cancer]][[gene]] <- list(p.value = stat$p.value)
          p <- ggplot(data = exps_gene)
          p <- p + geom_violin(mapping = aes(x = MSI_status, y = score, fill = MSI_status), alpha = 0.5, color = NA)
          p <- p + scale_fill_manual(values = c("MSI-H_tumors" = set1[1], "other_tumors" = "grey"))
          p <- p + geom_boxplot(mapping = aes(x = MSI_status, y = score), width=0.1, alpha = 0.8)
          p <- p + geom_point(mapping = aes(x = MSI_status, y = score), alpha = 0.5, stroke = 0, color = "black", position = position_jitter())
          p <- p + ylab("mRNA abundance (log2(FPKM+0.001))")
          p <- p + annotate("text", x = 1, y = max(exps_gene$score) , label = paste0("Mann-Whitney, p = ", signif(stat$p.value, digits = 2)))
          p <- p + ggtitle(label = paste0(gene, " mRNA abundance between MSI-H tumor(N=", length(exps_msih),") and other tumors(N=", length(exps_msil), ")"))
          p <- p + theme_nogrid()
          p <- p + theme(axis.title.x = element_blank()) + theme(legend.position = "none") + theme(plot.title = element_text(size=8))
          p
          fn = paste0(makeOutDir(), datatype, "/", cancer, "/msi_thres", msi_thres, "/", gene, '_', datatype, '_abundance_between_MSI-H_and_rest_in_', cancer, ".pdf")
          ggsave(file=fn, height=4, width=5, useDingbats=FALSE)
        } else {
          print(paste0(gene, " not enough data!"))
        }
      }
    }
  }
  for (cancer in c("CO", "UCEC")) {
    for (datatype in c( "protein", "phosphoprotein")) {
      exps <- readRDS(file = paste0(inputD, "analysis_results/expression_matrices/tables/intergrate_msi_expression/mmr_gene_msi_exps_list.RDS"))
      exps <- exps[[cancer]][[datatype]][, c("Gene", "score", "MSI_status")]
      exps <- exps[exps$Gene %in% genes_extra,]
      mann_results[[cancer]] <- vector("list")
      dir.create(paste0(makeOutDir(), datatype))
      dir.create(paste0(makeOutDir(), datatype, "/", cancer))
      dir.create(paste0(makeOutDir(), datatype, "/", cancer, "/msi_thres", msi_thres))
      
      for (gene in genes_extra) {
        exps_gene <- exps[exps$Gene == gene & !is.na(exps$score),]
        exps_msih <- exps_gene$score[exps_gene$MSI_status == "MSI-H_tumors"]
        exps_msil <- exps_gene$score[exps_gene$MSI_status == "other_tumors"]
        
        if (length(exps_msih) > 0) {
          stat <- wilcox.test(x = exps_msih, y = exps_msil)
          mann_results[[cancer]][[gene]] <- list(p.value = stat$p.value)
          p <- ggplot(data = exps_gene)
          p <- p + geom_violin(mapping = aes(x = MSI_status, y = score, fill = MSI_status), alpha = 0.5, color = NA)
          p <- p + scale_fill_manual(values = c("MSI-H_tumors" = set1[1], "other_tumors" = "grey"))
          p <- p + geom_boxplot(mapping = aes(x = MSI_status, y = score), width=0.1, alpha = 0.8)
          p <- p + geom_point(mapping = aes(x = MSI_status, y = score), alpha = 0.5, stroke = 0, color = "black", position = position_jitter())
          p <- p + ylab(paste0(datatype, " abundance (log2 ratio)"))
          p <- p + annotate("text", x = 1, y = max(exps_gene$score) , label = paste0("Mann-Whitney, p = ", signif(stat$p.value, digits = 2)))
          p <- p + ggtitle(label = paste0(gene, " ", datatype, " abundance between MSI-H tumor(N=", length(exps_msih),") and other tumors(N=", length(exps_msil), ")"))
          p <- p + theme_nogrid()
          p <- p + theme(axis.title.x = element_blank()) + theme(legend.position = "none") + theme(plot.title = element_text(size=8))
          p
          fn = paste0(makeOutDir(), datatype, "/", cancer, "/msi_thres", msi_thres, "/", gene, '_', datatype, '_abundance_between_MSI-H_and_rest_in_', cancer, ".pdf")
          ggsave(file=fn, height=4, width=5, useDingbats=FALSE)
        } else {
          print(paste0(gene, " not enough data!"))
        }
      }
    }
  }
  
  
}













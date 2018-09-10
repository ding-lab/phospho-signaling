# Yige Wu @ WashU 2018 Feb
# check for tumor-normal swap


# source ------------------------------------------------------------------
source('../phospho_network/phospho_network_shared.R')
source('cptac2p_analysis/preprocess_files/get_tumor_normal_pair.R')
source('cptac2p_analysis/phospho_network/phospho_network_plotting.R')
library(reshape)
library(ggplot2)
resultDnow <- makeOutDir()

# compare CDK1 for tumor-normal-matched pairs -----------------------------
genes <- c("CDK1", "TERT", "ERBB2", "EGFR")

for (cancer in c("BRCA","OV","CO")) {
  Pro.n_aveRep <- fread(paste(inputD, cancer, "/", prefix[cancer], "_PRO_formatted_normalized_replicate_averaged.txt",sep=""),
                        data.table = F)
  Pro.n_tumor <- get_tumor(Pro.n_aveRep, clinical.m = clinical)
  
  tumor_normal_match_tmp <- tumorSamp2normal(as.vector(colnames(Pro.n_tumor)), clinical = clinical)
  tumor_normal_match[[cancer]] <- tumor_normal_match_tmp
  
  for (gene in genes) {
    gene_pro <- Pro.n_aveRep[Pro.n_aveRep$Gene == gene,]
    if (nrow(gene_pro) > 0) {
      gene_pro_m <- melt(gene_pro)
      colnames(gene_pro_m)[2] <- "sampIDs"
      tn_match <- tumor_normal_match[[cancer]]
      tn_match_m <- melt.data.frame(tn_match, measure.vars = c("tumorSampIDs", "normalSampIDs"))
      tn_match_m$tumor_or_normal <- "tumor"
      tn_match_m$tumor_or_normal[tn_match_m$variable == "normalSampIDs"] <- "normal"
      gene_pro_m <- merge(gene_pro_m, tn_match_m, by.x = c("sampIDs"), by.y = c("value") )
      
      ## split patients into normal<tumor or normal>=tumor
      tn_exp <- merge(tn_match, gene_pro_m[,c("sampIDs", "value")], by.x = c("tumorSampIDs"), by.y = c("sampIDs"), all.x = T)
      colnames(tn_exp)[ncol(tn_exp)] <- "tumor_exp"
      tn_exp <- merge(tn_exp, gene_pro_m[,c("sampIDs", "value")], by.x = c("normalSampIDs"), by.y = c("sampIDs"), all.x = T)
      colnames(tn_exp)[ncol(tn_exp)] <- "normal_exp"
      tn_exp$tumor_minus_normal <- ifelse((tn_exp$tumor_exp > tn_exp$normal_exp), "normal<tumor", "normal>=tumor")
      gene_pro_m$tumor_minus_normal <- ifelse((gene_pro_m$partID %in% tn_exp$partID[tn_exp$tumor_minus_normal == "normal<tumor"]), "normal<tumor", "normal>=tumor")
      
      p <- ggplot(data = gene_pro_m, mapping = aes(x = tumor_or_normal, y = value, color = partID, group = partID))
      p <- p + geom_line()
      p <- p + theme_bw()
      p <- p + theme(legend.position="none",
                     axis.title.x = element_blank(),
                     axis.title.y = element_blank())
      p <- p + ggtitle(label = paste(gene, " protein expression"))
      p <- p + facet_grid(facets = . ~ tumor_minus_normal)
      p <- p + geom_text(aes(label= ifelse(tumor_or_normal == "normal", as.character(partID), NA), 
                             color = partID), size= 2, nudge_x = -0.2, alpha = 0.6)
      p
      ggsave(filename = paste0(resultDnow, gene, "_", cancer, "_tumor_vs_normal_Protein_expression.pdf"),
             height = 4, width = 6)
    }
  }
}




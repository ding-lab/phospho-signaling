# Yige Wu @ WashU 2018 Mar
## heatmaps showing prospective BRCA PRM peptide-level abundance and TMT protein abundance from CDAP
## match parent specimen labels

# souce -------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/preprocess_files/preprocess_files.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_plotting.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')

# inputs ------------------------------------------------------------------
subtype_map <- read_excel("~/Box Sync/cptac2p/cptac_shared/5_CPTAC2_Breast_Prospective_Collection_BI/proteome-data-v1.01011/data/20171106_CPTAC2_ProspectiveBreastCancer_Sample_Annotations_Table_v79.xlsx")
prm_parent_labeled <- fread(input = "./cptac2p/cptac_shared/analysis_results/validation/cor_prm_ms_match_specimen/BRCA_prm_trplicatemean_median_parent_labeled.txt", data.table = F)
ms_parent_labeled <- fread(input = "./cptac2p/cptac_shared/analysis_results/validation/cor_prm_ms_match_specimen/BRCA_tmt_match_prm_parent_labeled.txt", data.table = F)
cancer <- "BRCA"
cor_df <- fread(input = paste0("./cptac2p/cptac_shared/analysis_results/validation/tables/cor_prm_ms_match_specimen/", cancer, "_",  "spearman", "_correlation_prm_ms_match_parent_specimen.txt"))

# functions ---------------------------------------------------------------
ecdf_fun = function(x,perc) ecdf(x)(perc)

expression_effect = function(m){ 
  cat("##### EXPRESSION ANALYSIS #####\n")
  minNum = 5
  m = as.matrix(m)
  num = nrow(m)
  m2 = as.matrix(m[rowSums(!is.na(m)) >= minNum, ])
  num_NA= nrow(m2)
  cat(paste("Original number of markers:", num, "; NA filtered:", num_NA, "\n", sep=" "))
  
  # initiate tables
  outlier = matrix(data = NA,nrow=dim(m2)[1],ncol=dim(m2)[2])
  row.names(outlier) = row.names(m2)
  colnames(outlier) = colnames(m2)
  exp_score = outlier
  exp_quantile = outlier
  
  # gene-wise expression score and quantile score
  for (i in 1:nrow(m2)){
    #IQR = quantile(m2[i,], probs=0.75, na.rm=T) - quantile(m2[i,], probs=0.25, na.rm=T) 
    exp_score[i,] = m2[i,]#(m2[i,] - quantile(m2[i,], probs=0.50, na.rm=T))/IQR
    if (length(m2[i, m2[i,] != 0 & !is.na(m2[i,] != 0)]) > 0) {
      exp_quantile[i, m2[i,] != 0 & !is.na(m2[i,] != 0)] = ecdf_fun(m2[i, m2[i,] != 0 & !is.na(m2[i,] != 0)], m2[i, m2[i,] != 0 & !is.na(m2[i,] != 0)])
    }
  }
  return(list("exp_score"=exp_score, "exp_quantile"=exp_quantile))
}

# get super table ---------------------------------------------------------
## melt prm table
prm_qt <- data.frame(Gene = prm_parent_labeled$Gene)
tmp <- expression_effect(prm_parent_labeled[, !(colnames(prm_parent_labeled) %in% c("Gene"))])
prm_qt <- cbind(prm_qt, tmp$exp_quantile)
prm_qt_m <- melt(prm_qt, id.vars = c("Gene"))
colnames(prm_qt_m) <- c("Gene", "Specimen.Label", "qt.prm")
prm_value_m <- melt(prm_parent_labeled, id.vars = c("Gene"))
colnames(prm_value_m) <- c("Gene", "Specimen.Label", "value.prm")

## melt tmt table
tmt_prm_genes <- ms_parent_labeled[ms_parent_labeled$Gene %in% prm_parent_labeled$Gene,]
tmt_qt <- data.frame(Gene = tmt_prm_genes$Gene)
tmp <- expression_effect(tmt_prm_genes[, !(colnames(tmt_prm_genes) %in% c("Gene"))])
tmt_qt <- cbind(tmt_qt, tmp$exp_quantile)
tmt_qt_m <- melt(tmt_qt, id.vars = c("Gene"))
colnames(tmt_qt_m) <- c("Gene", "Specimen.Label", "qt.tmt")
tmt_value_m <- melt(tmt_prm_genes, id.vars = c("Gene"))
colnames(tmt_value_m) <- c("Gene", "Specimen.Label", "value.tmt")

prm_tmt_sup <- merge(prm_qt_m, tmt_qt_m, all = T)
prm_tmt_sup <- merge(prm_tmt_sup, prm_value_m, all = T)
prm_tmt_sup <- merge(prm_tmt_sup, tmt_value_m, all = T)
prm_tmt_sup$qt.tmt_prm = prm_tmt_sup$qt.tmt - prm_tmt_sup$qt.prm

samples.t <- unique(prm_tmt_sup$Specimen.Label)
samples.t.pam50 <- sampID2pam50(sampleID_vector = samples.t, subtype_map = subtype_map)
samples.t.pam50[is.na(samples.t.pam50)] <- "NA"
samples.t.pam50 <- data.frame(Specimen.Label = samples.t, pam50 = samples.t.pam50)
prm_tmt_sup <- merge(prm_tmt_sup, samples.t.pam50, all = T)

# heatmap for quantile offset (TMT-PRM) all genes ---------------------------------------------
df <- prm_tmt_sup[!is.na(prm_tmt_sup$qt.tmt_prm),]
p = ggplot(df)
p = p + geom_tile(aes(x=Specimen.Label, y=Gene, fill=qt.tmt_prm), color=NA)#, linetype="blank") 
p = p + facet_grid(.~pam50, drop=T, space = "free",scales = "free")#, space = "free", scales = "free")
p = p + scale_fill_gradientn(name= "protein abundance quantile (TMT-PRM)", na.value=NA, colours=RdBu1024, limit=c(-1,1))
p = p + theme_bw() + theme_nogrid()
p = p + theme(axis.title = element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(size = 5))
p = p + theme(axis.ticks=element_blank(),legend.position="bottom")
p = p + theme(panel.spacing.x=unit(0, "lines"), panel.spacing.y=unit(0, "lines"))
p
resultDnow <- makeOutDir()
fn = paste0(resultDnow, 'qt.tmt_prm_in_', cancer, ".pdf")
ggsave(file=fn, height=12, width=5, useDingbats=FALSE)

# heatmap for quantile offset (TMT-PRM) positive correlated (fdr<0.1) ---------------------------------------------
library(reshape2)
df <- prm_tmt_sup[!is.na(prm_tmt_sup$qt.tmt_prm) & (prm_tmt_sup$Gene %in% cor_df$Gene[cor_df$fdr<0.1]) ,]
df <- unique(df)
casted = dcast( data.frame(df[,c("Gene", "Specimen.Label", "qt.tmt_prm")]) , Gene~Specimen.Label, value.var = "qt.tmt_prm")
# tmp <- as.vector(cor_df$Gene)[order(cor_df$rho)]
# df$Gene <- factor(df$Gene, levels = tmp[tmp %in% df$Gene])
df$Gene <- factor(df$Gene, levels = as.vector(casted$Gene)[order(rowMeans(casted[,-1], na.rm = T))])
df$Specimen.Label <- factor(df$Specimen.Label, levels = as.vector(colnames(casted))[order(colMeans(casted[,-1], na.rm = T))])

p = ggplot(df)
p = p + geom_tile(aes(x=Specimen.Label, y=Gene, fill=qt.tmt_prm), color=NA)#, linetype="blank") 
p = p + facet_grid(.~pam50, drop=T, space = "free",scales = "free")#, space = "free", scales = "free")
p = p + scale_fill_gradientn(name= "protein abundance quantile (TMT-PRM)", na.value=NA, colours=RdBu1024, limit=c(-1,1))
p = p + theme_bw() + theme_nogrid()
p = p + theme(axis.title = element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(size = 5))
p = p + theme(axis.ticks=element_blank(),legend.position="bottom")
p = p + theme(panel.spacing.x=unit(0, "lines"), panel.spacing.y=unit(0, "lines"))
p
resultDnow <- makeOutDir()
fn = paste0(resultDnow, 'qt.tmt_prm_in_', cancer, "_sig_cor_fdr_0.1.pdf")
ggsave(file=fn, height=8, width=5, useDingbats=FALSE)

p = ggplot(df)
p = p + geom_tile(aes(x=Specimen.Label, y=Gene, fill=qt.prm), color=NA)#, linetype="blank") 
p = p + facet_grid(.~pam50, drop=T, space = "free",scales = "free")#, space = "free", scales = "free")
p = p + scale_fill_gradientn(name= "protein abundance quantile (PRM)", na.value=NA, colours=RdBu1024, limit=c(0,1))
p = p + theme_bw() + theme_nogrid()
p = p + theme(axis.title = element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(size = 5))
p = p + theme(axis.ticks=element_blank(),legend.position="bottom")
p = p + theme(panel.spacing.x=unit(0, "lines"), panel.spacing.y=unit(0, "lines"))
p
resultDnow <- makeOutDir()
fn = paste0(resultDnow, 'qt.prm_in_', cancer, "_sig_cor_fdr_0.1.pdf")
ggsave(file=fn, height=8, width=5, useDingbats=FALSE)



p = ggplot(df)
p = p + geom_tile(aes(x=Specimen.Label, y=Gene, fill=qt.tmt), color=NA)#, linetype="blank") 
p = p + facet_grid(.~pam50, drop=T, space = "free",scales = "free")#, space = "free", scales = "free")
p = p + scale_fill_gradientn(name= "protein abundance quantile (TMT)", na.value=NA, colours=RdBu1024, limit=c(0,1))
p = p + theme_bw() + theme_nogrid()
p = p + theme(axis.title = element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(size = 5))
p = p + theme(axis.ticks=element_blank(),legend.position="bottom")
p = p + theme(panel.spacing.x=unit(0, "lines"), panel.spacing.y=unit(0, "lines"))
p
resultDnow <- makeOutDir()
fn = paste0(resultDnow, 'qt.tmt_in_', cancer, "_sig_cor_fdr_0.1.pdf")
ggsave(file=fn, height=8, width=5, useDingbats=FALSE)


# prm not in tmt ----------------------------------------------------------
df <- prm_tmt_sup[!is.na(prm_tmt_sup$value.prm) & !(prm_tmt_sup$Gene %in% prm_tmt_sup$Gene[!is.na(prm_tmt_sup$value.tmt)]),]
df <- unique(df)
df$value.prm_capped <- df$value.prm
df$value.prm_capped[df$value.prm_capped > 1] <- 1
tmp <- as.vector(df$pam50)
tmp[is.na(tmp) | (df$pam50) == "NA" | (df$pam50) == "xx" ] <- "unknown_PAM50_subtype"
df$pam50 <- tmp
df$pam50 <- factor(df$pam50, levels = c("Basal", "Her2", "LumA", "LumB", "Normal", "unknown_PAM50_subtype"))
p = ggplot(df)
p = p + geom_tile(aes(x=Specimen.Label, y=Gene, fill=value.prm), color=NA)#, linetype="blank") 
p = p + facet_grid(.~pam50, drop=T, space = "free",scales = "free")#, space = "free", scales = "free")
p = p + scale_fill_gradientn(name= "absolute protein abundance by PRM (fmol/ul, capped at 1)", na.value=NA, colors = RdBu1024[513:1024], limit=c(0,1))
p = p + theme_bw() + theme_nogrid()
p = p + theme(axis.title = element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(size = 5))
p = p + theme(axis.ticks=element_blank(),legend.position="bottom")
p = p + theme(panel.spacing.x=unit(0, "lines"), panel.spacing.y=unit(0, "lines"))
p
resultDnow <- makeOutDir()
fn = paste0(resultDnow, 'protein.level.prm_in_', cancer, "_not_in_tmt.pdf")
ggsave(file=fn, height=3, width=5, useDingbats=FALSE)

write.table(unique(df$Gene), file = paste0(resultDnow, "proteins_in_prm_notin_tmt.txt"), row.names = F, quote = F, col.names = F)

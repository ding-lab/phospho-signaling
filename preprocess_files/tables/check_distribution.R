# Yige Wu @ WashU Jan 2018
# check the distribution of log ratios per sample (for deciding normalization method)
# load packages -----------------------------------------------------------
library(ggplot2)

# inputs ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/preprocess_files/preprocess_files.R')
clinical <- read_excel(paste0(baseD, "Specimen Data_20161005_Yige_20171212.xls"))
clinical <- data.frame(clinical)
clinical_anno <- clinical

# process proteome and phosphoproteome data ---------------------------------------------
for (cancer in c("BRCA", "OV", "CO")) {
#  cancer <- "BRCA"
  ## process proteome
  Pro <- fread(raw[cancer, "protein"], data.table=FALSE, verbose=TRUE)
  
  # keep shared peptide log ratio
  Pro.f = unshared_pro(Pro) 
  
  ## plot the distribution of un-normalized unshared log ratio
  Pro.f.m <- melt(Pro.f); colnames(Pro.f.m)[1] <- 'specimen_id'
  p <- ggplot(Pro.f.m, aes(x=value))
  p <- p + geom_density(aes(group=specimen_id, colour=specimen_id))
  p <- p + theme(legend.position="none")
  p
  ggsave(filename = paste0(outputD, cancer, '_Pro_unnormalized_log_ratios_by_sample.pdf'))
  rm(Pro.f.m)
  rm(p)
  
  # normalize by sample
  Pro.n = normalize_by_sample(Pro.f)
  
  ## plot the distribution of un-normalized unshared log ratio
  Pro.n.m <- melt(Pro.n); colnames(Pro.n.m) <- c('gene', 'specimen_id', 'value')
  p <- ggplot(Pro.n.m, aes(x=value))
  p <- p + geom_density(aes(group=specimen_id, colour=specimen_id))
  p <- p + theme(legend.position="none")
  p
  ggsave(filename = paste0(outputD, cancer, '_Pro_normalized_log_ratios_by_sample.pdf'))
  rm(Pro.n.m)
  rm(p)
  
  ## process phosphoproteome
  Pho <- fread(raw[cancer, "phosphoprotein"], data.table=FALSE)
  Pho.f = unshared_pho(Pho, Pho[,c("Gene","Phosphosite")])
  Pho.n = normalize_by_sample(Pho.f[,!(colnames(Pho.f) %in% c("Gene", "Phosphosite"))])
  
  ## plot the distribution of un-normalized unshared log ratio
  Pho.f.m <- melt(Pho.f); 
  colnames(Pho.f.m)[3] <- 'specimen_id'
  p <- ggplot(Pho.f.m, aes(x=value))
  p <- p + geom_density(aes(group=specimen_id, colour=specimen_id))
  p <- p + theme(legend.position="none")
  p
  ggsave(filename = paste0(outputD, cancer, '_Pho_unnormalized_log_ratios_by_sample.pdf'))
  rm(Pho.f.m)
  rm(p)
  
  ## plot the distribution of un-normalized unshared log ratio
  Pho.n.m <- melt(Pho.n); colnames(Pho.n.m) <- c('gene', 'specimen_id', 'value')
  p <- ggplot(Pho.n.m, aes(x=value))
  p <- p + geom_density(aes(group=specimen_id, colour=specimen_id))
  p <- p + theme(legend.position="none")
  p
  ggsave(filename = paste0(outputD, cancer, '_Pho_normalized_log_ratios_by_sample.pdf'))
  rm(Pho.n.m)
  rm(p)
  
}

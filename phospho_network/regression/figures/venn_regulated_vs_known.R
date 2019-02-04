# Yige Wu @ WashU 2018 Apr
# venn diagram showing whether regulated enzyme-substrate pairs only in breast are detected in colorectal and ovarian cancer

# source ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_plotting.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
library(eulerr)


# set variables -----------------------------------------------------------

# inputs ------------------------------------------------------------------
reg_nonNA <- 20
fdr_thres <- c(0.05, 0.1); names(fdr_thres) <- c("kinase", "phosphatase")
SELF = "trans"
num_top <- 10


# first plot site-level comparison with known ones ------------------------
for (enzyme_type in c("kinase")) {
  # inputs ------------------------------------------------------------------
  regression <- fread(input = paste0(ppnD, "regression/tables/", "change_regression_nonNA/", 
                                     "regression", "_", "cptac2p_cptac3", "_", "tumor", "_reg_nonNA", reg_nonNA, ".txt"), data.table = F)
  ## filtering
  regression <- regression[regression$enzyme_type == enzyme_type & regression$SELF == SELF,]
  regression <- markSigKS(regression = regression, sig_thres = fdr_thres[enzyme_type], enzyme_type = enzyme_type)
  
  
  dat <- data.frame(pair = unique(regression$pair[regression$pair_pro %in% regression$pair_pro[!is.na(regression$Source) & !(regression$Source %in% c("PhosphoNetworks", "MIMP", "NetKIN"))]]))
  dat$covariation <- (dat$pair %in% regression$pair[regression$regulated])
  dat$experimental <- (dat$pair %in% regression$pair[!is.na(regression$Source) & !(regression$Source %in% c("PhosphoNetworks", "MIMP", "NetKIN"))])
  dat$PhosphoNetworks <- (dat$pair %in% regression$pair[!is.na(regression$Source) & grepl(pattern = "PhosphoNetworks", x = regression$Source)])
  dat$MIMP <- (dat$pair %in% regression$pair[!is.na(regression$Source) & grepl(pattern = "MIMP", x = regression$Source)])
  dat$NetKIN <- (dat$pair %in% regression$pair[!is.na(regression$Source) & grepl(pattern = "NetKIN", x = regression$Source)])
  
  
  for (sequence_prediction in c("PhosphoNetworks", "MIMP", "NetKIN")) {
    fn = paste(makeOutDir(resultD = resultD), enzyme_type, '_substrate_site_level_pairs_covariation_vs_experimental_vs', sequence_prediction, '.pdf',sep ="")
    grid.newpage()
    pdf(fn, height = 12, width = 20, useDingbats = FALSE)
    fit <- euler(combinations = dat[, c("covariation", "experimental", sequence_prediction)], input = "disjoint", shape = 'circle')
    p <-plot(fit, quantities = list(fontsize = 30), legend = list(fontsize = 30))
    grid.draw(p)
    dev.off()
  }
  
  fn = paste(makeOutDir(resultD = resultD), enzyme_type, '_substrate_site_level_pairs_covariation_vs_experimental.pdf',sep ="")
  grid.newpage()
  pdf(fn, height = 12, width = 20, useDingbats = FALSE)
  fit <- euler(combinations = dat[, c("covariation", "experimental")], input = "disjoint", shape = 'circle')
  p <-plot(fit, quantities = list(fontsize = 80), legend = F)
  grid.draw(p)
  dev.off()
}

# first plot protein-level comparison with known ones ------------------------
for (enzyme_type in c("kinase")) {
  # inputs ------------------------------------------------------------------
  regression <- fread(input = paste0(ppnD, "regression/tables/", "change_regression_nonNA/", 
                                     "regression", "_", "cptac2p_cptac3", "_", "tumor", "_reg_nonNA", reg_nonNA, ".txt"), data.table = F)
  ## filtering
  regression <- regression[regression$enzyme_type == enzyme_type & regression$SELF == SELF,]
  regression <- markSigKS(regression = regression, sig_thres = fdr_thres[enzyme_type], enzyme_type = enzyme_type)
  
  
  dat <- data.frame(pair_pro = unique(regression$pair_pro[!is.na(regression$Source) & !(regression$Source %in% c("PhosphoNetworks", "MIMP", "NetKIN"))]))
  dat$covariation <- (dat$pair_pro %in% regression$pair_pro[regression$regulated])
  dat$experimental <- (dat$pair_pro %in% regression$pair_pro[!is.na(regression$Source) & !(regression$Source %in% c("PhosphoNetworks", "MIMP", "NetKIN"))])
  dat$PhosphoNetworks <- (dat$pair_pro %in% regression$pair_pro[!is.na(regression$Source) & grepl(pattern = "PhosphoNetworks", x = regression$Source)])
  dat$MIMP <- (dat$pair_pro %in% regression$pair_pro[!is.na(regression$Source) & grepl(pattern = "MIMP", x = regression$Source)])
  dat$NetKIN <- (dat$pair_pro %in% regression$pair_pro[!is.na(regression$Source) & grepl(pattern = "NetKIN", x = regression$Source)])
  
  
  for (sequence_prediction in c("PhosphoNetworks", "MIMP", "NetKIN")) {
    fn = paste(makeOutDir(resultD = resultD), enzyme_type, '_substrate_protein_level_pairs_covariation_vs_experimental_vs', sequence_prediction, '.pdf',sep ="")
    grid.newpage()
    pdf(fn, height = 12, width = 20, useDingbats = FALSE)
    fit <- euler(combinations = dat[, c("covariation", "experimental", sequence_prediction)], input = "disjoint", shape = 'circle')
    p <-plot(fit, quantities = list(fontsize = 30), legend = list(fontsize = 30))
    grid.draw(p)
    dev.off()
  }
  
  fn = paste(makeOutDir(resultD = resultD), enzyme_type, '_substrate_protein_level_pairs_covariation_vs_experimental.pdf',sep ="")
  grid.newpage()
  pdf(fn, height = 12, width = 20, useDingbats = FALSE)
  fit <- euler(combinations = dat[, c("covariation", "experimental")], input = "disjoint", shape = 'circle')
  p <-plot(fit, quantities = list(fontsize = 80), legend = list(fontsize = 80))
  grid.draw(p)
  dev.off()
}


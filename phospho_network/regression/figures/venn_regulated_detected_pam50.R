# Yige Wu @ WashU 2018 Apr
# venn diagram showing whether regulated enzyme-substrate pairs only in breast are detected in colorectal and ovarian cancer

# source ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_plotting.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
library(eulerr)
library(colorspace)

# set variables -----------------------------------------------------------
subtypes <- c("Her2", "LumA", "LumB", "Basal")


# inputs ------------------------------------------------------------------
## input enzyme-substrate pairs examined
for (enzyme_type in c("phosphatase", "kinase")) {
  sup_cans_tab_en <- fread(input = paste0("./cptac2p/analysis_results/phospho_network/regression/tables/regression_omnipath&newworkin&depod_protein_level/", 
                                          enzyme_type, "_substrate_regression_cptac2p_3can_tumor.txt"), data.table = F)
  sup_cans_tab_en$GENE <- as.vector(sup_cans_tab_en$KINASE)
  sup_cans_tab_en$SUB_GENE <- as.vector(sup_cans_tab_en$SUBSTRATE)
  sup_cans_tab_en <- markSigSiteCan(regression = sup_cans_tab_en, sig_thres = reg_sig[enzyme_type], enzyme_type = enzyme_type)
  sup_cans_tab_en$regulated <- (sup_cans_tab_en$fdr_sig & sup_cans_tab_en$coef_sig)
  
  detected_cans <- NULL
  regulated_cans <- NULL
  pairs_all <- NULL
  detected_cans <- list()
  regulated_cans <- list()
  for (cancer in cancers_sort) {
    detected_cans[[cancer]] <- as.vector(unique(sup_cans_tab_en$pair[sup_cans_tab_en$Cancer == cancer]))
    regulated_cans[[cancer]] <- as.vector(unique(sup_cans_tab_en$pair[sup_cans_tab_en$regulated & sup_cans_tab_en$Cancer == cancer]))
    pairs_all <- unique(c(pairs_all, detected_cans[[cancer]]))
  }
  for (cancer in c("BRCA")) {
    sup_cans_tab_en <- fread(input = paste0("./cptac2p/analysis_results/phospho_network/regression/tables/regression_omnipath&newworkin&depod_protein_level_pam50/", 
                                            enzyme_type, "_substrate_regression_cptac2p_", cancer, "_pam50_tumor.txt"), data.table = F)
    sup_cans_tab_en$GENE <- as.vector(sup_cans_tab_en$KINASE)
    sup_cans_tab_en$SUB_GENE <- as.vector(sup_cans_tab_en$SUBSTRATE)
    sup_cans_tab_en <- markSigKS(regression = sup_cans_tab_en, sig_thres = reg_sig[enzyme_type], enzyme_type = enzyme_type)
    sup_cans_tab_en$regulated <- (sup_cans_tab_en$fdr_sig & sup_cans_tab_en$coef_sig)
    
    for (subtype in subtypes) {
      detected_cans[[subtype]] <- as.vector(unique(sup_cans_tab_en$pair[sup_cans_tab_en$Cancer == cancer & sup_cans_tab_en$subtype == subtype]))
      regulated_cans[[subtype]] <- as.vector(unique(sup_cans_tab_en$pair[sup_cans_tab_en$regulated & sup_cans_tab_en$Cancer == cancer & sup_cans_tab_en$subtype == subtype]))
    }
    
  }
  
  dat <- data.frame(pair = pairs_all)
  for (cancer in cancers_sort) {
    dat[, paste0("detected_", cancer)] <- FALSE
    dat[dat$pair %in% detected_cans[[cancer]], paste0("detected_", cancer)]  <- TRUE
    
    dat[, paste0("regulated_", cancer)] <- FALSE
    dat[dat$pair %in% regulated_cans[[cancer]], paste0("regulated_", cancer)]  <- TRUE
  }
  for (subtype in subtypes) {
    dat[, paste0("detected_", subtype)] <- FALSE
    dat[dat$pair %in% detected_cans[[subtype]], paste0("detected_", subtype)]  <- TRUE
    
    dat[, paste0("regulated_", subtype)] <- FALSE
    dat[dat$pair %in% regulated_cans[[subtype]], paste0("regulated_", subtype)]  <- TRUE
  }
  
  dat <- data.frame(dat)
  
  for (cancer in cancers_sort) {
    fn = paste(makeOutDir(resultD = resultD), enzyme_type, '_substrate_pairs_regulated_pam50&', cancer, '_venn.pdf',sep ="")
    grid.newpage()
    pdf(fn, height = 12, width = 12, useDingbats = FALSE)
    fit <- euler(combinations = dat[, c("pair", paste0("regulated_", c(cancer, subtypes)))], input = "disjoint", shape = "ellipse")
    p <-plot(fit, quantities = list(fontsize = 30), legend = list(fontsize = 20), fills = c(color_cancers2[cancer], diverge_hcl(5)[rainbow_hcl(5) != color_cancers2]))
    grid.draw(p)
    dev.off()
  }
  as.vector(dat$pair[dat$regulated_Basal & !(dat$regulated_BRCA)])
  as.vector(dat$pair[dat$regulated_LumA & !(dat$regulated_BRCA)])
  as.vector(dat$pair[dat$regulated_Basal & (dat$regulated_OV)])
}

## get list of detected and regulated

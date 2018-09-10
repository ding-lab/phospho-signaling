# Yige Wu @ WashU 2018 Aug
## take the wilcox test p-values from compari the substrate phosphorylation level of kinase-mutated tumors, kinase-amplified tumors, kinase-deleted tumors against other tumors for each kinase-substrate pair
## adjust to FDR 

# source ------------------------------------------------------------------
source("./cptac2p_analysis/phospho_network/phospho_network_shared.R")

# set variables -----------------------------------------------------------
amp_log2cn <- 0.3
del_log2cn <- -0.3
num_control_thres <- 10
num_genoalt_bottom <- 5
num_genoalt_upper <- 30
networkin_thres <- 10

# inputs ------------------------------------------------------------------
## input enzyme-substrate table
# ptms_site_pairs_sup <- read_csv(paste0(ppnD, "compile_enzyme_substrate/tables/compile_omnipath_networkin_depod/omnipath_networkin_enzyme_substrate_site_level_union_source_cleaned_addDEPOD_extended.csv"))
ptms_site_pairs_sup <- read_csv(paste0(ppnD, "compile_enzyme_substrate/tables/compile_omnipath_networkin_depod_signor/omnipath_networkin_DEPOD_SignorNotSiteMapped.csv"))
# ptms_site_pairs_sup <- ptms_site_pairs_sup[ptms_site_pairs_sup$Source == "NetKIN" | (ptms_site_pairs_sup$Source != "NetKIN" & ptms_site_pairs_sup$networkin_score >= 5),]
ptms_site_pairs_sup <- ptms_site_pairs_sup[ptms_site_pairs_sup$Source != "NetKIN" | (ptms_site_pairs_sup$Source == "NetKIN" & ptms_site_pairs_sup$networkin_score >= networkin_thres),]
ptms_site_pairs_sup <- data.frame(ptms_site_pairs_sup)

# adjust p-value based on number of deleted/amplified/mutated samples & number of control samples ----------------------------------------------------------------------
for (cancer in c("OV", "BRCA")) {
  mut_cnv_tab <- fread(input = paste0(ppnD, "genoalt/tables/mut_cnv_impact/", cancer, "_mut_cnv_tab.txt"), data.table = F)
  mut_cnv_tab$pair_pro <- paste0(mut_cnv_tab$GENE, ":", mut_cnv_tab$SUB_GENE)
  mut_cnv_tab <- mut_cnv_tab[mut_cnv_tab$pair_pro %in% ptms_site_pairs_sup$pair_pro,]
  mut_cnv_tab$pair <- paste0(mut_cnv_tab$GENE, ":", mut_cnv_tab$SUB_GENE, ":", mut_cnv_tab$SUB_MOD_RSD)
  mut_cnv_tab$pair_know <- (mut_cnv_tab$pair %in% ptms_site_pairs_sup$pair)
  for (genoalt_type in c("mut", "amp", "del")) {
    col2adjust <- paste0("p_", genoalt_type)
    row2adjust <- which(mut_cnv_tab[, col2adjust] > 0 & mut_cnv_tab[, paste0("num_", genoalt_type)] >= num_genoalt_bottom & mut_cnv_tab[, paste0("num_", genoalt_type)] <= num_genoalt_upper & mut_cnv_tab$num_control >= num_control_thres)
    mut_cnv_tab[, paste0("fdr_", genoalt_type)] <- 0
    mut_cnv_tab[row2adjust, paste0("fdr_", genoalt_type)] <- p.adjust(mut_cnv_tab[row2adjust, paste0("p_", genoalt_type)], method = "fdr")
  }
  mut_cnv_mut_tab <- mut_cnv_tab[mut_cnv_tab$fdr_mut > 0,]; mut_cnv_mut_tab <- mut_cnv_mut_tab[order(mut_cnv_mut_tab$fdr_mut, mut_cnv_mut_tab$p_mut),]
  mut_cnv_amp_tab <- mut_cnv_tab[mut_cnv_tab$fdr_amp > 0,]; mut_cnv_amp_tab <- mut_cnv_amp_tab[order(mut_cnv_amp_tab$fdr_amp, mut_cnv_amp_tab$p_amp),]
  mut_cnv_del_tab <- mut_cnv_tab[mut_cnv_tab$fdr_del > 0,]; mut_cnv_del_tab <- mut_cnv_del_tab[order(mut_cnv_del_tab$fdr_del, mut_cnv_del_tab$p_del),]
}

# adjust p-value based on number of deleted/amplified/mutated samples & number of control samples & known enzyme-substrate pairs----------------------------------------------------------------------
for (cancer in c("OV", "BRCA")) {
  mut_cnv_tab <- fread(input = paste0(ppnD, "genoalt/tables/mut_cnv_impact/", cancer, "_mut_cnv_tab.txt"), data.table = F)
  mut_cnv_tab$pair_pro <- paste0(mut_cnv_tab$GENE, ":", mut_cnv_tab$SUB_GENE)
  mut_cnv_tab <- mut_cnv_tab[mut_cnv_tab$pair_pro %in% ptms_site_pairs_sup$pair_pro,]
  mut_cnv_tab$pair <- paste0(mut_cnv_tab$GENE, ":", mut_cnv_tab$SUB_GENE, ":", mut_cnv_tab$SUB_MOD_RSD)
  mut_cnv_tab$pair_know <- (mut_cnv_tab$pair %in% ptms_site_pairs_sup$pair)
  for (genoalt_type in c("mut", "amp", "del")) {
    col2adjust <- paste0("p_", genoalt_type)
    row2adjust <- which(mut_cnv_tab[, col2adjust] > 0 & mut_cnv_tab$pair_know & mut_cnv_tab[, paste0("num_", genoalt_type)] >= num_genoalt_bottom & mut_cnv_tab[, paste0("num_", genoalt_type)] <= num_genoalt_upper & mut_cnv_tab$num_control >= num_control_thres)
    mut_cnv_tab[, paste0("fdr_", genoalt_type)] <- 0
    mut_cnv_tab[row2adjust, paste0("fdr_", genoalt_type)] <- p.adjust(mut_cnv_tab[row2adjust, paste0("p_", genoalt_type)], method = "fdr")
  }
  mut_cnv_mut_tab <- mut_cnv_tab[mut_cnv_tab$fdr_mut > 0,]; mut_cnv_mut_tab <- mut_cnv_mut_tab[order(mut_cnv_mut_tab$fdr_mut, mut_cnv_mut_tab$p_mut),]
  mut_cnv_amp_tab <- mut_cnv_tab[mut_cnv_tab$fdr_amp > 0,]; mut_cnv_amp_tab <- mut_cnv_amp_tab[order(mut_cnv_amp_tab$fdr_amp, mut_cnv_amp_tab$p_amp),]
  mut_cnv_del_tab <- mut_cnv_tab[mut_cnv_tab$fdr_del > 0,]; mut_cnv_del_tab <- mut_cnv_del_tab[order(mut_cnv_del_tab$fdr_del, mut_cnv_del_tab$p_del),]
}

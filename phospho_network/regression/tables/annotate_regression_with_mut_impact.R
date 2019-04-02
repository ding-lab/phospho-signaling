# Yige Wu @ WashU 2019 Feb
## annotate the regression result with mutational impact


# source ------------------------------------------------------------------
baseD = "/Users/yigewu/Box\ Sync/"
wd <- getwd()
if (wd != baseD) {
  setwd(baseD)
}
source("./cptac2p_analysis/phospho_network/phospho_network_shared.R")

# set variables -----------------------------------------------------------
num_genoalt_thres <- 5
reg_nonNA <- 20
fdr_thres <- c(0.05, 0.05); names(fdr_thres) <- c("kinase", "phosphatase")
affected_exp_types2process <- c("PRO", "PHO")
cancers2process <- c("BRCA", "OV", "CO", "UCEC", "CCRCC")

# input mutational impact on proteome -------------------------------------
mut_impact_tab <- load_mut_impact_proteome()
mut_impact_tab %>% head()

## filter for protein and phospho
mut_impact_tab <- mut_impact_tab[mut_impact_tab$affected_exp_type %in% affected_exp_types2process,]
mut_impact_tab %>% head()

# input regression super table --------------------------------------------
regression <- fread(input = paste0(ppnD, "regression/tables/", "change_regression_nonNA/", 
                                "regression", "_", "cptac2p_cptac3", "_", "tumor", "_reg_nonNA", reg_nonNA, ".txt"), data.table = F)
regression %>% nrow()

# annotate with regression significance -----------------------------------
regression <- change_regression_nonNA(regression = regression, reg_nonNA = reg_nonNA, reg_sig = fdr_thres)
regression %>% nrow()

regression$pair_pro <- paste0(regression$GENE, ":", regression$SUB_GENE)
regression$pair_pro_cancer <- paste0(regression$pair_pro, ":", regression$Cancer)
regression$pair_pro_rev <- paste0(regression$SUB_GENE, ":", regression$GENE)
regression$pair_pro_rev_cancer <- paste0(regression$pair_pro_rev, ":", regression$Cancer)
regression$pair_cancer <- paste0(regression$pair, ":", regression$Cancer)

# add mutational impact ---------------------------------------------------
## columns: GENE.num_mut, SUB_GENE.num_mut, GENE.p, SUB_GENE.p, GENE.fdr, SUB_GENE.fdr
### add enzyme mutation impact
mut_impact_tmp <- mut_impact_tab[mut_impact_tab$pair_pro_cancer %in% regression$pair_pro_cancer,]
### take the most significant result, either protein or phosphoprotein
mut_impact_tmp <- mut_impact_tmp[order(mut_impact_tmp$fdr, mut_impact_tmp$pair_pro_cancer),]
mut_impact_tmp <- mut_impact_tmp[!duplicated(mut_impact_tmp$pair_pro_cancer),]
regression <- merge(regression, mut_impact_tmp[, c("pair_pro_cancer", "p", "fdr", "num")], 
                    by = c("pair_pro_cancer"), all.x = T)
### add substrate mutation impact
mut_impact_tmp <- mut_impact_tab[mut_impact_tab$pair_pro_cancer %in% regression$pair_pro_rev_cancer,]
### take the most significant result, either protein or phosphoprotein
mut_impact_tmp <- mut_impact_tmp[order(mut_impact_tmp$p, mut_impact_tmp$pair_pro_cancer),]
mut_impact_tmp <- mut_impact_tmp[!duplicated(mut_impact_tmp$pair_pro_cancer),]
regression <- merge(regression, mut_impact_tmp[, c("pair_pro_cancer", "p", "fdr", "num")], 
                    by.y = c("pair_pro_cancer"), by.x = c("pair_pro_rev_cancer"), 
                    suffixes = c(".GENE", ".SUB_GENE"),
                    all.x = T)
regression %>% nrow()


# annotate with cancer type specificity -----------------------------------
size <- 83
regression$is_detected_in_all_downsized <- F
regression$regulated_uniq <- F

for (enzyme_type in c("kinase", "phosphatase")) {
  ## input the regression table for only pairs detectable in all cancers
  file_path_tmp <- paste0("./cptac2p/analysis_results/phospho_network/regression/tables/generate_regression_regulated_uniq_marked/", "regression_", "size", size, "_FDR0.05_detected_in_", paste0(cancers2process, collapse = "_"),".txt")
  tab_tmp <- fread(input = file_path_tmp, data.table = F)
  tab_tmp <- adjust_regression_by_nonNA(regression = tab_tmp, reg_nonNA = 20, reg_sig = reg_sig)

  tab_tmp$pair_cancer <- paste0(tab_tmp$pair, ":", tab_tmp$Cancer)
  regression$is_detected_in_all_downsized[regression$pair %in% tab_tmp$pair] <- T
  
  regression$regulated_uniq[regression$pair_cancer %in% tab_tmp$pair_cancer[tab_tmp$regulated_uniq]] <- T
}

fn <- paste0(makeOutDir(resultD = resultD), 
             "regression", "_", "cptac2p_cptac3", "_", "tumor", "_reg_nonNA", reg_nonNA, "_mut_impact_cancer_specificity_annotated", ".txt")
write.table(x = regression, file = fn, append = F, sep = "\t", row.names = F)

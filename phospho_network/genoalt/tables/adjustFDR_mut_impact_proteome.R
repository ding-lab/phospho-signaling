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

# input mutational impact on proteome -------------------------------------
mut_impact_tab <- fread(input = paste0(ppnD, "genoalt/tables/test_mut_impact_proteome/mut_impact_proteome_RNA_cptac2p_cptac3_tab.txt"), data.table = F, sep = "\t")

## adjust for each cancer type
mut_impact_tab$fdr <- FDR_by_id_columns(p_vector = mut_impact_tab$p, id_columns = c("SUB_GENE.is_substrate",  "SUB_GENE.is_kinase", "SUB_GENE.is_phosphatase", "SUB_GENE.is_complex", "SELF", "affected_exp_type", "variant_class", "cancer"), df = mut_impact_tab)
mut_impact_tab$pair_pro_cancer <- paste0(mut_impact_tab$pair_pro, ":", mut_impact_tab$cancer)

write.table(x = mut_impact_tab, file = paste0(makeOutDir(resultD = resultD), "mut_impact_proteome_RNA_cptac2p_cptac3_tab_wFDR.txt"), quote = F, sep = "\t", row.names = F)
# Yige Wu @ WashU 2018 Feb
# overlap differentially phosphorylated sites and trans-regulated phosphosites, find common kinases


# input -------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
resultDnow <- makeOutDir()

# overlap DE phosphosites and trans-regulated phosphosites ----------------
Pho.diffexp3can <- fread(paste0(resultD, "diffexp/tables/differential_expression/", "Pho.diffexp3can.txt"),
                         data.table = F)
Pho.diffexp3can <- cbind(Pho.diffexp3can, formatPhosphosite(phosphosite_vector = Pho.diffexp3can$Phosphosite, Pho.diffexp3can$Gene))
Pho.diffexp3can <- Pho.diffexp3can[, !(colnames(Pho.diffexp3can) %in% c("Gene", "Phosphosite"))]
for (protein in c("kinase")) {
  tn = paste0(resultD,"regression/tables/",protein,"_substrate_regression_cptac2p_3can.txt")
  table_3can <- fread(tn, data.table = F)
  table_3can_diffexp <- merge(table_3can, Pho.diffexp3can, by = c("SUBSTRATE", "SUB_MOD_RSD", "Cancer", "transcript"), all.x = T)
  write.table(table_3can_diffexp, row.names = F, quote=F, sep = '\t', file=paste0(resultDnow, "/", protein,"_substrate_regression_overlap_cptac2p_3can_diffPhosphositeAnnotated.txt"))
  write.table(table_3can_diffexp[!is.na(table_3can_diffexp$`Fold Change`),], row.names = F, quote=F, sep = '\t', file=paste0(resultDnow, "/", protein,"_substrate_regression_overlap_cptac2p_3can_diffPhosphositeFiltered.txt"))
}


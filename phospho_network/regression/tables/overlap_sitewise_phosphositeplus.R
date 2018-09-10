# Yige Wu @ WashU 2018 Jan
# examine the how many exact phosphosites from Phosphositeplus our result validated
# good to use their annotation(experimental etc.) to say something about these sites

# input -------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
resultDnow <- makeOutDir()
sig = 0.05

for (protein in c("kinase")) {
  tn = paste0(resultD,"regression/tables/",protein,"_substrate_regression_cptac2p_3can.txt")
  table_3can <- fread(tn, data.table = F)
  table_3can <- markSigKS(table_3can, sig_thres = sig, protein_type = protein)
  
  examined_summary <- data.frame(table(table_3can[,c("SELF", "Cancer")]))
  colnames(examined_summary)[which(colnames(examined_summary) == "Freq")] <- "examined_Freq"
  phosphositeplus_vad <- table_3can[table_3can$fdr_sig & table_3can$coef_sig,]
  validated_summary <- data.frame(table(phosphositeplus_vad[,c("SELF", "Cancer")]))
  colnames(validated_summary)[which(colnames(validated_summary) == "Freq")] <- "validated_Freq"
  validated_summary <- merge(validated_summary, examined_summary)
  tn = paste0(resultDnow, protein,"_substrate_regression_validate.txt")
  write.table(validated_summary, file=tn, quote=F, sep = '\t', row.names = FALSE)
  
  ## summarize the ones overlapping phosphositeplus
  k_s_table <- load_ks_table(protein)
  phosphositeplus_anno <- merge(table_3can, 
                               k_s_table[,c("GENE", "SUB_GENE", "SUB_MOD_RSD", 
                                            "DOMAIN", "IN_VIVO_RXN", "IN_VITRO_RXN", "CST_CAT.")], 
                               by.x = c("KINASE", "SUBSTRATE", "SUB_MOD_RSD"),
                               by.y = c("GENE", "SUB_GENE", "SUB_MOD_RSD"),
                               all.x = T)
  phosphositeplus_anno <- phosphositeplus_anno[!is.na(phosphositeplus_anno$IN_VIVO_RXN),]
  phosphositeplus_vad <- phosphositeplus_anno[phosphositeplus_anno$fdr_sig & phosphositeplus_anno$coef_sig,]
  
  examined_summary <- data.frame(table(phosphositeplus_anno[,c("SELF", "Cancer")]))
  colnames(examined_summary)[which(colnames(examined_summary) == "Freq")] <- "examined_Freq"
  validated_summary <- data.frame(table(phosphositeplus_vad[,c("SELF", "Cancer")]))
  colnames(validated_summary)[which(colnames(examined_summary) == "Freq")] <- "validated_Freq"
  validated_summary <- merge(validated_summary, examined_summary)
  tn = paste0(resultDnow, protein,"_substrate_regression_validate_sitewise_PhosphositePlus.txt")
  write.table(validated_summary, file=tn, quote=F, sep = '\t', row.names = FALSE)
  
  examined_summary <- data.frame(table(phosphositeplus_anno[,c("SELF", "Cancer", "IN_VIVO_RXN", "IN_VITRO_RXN")]))
  colnames(examined_summary)[which(colnames(examined_summary) == "Freq")] <- "examined_Freq"
  validated_summary <- data.frame(table(phosphositeplus_vad[,c("SELF", "Cancer", "IN_VIVO_RXN", "IN_VITRO_RXN")]))
  colnames(validated_summary)[which(colnames(validated_summary) == "Freq")] <- "validated_Freq"
  validated_summary <- merge(validated_summary, examined_summary)
  validated_summary <- validated_summary[validated_summary$examined_Freq > 0,]
  tn = paste0(resultDnow, protein,"_substrate_regression_validate_sitewise_PhosphositePlus_in_vivo_vitro.txt")
  write.table(validated_summary, file=tn, quote=F, sep = '\t', row.names = FALSE)
  
  for (self in c('cis', 'trans')) {
    phosphositeplus_vad_self <-phosphositeplus_vad[phosphositeplus_vad$SELF == self,]
    tn = paste0(resultDnow, protein,"_substrate_regression_", self, "_validate_sitewise_PhosphositePlus.txt")
    write.table(phosphositeplus_vad_self, file=tn, quote=F, sep = '\t', row.names = FALSE)
  }
}
    
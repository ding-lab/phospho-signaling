# Yige Wu @ WashU 2016 Dec
# calculate the validation statistics for regression results (coefficient > 0, FDR < alpha)

# choose kinase/phosphotase, cancer , significance level -----------------------------------------------
protein <- "kinase"

# input -------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
source(paste0(baseD, "cptac2p_analysis/pan3can_aes.R")) # aes for general purposes; it should be one directory out of the working directory
tn = paste0(resultD,"regression/tables/",protein,"_substrate_regression_cptac2p_3can.txt")
table_3can <- fread(tn, data.table = F)
table_3can <- markVadRetro(regression = table_3can, sig_thres = sig, protein_type = protein)
k_s_table <- load_ks_table(protein)
resultDnow <- makeOutDir()


# process while looping around cancers ------------------------------------
for (cancer in c("BRCA","OV","CO")) {
  for (self in c("cis", "trans")) {
    table1 <- table_3can[table_3can$Cancer== cancer & table_3can$SELF == self,]
    unique_kinase <- as.vector(unique(table1$KINASE))
    # initiate for gene-level validation
    valid_stat <- data.frame(kinase = unique_kinase)
    rownames(valid_stat) <- unique_kinase
    for( kinase in unique_kinase) {
      temp <- table1[table1$KINASE==kinase,]
      valid_stat[kinase,"valid_count"] <- length(which(temp$fdr_sig & temp$coef_sig))
      valid_stat[kinase,"all_count"] <- nrow(temp)
      valid_stat[kinase,"vad_retro_sitewise_count"] <- length(which(temp$fdr_sig & temp$coef_sig & temp$vad_retro_sitewise))
      valid_stat[kinase,"vad_retro_pairwise_count"] <- length(which(temp$fdr_sig & temp$coef_sig & temp$vad_retro_pairwise))
    }
    valid_stat$valid_ratio <- valid_stat$valid_count/valid_stat$all_count
    tn2 = paste(resultDnow,"cptac2p_",cancer,"_",protein,"_", self,"_regression_validation_statistics.txt", sep="")
    write.table(valid_stat, file=tn2, quote=F, sep = '\t', row.names = FALSE)
  }
}





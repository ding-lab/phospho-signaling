# Yige Wu @ WashU 2018 Jul
# make table for phosphosites shared between top regulating kinases and phosphotase
# table with regulating relationships between kinases and phosphotases

# Yige Wu @ WashU 2018 Apr
# check phosphosites being regulated by more than 2 kinases
# check phosphosites being regulated by both kinase and phosphatase


# source ------------------------------------------------------------------
source("./cptac2p_analysis/phospho_network/phospho_network_shared.R")
library(dplyr)

# inputs ------------------------------------------------------------------
top_num <- 20
sig_pk <- 0.1
sig_pp <- 0.2

# reformat --------------------------------------------------------------------
for (cancer in c("BRCA", "OV", "CO")) {
  sup_tab <- fread(paste0(ppnD, "kinase_activity/tables/fisher_es_pairs/", 
                          cancer, "_KSEA_enzyme_diffexp_substrate.txt"), data.table = F)
  sup_tab <- data.frame(sup_tab)
  sup_tab <- sup_tab[sup_tab$aa_reported,]
  sup_tab$pair <- paste0(sup_tab$GENE, ":", sup_tab$SUB_GENE, ":", sup_tab$SUB_MOD_RSD)
  sup_tab$site <- paste0(sup_tab$SUB_GENE, ":", sup_tab$SUB_MOD_RSD)
  trans_tab <- sup_tab[!is.na(sup_tab$FDR_pho_kin.aa_reported),]
  trans_pk_tab <- trans_tab[trans_tab$enzyme_type == "kinase",]
  trans_pp_tab <- trans_tab[trans_tab$enzyme_type == "phosphotase",]
  
  trans_pk_sig_tab <- trans_pk_tab[trans_pk_tab$FDR_pho_kin.aa_reported < sig_pk & trans_pk_tab$coef_pho_kin > 0,]
  trans_pp_sig_tab <- trans_pp_tab[trans_pp_tab$FDR_pho_kin.aa_reported < sig_pp & trans_pp_tab$coef_pho_kin < 0,]
  
  trans_pk_sig_multi <- data.frame(table(trans_pk_sig_tab[, c("site")]))
  trans_pk_sig_multi <- trans_pk_sig_multi[trans_pk_sig_multi$Freq > 1,]
  
  ## choose top regulating kinases
  trans_pk_count <- data.frame(table(trans_pk_sig_tab[, "GENE"]))
  trans_pk_count <- trans_pk_count[order(trans_pk_count$Freq, decreasing = T),]
  top_trans_pks <- as.vector(trans_pk_count$Var1[1:top_num])
  
  count <- 0
  kinase1 <- NULL
  kinase2 <- NULL
  shared_num <- NULL
  regulated <- NULL
  for (i in 1:(top_num-1)) {
    ki <- top_trans_pks[i]
    for (j in (i+1):top_num) {
      kj<- top_trans_pks[j]
      count <- count + 1
      kinase1[count] <- ki
      kinase2[count] <- kj
      shared_sites <- intersect(trans_pk_sig_tab$site[trans_pk_sig_tab$GENE == ki], trans_pk_sig_tab$site[trans_pk_sig_tab$GENE == kj])
      shared_num[count] <- length(shared_sites)
      if (length(which((trans_pk_sig_tab$GENE == ki & trans_pk_sig_tab$SUB_GENE == kj) | (trans_pk_sig_tab$GENE == kj & trans_pk_sig_tab$SUB_GENE == ki))) > 0){
        regulated[count] <- TRUE
      } else {
        regulated[count] <- FALSE
      }
    }
  }
  
  
  ## annotate each kinase with how many phosphosites they regulate
  shared_tab <- data.frame(kinase1, kinase2, shared_num, regulated)
  colnames(trans_pk_count) <- c("kinase", "sig_num")
  shared_tab <- merge(shared_tab, trans_pk_count, by.x = c("kinase1"), by.y = c("kinase"), all.x = T)
  colnames(shared_tab)[ncol(shared_tab)] <- "sig_num1"
  shared_tab <- merge(shared_tab, trans_pk_count, by.x = c("kinase2"), by.y = c("kinase"), all.x = T)
  colnames(shared_tab)[ncol(shared_tab)] <- "sig_num2"
  
  ## annotate whether one kinase is regulating another
  
  
  write.table(x = shared_tab, file = paste0(makeOutDir(resultD = resultD), cancer, "_shared_site_sig", sig_pk, ".txt"), 
              quote = F, row.names = F, sep = "\t")
}

# Yige Wu @ WashU 2018 Mar
# annotate table f hot and cold kinase-substrate pairs by fisher's test

# source ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')

# inputs ------------------------------------------------------------------
## outlier expression is 1*IQR + median
iqr_thres <- 1

# fisher test the number of samples with outlier enzyme phospho and substrate phospho --------

for (cancer in c("CO")) {
# for (cancer in c("BRCA", "OV")) {
  ksea_diffexp <- fread(input = paste0(ppnD, "kinase_activity/tables/integrate_enzyme_KSEA_diffexp_sub_FC_plus_cistrans_regression_plus_genoalt/", 
                                       cancer, "_KSEA_enzyme_diffexp_substrate.txt"), data.table = F)
  kinases <- as.vector(ksea_diffexp$GENE)
  substrates <- as.vector(ksea_diffexp$SUB_GENE)
  rsds <- as.vector(ksea_diffexp$SUB_MOD_RSD)
  
  pho_data <- fread(input = paste(cptac_sharedD, cancer,"/",prefix[cancer], "_PHO_formatted_normalized_noControl.txt",sep=""), data.table = F)
  pho_gdata <- fread(input = paste(cptac_sharedD, cancer,"/",prefix[cancer], "_PHO_normalized_collapsed_noControl.txt",sep=""), data.table = F)
  pho_rsd_split <- formatPhosphosite(pho_data$Phosphosite, pho_data$Gene)
  samples <- colnames(pho_gdata)[!(colnames(pho_gdata) %in% c("Gene"))]
  
  ## initiate 4 vectors for kinase outlier/not * substrate outlier/not
  kin_up_sub_up <- vector(mode = "numeric", length = nrow(ksea_diffexp))
  kin_up_sub_norm <- kin_up_sub_up
  kin_norm_sub_up <- kin_up_sub_up
  kin_norm_sub_norm <- kin_up_sub_up
  fisher_pvalues <- kin_up_sub_up
  
  for (i in which(!is.na(ksea_diffexp$FDR_pho_kin.aa_reported))) {
    ## get kinase
    kinase <- kinases[i]
    ## get substrate
    substrate <- substrates[i]
    ## get substrate phosphosite
    rsd <- rsds[i]
    
    pho_kinase_g <- pho_gdata[pho_gdata$Gene == kinase,samples]
    pho_sub <- pho_data[pho_rsd_split$SUBSTRATE == substrate & pho_rsd_split$SUB_MOD_RSD == rsd, samples]
    pho_sub <- pho_sub[1,]
    
    ## bind the phospho tgt
    pho_tab <- data.frame(t(rbind(pho_sub,pho_kinase_g)))
    colnames(pho_tab) <- c("pho_sub","pho_kinase_g")
    pho_tab <- pho_tab[complete.cases(pho_tab),]
    
    pho_kin_out <- iqr_thres*IQR(pho_kinase_g, na.rm = T) + median(t(pho_kinase_g), na.rm = T)
    pho_sub_out <- iqr_thres*IQR(pho_sub, na.rm = T) + median(t(pho_sub), na.rm = T)
    
    out_tab <- table(pho_tab$pho_kinase_g > pho_kin_out, pho_tab$pho_sub > pho_sub_out)
    if (nrow(out_tab) == 2 && ncol(out_tab) == 2) {
      stat <- fisher.test(out_tab)
      ## store the values
      fisher_pvalues[i] <- stat$p.value
      kin_up_sub_up[i] <- out_tab["TRUE", "TRUE"]
      kin_up_sub_norm[i] <- out_tab["TRUE", "FALSE"]
      kin_norm_sub_up[i] <- out_tab["FALSE", "TRUE"]
      kin_norm_sub_norm[i] <- out_tab["FALSE", "FALSE"]
    }
  }
  ksea_diffexp <- cbind(ksea_diffexp, fisher_pvalues, kin_up_sub_up, kin_up_sub_norm, kin_norm_sub_up, kin_norm_sub_norm)
  # ksea_diffexp <- cbind(ksea_diffexp, sup_tab[,c('fisher_pvalues', 'kin_up_sub_up', 'kin_up_sub_norm', 'kin_norm_sub_up', 'kin_norm_sub_norm')])
  write.table(x = ksea_diffexp,
              file = paste0(makeOutDir(resultD = resultD),
                            cancer, "_KSEA_enzyme_diffexp_substrate.txt"),
              row.names = F, col.names = T, quote = F, sep = "\t")
}



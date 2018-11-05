# Yige Wu @ WashU 2018 Apr
# venn diagram showing whether regulated enzyme-substrate pairs only in breast are detected in colorectal and ovarian cancer

# source ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')


# set variables -----------------------------------------------------------
# cptac_phase2process <- "cptac2p_3can"
# cancers2process <- cancers_sort
cptac_phase2process <- "cptac3"
cancers2process <- "UCEC"

# split table for smaller supplementary tables ----------------------------
for (enzyme_type in c("kinase", "phosphatase")) {
  # sup_cans_tab_en <- fread(input = paste0("./cptac2p/analysis_results/phospho_network/regression/tables/regression_manualAdded_protein_level/", 
  #                                         enzyme_type, "_substrate_regression_", cptac_phase2process, "_tumor.txt"), data.table = F)
  sup_cans_tab_en <- fread(input = paste0("./cptac2p/analysis_results/phospho_network/regression/tables/regression_omnipath&newworkin&depod&signor&manual_protein_level/", 
                                          enzyme_type, "_substrate_regression_", cptac_phase2process, "_tumor.txt"), data.table = F)
  
  for (cancer in cancers2process) {
    for (self in c("trans", "cis")) {
      tab2w <- sup_cans_tab_en[sup_cans_tab_en$Cancer == cancer & sup_cans_tab_en$SELF == self,]
      tab2w$self <- NULL
      tab2w$pair <- NULL
      
      write.table(x = tab2w, file = paste0(makeOutDir(resultD = resultD), enzyme_type, "_substrate_regression_", cancer, "_", self , "_tumor.txt"), quote = F, row.names = F, sep = "\t")
      write.table(x = tab2w, file = paste0(makeOutDir(resultD = resultD), enzyme_type, "_substrate_regression_", cancer, "_", self , "_tumor.csv"), quote = F, row.names = F, sep = ",")
      
      if (self == "trans") {
        tab2w <- tab2w[tab2w$FDR_pho_kin < 0.1,]
        write.table(x = tab2w, file = paste0(makeOutDir(resultD = resultD), enzyme_type, "_substrate_regression_", cancer, "_", self , "_FDR0.1", "_tumor.csv"), quote = F, row.names = F, sep = ",")
        
        tab2w <- tab2w[tab2w$FDR_pho_kin < 0.05,]
        write.table(x = tab2w, file = paste0(makeOutDir(resultD = resultD), enzyme_type, "_substrate_regression_", cancer, "_", self , "_FDR0.05", "_tumor.csv"), quote = F, row.names = F, sep = ",")
    
      }
    }
  }
}


# inputs ------------------------------------------------------------------
## input enzyme-substrate pairs examined
regulated_ratio_across_thres <- NULL
for (reg_nonNA in seq(from = 5, to = 30, by = 5)) {
  for (enzyme_type in c("kinase", "phosphatase")) {
    # sup_cans_tab_en <- fread(input = paste0("./cptac2p/analysis_results/phospho_network/regression/tables/regression_manualAdded_protein_level/", 
    #                                         enzyme_type, "_substrate_regression_", cptac_phase2process, "_tumor.txt"), data.table = F)
    sup_cans_tab_en <- fread(input = paste0("./cptac2p/analysis_results/phospho_network/regression/tables/regression_omnipath&newworkin&depod&signor&manual_protein_level/", 
                                            enzyme_type, "_substrate_regression_", cptac_phase2process, "_tumor.txt"), data.table = F)
    sup_cans_tab_en <- sup_cans_tab_en[sup_cans_tab_en$Size >= reg_nonNA,]
    name = c("pro_kin","pro_sub","pho_kin")
    
    ## adjust p-values to FDR
    for (cancer in cancers2process) {
      for(self in c(TRUE,FALSE)) {
        for(coln in name) {#adjust pvalues for each variable
          row <- (sup_cans_tab_en$self==self) & (sup_cans_tab_en$Cancer==cancer)
          sup_cans_tab_en[row,paste("FDR_",coln,sep = "")] <-p.adjust(sup_cans_tab_en[row,paste("P_",coln,sep = "")],method = "fdr")
        }
      }
    }
    
    sup_cans_tab_en$GENE <- as.vector(sup_cans_tab_en$KINASE)
    sup_cans_tab_en$SUB_GENE <- as.vector(sup_cans_tab_en$SUBSTRATE)
    write.table(x = sup_cans_tab_en, file = paste0(makeOutDir(resultD = resultD), 
                                                   enzyme_type, "_substrate_regression_", cptac_phase2process, "_tumor",
                                                   "_reg_nonNA", reg_nonNA, ".txt"), row.names = F, quote = F, sep = "\t")
    
    ## get common detected phosphosites
    num_nonNA_common <- reg_nonNA
    pho_sites_cans <- list()
    pho_sites_all <- NULL
    for (cancer in cancers_sort) {
      pho_data <- loadPhosphositeNormalizedTumor(cancer = cancer)
      samples <- colnames(pho_data)[!(colnames(pho_data) %in% c("Gene", "Phosphosite"))]
      pho_data <- pho_data[rowSums(!is.na(pho_data[, samples])) >= num_nonNA_common,]
      pho_head <- formatPhosphosite(phosphosite_vector = pho_data$Phosphosite, gene_vector = pho_data$Gene)
      pho_sites <- unique(pho_head[, c("SUBSTRATE", "SUB_MOD_RSD")])
      pho_sites$site <- paste0(pho_sites$SUBSTRATE, "_", pho_sites$SUB_MOD_RSD)
      pho_sites_cans[[cancer]] <- unique(as.vector(pho_sites$site))
      pho_sites_all <- unique(c(pho_sites_all, as.vector(pho_sites$site)))
    }
    dat <- data.frame(site = pho_sites_all)
    for (cancer in cancers_sort) {
      dat[, cancer] <- FALSE
      dat[dat$site %in% pho_sites_cans[[cancer]], cancer]  <- TRUE
    }
    pho_site_common <- dat[dat$BRCA & dat$OV & dat$CO,]

    ## filter by phosphosite
    sup_cans_tab_en$site <- paste0(sup_cans_tab_en$SUB_GENE, "_", sup_cans_tab_en$SUB_MOD_RSD)
    sup_cans_tab_en <- sup_cans_tab_en[sup_cans_tab_en$site %in% pho_site_common$site,]
    
    write.table(x = sup_cans_tab_en, file = paste0(makeOutDir(resultD = resultD), 
                                                   enzyme_type, "_substrate_regression_", cptac_phase2process, "_tumor",
                                                   "_reg_nonNA", reg_nonNA, "_commonsite_nonNA", num_nonNA_common,".txt"), row.names = F, quote = F, sep = "\t")
  }
}



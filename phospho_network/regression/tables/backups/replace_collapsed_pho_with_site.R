# Yige Wu @ WashU 2018 Oct
# replace the collapsed phosphorylation with particular site on the kinase, see if correlated pair stay significant

# source ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')

# set variables -----------------------------------------------------------
cptac_phase2process <- "cptac2p_3can"
cancers2process <- cancers_sort
reg_nonNA2test <- c(20)
fdr_thres2process <- c(0.05)
sample_type <- "tumor"

# inputs ------------------------------------------------------------------
## input the pairs
tab_pairs <- data.frame(GENE = c("RPS6KB1", "MAPK1", "MAPK3"),
                            SUB_GENE = c("GSK3B", "RPS6KB1", "RPS6KB1"),
                            SUB_MOD_RSD = c("S9", "T421S424", "T421S424"))


## input the already analyzed regression model
tab_pairs_replaced <- NULL
for (reg_nonNA in reg_nonNA2test) {
  for (enzyme_type in c("kinase")) {
    regression <- fread(input = paste0(ppnD, "regression/tables/", "change_regression_nonNA/", 
                                       enzyme_type, "_substrate_regression_", cptac_phase2process, "_tumor",
                                       "_reg_nonNA", reg_nonNA, ".txt"), data.table = F)
    regression$site <- paste0(regression$SUB_GENE, ":", regression$SUB_MOD_RSD)

    least_samples <- reg_nonNA
    ## calculate the per pair per cancer per kinase site
    for (cancer in cancers2process) {
      ## input protein abundance
      if (cptac_phase2process == "cptac2p_3can") {
        pro_data <- loadProteinNormalizedTumor(cancer = cancer)
        pho_data <- loadPhosphositeNormalizedTumor(cancer = cancer)
        pho_gdata <- loadPhosphoproteinNormalizedTumor(cancer = cancer)
        pho_rsd_split <- formatPhosphosite(pho_data$Phosphosite, pho_data$Gene)
        
        transcripts <- as.vector(pho_rsd_split$transcript)
        sub_mod_rsd <- as.vector(pho_rsd_split$SUB_MOD_RSD)
        samples <- colnames(pro_data)[!(colnames(pro_data) %in% c("Gene"))]
      }

      for (i in 1:nrow(tab_pairs)) {
        kinase <- as.vector(tab_pairs[i, "GENE"])
        substrate <- as.vector(tab_pairs[i, "SUB_GENE"])
        rsd <- as.vector(tab_pairs[i, "SUB_MOD_RSD"])
        pair <- paste0(kinase, ":", substrate, ":", rsd)
        pho_sub <- pho_data[pho_rsd_split$SUBSTRATE == substrate & pho_rsd_split$SUB_MOD_RSD == rsd,samples]
        
        # go through all the phosphosites
        k_pho_table <- which(pho_rsd_split$SUBSTRATE == kinase)
        
        for (j in k_pho_table) {
          # find phosphorylation level
          pho_kinase_g <- pho_data[j, samples]
          
          # find substrate expressio level and normaize
          pro_sub <- pro_data[pro_data$Gene == substrate,samples]
          
          if(nrow(pro_sub) != 0 & nrow(pho_sub) != 0){
            #prepare regression data for model2
            data2 <- data.frame(t(rbind(pho_sub,pro_sub,pho_kinase_g)))
            data_complete <- data2[complete.cases(data2),]
            size <- nrow(data_complete)
            if(size > least_samples ){
              colnames(data2) <- c("pho_sub","pro_sub","pho_kinase_g")
              fit2 <- lm(pho_sub ~ pro_sub + pho_kinase_g, data = data2)
              
              ## initiate
              table_trans <- regression[regression$pair == pair & regression$Cancer == cancer,]
              table_trans$Size <- size
              table_trans$P_pro_sub <- c(coef(summary(fit2))[2,4])
              table_trans$coef_pro_sub <- fit2$coefficients[2]
              table_trans$P_pho_kin <- c(coef(summary(fit2))[3,4])
              table_trans$coef_pho_kin <- fit2$coefficients[3]
              
              regression_new <- rbind(table_trans, regression[regression$pair != pair & regression$Cancer == cancer & regression$SELF == "trans",])
              
              ## adjust FDR
              name = c("pro_kin","pro_sub","pho_kin")
              for(coln in name) {
                regression_new[, paste("FDR_",coln,sep = "")] <- p.adjust(regression_new[, paste("P_",coln,sep = "")], method = "fdr")
              }
              ## add the new value
              tab2add <- regression_new[regression_new$pair == pair,]
              tab2add$ENZ_MOD_RSD <- pho_rsd_split[j, "SUB_MOD_RSD"]
              tab_pairs_replaced <- rbind(tab_pairs_replaced, tab2add)
              
              # if (cancer == "BRCA") {
              #   stop("")
              # }
              
            }
          }
        }
      }
    }
  }
}

## write table

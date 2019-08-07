# Yige Wu @ WashU 2018 Jan
# look at correlations of kinase and downstream substrates phosphorylation status
# pho_sub~pro_sub+pho_kin(callapsed)

# source-------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')


# set variables -----------------------------------------------------------
least_partIDs <- 5# least number of partIDs with complete data for each model
inputnames <- c("noControl", "Control")
names(inputnames) <- c("tumor", "normal")
cptac_phase2process <- "cptac2p_3can"
# cptac_phase2process <- "cptac3"
cancers2process <- cancers_sort
# cancers2process <- "UCEC"

# inputs ------------------------------------------------------------------
ptms_site_pairs_sup <- read_csv(paste0(ppnD, "compile_enzyme_substrate/tables/compile_omnipath_networkin_depod_signor_manual/Omnipath_NetworKin_DEPOD_SignorNotSiteMapped_manualAdded.csv"))
ptms_site_pairs_sup <- ptms_site_pairs_sup[ptms_site_pairs_sup$Source != "NetKIN" | (ptms_site_pairs_sup$Source == "NetKIN" & ptms_site_pairs_sup$networkin_score >= 5),]
ptms_site_pairs_sup <- data.frame(ptms_site_pairs_sup)
ptms_site_pairs_sup$pair <- paste0(ptms_site_pairs_sup$pair_pro, ":", ptms_site_pairs_sup$SUB_MOD_RSD)

# loop -----------------------------------------------------------------
for (sample_type in c("tumor")) {
  for( enzyme_type in c("kinase")) {
  # for( enzyme_type in c("phosphatase", "kinase")) {
    k_s_table <- ptms_site_pairs_sup[ptms_site_pairs_sup$enzyme_type == enzyme_type,]
    genes4cis <- unique(c(as.vector(k_s_table$GENE), as.vector(k_s_table$SUB_GENE)))

    table_3can <- NULL
    for (cancer in cancers2process) {
    # for (cancer in "BRCA") {
        
      ## input protein abundance
      pro_data <- loadParseProteomicsData(data_type = "PRO", cancer = cancer, sample_type = sample_type)

      ## input phosphosite abundance
      pho_data <- loadParseProteomicsData(data_type = "PHO", cancer = cancer, sample_type = sample_type)
      
      ## input phosphoprotein abundance
      pho_gdata <- loadParseProteomicsData(data_type = "PHO_collapsed", cancer = cancer, sample_type = sample_type)
      
      ## input somatic mutation matrix
      mut_mat <- fread(input = paste0(ppnD, "genoalt/tables/generate_somatic_mutation_matrix/", cancer, "_somatic_mutation.txt"), data.table = F)
      rownames(mut_mat) <- as.vector(mut_mat$Hugo_Symbol)
      mut_count <- rowSums(x = ((mut_mat[, -1] != "") & (mut_mat[, -1] != "Silent")))
      names(mut_count) <- as.vector(mut_mat$Hugo_Symbol)

      ## the participant IDS
      partIDs <- colnames(pro_data)[-1]
      partIDs <- intersect(partIDs, colnames(mut_mat))
      
      # initiate ----------------------------------------------------------------
      # calculate the length of trans table
      k_s_siteswdata <- merge(unique(rbind(k_s_table[, c("GENE", "SUB_GENE")], data.frame(GENE = genes4cis, SUB_GENE = genes4cis))), 
                         pho_data, by.x = c("SUB_GENE"), by.y = c("Gene"))
      k_s_siteswdata <- k_s_siteswdata[k_s_siteswdata$GENE %in% names(mut_count >= least_partIDs) & k_s_siteswdata$SUB_GENE %in% pro_data$Gene,]
      
      mut_mat2overlap <- as.matrix(mut_mat[as.vector(k_s_siteswdata$GENE), partIDs])
      mut_mat2overlap <- (mut_mat2overlap != "" & mut_mat2overlap != "Silent")
      rownames(mut_mat2overlap) <- NULL
      pho_data2overlap <- as.matrix(k_s_siteswdata[, partIDs])
      pho_data2overlap <- as.matrix(!is.na(pho_data2overlap))
      sizes_mut <- rowSums(mut_mat2overlap & pho_data2overlap)
      
      k_s_sites <- k_s_siteswdata[sizes_mut > 0, c("GENE", "SUB_GENE", "Phosphosite", "Peptide_ID")]
      if (nrow(k_s_sites) != nrow(unique(k_s_sites))) {
        stop("check k_s_sites!")
      }
      ntrans <- nrow(k_s_sites)
      
      # looping over kinases for trans pairs -----------------------------------------------------------------
      # initiating the table for trans
      vec_char <- vector(mode = "character", length = ntrans)
      vec_num <- vector(mode = "numeric", length = ntrans) + NA
      KINASE <- vec_char;SUBSTRATE <- vec_char; SUB_MOD_RSD <- vec_char;
      FDR_pro_kin <- vec_num;FDR_pro_sub <- vec_num;FDR_mut_kin <- vec_num;
      coef_pro_kin <- vec_num;coef_pro_sub <- vec_num;coef_mut_kin <- vec_num;
      Cancer <- vec_char; Peptide_ID <- vec_char; model <- vec_char;
      Size <- vec_num;
      P_pro_kin <- vec_num;P_pro_sub <- vec_num;P_mut_kin <- vec_num;
      sd_pro_kin <- vec_num;sd_pro_sub <- vec_num;sd_pho_sub <- vec_num; num_mut_kin <- vec_num;
      
      count <- 0
      for (kinase in unique(k_s_sites$GENE)){
        mut_kin <- mut_mat[mut_mat$Hugo_Symbol == kinase, partIDs]
        
        if ( nrow(mut_kin) > 0 ){
          k_sub <- unique(k_s_sites$SUB_GENE[k_s_sites$GENE == kinase])
          
          for (substrate in k_sub){# for each substrate for one kinase
            
            for (j in which(k_s_sites$SUB_GENE == substrate & k_s_sites$GENE == kinase)) {
              rsd <- k_s_sites$Phosphosite[j]
              peptide_id <- k_s_sites$Peptide_ID[j]
              i <- which(pho_data$Phosphosite == rsd & pho_data$Peptide_ID == peptide_id)
              
              # find phosphorylation level
              pho_sub <- pho_data[i,partIDs]
              
              # find substrate expressio level and normaize
              pro_sub <- pro_data[pro_data$Gene == substrate, partIDs]
              
              if(nrow(pro_sub) != 0){
                #prepare regression data for model2
                data2 <- data.frame(t(rbind(pho_sub,pro_sub,mut_kin)))
                colnames(data2) <- c("pho_sub","pro_sub","mut_kin")
                data2$pho_sub <- as.numeric(as.vector(data2$pho_sub))
                data2$pro_sub <- as.numeric(as.vector(data2$pro_sub))
                data2$mut_kin <- as.numeric(!(data2$mut_kin %in% c("Silent", "")))
                data_complete <- data2[complete.cases(data2),]
                size <- nrow(data_complete)
                if(size > least_partIDs & length(unique(data_complete$mut_kin)) > 1){
                  fit2 <- lm(pho_sub ~ pro_sub + mut_kin, data = data2)
                  
                  count <- count + 1
                  KINASE[count] <- kinase
                  SUBSTRATE[count] <- substrate
                  SUB_MOD_RSD[count] <- pho_data$Phosphosite[i]
                  P_pro_sub[count] <- c(coef(summary(fit2))[2,4])
                  coef_pro_sub[count] <- fit2$coefficients[2]
                  P_mut_kin[count] <- c(coef(summary(fit2))[3,4])
                  coef_mut_kin[count] <- fit2$coefficients[3]
                  Size[count] <- size
                  Peptide_ID[count] <- pho_data$Peptide_ID[i]
                  sd_pro_sub[count] <- sd(data_complete$pro_sub)
                  sd_pho_sub[count] <- sd(data_complete$pho_sub)
                  num_mut_kin[count] <- length(which(data_complete$mut_kin == 1))
                }
              }
              
              # go through all the phosphosites
              for (i in s_pho_table) {
                

              }
            }
 
          }
        }
      }
      table_trans <- data.frame(KINASE,SUBSTRATE,SUB_MOD_RSD,
                                FDR_pro_sub,FDR_mut_kin,
                                coef_pro_sub,coef_mut_kin,
                                Cancer,model,Size, Peptide_ID,
                                P_pro_sub,P_mut_kin,
                                sd_pro_kin,sd_pho_sub, num_mut_kin)
      
      # integrate table from all the models(need to repetite again for another cancer dataset) --------------------------------------------------------
      table_trans$model <- "pho_sub~pro_sub+mut_kin"
      tabletrans <- table_trans[!is.na(table_trans$P_pro_sub),]
      
      # combine table
      table <- tabletrans
      
      # mark cancer and self-regulation
      table$Cancer <- cancer
      
      table_3can <- rbind(table_3can, table)
    }
    
    # combine table from 3 cancers and adjust P-values -------------------------------------------------------
    ## combine table from BRCA and OV
    table_3can$SELF <- ifelse(test = as.vector(table_3can$KINASE) == as.vector(table_3can$SUBSTRATE), "cis", "trans")
    name = c("pro_sub","mut_kin")
    
    ## adjust p-values to FDR
    for (cancer in cancers2process) {
      for(SELF in c("cis","trans")) {
        for(coln in name) {#adjust pvalues for each variable
          row <- (table_3can$SELF==SELF) & (table_3can$Cancer==cancer)
          table_3can[row,paste("FDR_",coln,sep = "")] <-p.adjust(table_3can[row,paste("P_",coln,sep = "")],method = "fdr")
        }
      }
    }
    
    table_3can$pair <- paste(table_3can$KINASE,table_3can$SUBSTRATE,table_3can$SUB_MOD_RSD,sep = ":")
    table_3can <- merge(table_3can, unique(ptms_site_pairs_sup[, c("pair", "Source")]), all.x = T)
    
    # write out tables --------------------------------------------------------
    tn = paste0(makeOutDir(resultD = resultD), enzyme_type,"_substrate_regression_", cptac_phase2process , "_", sample_type, ".txt")
    write.table(table_3can, file=tn, quote=F, sep = '\t', row.names = FALSE)
  }
}

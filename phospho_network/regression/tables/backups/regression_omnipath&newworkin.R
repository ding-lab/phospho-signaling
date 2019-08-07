# Yige Wu @ WashU 2018 Jan
# look at correlations of kinase and downstream substrates phosphorylation status
# pho_sub~pro_sub+pho_kin(callapsed)

# edit parameters-------------------------
## choose kinase/phosphotase, significance level, outlier threshold and least sample numbe
least_samples <- 5# least number of samples with complete data for each model
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
resultDnow <- makeOutDir()
inputnames <- c("noControl", "Control")
names(inputnames) <- c("tumor", "normal")

# input -----------------------------------------------------------------
for (sample_type in c("tumor", "normal")) {
  for( protein in c("kinase", "phosphotase")) {
    k_s_table <- load_ks_table(protein)
    kinase_trans <- as.vector(unique(k_s_table$GENE[as.vector(k_s_table$GENE)!=as.vector(k_s_table$SUB_GENE)]))
    kinase_cis <- as.vector(unique(k_s_table$GENE[as.vector(k_s_table$GENE)==as.vector(k_s_table$SUB_GENE)]))
    table_3can <- NULL
    for (cancer in c("BRCA","OV","CO")) { 
      pro_data <- read.table(header=TRUE, sep="\t", check.names=F, file=paste(inputD, cancer,"/",prefix[cancer], "_PRO_formatted_normalized_", inputnames[sample_type], ".txt",sep=""))
      pho_data <- read.table(header=TRUE, sep="\t", check.names=F, file=paste(inputD, cancer,"/",prefix[cancer], "_PHO_formatted_normalized_", inputnames[sample_type], ".txt",sep=""))
      pho_gdata <- read.table(header=TRUE, sep="\t", check.names=F, file=paste(inputD, cancer,"/",prefix[cancer], "_PHO_normalized_collapsed_", inputnames[sample_type], ".txt",sep=""))
      pho_rsd_split <- formatPhosphosite(pho_data$Phosphosite, pho_data$Gene)
      
      transcripts <- as.vector(pho_rsd_split$transcript)
      sub_mod_rsd <- as.vector(pho_rsd_split$SUB_MOD_RSD)
      
      samples <- colnames(pro_data)[!(colnames(pro_data) %in% c("Gene"))]
      # initiate ----------------------------------------------------------------
      # calculate the length of cis table
      ncis <- 0
      for (gene in kinase_cis) {
        ncis <- ncis + length(which(pho_rsd_split$SUBSTRATE==gene))
      }
      
      # calculate the length of trans table
      ntrans <- 0
      for (gene in kinase_trans) {
        subs <- k_s_table$SUB_GENE[k_s_table$GENE==gene & k_s_table$SUB_GENE!=gene]
        for ( sub in unique(subs)) {
          ntrans <- ntrans + length(which(pho_rsd_split$SUBSTRATE==sub))
        }
      }
      
      # looping over kinases for cis pairs -----------------------------------------------------------------
      # initiating the table for cis
      vec_char <- vector(mode = "character", length = ncis)
      vec_num <- vector(mode = "numeric", length = ncis) + NA
      KINASE <- vec_char;SUBSTRATE <- vec_char; SUB_MOD_RSD <- vec_char;
      FDR_pro_kin <- vec_num;FDR_pro_sub <- vec_num;FDR_pho_kin <- vec_num;
      coef_pro_kin <- vec_num;coef_pro_sub <- vec_num;coef_pho_kin <- vec_num;
      Cancer <- vec_char;transcript <- vec_char;model <- vec_char;
      Size <- vec_num;P_pro_kin <- vec_num;P_pro_sub <- vec_num;P_pho_kin <- vec_num;
      
      count <- 0
      for (kinase in kinase_cis){
        #for (kinase in "ERBB2"){#test
        # find protein expression level for the kinase
        pro_kin <- pro_data[pro_data$Gene == kinase,samples]
        
        substrate <- kinase
        if(nrow(pro_kin) != 0){
          s_pho_table <- which(pho_rsd_split$SUBSTRATE==substrate)
          
          # go through all the phosphosites
          for (i in s_pho_table) {
            
            # find phosphorylation level
            pho_sub <- pho_data[i,samples]
            
            #prepare regression data for model1
            data1 <- data.frame(t(rbind(pho_sub,pro_kin)))
            colnames(data1) <- c("pho_sub","pro_kin")
            
            size <- nrow(data1[complete.cases(data1),])
            if( size > least_samples ){
              fit1 <- lm(pho_sub ~ pro_kin,data = data1)
              
              count <- count + 1
              KINASE[count] <- kinase
              SUBSTRATE[count] <- substrate
              SUB_MOD_RSD[count] <- sub_mod_rsd[i]
              transcript[count] <- transcripts[i]
              P_pro_kin[count] <- c(coef(summary(fit1))[2,4])
              coef_pro_kin[count] <- fit1$coefficients[2]
              Size[count] <- size
            }
          }
        }
      }
      
      table_cis <- data.frame(KINASE,SUBSTRATE,SUB_MOD_RSD,
                              FDR_pro_kin,FDR_pro_sub,FDR_pho_kin,
                              coef_pro_kin,coef_pro_sub,coef_pho_kin,
                              Cancer,transcript,model,Size,
                              P_pro_kin,P_pro_sub,P_pho_kin)
      
      # looping over kinases for trans pairs -----------------------------------------------------------------
      # initiating the table for trans
      vec_char <- vector(mode = "character", length = ntrans)
      vec_num <- vector(mode = "numeric", length = ntrans) + NA
      KINASE <- vec_char;SUBSTRATE <- vec_char; SUB_MOD_RSD <- vec_char;
      FDR_pro_kin <- vec_num;FDR_pro_sub <- vec_num;FDR_pho_kin <- vec_num;
      coef_pro_kin <- vec_num;coef_pro_sub <- vec_num;coef_pho_kin <- vec_num;
      Cancer <- vec_char;transcript <- vec_char;model <- vec_char;
      Size <- vec_num;P_pro_kin <- vec_num;P_pro_sub <- vec_num;P_pho_kin <- vec_num;
      
      count <- 0
      for (kinase in kinase_trans){
        pho_kinase_g <- pho_gdata[pho_gdata$Gene == kinase,samples]
        
        if ( nrow(pho_kinase_g) > 0 ){
          k_sub <- unique(k_s_table$SUB_GENE[k_s_table$GENE == kinase & k_s_table$SUB_GENE != kinase])
          
          for (substrate in k_sub){# for each substrate for one kinase
            # find its phosphosites-row numbers
            s_pho_table <- which(pho_rsd_split$SUBSTRATE==substrate)
            
            # go through all the phosphosites
            for (i in s_pho_table) {
              
              # find phosphorylation level
              pho_sub <- pho_data[i,samples]
              
              # find substrate expressio level and normaize
              pro_sub <- pro_data[pro_data$Gene == substrate,samples]
              
              if(nrow(pro_sub) != 0){
                #prepare regression data for model2
                data2 <- data.frame(t(rbind(pho_sub,pro_sub,pho_kinase_g)))
                colnames(data2) <- c("pho_sub","pro_sub","pho_kinase_g")
                
                size <- nrow(data2[complete.cases(data2),])
                if(size > least_samples ){
                  fit2 <- lm(pho_sub ~ pro_sub + pho_kinase_g, data = data2)
                  
                  count <- count + 1
                  KINASE[count] <- kinase
                  SUBSTRATE[count] <- substrate
                  SUB_MOD_RSD[count] <- sub_mod_rsd[i]
                  transcript[count] <- transcripts[i]
                  P_pro_sub[count] <- c(coef(summary(fit2))[2,4])
                  coef_pro_sub[count] <- fit2$coefficients[2]
                  P_pho_kin[count] <- c(coef(summary(fit2))[3,4])
                  coef_pho_kin[count] <- fit2$coefficients[3]
                  Size[count] <- size
                }
              }
            }
          }
        }
      }
      table_trans <- data.frame(KINASE,SUBSTRATE,SUB_MOD_RSD,
                                FDR_pro_kin,FDR_pro_sub,FDR_pho_kin,
                                coef_pro_kin,coef_pro_sub,coef_pho_kin,
                                Cancer,transcript,model,Size,
                                P_pro_kin,P_pro_sub,P_pho_kin)
      
      # integrate table from all the models(need to repetite again for another cancer dataset) --------------------------------------------------------
      table_cis$model <- "pho_sub~pro_kin"
      tablecis <- table_cis[!is.na(table_cis$P_pro_kin),]
      
      table_trans$model <- "pho_sub~pro_sub+pho_kin"
      tabletrans <- table_trans[!is.na(table_trans$P_pro_sub),]
      
      # combine table
      table <- rbind(tablecis,tabletrans)
      
      # mark cancer and self-regulation
      table$Cancer <- cancer
      
      table_3can <- rbind(table_3can, table)
    }
    
    # combine table from 3 cancers and adjust P-values -------------------------------------------------------
    ## combine table from BRCA and OV
    table_3can$self <- as.character(table_3can$KINASE) == as.character(table_3can$SUBSTRATE)
    table_3can$SELF <- "trans"; table_3can$SELF[table_3can$self] <- "cis"
    name = c("pro_kin","pro_sub","pho_kin")
    
    ## adjust p-values to FDR
    for (cancer in c("BRCA","OV","CO")) {
      for(self in c(TRUE,FALSE)) {
        for(coln in name) {#adjust pvalues for each variable
          row <- (table_3can$self==self) & (table_3can$Cancer==cancer)
          table_3can[row,paste("FDR_",coln,sep = "")] <-p.adjust(table_3can[row,paste("P_",coln,sep = "")],method = "fdr")
        }
      }
    }
    
    table_3can$pair <- paste(table_3can$KINASE,table_3can$SUBSTRATE,table_3can$SUB_MOD_RSD,sep = ":")
    # write out tables --------------------------------------------------------
    
    tn = paste0(resultDnow, protein,"_substrate_regression_cptac2p_3can_", sample_type, ".txt")
    write.table(table_3can, file=tn, quote=F, sep = '\t', row.names = FALSE)
  }
}

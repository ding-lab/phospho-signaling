# Yige Wu @ WashU 2018 Jan
# look at correlations of kinase and downstream substrates phosphorylation status within PAM50 subtypes
# pho_sub~pro_sub+pho_kin(callapsed)

# source-------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')


# set variables -----------------------------------------------------------
least_samples <- 5# least number of samples with complete data for each model
inputnames <- c("noControl", "Control")
names(inputnames) <- c("tumor", "normal")

# inputs ------------------------------------------------------------------
ptms_site_pairs_sup <- read_csv(paste0(ppnD, "compile_enzyme_substrate/tables/compile_omnipath_networkin_depod/omnipath_networkin_enzyme_substrate_site_level_union_source_cleaned_addDEPOD_extended.csv"))
ptms_site_pairs_sup <- ptms_site_pairs_sup[ptms_site_pairs_sup$Source == "NetKIN" | (ptms_site_pairs_sup$Source != "NetKIN" & ptms_site_pairs_sup$networkin_score >= 5),]
ptms_site_pairs_sup <- data.frame(ptms_site_pairs_sup)

pam50_map <- loadPAM50Map()
subtypes <- c("Her2", "LumA", "LumB", "Basal")

# loop by 2 subtype2-----------------------------------------------------------------
## get combinations of 2 subtypes
subtype2_comb <- expand.grid(subtypes, subtypes)
subtype2_comb <- subtype2_comb[as.numeric(subtype2_comb$Var1) > as.numeric(subtype2_comb$Var2),]

for (sample_type in c("tumor")) {
  # for( enzyme_type in c("kinase")) {
  for( enzyme_type in c("phosphatase")) {
    k_s_table <- ptms_site_pairs_sup[ptms_site_pairs_sup$enzyme_type == enzyme_type,]
    k_s_psp <- load_ks_table(protein = enzyme_type)
    kinase_trans <- as.vector(unique(k_s_table$GENE[as.vector(k_s_table$GENE)!=as.vector(k_s_table$SUB_GENE)]))
    kinase_cis <- as.vector(unique(k_s_psp$GENE[as.vector(k_s_psp$GENE)==as.vector(k_s_psp$SUB_GENE)]))
    
    table_subtypes <- NULL
    for (cancer in c("BRCA")) { 
      pro_all <- fread(input = paste0(cptac_sharedD, cancer,"/",prefix[cancer], "_PRO_formatted_normalized_", inputnames[sample_type], ".txt"), data.table = F)
      pho_all <- fread(input = paste0(cptac_sharedD, cancer,"/",prefix[cancer], "_PHO_formatted_normalized_", inputnames[sample_type], ".txt"), data.table = F)
      pho_gall <- fread(input = paste0(cptac_sharedD, cancer,"/",prefix[cancer], "_collapsed_PHO_formatted_normalized_replicate_averaged_Tumor.txt"), data.table = F)
      pho_rsd_split <- formatPhosphosite(pho_all$Phosphosite, pho_all$Gene)
      
      transcripts <- as.vector(pho_rsd_split$transcript)
      sub_mod_rsd <- as.vector(pho_rsd_split$SUB_MOD_RSD)
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
      
      samples <- colnames(pro_all)[!(colnames(pro_all) %in% c("Gene"))]
      pam50 <- sampID2pam50(sampleID_vector = samples, sample_map = loadSampMap(), pam50_map = pam50_map)
      for (i in 1:nrow(subtype2_comb)) { 
        subtype1 <- as.character(subtype2_comb[i, "Var1"])
        subtype2 <- as.character(subtype2_comb[i, "Var2"])
        
        samples_subtype <- samples[pam50 %in% c(subtype1, subtype2)]
        pro_data <- pro_all[, c("Gene", samples_subtype)]
        pho_data <- pho_all[, c("Gene", "Phosphosite", samples_subtype)]
        pho_gdata <- pho_gall[, c("Gene", samples_subtype)]
        
        # looping over kinases for cis pairs -----------------------------------------------------------------
        # initiating the table for cis
        vec_char <- vector(mode = "character", length = ncis)
        vec_num <- vector(mode = "numeric", length = ncis) + NA
        KINASE <- vec_char;SUBSTRATE <- vec_char; SUB_MOD_RSD <- vec_char;
        FDR_pro_kin <- vec_num;FDR_pro_sub <- vec_num;FDR_pho_kin <- vec_num;
        coef_pro_kin <- vec_num;coef_pro_sub <- vec_num;coef_pho_kin <- vec_num;
        Cancer <- vec_char;transcript <- vec_char;model <- vec_char;
        Size <- vec_num;P_pro_kin <- vec_num;P_pro_sub <- vec_num;P_pho_kin <- vec_num;
        sd_pro_kin <- vec_num;sd_pro_sub <- vec_num;sd_pho_kin <- vec_num; sd_pho_sub <- vec_num;
        
        count <- 0
        for (kinase in kinase_cis){
          #for (kinase in "ERBB2"){#test
          # find enzyme_type expression level for the kinase
          pro_kin <- pro_data[pro_data$Gene == kinase,samples_subtype]
          
          substrate <- kinase
          if(nrow(pro_kin) != 0){
            s_pho_table <- which(pho_rsd_split$SUBSTRATE==substrate)
            
            # go through all the phosphosites
            for (i in s_pho_table) {
              
              # find phosphorylation level
              pho_sub <- pho_data[i,samples_subtype]
              
              #prepare regression data for model1
              data1 <- data.frame(t(rbind(pho_sub,pro_kin)))
              colnames(data1) <- c("pho_sub","pro_kin")
              data_complete <- data1[complete.cases(data1),]
              size <- nrow(data_complete)
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
                sd_pro_kin[count] <- sd(data_complete$pro_kin)
                sd_pho_sub[count] <- sd(data_complete$pho_sub)
              }
            }
          }
        }
        
        table_cis <- data.frame(KINASE,SUBSTRATE,SUB_MOD_RSD,
                                FDR_pro_kin,FDR_pro_sub,FDR_pho_kin,
                                coef_pro_kin,coef_pro_sub,coef_pho_kin,
                                Cancer,transcript,model,Size,
                                P_pro_kin,P_pro_sub,P_pho_kin,
                                sd_pro_kin,sd_pro_sub,sd_pho_kin,sd_pho_sub)
        
        # looping over kinases for trans pairs -----------------------------------------------------------------
        # initiating the table for trans
        vec_char <- vector(mode = "character", length = ntrans)
        vec_num <- vector(mode = "numeric", length = ntrans) + NA
        KINASE <- vec_char;SUBSTRATE <- vec_char; SUB_MOD_RSD <- vec_char;
        FDR_pro_kin <- vec_num;FDR_pro_sub <- vec_num;FDR_pho_kin <- vec_num;
        coef_pro_kin <- vec_num;coef_pro_sub <- vec_num;coef_pho_kin <- vec_num;
        Cancer <- vec_char;transcript <- vec_char;model <- vec_char;
        Size <- vec_num;
        P_pro_kin <- vec_num;P_pro_sub <- vec_num;P_pho_kin <- vec_num;
        sd_pro_kin <- vec_num;sd_pro_sub <- vec_num;sd_pho_kin <- vec_num; sd_pho_sub <- vec_num;
        
        count <- 0
        for (kinase in kinase_trans){
          pho_kinase_g <- pho_gdata[pho_gdata$Gene == kinase,samples_subtype]
          
          if ( nrow(pho_kinase_g) > 0 ){
            k_sub <- unique(k_s_table$SUB_GENE[k_s_table$GENE == kinase & k_s_table$SUB_GENE != kinase])
            
            for (substrate in k_sub){# for each substrate for one kinase
              # find its phosphosites-row numbers
              s_pho_table <- which(pho_rsd_split$SUBSTRATE==substrate)
              
              # go through all the phosphosites
              for (i in s_pho_table) {
                
                # find phosphorylation level
                pho_sub <- pho_data[i,samples_subtype]
                
                # find substrate expressio level and normaize
                pro_sub <- pro_data[pro_data$Gene == substrate,samples_subtype]
                
                if(nrow(pro_sub) != 0){
                  #prepare regression data for model2
                  data2 <- data.frame(t(rbind(pho_sub,pro_sub,pho_kinase_g)))
                  colnames(data2) <- c("pho_sub","pro_sub","pho_kinase_g")
                  data_complete <- data2[complete.cases(data2),]
                  size <- nrow(data_complete)
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
                    sd_pro_sub[count] <- sd(data_complete$pro_sub)
                    sd_pho_kin[count] <- sd(data_complete$pho_kinase_g)
                    sd_pho_sub[count] <- sd(data_complete$pho_sub)
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
                                  P_pro_kin,P_pro_sub,P_pho_kin,
                                  sd_pro_kin,sd_pro_sub,sd_pho_kin,sd_pho_sub)
        
        # integrate table from all the models(need to repetite again for another cancer dataset) --------------------------------------------------------
        table_cis$model <- "pho_sub~pro_kin"
        tablecis <- table_cis[!is.na(table_cis$P_pro_kin),]
        
        table_trans$model <- "pho_sub~pro_sub+pho_kin"
        tabletrans <- table_trans[!is.na(table_trans$P_pro_sub),]
        
        # combine table
        table <- rbind(tablecis,tabletrans)
        
        # mark cancer and self-regulation
        table$Cancer <- cancer
        table$subtype <- paste0(subtype1, "_", subtype2)
        table_subtypes <- rbind(table_subtypes, table)
      }
    }
    
    # combine table from 3 cancers and adjust P-values -------------------------------------------------------
    ## combine table from BRCA and OV
    table_subtypes$self <- as.character(table_subtypes$KINASE) == as.character(table_subtypes$SUBSTRATE)
    table_subtypes$SELF <- "trans"; table_subtypes$SELF[table_subtypes$self] <- "cis"
    name = c("pro_kin","pro_sub","pho_kin")
    
    ## adjust p-values to FDR
    for (cancer in unique(table_subtypes$Cancer)) {
      for (subtype in unique(table_subtypes$subtype)) {
        for(self in c(TRUE,FALSE)) {
          for(coln in name) {#adjust pvalues for each variable
            row <- (table_subtypes$self==self) & (table_subtypes$Cancer==cancer) & (table_subtypes$subtype == subtype)
            table_subtypes[row,paste("FDR_",coln,sep = "")] <-p.adjust(table_subtypes[row,paste("P_",coln,sep = "")],method = "fdr")
          }
        }
      }
    }
    
    table_subtypes$pair <- paste(table_subtypes$KINASE,table_subtypes$SUBSTRATE,table_subtypes$SUB_MOD_RSD,sep = ":")
    table_subtypes <- merge(table_subtypes, unique(ptms_site_pairs_sup[, c("pair", "Source")]), all.x = T)
    
    # write out tables --------------------------------------------------------
    tn = paste0(makeOutDir(resultD = resultD), enzyme_type,"_substrate_regression_cptac2p_", cancer, "_pam50_2subtypes_", sample_type, ".txt")
    write.table(table_subtypes, file=tn, quote=F, sep = '\t', row.names = FALSE)
  }
}


# loop by single subtype-----------------------------------------------------------------
for (sample_type in c("tumor")) {
  # for( enzyme_type in c("kinase")) {
  for( enzyme_type in c("phosphatase")) {
    k_s_table <- ptms_site_pairs_sup[ptms_site_pairs_sup$enzyme_type == enzyme_type,]
    k_s_psp <- load_ks_table(protein = enzyme_type)
    kinase_trans <- as.vector(unique(k_s_table$GENE[as.vector(k_s_table$GENE)!=as.vector(k_s_table$SUB_GENE)]))
    kinase_cis <- as.vector(unique(k_s_psp$GENE[as.vector(k_s_psp$GENE)==as.vector(k_s_psp$SUB_GENE)]))
    
    table_subtypes <- NULL
    for (cancer in c("BRCA")) { 
      pro_all <- fread(input = paste0(cptac_sharedD, cancer,"/",prefix[cancer], "_PRO_formatted_normalized_", inputnames[sample_type], ".txt"), data.table = F)
      pho_all <- fread(input = paste0(cptac_sharedD, cancer,"/",prefix[cancer], "_PHO_formatted_normalized_", inputnames[sample_type], ".txt"), data.table = F)
      pho_gall <- fread(input = paste0(cptac_sharedD, cancer,"/",prefix[cancer], "_collapsed_PHO_formatted_normalized_replicate_averaged_Tumor.txt"), data.table = F)
      pho_rsd_split <- formatPhosphosite(pho_all$Phosphosite, pho_all$Gene)
      
      transcripts <- as.vector(pho_rsd_split$transcript)
      sub_mod_rsd <- as.vector(pho_rsd_split$SUB_MOD_RSD)
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
      
      samples <- colnames(pro_all)[!(colnames(pro_all) %in% c("Gene"))]
      pam50 <- sampID2pam50(sampleID_vector = samples, sample_map = loadSampMap(), pam50_map = pam50_map)
      subtypes <- c("Her2", "LumA", "LumB", "Basal")
      for (subtype in subtypes) { 
        samples_subtype <- samples[pam50 == subtype]
        pro_data <- pro_all[, c("Gene", samples_subtype)]
        pho_data <- pho_all[, c("Gene", "Phosphosite", samples_subtype)]
        pho_gdata <- pho_gall[, c("Gene", samples_subtype)]
        
        # looping over kinases for cis pairs -----------------------------------------------------------------
        # initiating the table for cis
        vec_char <- vector(mode = "character", length = ncis)
        vec_num <- vector(mode = "numeric", length = ncis) + NA
        KINASE <- vec_char;SUBSTRATE <- vec_char; SUB_MOD_RSD <- vec_char;
        FDR_pro_kin <- vec_num;FDR_pro_sub <- vec_num;FDR_pho_kin <- vec_num;
        coef_pro_kin <- vec_num;coef_pro_sub <- vec_num;coef_pho_kin <- vec_num;
        Cancer <- vec_char;transcript <- vec_char;model <- vec_char;
        Size <- vec_num;P_pro_kin <- vec_num;P_pro_sub <- vec_num;P_pho_kin <- vec_num;
        sd_pro_kin <- vec_num;sd_pro_sub <- vec_num;sd_pho_kin <- vec_num; sd_pho_sub <- vec_num;
        
        count <- 0
        for (kinase in kinase_cis){
          #for (kinase in "ERBB2"){#test
          # find enzyme_type expression level for the kinase
          pro_kin <- pro_data[pro_data$Gene == kinase,samples_subtype]
          
          substrate <- kinase
          if(nrow(pro_kin) != 0){
            s_pho_table <- which(pho_rsd_split$SUBSTRATE==substrate)
            
            # go through all the phosphosites
            for (i in s_pho_table) {
              
              # find phosphorylation level
              pho_sub <- pho_data[i,samples_subtype]
              
              #prepare regression data for model1
              data1 <- data.frame(t(rbind(pho_sub,pro_kin)))
              colnames(data1) <- c("pho_sub","pro_kin")
              data_complete <- data1[complete.cases(data1),]
              size <- nrow(data_complete)
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
                sd_pro_kin[count] <- sd(data_complete$pro_kin)
                sd_pho_sub[count] <- sd(data_complete$pho_sub)
              }
            }
          }
        }
        
        table_cis <- data.frame(KINASE,SUBSTRATE,SUB_MOD_RSD,
                                FDR_pro_kin,FDR_pro_sub,FDR_pho_kin,
                                coef_pro_kin,coef_pro_sub,coef_pho_kin,
                                Cancer,transcript,model,Size,
                                P_pro_kin,P_pro_sub,P_pho_kin,
                                sd_pro_kin,sd_pro_sub,sd_pho_kin,sd_pho_sub)
        
        # looping over kinases for trans pairs -----------------------------------------------------------------
        # initiating the table for trans
        vec_char <- vector(mode = "character", length = ntrans)
        vec_num <- vector(mode = "numeric", length = ntrans) + NA
        KINASE <- vec_char;SUBSTRATE <- vec_char; SUB_MOD_RSD <- vec_char;
        FDR_pro_kin <- vec_num;FDR_pro_sub <- vec_num;FDR_pho_kin <- vec_num;
        coef_pro_kin <- vec_num;coef_pro_sub <- vec_num;coef_pho_kin <- vec_num;
        Cancer <- vec_char;transcript <- vec_char;model <- vec_char;
        Size <- vec_num;
        P_pro_kin <- vec_num;P_pro_sub <- vec_num;P_pho_kin <- vec_num;
        sd_pro_kin <- vec_num;sd_pro_sub <- vec_num;sd_pho_kin <- vec_num; sd_pho_sub <- vec_num;
        
        count <- 0
        for (kinase in kinase_trans){
          pho_kinase_g <- pho_gdata[pho_gdata$Gene == kinase,samples_subtype]
          
          if ( nrow(pho_kinase_g) > 0 ){
            k_sub <- unique(k_s_table$SUB_GENE[k_s_table$GENE == kinase & k_s_table$SUB_GENE != kinase])
            
            for (substrate in k_sub){# for each substrate for one kinase
              # find its phosphosites-row numbers
              s_pho_table <- which(pho_rsd_split$SUBSTRATE==substrate)
              
              # go through all the phosphosites
              for (i in s_pho_table) {
                
                # find phosphorylation level
                pho_sub <- pho_data[i,samples_subtype]
                
                # find substrate expressio level and normaize
                pro_sub <- pro_data[pro_data$Gene == substrate,samples_subtype]
                
                if(nrow(pro_sub) != 0){
                  #prepare regression data for model2
                  data2 <- data.frame(t(rbind(pho_sub,pro_sub,pho_kinase_g)))
                  colnames(data2) <- c("pho_sub","pro_sub","pho_kinase_g")
                  data_complete <- data2[complete.cases(data2),]
                  size <- nrow(data_complete)
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
                    sd_pro_sub[count] <- sd(data_complete$pro_sub)
                    sd_pho_kin[count] <- sd(data_complete$pho_kinase_g)
                    sd_pho_sub[count] <- sd(data_complete$pho_sub)
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
                                  P_pro_kin,P_pro_sub,P_pho_kin,
                                  sd_pro_kin,sd_pro_sub,sd_pho_kin,sd_pho_sub)
        
        # integrate table from all the models(need to repetite again for another cancer dataset) --------------------------------------------------------
        table_cis$model <- "pho_sub~pro_kin"
        tablecis <- table_cis[!is.na(table_cis$P_pro_kin),]
        
        table_trans$model <- "pho_sub~pro_sub+pho_kin"
        tabletrans <- table_trans[!is.na(table_trans$P_pro_sub),]
        
        # combine table
        table <- rbind(tablecis,tabletrans)
        
        # mark cancer and self-regulation
        table$Cancer <- cancer
        table$subtype <- subtype
        table_subtypes <- rbind(table_subtypes, table)
      }
      
    }
    
    # combine table from 3 cancers and adjust P-values -------------------------------------------------------
    ## combine table from BRCA and OV
    table_subtypes$self <- as.character(table_subtypes$KINASE) == as.character(table_subtypes$SUBSTRATE)
    table_subtypes$SELF <- "trans"; table_subtypes$SELF[table_subtypes$self] <- "cis"
    name = c("pro_kin","pro_sub","pho_kin")
    
    ## adjust p-values to FDR
    for (cancer in cancers_sort) {
      for (subtype in subtypes) {
        for(self in c(TRUE,FALSE)) {
          for(coln in name) {#adjust pvalues for each variable
            row <- (table_subtypes$self==self) & (table_subtypes$Cancer==cancer) & (table_subtypes$subtype == subtype)
            table_subtypes[row,paste("FDR_",coln,sep = "")] <-p.adjust(table_subtypes[row,paste("P_",coln,sep = "")],method = "fdr")
          }
        }
      }
    }
    
    table_subtypes$pair <- paste(table_subtypes$KINASE,table_subtypes$SUBSTRATE,table_subtypes$SUB_MOD_RSD,sep = ":")
    table_subtypes <- merge(table_subtypes, unique(ptms_site_pairs_sup[, c("pair", "Source")]), all.x = T)

    # write out tables --------------------------------------------------------
    tn = paste0(makeOutDir(resultD = resultD), enzyme_type,"_substrate_regression_cptac2p_", cancer, "_pam50_", sample_type, ".txt")
    write.table(table_subtypes, file=tn, quote=F, sep = '\t', row.names = FALSE)
  }
}


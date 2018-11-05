# Yige Wu @ WashU 2018 Jan
# look at correlations of kinase and downstream substrates phosphorylation status
# pho_sub~pro_sub+pho_kin(callapsed)

# source-------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')


# set variables -----------------------------------------------------------
least_samples <- 5# least number of samples with complete data for each model
inputnames <- c("noControl", "Control")
names(inputnames) <- c("tumor", "normal")
cptac_phase2process <- "cptac2p_3can"
cancers2process <- cancers_sort

# inputs ------------------------------------------------------------------
## input the new enzyme-substrate table
es_known_new <- read_csv(paste0(ppnD, "compile_enzyme_substrate/tables/compile_omnipath_networkin_depod_signor_manual/Omnipath_NetworKin_DEPOD_SignorNotSiteMapped_manualAdded.csv"))
es_known_new <- es_known_new[es_known_new$Source != "NetKIN" | (es_known_new$Source == "NetKIN" & es_known_new$networkin_score >= 5),]
es_known_new <- data.frame(es_known_new)

## input the old enzyme-substrate table
es_known_old <- read_csv(paste0(ppnD, "compile_enzyme_substrate/tables/compile_omnipath_networkin_depod_signor/omnipath_networkin_DEPOD_SignorNotSiteMapped.csv"))
es_known_old <- es_known_old[es_known_old$Source != "NetKIN" | (es_known_old$Source == "NetKIN" & es_known_old$networkin_score >= 5),]
es_known_old <- data.frame(es_known_old)

es_known_new <- es_known_new[!(es_known_new$pair_pro %in% es_known_old$pair_pro),]
es_known_new$pair <- paste0(es_known_new$GENE, ":", es_known_new$SUB_GENE, "_", es_known_new$SUB_MOD_RSD)

# loop -----------------------------------------------------------------
for (sample_type in c("tumor")) {
  for( enzyme_type in c("phosphatase", "kinase")) {
    k_s_table <- es_known_new[es_known_new$enzyme_type == enzyme_type,]
    kinase_trans <- as.vector(unique(k_s_table$GENE[as.vector(k_s_table$GENE)!=as.vector(k_s_table$SUB_GENE)]))
    kinase_cis <- unique(c(as.vector(k_s_table$GENE[as.vector(k_s_table$GENE) == as.vector(k_s_table$SUB_GENE)])))

    table_3can <- NULL
    for (cancer in cancers2process) { 
      pro_data <- loadProteinNormalizedTumor(cancer = cancer)
      pho_data <- loadPhosphositeNormalizedTumor(cancer = cancer)
      pho_gdata <- loadPhosphoproteinNormalizedTumor(cancer = cancer)
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
      tablecis <- NULL
      if (ncis > 0) {
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
        table_cis$model <- "pho_sub~pro_kin"
        tablecis <- table_cis[!is.na(table_cis$P_pro_kin),]
      }

      tabletrans <- NULL
      if (ntrans > 0){
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
        table_trans$model <- "pho_sub~pro_sub+pho_kin"
        tabletrans <- table_trans[!is.na(table_trans$P_pro_sub),]
      }

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
    
    ## input previous processed table
    table_old <- fread(input = paste0("./cptac2p/analysis_results/phospho_network/regression/tables/regression_omnipath&newworkin&depod&signor_protein_level/", 
                                            enzyme_type, "_substrate_regression_", cptac_phase2process, "_tumor.txt"), data.table = F)
    colnames2bind <- intersect(colnames(table_3can), colnames(table_old))
    table_3can <- rbind(table_3can[, colnames2bind], table_old[, colnames2bind])
    
    ## adjust p-values to FDR
    name = c("pro_kin","pro_sub","pho_kin")
    for (cancer in cancers2process) {
      for(self in c(TRUE,FALSE)) {
        for(coln in name) {#adjust pvalues for each variable
          row <- (table_3can$self==self) & (table_3can$Cancer==cancer)
          table_3can[row,paste("FDR_",coln,sep = "")] <-p.adjust(table_3can[row,paste("P_",coln,sep = "")],method = "fdr")
        }
      }
    }
    
    table_3can$pair <- paste(table_3can$KINASE,table_3can$SUBSTRATE,table_3can$SUB_MOD_RSD,sep = ":")
    table_3can <- merge(table_3can, unique(es_known_new[, c("pair", "Source")]), all.x = T)

    # write out tables --------------------------------------------------------
    tn = paste0(makeOutDir(resultD = resultD), enzyme_type,"_substrate_regression_", cptac_phase2process , "_", sample_type, ".txt")
    write.table(table_3can, file=tn, quote=F, sep = '\t', row.names = FALSE)
  }
}

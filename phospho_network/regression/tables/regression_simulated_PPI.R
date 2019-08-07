# Yige Wu @ WashU 2018 Jan
# look at correlations of kinase and downstream substrates phosphorylation status
# pho_sub~pro_sub+pho_kin(callapsed)

# source-------------------------
setwd("/Users/yigewu/Box Sync")
source('./cptac2p_analysis/phospho_network/phospho_network_shared.R')

# set variables -----------------------------------------------------------
least_samples <- 20 # least number of samples with complete data for each model
inputnames <- c("noControl", "Control")
names(inputnames) <- c("tumor", "normal")
data2process <- matrix(data = c("BRCA", "CDAP", "tumor", "scaled", "cptac2p",
                                "CO", "CDAP", "tumor", "scaled", "cptac2p",
                                "UCEC", "PGDAC", "tumor", "median_polishing", "cptac3",
                                "OV", "CDAP", "tumor", "scaled", "cptac2p",
                                "CCRCC", "PGDAC", "tumor", "MD_MAD", "cptac3",
                                "LIHC", "PGDAC", "tumor", "MD", "cptac3"), ncol = 5, byrow = T)
num_iterations <- 1
num_iterations <- 2

# inputs ------------------------------------------------------------------
omnipath_tab <- load_omnipath()
# omnipath_tab <- omnipath_tab[omnipath_tab$Source != "NetKIN" | (omnipath_tab$Source == "NetKIN" & omnipath_tab$networkin_score >= 5),]
omnipath_tab <- annotate_ks_source(omnipath_tab)
omnipath_tab <- omnipath_tab %>%
  filter(is.direct == T)
omnipath_tab <- data.frame(omnipath_tab)
k_s_table <- omnipath_tab
psp_tab <- load_psp()



# loop -----------------------------------------------------------------
for (iteration in 1:num_iterations) {
  subdir1 <- paste0(makeOutDir(resultD = resultD), "iteration", iteration, "/")
  dir.create(subdir1)
  
  # Use protein detected in CCRCC phosphorylation data to simulate----------------------
  i <- 5
  cancer <- data2process[i,1]
  pipeline_type <- data2process[i,2]
  sample_type <- data2process[i,3]
  norm_type <- data2process[i,4]
  cptac_phase <- data2process[i,5]
  ## simulate a set of proteins to substitute for the substrates
  pho_data <- loadParseProteomicsData(cancer = cancer, expression_type = "PHO", sample_type = sample_type, pipeline_type = pipeline_type, norm_type = norm_type)
  all_proteins <- unique(pho_data$Gene)
  k_s_table_trans <- k_s_table[as.vector(k_s_table$GENE) != as.vector(k_s_table$SUB_GENE),]
  
  ## simulate substrates
  substrates_trans  <- as.vector(k_s_table_trans$SUB_GENE)
  substrates_trans_uniq <- unique(substrates_trans)
  substrates_trans_uniq_replace <- sample(x = all_proteins, size = length(substrates_trans_uniq), replace = F)
  names(substrates_trans_uniq_replace) <- substrates_trans_uniq
  k_s_table_trans$SUB_GENE <- substrates_trans_uniq_replace[substrates_trans]
  
  ## simulate kinases
  enzymes_trans <- as.vector(k_s_table_trans$GENE)
  enzymes_trans_uniq <- unique(enzymes_trans)
  enzymes_trans_uniq_replace <- sample(x = all_proteins, size = length(enzymes_trans_uniq), replace = F)
  names(enzymes_trans_uniq_replace) <- enzymes_trans_uniq
  k_s_table_trans$GENE <- enzymes_trans_uniq_replace[enzymes_trans]
  kinase_trans <- as.vector(unique(k_s_table$GENE[as.vector(k_s_table$GENE)!=as.vector(k_s_table$SUB_GENE)]))
  
  
  kinase_cis <- unique(c(as.vector(k_s_table$GENE[as.vector(k_s_table$GENE) == as.vector(k_s_table$SUB_GENE)]), 
                         as.vector(psp_tab$GENE[as.vector(psp_tab$GENE)==as.vector(psp_tab$SUB_GENE)])))
  kinase_cis <- sample(x = all_proteins, size = length(kinase_cis), replace = F)
  
  write.table(x = k_s_table_trans, file = paste0(subdir1, "k_s_table_trans.txt"), quote = F, sep = "\t", row.names = F)
  write.table(x = kinase_cis, file = paste0(subdir1, "kinase_cis.txt"), quote = F, sep = "\t", row.names = F)
  
  for (i in 1:nrow(data2process)) {
    cancer <- data2process[i,1]
    pipeline_type <- data2process[i,2]
    sample_type <- data2process[i,3]
    norm_type <- data2process[i,4]
    cptac_phase <- data2process[i,5]
    
    subdir2 <- paste0(subdir1, cancer, "/")
    dir.create(subdir2)
    tn = paste0(subdir2, "regression_", cptac_phase , "_", cancer, "_", sample_type, "_", pipeline_type, "_", norm_type, "_iteration", iteration, ".txt")
    
    if (!file.exists(tn)) {
      ## input protein abundance
      if (cancer == "LIHC") {
        pro_data <- fread("./cptac2p/analysis_results/preprocess_files/tables/parse_China_Liver_log2_noimputation_NA50_median/LIHC_PRO_tumor_PGDAC_MD_partID.txt", data.table = F)
      } else {
        pro_data <- loadParseProteomicsData(cancer = cancer, expression_type = "PRO", sample_type = sample_type, pipeline_type = pipeline_type, norm_type = norm_type)
        pro_data <- data.frame(pro_data)
        colnames(pro_data)[1] <- "Gene"
      }
      
      pho_data <- loadParseProteomicsData(cancer = cancer, expression_type = "PHO", sample_type = sample_type, pipeline_type = pipeline_type, norm_type = norm_type)
      pho_data <- data.frame(pho_data)
      
      pho_gdata <- loadParseProteomicsData(cancer = cancer, expression_type = "collapsed_PHO", sample_type = sample_type, pipeline_type = pipeline_type, norm_type = norm_type)
      pho_gdata <- data.frame(pho_gdata)
      colnames(pho_gdata)[1] <- "Gene"
      
      pho_rsd_split <- data.frame(SUBSTRATE = pho_data$Gene, SUB_MOD_RSD = pho_data$Phosphosite)
      colnames(pho_rsd_split) <- c("SUBSTRATE", "SUB_MOD_RSD")
      
      if (any(is.na(pho_rsd_split$SUBSTRATE) | is.na(pho_rsd_split$SUB_MOD_RSD))) {
        stop()
      }
      
      transcripts <- NA
      sub_mod_rsd <- as.vector(pho_rsd_split$SUB_MOD_RSD)
      samples <- colnames(pro_data)[!(colnames(pro_data) %in% c("Gene"))]
      samples <- intersect(samples, colnames(pho_data))
      samples <- intersect(samples, colnames(pho_gdata))
      
      
      ## substrate need to be different from kinases
      
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
      sd_pro_kin <- vec_num;sd_pro_sub <- vec_num;sd_pho_kin <- vec_num; sd_pho_sub <- vec_num;
      
      count <- 0
      for (kinase in kinase_cis){
        #for (kinase in "ERBB2"){#test
        # find enzyme_type expression level for the kinase
        pro_kin <- pro_data[pro_data$Gene == kinase, samples]
        
        substrate <- kinase
        if(nrow(pro_kin) != 0){
          s_pho_table <- which(pho_rsd_split$SUBSTRATE==substrate)
          
          # go through all the phosphosites
          for (i in s_pho_table) {
            
            # find phosphorylation level
            pho_sub <- pho_data[i, samples]
            
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
        pho_kinase_g <- pho_gdata[pho_gdata$Gene == kinase,samples]
        
        if ( nrow(pho_kinase_g) > 0 ){
          k_sub <- unique(k_s_table_trans$SUB_GENE[k_s_table_trans$GENE == kinase])
          
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
      
      # integrate table from all the models(need to repetite again for another cancer dataset) --------------------------------------------------------
      table_cis$model <- "pho_sub~pro_kin"
      tablecis <- table_cis[!is.na(table_cis$P_pro_kin),]
      
      table_trans$model <- "pho_sub~pro_sub+pho_kin"
      tabletrans <- table_trans[!is.na(table_trans$P_pro_sub),]
      
      # combine table
      table <- rbind(tablecis,tabletrans)
      
      # mark cancer and self-regulation
      table$Cancer <- cancer
      table$self <- as.character(table$KINASE) == as.character(table$SUBSTRATE)
      table$SELF <- "trans"; table$SELF[table$self] <- "cis"
      table$pair <- paste(table$KINASE,table$SUBSTRATE,table$SUB_MOD_RSD,sep = ":")
      
      ## write out
      write.table(table, file=tn, quote=F, sep = '\t', row.names = FALSE)
    }
  }
  
}


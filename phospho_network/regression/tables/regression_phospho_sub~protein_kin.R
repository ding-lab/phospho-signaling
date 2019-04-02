# Yige Wu @ WashU 2018 Jan
# look at correlations of kinase and downstream substrates phosphorylation status
# pho_sub~pro_sub+pho_kin(callapsed)

# source-------------------------
setwd("/Users/yigewu/Box Sync")
source('./cptac2p_analysis/phospho_network/phospho_network_shared.R')

# set variables -----------------------------------------------------------
least_samples <- 20# least number of samples with complete data for each model
inputnames <- c("noControl", "Control")
names(inputnames) <- c("tumor", "normal")
data2process <- matrix(data = c("BRCA", "CDAP", "tumor", "scaled", "cptac2p",
                                "CO", "CDAP", "tumor", "scaled", "cptac2p",
                                "UCEC", "PGDAC", "tumor", "median_polishing", "cptac3",
                                "OV", "CDAP", "tumor", "scaled", "cptac2p",
                                "CCRCC", "PGDAC", "tumor", "MD_MAD", "cptac3",
                                "LIHC", "PGDAC", "tumor", "MD", "cptac3"), ncol = 5, byrow = T)

# inputs ------------------------------------------------------------------
omnipath_tab <- load_omnipath()
omnipath_tab <- omnipath_tab[!(omnipath_tab$Source %in% c("NetKIN", "PhosphoNetworks", "MIMP")),]
omnipath_tab <- data.frame(omnipath_tab)

k_s_table <- omnipath_tab
psp_tab <- load_psp()
kinase_trans <- as.vector(unique(k_s_table$GENE[as.vector(k_s_table$GENE)!=as.vector(k_s_table$SUB_GENE)]))
kinase_cis <- unique(c(as.vector(k_s_table$GENE[as.vector(k_s_table$GENE) == as.vector(k_s_table$SUB_GENE)]), 
                       as.vector(psp_tab$GENE[as.vector(psp_tab$GENE)==as.vector(psp_tab$SUB_GENE)])))

# loop -----------------------------------------------------------------
for (i in 1:nrow(data2process)) {
  cancer <- data2process[i,1]
  pipeline_type <- data2process[i,2]
  sample_type <- data2process[i,3]
  norm_type <- data2process[i,4]
  cptac_phase <- data2process[i,5]
  
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
  
  # initiate ----------------------------------------------------------------
  # calculate the length of trans table
  ntrans <- 0
  for (gene in kinase_trans) {
    subs <- k_s_table$SUB_GENE[k_s_table$GENE==gene & k_s_table$SUB_GENE!=gene]
    for ( sub in unique(subs)) {
      ntrans <- ntrans + length(which(pho_rsd_split$SUBSTRATE==sub))
    }
  }

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
    pho_kinase_g <- pro_data[pro_data$Gene == kinase,samples]
    
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
  
  # integrate table from all the models(need to repetite again for another cancer dataset) --------------------------------------------------------
  # table_cis$model <- "pho_sub~pro_kin"
  # tablecis <- table_cis[!is.na(table_cis$P_pro_kin),]
  tablecis <- NULL
  
  table_trans$model <- "pho_sub~pro_sub+pro_kin"
  tabletrans <- table_trans[!is.na(table_trans$P_pro_sub),]
  
  # combine table
  table <- rbind(tablecis,tabletrans)
  
  # mark cancer and self-regulation
  table$Cancer <- cancer
  table$self <- as.character(table$KINASE) == as.character(table$SUBSTRATE)
  table$SELF <- "trans"; table$SELF[table$self] <- "cis"
  table$pair <- paste(table$KINASE,table$SUBSTRATE,table$SUB_MOD_RSD,sep = ":")
  table <- merge(table, unique(omnipath_tab[, c("pair", "Source", "enzyme_type")]), all.x = T)
  
  ## write out
  tn = paste0(makeOutDir(resultD = resultD), "regression_", cptac_phase , "_", cancer, "_", sample_type, "_", pipeline_type, "_", norm_type, ".txt")
  write.table(table, file=tn, quote=F, sep = '\t', row.names = FALSE)
}


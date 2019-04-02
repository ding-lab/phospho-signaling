# Yige Wu @ WashU 2018 Jan
# look at correlations of kinase and downstream substrates phosphorylation status
# pho_sub~pro_sub+pho_kin(callapsed)

# source-------------------------
wd <- getwd()
if (wd != "/Users/yigewu/Box Sync") {
  setwd("/Users/yigewu/Box Sync")
}
source('./cptac2p_analysis/phospho_network/phospho_network_shared.R')

# set variables -----------------------------------------------------------
least_samples <- 5# least number of samples with complete data for each model
inputnames <- c("noControl", "Control")
names(inputnames) <- c("tumor", "normal")
data2process <- matrix(data = c("BRCA", "CDAP", "tumor", "scaled", "cptac2p",
                                "CO", "CDAP", "tumor", "scaled", "cptac2p",
                                "UCEC", "PGDAC", "tumor", "median_polishing", "cptac3",
                                "OV", "CDAP", "tumor", "scaled", "cptac2p",
                                "CCRCC", "PGDAC", "tumor", "MD_MAD", "cptac3"), ncol = 5, byrow = T)
data2process

# input regression model --------------------------------------------------
regression <- fread(input = paste0(ppnD, "regression/tables/annotate_regression_with_mut_impact/", 
                                   "regression_cptac2p_cptac3_tumor_reg_nonNA20_mut_impact_cancer_specificity_annotated.txt"),
                    data.table = F)
regression <- regression %>%
  filter(pair_pro %in% omnipath_tab$pair_pro[!(omnipath_tab$Source %in% c("NetKIN", "PhosphoNetworks", "MIMP"))] | pair_pro %in% psp_tab$pair_pro)

regression <- adjust_regression_by_nonNA(regression = regression, reg_nonNA = 20, reg_sig = reg_sig)
regression <- annotate_ks_source(regression = regression)

# overlap phosphatase and kinase ------------------------------------------
enzyme_type <- "phosphatase"
tab_pp <- regression[regression$enzyme_type == enzyme_type & regression$SELF == SELF,]
tab_pp_reg <- tab_pp[tab_pp$regulated,] 

enzyme_type <- "kinase"
tab_pk <- regression[regression$enzyme_type == enzyme_type & regression$SELF == SELF,]
tab_pk_reg <- tab_pk[tab_pk$regulated,] 

tab_pk_pp <- merge(tab_pk_reg[tab_pk_reg$SELF == "trans", c("GENE", "SUB_GENE", "SUB_MOD_RSD", "Cancer", "pair", "Source")], 
                   tab_pp_reg[tab_pp_reg$SELF == "trans", c("GENE", "SUB_GENE", "SUB_MOD_RSD", "Cancer", "pair", "Source")], 
                   by = c("SUB_GENE", "SUB_MOD_RSD", "Cancer"), suffixes = c(".pk", ".pp"))

# loop -----------------------------------------------------------------
# for (i in 5) {
regression_sup <- NULL
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
  
  # get the substrate phosphosites to work on -------------------------------
  tab_sub <- tab_pk_pp %>%
    filter(Cancer == cancer) %>%
    select(SUB_GENE, SUB_MOD_RSD) %>%
    unique()
  
  # initiate ----------------------------------------------------------------
  # calculate the length of trans table
  tab_trans <- regression %>%
    filter(Cancer == cancer) %>%
    filter((pair %in% tab_pk_pp$pair.pk[tab_pk_pp$Cancer == cancer]) | (pair %in% tab_pk_pp$pair.pp[tab_pk_pp$Cancer == cancer]))
  
  ntrans <- nrow(tab_trans)
  
  # looping over kinases for trans pairs -----------------------------------------------------------------
  # initiating the table for trans
  vec_char <- vector(mode = "character", length = ntrans)
  vec_num <- vector(mode = "numeric", length = ntrans) + NA
  KINASE <- vec_char;SUBSTRATE <- vec_char; SUB_MOD_RSD <- vec_char;
  FDR_pro_kin <- vec_num;FDR_pro_sub <- vec_num;FDR_pho_kin <- vec_num;
  coef_pro_kin <- vec_num;coef_pro_sub <- vec_num;coef_pho_kin <- vec_num;
  Cancer <- vec_char;transcript <- vec_char;model <- vec_char;
  Size <- vec_num;
  model <- vec_char
  P_pro_kin <- vec_num;P_pro_sub <- vec_num;P_pho_kin <- vec_num;
  sd_pro_kin <- vec_num;sd_pro_sub <- vec_num;sd_pho_kin <- vec_num; sd_pho_sub <- vec_num;
  
  for (j in 1:nrow(tab_sub)) {
    substrate_gene <- tab_sub$SUB_GENE[j]
    substrate_rsd <- tab_sub$SUB_MOD_RSD[j]
    
    pho_sub <- pho_data %>%
      filter(Gene == substrate_gene & Phosphosite == substrate_rsd) %>%
      select(samples)
    pho_sub <- pho_sub[1,]
    
    pro_sub <- pro_data %>%
      filter(Gene == substrate_gene) %>%
      select(samples)
    
    kinase_genes2process <- tab_pk_pp %>%
      filter(SUB_GENE == substrate_gene & SUB_MOD_RSD == substrate_rsd) %>%
      filter(Cancer == cancer) %>%
      select(GENE.pk) %>%
      as.vector() %>%
      unique
    kinase_genes2process
    
    phosphatase_genes2process <- tab_pk_pp %>%
      filter(SUB_GENE == substrate_gene & SUB_MOD_RSD == substrate_rsd) %>%
      filter(Cancer == cancer) %>%
      select(GENE.pp) %>%
      as.vector() %>%
      unique
    phosphatase_genes2process
    
    df4regression <- data.frame(t(rbind(pho_sub,pro_sub)))
    colnames(df4regression) <- c("pho_sub","pro_sub")
    
    enzyme_genes2process <- c(kinase_genes2process$GENE.pk, phosphatase_genes2process$GENE.pp)
    for (enzyme_gene in enzyme_genes2process) {
      pho_kinase_g <- pho_gdata %>%
        filter(Gene == enzyme_gene) %>%
        select(samples)
      
      df4regression[, paste0("pho_", enzyme_gene)] <- t(pho_kinase_g)
    }

    data_complete <- df4regression[complete.cases(df4regression),]
    size <- nrow(data_complete)
    if(size > least_samples ){
      
      model_string_tmp <- paste0("pho_sub ~ ", paste0(paste0("pho_", enzyme_genes2process, " + "), collapse = ""), "pro_sub")
      formula_tmp <- as.formula(model_string_tmp)
      
      fit2 <- lm(formula = formula_tmp, data = data_complete)
      
      count <- count + 1
      for (enzyme_gene in enzyme_genes2process) {
        row_tmp <- which(tab_trans$SUBSTRATE == substrate_gene & tab_trans$SUB_MOD_RSD == substrate_rsd & tab_trans$GENE == enzyme_gene)
        var_tmp <- paste0("pho_", enzyme_gene)
        if (var_tmp %in% rownames(coef(summary(fit2)))) {
          P_pho_kin[row_tmp] <- c(coef(summary(fit2))[paste0("pho_", enzyme_gene),4])
          coef_pho_kin[row_tmp] <- fit2$coefficients[paste0("pho_", enzyme_gene)]
        } else {
          P_pho_kin[row_tmp] <- NA
          coef_pho_kin[row_tmp] <- NA
        }
        P_pro_sub[row_tmp] <- c(coef(summary(fit2))["pro_sub",4])
        coef_pro_sub[row_tmp] <- fit2$coefficients["pro_sub"]
        
        Size[row_tmp] <- size
        model[row_tmp] <- model_string_tmp
      }
    }
  }
  tab_trans$P_pho_kin <- P_pho_kin
  tab_trans$coef_pho_kin <- coef_pho_kin
  tab_trans$P_pro_sub <- P_pro_sub
  tab_trans$coef_pro_sub <- coef_pro_sub
  tab_trans$Size <- Size
  tab_trans$model <- model
  
  ## write out
  tn = paste0(makeOutDir(resultD = resultD), "regression_", cptac_phase , "_", cancer, "_", sample_type, "_", pipeline_type, "_", norm_type, ".txt")
  write.table(tab_trans, file=tn, quote=F, sep = '\t', row.names = FALSE)
  
  regression_sup <- rbind(tab_trans, regression_sup)
}


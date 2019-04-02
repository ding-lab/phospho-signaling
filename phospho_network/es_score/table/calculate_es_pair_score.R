# Yige Wu @ WashU 2019 Jan
# calculate enzyme-substrate pair score for each sample for regulated pairs
## to calculate the kinase-substrate score, determine outlier samples and do fisher's test on determining subtype enrichment

# source ------------------------------------------------------------------
baseD = "/Users/yigewu/Box\ Sync/"
wd <- getwd()
if (wd != baseD) {
  setwd(baseD)
}
source('./cptac2p_analysis/preprocess_files/preprocess_files_shared.R')
source('./cptac2p_analysis/phospho_network/phospho_network_shared.R')


# set variables -----------------------------------------------------------
reg_nonNA <- 20
## just filter for kinase now
enzyme_type <- "kinase"
fdr_pp <- 0.1
fdr_pk <- 0.05
cancers2process <- c("BRCA", "OV", "CO", "UCEC", "CCRCC")
## the parameters for the input proteome data files
data2process <- matrix(data = c("BRCA", "CDAP", "tumor", "scaled", "cptac2p",
                                "CO", "CDAP", "tumor", "scaled", "cptac2p",
                                "UCEC", "PGDAC", "tumor", "median_polishing", "cptac3",
                                "OV", "CDAP", "tumor", "scaled", "cptac2p",
                                "CCRCC", "PGDAC", "tumor", "MD_MAD", "cptac3",
                                "LIHC", "PGDAC", "tumor", "MD", "cptac3"), ncol = 5, byrow = T)
## outlier threshold
outlier_sd <- 1.5


# input regression result table -------------------------------------------
# regression <- fread(input = paste0(ppnD, "regression/tables/", "change_regression_nonNA/", 
#                                    "regression", "_", "cptac2p_cptac3", "_", "tumor", "_reg_nonNA", reg_nonNA, ".txt"), data.table = F)
regression <- fread(input = paste0(ppnD, "regression/tables/annotate_regression_with_mut_impact/", 
                                   "regression_cptac2p_cptac3_tumor_reg_nonNA20_mut_impact_cancer_specificity_annotated.txt"),
                    data.table = F)
## already annotated with the supporting evidence for each pair

## divide into kinase and phosphatase related
regression <- regression[regression$enzyme_type == enzyme_type,]
## get the significant pairs
regression <- markSigSiteCan(regression = regression, sig_thres = ifelse(enzyme_type == "kinase" ,fdr_pk, fdr_pp), enzyme_type = enzyme_type)
regression$regulated <- (regression$fdr_sig & regression$coef_sig)

# decide on the pairs to generate outlier score for -----------------------
## better if generater a range or pairs encompassing as much as we will possibly cover in this papaer
## but for test run sake, I'll just focus on known pairs (experimentally supported) no matter it's correlated or not
## the scores for the same pairs are generated for all cancers
pairs2generate <- unique(regression$pair[!is.na(regression$Source) & !(regression$Source %in% c("MIMP", "PhosphoNetworks", "NetKIN"))])

## add mutational associated pairs
pairs2generate <- unique(c(pairs2generate, 
                           unique(regression$pair[regression$SUB_GENE %in% unlist(SMGs) & regression$p.SUB_GENE < 0.05 & !is.na(regression$p.SUB_GENE)]),
                           unique(regression$pair[regression$GENE %in% unlist(SMGs) & regression$p.GENE < 0.05 & !is.na(regression$p.GENE)])))

## since the cis pairs are most likley going to show effect however the site-levelr annotation difference impede the exact match, here I'll just include all regulated ones
pairs2generate <- unique(c(pairs2generate, unique(regression$pair[regression$regulated & regression$SELF == "cis"])))
pairs2generate %>% length()
## divide the pairs to process into cis and trans
pairs2generate_head <- unique(regression[regression$pair %in% pairs2generate,  c("GENE", "SUB_GENE", "SUB_MOD_RSD", "SELF", "pair")])
rownames(pairs2generate_head) <- pairs2generate_head$pair
pairs2generate_head$site_id <- paste0(pairs2generate_head$SUB_GENE, ".", pairs2generate_head$SUB_MOD_RSD)

# kinases -----------------------------------------------------------------
subdir1 <- paste0(makeOutDir(resultD = resultD), enzyme_type, "/")
dir.create(subdir1)

# for (i in c(2))) {
for (i in 1:5) {
  cancer <- data2process[i,1]
  pipeline_type <- data2process[i,2]
  sample_type <- data2process[i,3]
  norm_type <- data2process[i,4]
  cptac_phase <- data2process[i,5]
  
  if (!file.exists(paste0(subdir1, "esscore_tab", "_", cancer, "_", enzyme_type, "_reg_nonNA", reg_nonNA, ".txt"))) {
    ## input gene-level and site-level phosphorylation and protein levels
    pro_data <- loadParseProteomicsData(cancer = cancer, expression_type = "PRO", sample_type = sample_type, pipeline_type = pipeline_type, norm_type = norm_type)
    pro_data <- data.frame(pro_data)
    
    pro_score <- as.matrix(pro_data[,-1])
    pro_score <- scale_by_row(pro_score)
    rownames(pro_score) <- pro_data$Gene
    pro_score <- data.frame(pro_score)
    
    phog_data <- loadParseProteomicsData(cancer = cancer, expression_type = "collapsed_PHO", sample_type = sample_type, pipeline_type = pipeline_type, norm_type = norm_type)
    phog_data <- data.frame(phog_data)
    rownames(phog_data) <- phog_data$Gene
    
    phog_score <- as.matrix(phog_data[,-1])
    phog_score <- scale_by_row(phog_score)
    rownames(phog_score) <- phog_data$Gene
    phog_score <- data.frame(phog_score)
    
    pho_data <- loadParseProteomicsData(cancer = cancer, expression_type = "PHO", sample_type = sample_type, pipeline_type = pipeline_type, norm_type = norm_type)
    partIDs <- colnames(pho_data); partIDs <- partIDs[!(partIDs %in% c("Gene", "Phosphosite", "Peptide_ID"))]
    pho_data$site_id <- paste0(pho_data$Gene, ".", pho_data$Phosphosite)
    
    pho_score <- as.matrix(pho_data[, partIDs])
    pho_score <- scale_by_row(pho_score)
    rownames(pho_score) <- pho_data$site_id
    pho_score <- data.frame(pho_score)
    rm(pho_data)
    
    ## generate kinase score table
    ### for cis pairs the kinase score will use protein level
    escore_tab_cis <- pairs2generate_head[pairs2generate_head$SELF == "cis",]
    escore_tab_cis <- cbind(escore_tab_cis, pro_score[as.vector(escore_tab_cis$GENE),])
    
    ### for trans pairs the kinase score will use phosphoprotein level
    escore_tab_trans <- pairs2generate_head[pairs2generate_head$SELF == "trans",]
    escore_tab_trans <- cbind(escore_tab_trans, phog_score[as.vector(escore_tab_trans$GENE),])
    
    ## bind cis and trans kinase score
    escore_tab <- rbind(escore_tab_cis, escore_tab_trans)
    escore_tab <- escore_tab[as.vector(pairs2generate_head$pair),]
    colnames(escore_tab) <- c(colnames(pairs2generate_head), partIDs)
    write.table(x = escore_tab, file = paste0(subdir1, "escore_tab", "_", cancer, "_", enzyme_type, "_reg_nonNA", reg_nonNA, ".txt"), row.names = F, quote = F, sep = "\t")
    
    ## generate substrate score table
    sscore_tab <- pairs2generate_head
    sscore_tab <- cbind(sscore_tab, pho_score[as.vector(sscore_tab$site_id),])
    colnames(sscore_tab) <- c(colnames(pairs2generate_head), partIDs)
    write.table(x = sscore_tab, file = paste0(subdir1, "sscore_tab", "_", cancer, "_", enzyme_type, "_reg_nonNA", reg_nonNA, ".txt"), row.names = F, quote = F, sep = "\t")
    
    ## add above 2 tables to generate k-s score table
    esscore_tab <- pairs2generate_head
    esscore_tab <- cbind(esscore_tab, (escore_tab[, partIDs] + sscore_tab[, partIDs]))
    colnames(esscore_tab) <- c(colnames(pairs2generate_head), partIDs)
    write.table(x = esscore_tab, file = paste0(subdir1, "esscore_tab", "_", cancer, "_", enzyme_type, "_reg_nonNA", reg_nonNA, ".txt"), row.names = F, quote = F, sep = "\t")
  } else {
    esscore_tab <- fread(input = paste0(subdir1, "esscore_tab", "_", cancer, "_", enzyme_type, "_reg_nonNA", reg_nonNA, ".txt"), data.table = F)
    partIDs <- colnames(esscore_tab); partIDs <- partIDs[!(partIDs %in% c("GENE", "SUB_GENE", "SUB_MOD_RSD", "SELF", "pair", "site_id"))]
  }
  fn <- paste0("./cptac2p/analysis_results/phospho_network/es_score/table/calculate_es_pair_score/fisher_stat_tab", "_outlier", outlier_sd, "SD_", cancer, "_", enzyme_type, "_reg_nonNA", reg_nonNA, ".txt")
  if (!file.exists(fn)) {
    ## calculate per row which samples are high outliers to generate outlier status table
    esscore_tab_scaled <- esscore_tab[, partIDs]
    esscore_tab_scaled <- scale_by_row(esscore_tab_scaled)
    esscore_tab_outlier <- cbind(pairs2generate_head, (esscore_tab_scaled > outlier_sd))
    colnames(esscore_tab_outlier) <- c(colnames(pairs2generate_head), partIDs)
    write.table(x = esscore_tab_outlier, file = paste0(subdir1, "esscore_tab_outlier", outlier_sd, "SD_", cancer, "_", enzyme_type, "_reg_nonNA", reg_nonNA, ".txt"), row.names = F, quote = F, sep = "\t")
    next()
    
    ## input subtype info and get the subtypes for degsinated order of patient IDs
    if (cancer == "BRCA") {
      subtypes <- partID2pam50(patientID_vector = partIDs, pam50_map = loadPAM50Map())
    } else if (cancer == "CO") {
      subtypes <- partID2MSI(patientID_vector = partIDs, subtype_map = loadMSIMap())
    } else if (cancer == "UCEC") {
      subtypes <- partID2UCECsubtype(patientID_vector = partIDs)
    } else {
      next()
    }
    
    ## do fishers's table per subtype to see subtype enrichment
    fisher_stat_tab <- NULL
    for (subtype in names(table(subtypes))[table(subtypes) > 4]) {
      is.subtype <- (subtypes == subtype)
      fisher_stat_subtype <- sapply(rownames(esscore_tab_outlier), FUN = function(pair, is.subtype, esscore_tab_outlier) {
        partIDs_nonNA <- colnames(esscore_tab_outlier)[which(esscore_tab_outlier[pair,] == TRUE | esscore_tab_outlier[pair,] == FALSE)]
        tab4fisher <- table(esscore_tab_outlier[pair, partIDs_nonNA], is.subtype[partIDs_nonNA])
        if (all(dim(tab4fisher) == c(2,2))) {
          fisher_stat <- fisher.test(tab4fisher, alternative = "greater")
          fisher_p.value <- fisher_stat$p.value
          fisher_or_bottom <- as.numeric(fisher_stat$conf.int[1])
          fisher_or_upper <- as.numeric(fisher_stat$conf.int[2])
          
          num_is.outlier_is.subtype = as.vector(tab4fisher)[4]
          num_not.outlier_not.subtype = as.vector(tab4fisher)[1]
          num_is.outlier_not.subtype = as.vector(tab4fisher)[2]
          num_not.outlier_is.subtype = as.vector(tab4fisher)[3]
        } else {
          fisher_p.value <- NA; fisher_or_bottom <- NA; fisher_or_upper <- NA
          num_is.outlier_is.subtype <- NA
          num_not.outlier_not.subtype <- NA
          num_is.outlier_not.subtype <- NA
          num_not.outlier_is.subtype <- NA
        }
        fisher_stat <- c(fisher_p.value, fisher_or_bottom, fisher_or_upper,
                         num_is.outlier_is.subtype, num_not.outlier_not.subtype,
                         num_is.outlier_not.subtype, num_not.outlier_is.subtype)
        return(fisher_stat)
      }, is.subtype = is.subtype, esscore_tab_outlier = esscore_tab_outlier)
      
      fisher_stat_subtype <- t(fisher_stat_subtype)
      colnames(fisher_stat_subtype) <- c("p.value", "or_bottom", "or_upper", 
                                         "num_is.outlier_is.subtype", "num_not.outlier_not.subtype",
                                         "num_is.outlier_not.subtype", "num_not.outlier_is.subtype")
      fisher_stat_subtype <- cbind(pairs2generate_head, fisher_stat_subtype)
      fisher_stat_subtype$subtype <- subtype
      fisher_stat_tab <- rbind(fisher_stat_tab, fisher_stat_subtype)
    }
    fisher_stat_tab$cancer <- cancer
    write.table(x = fisher_stat_tab, file = fn, 
                row.names = F, quote = F, sep = "\t")
  }
}


# clean up ----------------------------------------------------------------
enzyme_type <- "kinase"
fisher_stat_cancers <- NULL
for (cancer in c("BRCA", "CO", "UCEC")) {
  fisher_stat_cancer <- fread(input = paste0(makeOutDir(resultD = resultD),
                                             "fisher_stat_tab_outlier", outlier_sd, "SD_", cancer, "_", enzyme_type, "_reg_nonNA", reg_nonNA, ".txt"))
  fisher_stat_cancer <- fisher_stat_cancer[!is.na(fisher_stat_cancer$p.value) & fisher_stat_cancer$p.value < 1,]
  fisher_stat_cancer <- unique(fisher_stat_cancer)
  for (SELF in c("cis", "trans")) {
    fisher_stat_cancer_self <- fisher_stat_cancer[fisher_stat_cancer$SELF == SELF,]
    for (subtype in unique(fisher_stat_cancer_self$subtype)) {
      fisher_stat_cancer_self_subtype <- fisher_stat_cancer_self[fisher_stat_cancer_self$subtype == subtype,]
      fisher_stat_cancer_self_subtype$fdr <- p.adjust(p = fisher_stat_cancer_self_subtype$p.value, method = "fdr")
      fisher_stat_cancers <- rbind(fisher_stat_cancers, fisher_stat_cancer_self_subtype)
    }
  }
}
fisher_stat_cancers <- unique(fisher_stat_cancers)
write.table(x = fisher_stat_cancers, file = paste0(makeOutDir(resultD = resultD),
                                               "fisher_stat_tab_outlier", outlier_sd, "SD_all_cancers_", enzyme_type, "_reg_nonNA", reg_nonNA, ".txt"),
            row.names = F, quote = F, sep = "\t")

fisher_stat_cancers_fdr_sig <- fisher_stat_cancers[fisher_stat_cancers$fdr < 0.05,]
fisher_stat_cancers_sig <- fisher_stat_cancers[fisher_stat_cancers$p.value < 0.05 & fisher_stat_cancers$SELF == "trans" & fisher_stat_cancers$num_is.outlier_is.subtype > 4,]

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


# input ccRCC clinical info -----------------------------------------------
ccRCC_clinical <- readxl::read_xlsx(path = "./Ding_Lab/Projects_Current/CPTAC/CPTACIII/CCRCC_AWG/CCRCC_shared_data/manuscripts/ccRCC Manuscript/CPTAC3-ccRCC-SupplementaryTables_Final/Table S1.xlsx",
                                    sheet = "ccrcc_clinical_characteristics")


# loop -----------------------------------------------------------------

# for (i in c(2))) {
for (i in 5) {
  cancer <- data2process[i,1]
  pipeline_type <- data2process[i,2]
  sample_type <- data2process[i,3]
  norm_type <- data2process[i,4]
  cptac_phase <- data2process[i,5]
  
  esscore_tab_fn <- paste0("./cptac2p/analysis_results/phospho_network/es_score/table/calculate_es_pair_score/kinase/escore_tab", "_", cancer, "_", enzyme_type, "_reg_nonNA", reg_nonNA, ".txt")
  esscore_tab <- fread(input = esscore_tab_fn, data.table = F)

  partIDs_exp <- colnames(esscore_tab)
  partIDs_exp <- partIDs_exp[!(partIDs_exp %in% c("GENE", "SUB_GENE", "SUB_MOD_RSD", "SELF", "pair", "site_id"))]
  
  esscore_tab_score <- esscore_tab[,partIDs_exp]
  esscore_tab_outlier_mat <- (esscore_tab_score > outlier_sd)
  esscore_tab_outlier <- cbind(esscore_tab[,c("GENE", "SUB_GENE", "SUB_MOD_RSD", "SELF", "pair", "site_id")], esscore_tab_outlier_mat)
  rownames(esscore_tab_outlier) <- as.vector(esscore_tab_outlier$pair)
  esscore_tab_outlier <- esscore_tab_outlier[regression$pair[regression$Cancer == cancer & regression$regulated == T],]
  esscore_tab_outlier <- esscore_tab_outlier[rowSums(esscore_tab_outlier[, partIDs_exp] == T, na.rm = T) >= 4,]
  
  is.highgrade <- matrix(ccRCC_clinical$Grade %in% c("G3", "G4"), nrow = 1)
  colnames(is.highgrade) <- ccRCC_clinical$Case_ID
  
  fisher_stat_tmp <- sapply(rownames(esscore_tab_outlier), FUN = function(pair, is.highgrade, outlier_df) {
    partIDs_nonNA <- colnames(outlier_df)[which(outlier_df[pair,] == TRUE | outlier_df[pair,] == FALSE)]
    partIDs_nonNA <- intersect(partIDs_nonNA, colnames(is.highgrade))
    tab4fisher <- table(outlier_df[pair, partIDs_nonNA], is.highgrade[,partIDs_nonNA])
    
    fisher_p.value <- NA; fisher_or_bottom <- NA; fisher_or_upper <- NA
    num_is.outlier_is.highgrade <- NA
    num_not.outlier_not.subtype <- NA
    num_is.outlier_not.subtype <- NA
    num_not.outlier_is.highgrade <- NA
    
    if (any(outlier_df[pair, partIDs_nonNA] & is.highgrade[,partIDs_nonNA])) {
      if (length(as.vector(tab4fisher)) == 4) {
        fisher_stat <- fisher.test(tab4fisher)
        fisher_p.value <- fisher_stat$p.value
        fisher_or_bottom <- as.numeric(fisher_stat$conf.int[1])
        fisher_or_upper <- as.numeric(fisher_stat$conf.int[2])
        
        num_is.outlier_is.highgrade = as.vector(tab4fisher)[4]
        num_not.outlier_not.subtype = as.vector(tab4fisher)[1]
        num_is.outlier_not.subtype = as.vector(tab4fisher)[2]
        num_not.outlier_is.highgrade = as.vector(tab4fisher)[3]
      } 
    }
    fisher_stat <- c(fisher_p.value, fisher_or_bottom, fisher_or_upper,
                     num_is.outlier_is.highgrade, num_not.outlier_not.subtype,
                     num_is.outlier_not.subtype, num_not.outlier_is.highgrade)
    return(fisher_stat)
  }, is.highgrade = is.highgrade, outlier_df = esscore_tab_outlier)
  
  fisher_stat_tab <- t(fisher_stat_tmp)
  colnames(fisher_stat_tab) <- c("p.value", "or_bottom", "or_upper", 
                                 "num_is.outlier_is.highgrade", "num_not.outlier_not.highgrade",
                                 "num_is.outlier_not.highgrade", "num_not.outlier_is.highgrade")
  fisher_stat_tab <- cbind(esscore_tab_outlier[,c("GENE", "SUB_GENE", "SUB_MOD_RSD", "SELF", "pair", "site_id")], fisher_stat_tab)
  
  fisher_stat_tab <- data.frame(fisher_stat_tab)
  fisher_stat_tab$cancer <- cancer
  fisher_stat_tab$clinical_info <- "highgrade"
  fisher_stat_tab$fdr_by_GENE <- FDR_by_id_columns(p_vector = fisher_stat_tab$p.value, id_columns = c("SELF", "GENE"), df = fisher_stat_tab)
  fisher_stat_tab$fdr <- FDR_by_id_columns(p_vector = fisher_stat_tab$p.value, id_columns = c("SELF", "cancer"), df = fisher_stat_tab)
  
  fn <- paste0(makeOutDir(resultD = resultD),
               "fisher_stat_tab_escore_outlier", outlier_sd, "SD_highgrade", cancer, "_", enzyme_type, "_reg_nonNA", reg_nonNA, ".txt")
  write.table(x = fisher_stat_tab, file = fn, 
              row.names = F, quote = F, sep = "\t")
}

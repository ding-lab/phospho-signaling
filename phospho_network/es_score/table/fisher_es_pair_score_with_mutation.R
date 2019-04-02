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
  
  esscore_tab_outlier_fn <- paste0("./cptac2p/analysis_results/phospho_network/es_score/table/calculate_es_pair_score/kinase/esscore_tab", "_outlier", outlier_sd, "SD_", cancer, "_", enzyme_type, "_reg_nonNA", reg_nonNA, ".txt")
  esscore_tab_outlier <- fread(input = esscore_tab_outlier_fn, data.table = F)
  partIDs_exp <- colnames(esscore_tab_outlier); partIDs_exp <- partIDs_exp[!(partIDs_exp %in% c("GENE", "SUB_GENE", "SUB_MOD_RSD", "SELF", "pair", "site_id"))]
  rownames(esscore_tab_outlier) <- as.vector(esscore_tab_outlier$pair)
  
  maf <- loadMaf(cancer = cancer, maf_files = maf_files)
  pair_tab <- unique(c(unlist(SMGs), esscore_tab_outlier$GENE, esscore_tab_outlier$SUB_GENE))
  mut_mat <- generate_somatic_mutation_matrix(pair_tab = pair_tab, maf = maf)
  genes2test <- unique(mut_mat$Hugo_Symbol)
    
  
  ## do fishers's table per mutated SMG gene to see subtype enrichment
  fisher_stat_tab <- NULL
  for (gene2test in genes2test) {
    is.genealtered <- mut_mat[mut_mat$Hugo_Symbol == gene2test, -1]
    is.genealtered <- (is.genealtered != "") & (is.genealtered != "Silent")
    partIDs_overlap <- intersect(partIDs_exp, colnames(is.genealtered))
    partIDs_genoaltered <- partIDs_overlap[is.genealtered[,partIDs_overlap]]
    if (length(partIDs_genoaltered) < 4) {
      next()
    }
    
    esscore_tab_outlier_tmp <- esscore_tab_outlier[rowSums(esscore_tab_outlier[, partIDs_genoaltered] == T, na.rm = T) >= 4,]
    if (nrow(esscore_tab_outlier_tmp) == 0) {
      next()
    }
    fisher_stat_tmp <- sapply(rownames(esscore_tab_outlier_tmp), FUN = function(pair, is.genealtered, outlier_df) {
      partIDs_nonNA <- colnames(outlier_df)[which(outlier_df[pair,] == TRUE | outlier_df[pair,] == FALSE)]
      partIDs_nonNA <- intersect(partIDs_nonNA, colnames(is.genealtered))
      tab4fisher <- table(outlier_df[pair, partIDs_nonNA], is.genealtered[,partIDs_nonNA])
      
      fisher_p.value <- NA; fisher_or_bottom <- NA; fisher_or_upper <- NA
      num_is.outlier_is.genealtered <- NA
      num_not.outlier_not.subtype <- NA
      num_is.outlier_not.subtype <- NA
      num_not.outlier_is.genealtered <- NA
      
      if (any(outlier_df[pair, partIDs_nonNA] & is.genealtered[,partIDs_nonNA])) {
        if (length(as.vector(tab4fisher)) == 4) {
          fisher_stat <- fisher.test(tab4fisher, alternative = "greater")
          fisher_p.value <- fisher_stat$p.value
          fisher_or_bottom <- as.numeric(fisher_stat$conf.int[1])
          fisher_or_upper <- as.numeric(fisher_stat$conf.int[2])
          
          num_is.outlier_is.genealtered = as.vector(tab4fisher)[4]
          num_not.outlier_not.subtype = as.vector(tab4fisher)[1]
          num_is.outlier_not.subtype = as.vector(tab4fisher)[2]
          num_not.outlier_is.genealtered = as.vector(tab4fisher)[3]
        } 
      }
      fisher_stat <- c(fisher_p.value, fisher_or_bottom, fisher_or_upper,
                       num_is.outlier_is.genealtered, num_not.outlier_not.subtype,
                       num_is.outlier_not.subtype, num_not.outlier_is.genealtered)
      return(fisher_stat)
    }, is.genealtered = is.genealtered, outlier_df = esscore_tab_outlier_tmp)
    
    fisher_stat_tmp <- t(fisher_stat_tmp)
    colnames(fisher_stat_tmp) <- c("p.value", "or_bottom", "or_upper", 
                                       "num_is.outlier_is.genealtered", "num_not.outlier_not.subtype",
                                       "num_is.outlier_not.subtype", "num_not.outlier_is.genealtered")
    fisher_stat_tmp <- cbind(esscore_tab_outlier_tmp[,c("GENE", "SUB_GENE", "SUB_MOD_RSD", "SELF", "pair", "site_id")], fisher_stat_tmp)
    fisher_stat_tmp$gene_altered <- gene2test
    fisher_stat_tab <- rbind(fisher_stat_tab, fisher_stat_tmp)
  }
  fisher_stat_tab <- data.frame(fisher_stat_tab)
  fisher_stat_tab$cancer <- cancer
  fisher_stat_tab$genoalt_type <- "mut"
  
  fn <- paste0(makeOutDir(resultD = resultD),
               "fisher_stat_tab_outlier", outlier_sd, "SD_", cancer, "_", enzyme_type, "_reg_nonNA", reg_nonNA, ".txt")
  write.table(x = fisher_stat_tab, file = fn, 
              row.names = F, quote = F, sep = "\t")
}


# clean up ----------------------------------------------------------------
enzyme_type <- "kinase"
fisher_stat_cancers <- NULL
for (cancer in c("BRCA", "CO", "UCEC", "OV", "CCRCC")) {
  fisher_stat_cancer <- fread(input = paste0(makeOutDir(resultD = resultD),
                                             "fisher_stat_tab_outlier", outlier_sd, "SD_", cancer, "_", enzyme_type, "_reg_nonNA", reg_nonNA, ".txt"))
  fisher_stat_cancer <- fisher_stat_cancer[!is.na(fisher_stat_cancer$p.value) & fisher_stat_cancer$p.value < 1,]
  fisher_stat_cancer <- unique(fisher_stat_cancer)
  for (SELF in c("cis", "trans")) {
    fisher_stat_cancer_self <- fisher_stat_cancer[fisher_stat_cancer$SELF == SELF,]
    for (gene_altered in unique(fisher_stat_cancer_self$gene_altered)) {
      fisher_stat_cancer_self_subtype <- fisher_stat_cancer_self[fisher_stat_cancer_self$gene_altered == gene_altered,]
      fisher_stat_cancer_self_subtype$fdr <- p.adjust(p = fisher_stat_cancer_self_subtype$p.value, method = "fdr")
      fisher_stat_cancers <- rbind(fisher_stat_cancers, fisher_stat_cancer_self_subtype)
    }
  }
}
fisher_stat_cancers <- unique(fisher_stat_cancers)
write.table(x = fisher_stat_cancers, file = paste0(makeOutDir(resultD = resultD),
                                                   "fisher_stat_tab_outlier", outlier_sd, "SD_all_cancers_", enzyme_type, "_reg_nonNA", reg_nonNA, ".txt"),
            row.names = F, quote = F, sep = "\t")

fisher_stat_cancers_fdr_sig <- fisher_stat_cancers[fisher_stat_cancers$fdr < 0.1 & fisher_stat_cancers$num_is.outlier_is.genealtered > 4,]
fisher_stat_cancers_sig <- fisher_stat_cancers[fisher_stat_cancers$p.value < 0.05 & fisher_stat_cancers$SELF == "trans" & fisher_stat_cancers$num_is.outlier_is.genealtered > 4,]

# Yige Wu @ WashU 2018 Aug
# calculate enzyme-substrate pair score for each sample for regulated pairs

# source ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
# inputs ------------------------------------------------------------------
## input enzyme substrate table
ptms_site_pairs_sup <- read_csv(paste0(ppnD, "compile_enzyme_substrate/tables/compile_omnipath_networkin_depod/omnipath_networkin_enzyme_substrate_site_level_union_source_cleaned_addDEPOD_extended.csv"))
sampmap <- loadSampMap()

# set variables -----------------------------------------------------------
fdr_pp <- 0.1
fdr_pk <- 0.05

# kinases -----------------------------------------------------------------
for( enzyme_type in c("kinase")) {
  regression <- fread(input = paste0("./cptac2p/analysis_results/phospho_network/regression/tables/regression_omnipath&newworkin&depod_protein_level/", enzyme_type, "_substrate_regression_cptac2p_3can_tumor.txt"), data.table = F)
  regression$GENE <- as.vector(regression$KINASE)
  regression$SUB_GENE <- as.vector(regression$SUBSTRATE)
  regression <- markSigSiteCan(regression = regression, sig_thres = fdr_pk, enzyme_type = enzyme_type)
  regression$regulated <- (regression$fdr_sig & regression$coef_sig)
  subdir1 <- paste0(makeOutDir(resultD = resultD), enzyme_type, "/")
  dir.create(subdir1)
  
  regression <- regression[regression$regulated,]
  for (cancer in cancers_sort) {
    regression_can <- regression[regression$Cancer == cancer,]
    phog <- fread(paste0(cptac_sharedD, cancer, "/", prefix[cancer], "_", "collapsed_PHO", "_formatted_normalized_replicate_averaged_Tumor.txt"),
                  data.table = F)
    pho <- fread(paste0(cptac_sharedD, cancer, "/", prefix[cancer], "_", "PHO", "_formatted_normalized_noControl.txt"),
                 data.table = F)
    
    pho_head  <- formatPhosphosite(phosphosite_vector = as.vector(pho$Phosphosite), gene_vector = as.vector(pho$Gene))
    pro <- fread(paste0(cptac_sharedD, cancer, "/", prefix[cancer], "_", "PRO", "_formatted_normalized_noControl.txt"),
                 data.table = F)
    samples <- colnames(phog); samples <- samples[!(samples == "Gene")]

    for (self in c("trans")) {
      regression_can_self <- regression_can[regression_can$SELF == self,]
      
      es_score_can <- NULL
      es_out_can <- NULL
      for (i in c(1:nrow(regression_can_self))) {
        enzyme <- regression_can_self$GENE[i]
        sub <- regression_can_self$SUB_GENE[i]
        rsd <- regression_can_self$SUB_MOD_RSD[i]
        
        pho_en <- phog[phog$Gene == enzyme, samples]
        pho_sub <- pho[pho_head$SUBSTRATE == sub & pho_head$SUB_MOD_RSD == rsd, samples]
        pho_sub <- pho_sub[1,]
        
        es_score <- pho_en+pho_sub
        es_score_med <- median(x = as.numeric(es_score), na.rm = T)
        es_score_sd <- sd(x = as.numeric(es_score), na.rm = T)
        
        es_score_can <- rbind(es_score_can, es_score)
        es_out_can <- rbind(es_out_can, (es_score > (es_score_med + 2*es_score_sd)))
      }
      tab_head <- data.frame(pair = as.vector(regression_can_self$pair))
      
      write.table(x = cbind(tab_head, es_score_can), file = paste0(subdir1, cancer, "_", self,  "_es_score_sampID.txt"), row.names = F, col.names = T, quote = F)
      write.table(x = cbind(tab_head, es_out_can), file = paste0(subdir1, cancer, "_", self,  "_es_score_outlier_sampID.txt"), row.names = F, col.names = T, quote = F)
      
      partIDs <- sampID2partID(sampleID_vector = colnames(es_score_can), sample_map = sampmap)
      colnames(es_score_can) <- partIDs
      write.table(x = cbind(tab_head, es_score_can), file = paste0(subdir1, cancer, "_", self,  "_es_score_partID.txt"), row.names = F, col.names = T, quote = F, sep = "\t")
      colnames(es_out_can) <- partIDs
      write.table(x = cbind(tab_head, es_out_can), file = paste0(subdir1, cancer, "_", self,  "_es_score_outlier_partID.txt"), row.names = F, col.names = T, quote = F, sep = "\t")
    }

    for (self in c("cis")) {
      regression_can_self <- regression_can[regression_can$SELF == self,]
      
      es_score_can <- NULL
      es_out_can <- NULL
      for (i in c(1:nrow(regression_can_self))) {
        enzyme <- regression_can_self$GENE[i]
        sub <- regression_can_self$SUB_GENE[i]
        rsd <- regression_can_self$SUB_MOD_RSD[i]
        
        pro_en <- pro[pro$Gene == enzyme, samples]
        pho_sub <- pho[pho_head$SUBSTRATE == sub & pho_head$SUB_MOD_RSD == rsd, samples]
        pho_sub <- pho_sub[1,]
        
        es_score <- pro_en+pho_sub
        es_score_med <- median(x = as.numeric(es_score), na.rm = T)
        es_score_sd <- sd(x = as.numeric(es_score), na.rm = T)
        
        es_score_can <- rbind(es_score_can, es_score)
        es_out_can <- rbind(es_out_can, (es_score > (es_score_med + 2*es_score_sd)))
      }
      tab_head <- data.frame(pair = as.vector(regression_can_self$pair))
      
      write.table(x = cbind(tab_head, es_score_can), file = paste0(subdir1, cancer, "_", self, "_es_score_sampID.txt"), row.names = F, col.names = T, quote = F)
      write.table(x = cbind(tab_head, es_out_can), file = paste0(subdir1, cancer, "_", self, "_es_score_outlier_sampID.txt"), row.names = F, col.names = T, quote = F)
      
      partIDs <- sampID2partID(sampleID_vector = colnames(es_score_can), sample_map = sampmap)
      colnames(es_score_can) <- partIDs
      write.table(x = cbind(tab_head, es_score_can), file = paste0(subdir1, cancer, "_", self,  "_es_score_partID.txt"), row.names = F, col.names = T, quote = F, sep = "\t")
      colnames(es_out_can) <- partIDs
      write.table(x = cbind(tab_head, es_out_can), file = paste0(subdir1, cancer, "_", self,  "_es_score_outlier_partID.txt"), row.names = F, col.names = T, quote = F, sep = "\t")
    }
    
  }
}

library(readxl)

subtype_map <- read_excel(paste0(cptac_sharedD, "5_CPTAC2_Breast_Prospective_Collection_BI/proteome-data-v1.01011/data/20171106_CPTAC2_ProspectiveBreastCancer_Sample_Annotations_Table_v79.xlsx"))

cancer <- "BRCA"


# Yige Wu @ WashU 2021 Sep
# test the trans correlation between individual kinase/phosphatase phosphosite with substrate phosphsoites

# set up libraries and output directory -----------------------------------
dir_base = "~/Box/Ding_Lab/Projects_Current/PanCan_Phospho-signaling/"
setwd(dir_base)
library(data.table)
library(stringr)
library(plyr)
library(dplyr)
library(doParallel)

# path2phospho_network_shared <- paste0(code_top_dir, "phospho_network/phospho_network_shared.R")
# source(path2phospho_network_shared)

# input dependencies ------------------------------------------------------
phosphosite_df <- fread(data.table = F, input = "~/Box/CPTAC_ccRCC/Data_Freeze_1.1/CPTAC_ccRCC_Discovery/Phosphoproteome/TMT/6_CPTAC3_CCRCC_Phospho_abundance_phosphosite_protNorm=2_CB_imputed.tsv")
protein_df <- fread(data.table = F, input = "~/Box/CPTAC_ccRCC/Data_Freeze_1.1/CPTAC_ccRCC_Discovery/Proteome/TMT/6_CPTAC3_CCRCC_Whole_abundance_gene_protNorm=2_CB.tsv")
ks_pair_df <- fread(data.table = F, input = "./Resources/Protein_Protein_Interactions/Omnipath_NetworKin_DEPOD_SignorNotSiteMapped_manualAdded.csv")
metadata_df <- fread(data.table = F, input = "~/Box/CPTAC_ccRCC/Data_Freeze_1.1/CPTAC_ccRCC_Discovery/Case_ID/CPTAC_ccRCC_discovery_caseID_v1.0.tsv")

# set variables -----------------------------------------------------------
least_samples <- 20# least number of samples with complete data for each model

# preprocess --------------------------------------------------------------
## get samples
samples <- metadata_df$Specimen.Label.tumor[metadata_df$Histologic_Type == "Clear cell renal cell carcinoma"]
samples <- intersect(samples, colnames(protein_df))
samples <- intersect(samples, colnames(phosphosite_df)) ## 103 tumors

## prepare phosphosite
phosphosite_idx_df <- data.frame(SUBSTRATE = phosphosite_df$Gene, id_phosphosite = phosphosite_df$Index)
phosphosite_idx_df <- phosphosite_idx_df %>%
  mutate(SUB_MOD_RSD = str_split_fixed(string = id_phosphosite, pattern = "_", n = 7)[,7])
phosphosite_value_df <- phosphosite_df[,samples] - phosphosite_df$ReferenceIntensity
phosphosite_idx_df$Number_nonna <- rowSums(!is.na(phosphosite_value_df))
## process protein data
protein_idx_df <- data.frame(Gene = protein_df$Index, id_protein = protein_df$Proteins)
protein_value_df <- protein_df[,samples] - protein_df$ReferenceIntensity
protein_idx_df$Number_nonna <- rowSums(!is.na(protein_value_df))

## get enzyme-substrate pairs with experimental evidence
pairs2test_df <- ks_pair_df %>%
  filter(Source != "NetKIN") %>%
  select(GENE, SUB_GENE, enzyme_type) %>%
  unique()
pairs2test_df <- merge(pairs2test_df, protein_idx_df, by.x = c("SUB_GENE"), by.y = c("Gene"))
pairs2test_df <- merge(pairs2test_df, phosphosite_idx_df, by.x = c("SUB_GENE"), by.y = c("SUBSTRATE"), suffixes = c(".pro_sub", ".pho_sub"))
pairs2test_df <- merge(pairs2test_df, 
                       phosphosite_idx_df %>%
                         mutate(ENZ_MOD_RSD = SUB_MOD_RSD) %>%
                         select(SUBSTRATE, ENZ_MOD_RSD, id_phosphosite), 
                       by.x = c("GENE"), by.y = c("SUBSTRATE"), suffixes = c(".pho_sub", ".pho_enz"))
pairs2test_df <- pairs2test_df %>%
  filter(SUB_MOD_RSD != "") %>%
  filter(ENZ_MOD_RSD != "") %>%
  filter(Number_nonna.pro_sub >= least_samples) %>%
  filter(SUB_GENE != GENE) %>%
  mutate(ENZ_phosphosite = paste0(GENE, "_", ENZ_MOD_RSD))
length(unique(pairs2test_df$GENE[pairs2test_df$enzyme_type == "kinase"])) ## 267 kinases
length(unique(pairs2test_df$GENE[pairs2test_df$enzyme_type == "phosphatase"])) ## 48 phosphatases

# initiate ----------------------------------------------------------------
# calculate the length of trans table
npairs2test <- nrow(pairs2test_df)

# initiating the table for trans
vec_char <- vector(mode = "character", length = npairs2test)
vec_num <- vector(mode = "numeric", length = npairs2test) + NA

FDR_pro_kin <- vec_num;FDR_pro_sub <- vec_num;FDR_pho_kin <- vec_num;
coef_pro_kin <- vec_num;coef_pro_sub <- vec_num;coef_pho_kin <- vec_num;
Cancer <- vec_char;
Transcript_sub <- vec_char; Transcript_enz <- vec_char;
model <- vec_char;
Size <- vec_num;
P_pro_kin <- vec_num;P_pro_sub <- vec_num;P_pho_kin <- vec_num;
sd_pro_kin <- vec_num;sd_pro_sub <- vec_num;sd_pho_kin <- vec_num; sd_pho_sub <- vec_num;

# foreach -----------------------------------------------------------------
idxs_proces <- 1:1000
registerDoParallel(cores = 4)
start_time <- Sys.time()
# test_list<-foreach(j=idxs_proces) %dopar% {
  test_list<-foreach(j=1:npairs2test) %dopar% {
  enzyme <- pairs2test_df[j, "GENE"]
  enz_mod_rsd <- pairs2test_df[j, "ENZ_MOD_RSD"]
  enz_phosphosite_id <- pairs2test_df[j, "id_phosphosite.pho_enz"]
  
  substrate <- pairs2test_df[j, "SUB_GENE"]
  sub_mod_rsd <- pairs2test_df[j, "SUB_MOD_RSD"]
  sub_phosphosite_id <- pairs2test_df[j, "id_phosphosite.pho_sub"]
  
  sub_protein_id <- pairs2test_df[j, "id_protein"]
  
  pho_enz <- phosphosite_value_df[phosphosite_idx_df$id_phosphosite == enz_phosphosite_id, ]
  pho_sub <- phosphosite_value_df[phosphosite_idx_df$id_phosphosite == sub_phosphosite_id, ]
  pro_sub <- protein_value_df[protein_idx_df$id_protein == sub_protein_id,]
  
  #prepare regression data for model2
  merged_data_df <- data.frame(t(rbind(pho_sub,pro_sub,pho_enz)))
  colnames(merged_data_df) <- c("pho_sub","pro_sub","pho_enz")
  data_complete <- merged_data_df[complete.cases(merged_data_df),]
  size <- nrow(data_complete)
  if(size > least_samples ){
    fit2 <- lm(pho_sub ~ pro_sub + pho_enz, data = merged_data_df)
    result_list <- c(coef(summary(fit2))[2,4], ##P_pro_sub
                     fit2$coefficients[2], ## coef_pro_sub
                     coef(summary(fit2))[3,4], ## P_pho_kin
                     fit2$coefficients[3], ## coef_pho_kin
                     size)
    return(result_list)
  } else {
    return(c(NA))
  }
}
end_time <- Sys.time()
end_time - start_time 
test_result_df <- do.call(rbind.data.frame, test_list)
colnames(test_result_df) <- c("P_pro_sub","coef_pro_sub", "P_pho_kin", "coef_pho_kin", "size")
test_result_df <- cbind(pairs2test_df[idxs_proces,], test_result_df)
test_result_df$FDR_pro_sub <- p.adjust(p = test_result_df$P_pro_sub, method = "fdr")
test_result_df$FDR_pho_kin <- p.adjust(p = test_result_df$P_pho_kin, method = "fdr")

test_result_df <- test_result_df[, c("GENE", "SUB_GENE", "SUB_MOD_RSD", 
                                     "FDR_pho_kin", "coef_pho_kin",
                                     "FDR_pro_sub", "coef_pro_sub", 
                                     "size", "P_pro_sub", "P_pho_kin", 
                                     "enzyme_type", "id_protein", "id_phosphosite.pho_sub", "id_phosphosite.pho_enz")]

# write output ------------------------------------------------------------
dir_out <- "./analysis_results/phospho_network/regression/tables/run_ks_regression_ccRCC_discovery_tumors/"
dir.create(dir_out)
file2write <- paste0(dir_out, "ccRCC_ks_trans.", "091521.", "tsv")
write.table(x = test_result_df, file = file2write, sep = "\t", row.names = F, quote = F)



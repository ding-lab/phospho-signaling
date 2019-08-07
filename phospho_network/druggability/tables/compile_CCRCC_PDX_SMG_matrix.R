# Yige Wu @ WashU 2019 July
## compile a table of SMG matrix for ccRCC PDXs

# source ------------------------------------------------------------------
source("./cptac2p_analysis/phospho_network/phospho_network_shared.R")


# set variables -----------------------------------------------------------
disease <- "CCRCC"
smg2search <- SMGs[[disease]]
smg2search
tissue_type2filter <- "PDX tumor"
cancer_subtype2filter <- "KIRC"

# input the sample info sheet to get the new IDs to search in maf files--------------------------
## read file
sample_info_tab <- read_excel("Ding_Lab/Projects_Current/ccRCC_drug/resources/PDXNet/summary_data_info/all_sample_info.v20190701.xlsx")
## filter by Cancer_Subtype and Tissue.v2
sample_info_filtered <- sample_info_tab %>%
  filter(Cancer_Subtype == cancer_subtype2filter) %>%
  filter(Tissue.v2 == tissue_type2filter)
  
# input the maf file ------------------------------------------------------
## read file
maf_tab <- fread(input = "./Ding_Lab/Projects_Current/ccRCC_drug/resources/PDXNet/summary_results/wxs.somaticMut/somaticMut_merged.20190624/u54.wxs_somaticMut_merged.maf.rc.caller", data.table = F)
## filter by SMG
maf_filtered <- maf_tab %>%
  filter(Hugo_Symbol %in% smg2search)
## add column for New_ID.v1
## filter by New_ID.v1

maf_filtered <- maf_filtered %>%
  mutate(New_ID.v1 = str_split_fixed(string = Tumor_Sample_Barcode, pattern = "_T", 2)[,1]) %>%
  filter(New_ID.v1 %in% sample_info_filtered$New_ID.v1) %>%
  select(New_ID.v1, Hugo_Symbol, HGVSp_Short)
maf_filtered_mat <- dcast(data = maf_filtered, formula = New_ID.v1 ~ Hugo_Symbol)
maf_filtered_df <- as.data.frame(maf_filtered_mat)

# compile the matrix ------------------------------------------------------
## New_ID.v1, ModelID, SMGs
pdx_smg_tab <- data.frame(New_ID.v1 = sample_info_filtered$New_ID.v1, ModelID = sample_info_filtered$ModelID)
pdx_smg_tab <- merge(pdx_smg_tab, maf_filtered_df, all.x = T, by = c("New_ID.v1"))
write.table(x = pdx_smg_tab, file = paste0(makeOutDir(resultD = resultD), "PDX_tumor_KIRC_SMG.tsv"), quote = F, row.names = F, sep = "\t")                     

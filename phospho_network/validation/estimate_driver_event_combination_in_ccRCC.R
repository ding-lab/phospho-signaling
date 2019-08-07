# Yige Wu @ WashU 2019 July
## estimate the frequency of driver events combination from CPTAC3 ccRCC cohort


# source ------------------------------------------------------------------
code_top_dir <- "~/Box/Ding_Lab/Projects_Current/PanCan_Phospho-signaling/phospho-signaling_analysis/"
path2phospho_network_shared <- paste0(code_top_dir, "phospho_network/phospho_network_shared.R")
source(path2phospho_network_shared)


# input maf file for ccRCC ------------------------------------------------
maf_tab <- loadMaf(cancer = "CCRCC", maf_files = maf_files)
maf_tab_smg <- maf_tab %>%
  filter(Hugo_Symbol %in% SMGs[["CCRCC"]]) %>%
  filter(Variant_Classification != "Silent") %>%
  select(Tumor_Sample_Barcode, Hugo_Symbol)

smg_mat <- dcast(data = maf_tab_smg, formula = Tumor_Sample_Barcode ~ Hugo_Symbol)
smg_mat <- data.frame(smg_mat)
smg_mat$driver_combination <- sapply(X = 1:nrow(smg_mat), FUN = function(i, df) {
  genotypes <- as.vector(df[i, 2:9])
  genotype_string <- paste0(genotypes, collapse = "_")
  return(genotype_string)
}, df = smg_mat)
table(smg_mat$driver_combination)
nrow(smg_mat)
maf_tab %>%
  select(Tumor_Sample_Barcode) %>%
  unique() %>%
  nrow()

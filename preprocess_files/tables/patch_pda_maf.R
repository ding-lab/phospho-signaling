# Yige Wu @ WashU 2020 Nov
## reference: https://www.nature.com/articles/bjc2014215

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/PanCan_Phospho-signaling/"
setwd(dir_base)
source("./phospho-signaling_analysis/load_pkgs.R")
source("./phospho-signaling_analysis/functions.R")
packages = c(
  "ggplot2",
  "ggrepel",
  "ComplexHeatmap",
  "circlize",
  "RColorBrewer"
)

for (pkg_name_tmp in packages) {
  library(package = pkg_name_tmp, character.only = T)
}

## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies --------------------------------------------
## input the mutation data
maf_df <- fread(data.table = F, input = "./Resources/PDA/Data/HTAN_PDAC_somatic_dnp.maf")
## input addtional mutation
mut_add_df <- fread(data.table = F, input = "./Resources/PDA/Data/Mutation_Additional.csv")

# combine -----------------------------------------------------------------
maf_sim_df <- maf_df %>%
  select(Hugo_Symbol, Tumor_Sample_Barcode, HGVSp_Short)
mut_add_df2 <- mut_add_df %>%
  rename(Hugo_Symbol = Gene) %>%
  mutate(Tumor_Sample_Barcode = paste0(Sample, "_T")) %>%
  rename(HGVSp_Short = AA_change) %>%
  dplyr::select(Hugo_Symbol, Tumor_Sample_Barcode, HGVSp_Short)

mut_combined_df <- rbind(maf_sim_df, mut_add_df2)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "HTAN_PDAC_somatic_dnp.patched.maf")
write.table(x = mut_combined_df , file = file2write, sep = "\t", quote = F, row.names = F)


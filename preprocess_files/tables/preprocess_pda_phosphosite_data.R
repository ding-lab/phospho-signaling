# Yige Wu @ WashU 2020 Nov

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/PanCan_Phospho-signaling/"
setwd(dir_base)
source("./phospho-signaling_analysis/load_pkgs.R")
source("./phospho-signaling_analysis/functions.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies --------------------------------------------
## input the mutation data
maf_df <- fread(data.table = F, input = "./Resources/PDA/Data/HTAN_PDAC_somatic_dnp.maf")
## input phosphosite data
phosite_index_df <- fread(data.table = F, input = "./Resources/PDA/Data/20200712-WU_HTAN-phosphoproteome-expression_matrices/WU_HTAN-phosphosite_matrix-log2_ratios-MD_norm_MAD_scaling.tsv")
phosite_df <- fread(data.table = F, input = "./Resources/PDA/Data/HTAN_PDA_phosphosite_matrix-log2_ratios-MD_norm_MAD_scaling_Formatted_NA20.txt")

# make unique id ----------------------------------------------------------
phosite_anno_df <- merge(x = phosite_df, 
                    y = phosite_index_df %>%
                      select(Flanking.Sequence.Phosphosites.6mer, Modifications, Phosphosite.Index),
                    by = c("Flanking.Sequence.Phosphosites.6mer", "Modifications"), all.x = T)
phosite_anno_df <- phosite_anno_df %>%
  mutate(Phosphosite.Index.short = str_split_fixed(string = Phosphosite.Index, pattern = "_", n = 7)[,7])
## reorder columns
colnames(phosite_anno_df)
colnames_dataid <- c("Gene", "Phosphosite.Index", "Flanking.Sequence.Phosphosites.6mer", "Modifications", "Phosphosite.Index.short")
colnames_sampleid <- colnames(phosite_anno_df)[!(colnames(phosite_anno_df) %in% colnames_dataid)]
colnames_sampleid
phosite_anno_df <- phosite_anno_df[, c(colnames_dataid, colnames_sampleid)]

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "HTAN_PDA_phosphosite_matrix-log2_ratios-MD_norm_MAD_scaling_Formatted_NA20.tsv")
write.table(x = phosite_anno_df, file = file2write, quote = F, sep = "\t", row.names = F)


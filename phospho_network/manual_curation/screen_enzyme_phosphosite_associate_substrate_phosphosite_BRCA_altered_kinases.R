# Yige Wu @ WashU 2019 July
## annotate the functions of substrate phosphosite for pairs associated with survival

# source ------------------------------------------------------------------
code_top_dir <- "~/Box/Ding_Lab/Projects_Current/PanCan_Phospho-signaling/phospho-signaling_analysis/"
path2phospho_network_shared <- paste0(code_top_dir, "phospho_network/phospho_network_shared.R")
source(path2phospho_network_shared)

# input regression result table --------------------------------------------------------
# regression <- fread(input = paste0(ppnD, "regression/tables/run_regression_with_enzyme_phosphosites_BRCA/", 
#                                    "regression_cptac3_BRCA_tumor_PGDAC_MD_MAD_nonNA20.txt"),
#                     data.table = F)

# regression <- fread(input = paste0(ppnD, "regression/tables/run_regression_with_enzyme_phosphosites_with_predicted_es_BRCA/", 
#                                    "regression_cptac3_BRCA_tumor_PGDAC_MD_MAD_nonNA20.txt"),
#                     data.table = F)

regression <- fread(input = paste0(ppnD, "regression/tables/run_regression_with_enzyme_phosphosites_BRCA_altered_kinases/", "regression_cptac2p_BRCA_tumor_CDAP_scaled_nonNA20.txt"),
                    data.table = F)

regression %>% nrow()
# regression <- regression %>%
#   mutate(GENE = KINASE) %>%
#   mutate(SUB_GENE = SUBSTRATE)
reg_sig["phosphatase"] <- 0.05
regression <- adjust_regression_by_nonNA(regression = regression, reg_nonNA = 20, reg_sig = reg_sig)
regression %>% nrow()
regression <- annotate_ks_source(regression)

# set variables -----------------------------------------------------------
gene2check <- "MAP2K4"
gene2check <- "AKT1"
gene2check <- "BRAF"
gene2check <- "ERBB2"
gene2check <- "EGFR"
gene2check <- "ATM"
gene2check <- "CDK12"
gene2check <- "MTOR"
gene2check <- "CHEK2"
gene2check <- "IGF1R"
gene2check <- "NTRK3"
gene2check <- "AKT3"
gene2check <- "CDK4"
gene2check <- "CDK6"


# select VHL related pairs ------------------------------------------------
regression2check <- regression %>%
  # filter(GENE == gene2check | SUB_GENE == gene2check) %>%
  filter(GENE == gene2check) %>%
  filter(Cancer == "BRCA") %>%
  filter(regulated == T) %>%
  mutate(residue = str_split_fixed(string = SUB_MOD_RSD, pattern = "[1-9]", n = 2)[,1]) %>%
  # filter(residue == "Y") %>%
  filter(residue != "Y") %>%
  # filter(ENZ_phosphosite.is_regulatory == T) %>%
  # filter(SUB_phosphosite.is_regulatory == T) %>%
  # filter(ENZ_phosphosite.is_regulatory == T | SUB_phosphosite.is_regulatory == T) %>%
  mutate(phosphosite = paste0(SUB_GENE, ":", SUB_MOD_RSD))


# Yige Wu @ WashU 2019 July
## annotate the functions of substrate phosphosite for pairs associated with survival


# source ------------------------------------------------------------------
code_top_dir <- "~/Box/Ding_Lab/Projects_Current/PanCan_Phospho-signaling/phospho-signaling_analysis/"
path2phospho_network_shared <- paste0(code_top_dir, "phospho_network/phospho_network_shared.R")
source(path2phospho_network_shared)

# input regression result table --------------------------------------------------------
regression <- fread(input = paste0("./Ding_Lab/Projects_Current/PanCan_Phospho-signaling/analysis_results/phospho_network/regression/tables/run_regression_with_enzyme_phosphosites_OV_altered_kinases/regression_cptac2p_OV_tumor_CDAP_scaled_nonNA20.20190820.v1.txt"),
                    data.table = F)

regression %>% nrow()
reg_sig["phosphatase"] <- 0.05
regression <- adjust_regression_by_nonNA(regression = regression, reg_nonNA = 20, reg_sig = reg_sig)
regression %>% nrow()
regression <- annotate_ks_source(regression)

# set variables -----------------------------------------------------------
gene2check <- "BRAF"
gene2check <- "CDK12"
gene2check <- "ATR"
gene2check <- "AKT1"
gene2check <- "AKT2"
gene2check <- "ATM"


# select VHL related pairs ------------------------------------------------
regression2check <- regression %>%
  # filter(GENE == gene2check | SUB_GENE == gene2check) %>%
  filter(GENE == gene2check) %>%
  filter(regulated == T) %>%
  mutate(residue = str_split_fixed(string = SUB_MOD_RSD, pattern = "[1-9]", n = 2)[,1]) %>%
  # filter(residue == "Y") %>%
  filter(residue != "Y") %>%
  # filter(ENZ_phosphosite.is_regulatory == T) %>%
  # filter(SUB_phosphosite.is_regulatory == T) %>%
  # filter(ENZ_phosphosite.is_regulatory == T | SUB_phosphosite.is_regulatory == T) %>%
  mutate(phosphosite = paste0(SUB_GENE, ":", SUB_MOD_RSD))

regression %>%
  # filter(GENE == gene2check | SUB_GENE == gene2check) %>%
  filter(GENE == gene2check) %>%
  filter(regulated == T) %>%
  mutate(residue = str_split_fixed(string = SUB_MOD_RSD, pattern = "[1-9]", n = 2)[,1]) %>%
  # filter(residue == "Y") %>%
  filter(residue != "Y") %>%
  # filter(ENZ_phosphosite.is_regulatory == T) %>%
  # filter(SUB_phosphosite.is_regulatory == T) %>%
  # filter(ENZ_phosphosite.is_regulatory == T | SUB_phosphosite.is_regulatory == T) %>%
  mutate(phosphosite = paste0(SUB_GENE, ":", SUB_MOD_RSD)) %>%
  select(ENZ_MOD_RSD, ENZ_phosphosite.is_regulatory) %>%
  unique()

regression %>%
  # filter(GENE == gene2check | SUB_GENE == gene2check) %>%
  filter(GENE == gene2check) %>%
  filter(regulated == T) %>%
  mutate(residue = str_split_fixed(string = SUB_MOD_RSD, pattern = "[1-9]", n = 2)[,1]) %>%
  # filter(residue == "Y") %>%
  filter(residue != "Y") %>%
  # filter(ENZ_phosphosite.is_regulatory == T) %>%
  # filter(SUB_phosphosite.is_regulatory == T) %>%
  # filter(ENZ_phosphosite.is_regulatory == T | SUB_phosphosite.is_regulatory == T) %>%
  mutate(phosphosite = paste0(SUB_GENE, ":", SUB_MOD_RSD)) %>%
  filter(ENZ_MOD_RSD == "S1981") %>%
  select(SUB_GENE, SUB_MOD_RSD, SUB_phosphosite.is_regulatory) %>%
  unique()


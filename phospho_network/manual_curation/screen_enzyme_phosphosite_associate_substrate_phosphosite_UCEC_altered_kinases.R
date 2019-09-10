# Yige Wu @ WashU 2019 July
## annotate the functions of substrate phosphosite for pairs associated with survival


# source ------------------------------------------------------------------
code_top_dir <- "~/Box/Ding_Lab/Projects_Current/PanCan_Phospho-signaling/phospho-signaling_analysis/"
path2phospho_network_shared <- paste0(code_top_dir, "phospho_network/phospho_network_shared.R")
source(path2phospho_network_shared)

# input regression result table --------------------------------------------------------
regression <- fread(input = paste0("./Ding_Lab/Projects_Current/PanCan_Phospho-signaling/analysis_results/phospho_network/regression/tables/run_regression_with_enzyme_phosphosites_UCEC_altered_kinases/regression_cptac3_UCEC_tumor_PGDAC_median_polishing_nonNA20.20190814.v1.txt"),
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
gene2check <- "FGFR1"
gene2check <- "FGFR2"
gene2check <- "FGFR3"
gene2check <- "FLT3"
gene2check <- "EGFR"
gene2check <- "ERBB2"
gene2check <- "ERBB3"
gene2check <- "JAK1"
gene2check <- "JAK2"
gene2check <- "MAP3K4"
gene2check <- "AKT1"
gene2check <- "MTOR"
gene2check <- "ATM"
gene2check <- "ATR"
gene2check <- "CHEK2"
gene2check <- "TGFBR2"

# select VHL related pairs ------------------------------------------------
regression2check <- regression %>%
  filter(GENE == gene2check | SUB_GENE == gene2check) %>%
  # filter(GENE == gene2check) %>%
  filter(Cancer == "UCEC") %>%
  filter(regulated == T) %>%
  mutate(residue = str_split_fixed(string = SUB_MOD_RSD, pattern = "[1-9]", n = 2)[,1]) %>%
  # filter(residue == "Y") %>%
  # filter(residue != "Y") %>%
  # filter(ENZ_phosphosite.is_regulatory == T) %>%
  # filter(SUB_phosphosite.is_regulatory == T) %>%
  # filter(ENZ_phosphosite.is_regulatory == T | SUB_phosphosite.is_regulatory == T) %>%
  mutate(phosphosite = paste0(SUB_GENE, ":", SUB_MOD_RSD))


regression %>%
  # filter(GENE == gene2check | SUB_GENE == gene2check) %>%
  filter(GENE == gene2check) %>%
  filter(Cancer == "UCEC") %>%
  filter(regulated == T) %>%
  mutate(residue = str_split_fixed(string = SUB_MOD_RSD, pattern = "[1-9]", n = 2)[,1]) %>%
  filter(residue == "Y") %>%
  # filter(residue != "Y") %>%
  # filter(ENZ_phosphosite.is_regulatory == T) %>%
  # filter(SUB_phosphosite.is_regulatory == T) %>%
  # filter(ENZ_phosphosite.is_regulatory == T | SUB_phosphosite.is_regulatory == T) %>%
  mutate(phosphosite = paste0(SUB_GENE, ":", SUB_MOD_RSD)) %>%
  select(ENZ_MOD_RSD, ENZ_phosphosite.is_regulatory) %>%
  unique()

regression %>%
  # filter(GENE == gene2check | SUB_GENE == gene2check) %>%
  filter(GENE == gene2check) %>%
  filter(Cancer == "UCEC") %>%
  filter(regulated == T) %>%
  mutate(residue = str_split_fixed(string = SUB_MOD_RSD, pattern = "[1-9]", n = 2)[,1]) %>%
  # filter(residue == "Y") %>%
  filter(residue != "Y") %>%
  # filter(ENZ_phosphosite.is_regulatory == T) %>%
  # filter(SUB_phosphosite.is_regulatory == T) %>%
  # filter(ENZ_phosphosite.is_regulatory == T | SUB_phosphosite.is_regulatory == T) %>%
  mutate(phosphosite = paste0(SUB_GENE, ":", SUB_MOD_RSD)) %>%
  filter(ENZ_MOD_RSD == "T1989") %>%
  select(SUB_GENE, SUB_MOD_RSD, SUB_phosphosite.is_regulatory) %>%
  unique()


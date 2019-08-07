# Yige Wu @ WashU 2019 July
## annotate the functions of substrate phosphosite for pairs associated with survival


# source ------------------------------------------------------------------
source("./cptac2p_analysis/phospho_network/phospho_network_shared.R")

# input regression result table --------------------------------------------------------
regression <- fread(input = paste0(ppnD, "regression/tables/annotate_regression_with_mut_impact/", 
                                   "regression_cptac2p_cptac3_tumor_reg_nonNA20_mut_impact_cancer_specificity_annotated.txt"),
                    data.table = F)
regression %>% nrow()
regression <- regression %>%
  filter(pair_pro %in% omnipath_tab$pair_pro[omnipath_tab$is.direct])
regression %>% nrow()
regression <- adjust_regression_by_nonNA(regression = regression, reg_nonNA = 20, reg_sig = reg_sig)
regression %>% nrow()
regression2merge <- regression %>%
  filter(regulated == T) %>%
  filter(GENE %in% SMGs[["CCRCC"]] | SUB_GENE %in% SMGs[["CCRCC"]]) %>%
  filter(Cancer == "CCRCC") %>%
  mutate(phosphosite = paste0(SUB_GENE, ":", SUB_MOD_RSD))

regression2merge <- regression %>%
  filter(regulated == T) %>%
  filter(GENE %in% c("KDR", "FLT1", "FLT4") | SUB_GENE %in% c("KDR", "FLT1", "FLT4")) %>%
  filter(Cancer == "CCRCC") %>%
  mutate(phosphosite = paste0(SUB_GENE, ":", SUB_MOD_RSD))
regression2merge$pair

# input regulatory phosphosites -------------------------------------------
regulatory_sites <- fread(input = "./Ding_Lab/Projects_Current/CPTAC/pan3can_shared_data/Phospho_databases/PhosphositePlus/May_28_2019/Regulatory_sites", data.table = F)
regulatory_sites <- regulatory_sites %>%
  mutate(SUB_MOD_RSD = str_split_fixed(string = MOD_RSD, pattern = "-", n = 2)[,1])
regulatory_sites2merge <- regulatory_sites %>%
  select(GENE, SUB_MOD_RSD, ON_FUNCTION, ON_PROCESS, ON_PROT_INTERACT, ON_OTHER_INTERACT) %>%
  mutate(phosphosite = paste0(GENE, ":", SUB_MOD_RSD)) %>%
  unique

regression2merged <- merge(regression2merge, regulatory_sites2merge, by = c("phosphosite"), all.x = T)

# input disease phosphosites ----------------------------------------------
disease_sites <- fread(input = "./Ding_Lab/Projects_Current/CPTAC/pan3can_shared_data/Phospho_databases/PhosphositePlus/Mar_04_2019/Disease-associated_sites", data.table = F)
disease_sites2merge <- disease_sites %>%
  mutate(SUB_MOD_RSD = str_split_fixed(string = MOD_RSD, pattern = "-", n = 2)[,1]) %>%
  mutate(phosphosite = paste0(GENE, ":", SUB_MOD_RSD)) %>%
  select(phosphosite, DISEASE, ALTERATION, ORGANISM, PMIDs)
regression2merged <- merge(regression2merged, disease_sites2merge, by = c("phosphosite"), all.x = T)

regression2merged %>%
  select(pair) %>%
  unique()

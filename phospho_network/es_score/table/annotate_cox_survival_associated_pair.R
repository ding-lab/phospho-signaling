# Yige Wu @ WashU 2019 Jun
## annotate the functions of substrate phosphosite for pairs associated with survival


# source ------------------------------------------------------------------
source("./cptac2p_analysis/phospho_network/phospho_network_shared.R")

# input table for associations -------------------------------------------------------------------
cox_stat_summary <- fread(input = "./cptac2p/analysis_results/phospho_network/es_score/table/run_survival_cox_with_esscore_demo/survival_cox_with_esscore_CCRCC_reg_nonNA20.txt", data.table = F)
cox_stat_summary2merge <- cox_stat_summary %>%
  filter(esscore_logrank_pvalue < 0.05) %>%
  mutate(phosphosite = str_split_fixed(string = pair, pattern = ":", n = 2)[,2])
  
# input regulatory phosphosites -------------------------------------------
regulatory_sites <- fread(input = "./Ding_Lab/Projects_Current/CPTAC/pan3can_shared_data/Phospho_databases/PhosphositePlus/May_28_2019/Regulatory_sites", data.table = F)
regulatory_sites <- regulatory_sites %>%
  mutate(SUB_MOD_RSD = str_split_fixed(string = MOD_RSD, pattern = "-", n = 2)[,1])
regulatory_sites2merge <- regulatory_sites %>%
  select(GENE, SUB_MOD_RSD, ON_FUNCTION, ON_PROCESS, ON_PROT_INTERACT, ON_OTHER_INTERACT) %>%
  mutate(phosphosite = paste0(GENE, ":", SUB_MOD_RSD)) %>%
  unique

cox_stat_summary2merged <- merge(cox_stat_summary2merge, regulatory_sites2merge, by = c("phosphosite"), all.x = T)

# input disease phosphosites ----------------------------------------------
disease_sites <- fread(input = "./Ding_Lab/Projects_Current/CPTAC/pan3can_shared_data/Phospho_databases/PhosphositePlus/Mar_04_2019/Disease-associated_sites", data.table = F)
disease_sites2merge <- disease_sites %>%
  mutate(SUB_MOD_RSD = str_split_fixed(string = MOD_RSD, pattern = "-", n = 2)[,1]) %>%
  mutate(phosphosite = paste0(GENE, ":", SUB_MOD_RSD)) %>%
  select(phosphosite, DISEASE, ALTERATION, ORGANISM, PMIDs)
cox_stat_summary2merged <- merge(cox_stat_summary2merged, disease_sites2merge, by = c("phosphosite"), all.x = T)


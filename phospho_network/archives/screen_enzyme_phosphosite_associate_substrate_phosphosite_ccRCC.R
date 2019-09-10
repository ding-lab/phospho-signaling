# Yige Wu @ WashU 2019 July
## annotate the functions of substrate phosphosite for pairs associated with survival


# source ------------------------------------------------------------------
code_top_dir <- "~/Box/Ding_Lab/Projects_Current/PanCan_Phospho-signaling/phospho-signaling_analysis/"
path2phospho_network_shared <- paste0(code_top_dir, "phospho_network/phospho_network_shared.R")
source(path2phospho_network_shared)

# input regression result table --------------------------------------------------------
# regression <- fread(input = paste0(ppnD, "regression/tables/run_regression_with_enzyme_phosphosites_ccRCC/", 
#                                    "regression_cptac3_CCRCC_tumor_PGDAC_MD_MAD_nonNA20.txt"),
#                     data.table = F)

# regression <- fread(input = paste0(ppnD, "regression/tables/run_regression_with_enzyme_phosphosites_with_predicted_es_ccRCC/", 
#                                    "regression_cptac3_CCRCC_tumor_PGDAC_MD_MAD_nonNA20.txt"),
#                     data.table = F)

regression <- fread(input = paste0(ppnD, "regression/tables/run_regression_with_enzyme_phosphosites_ccRCC_altered_kinases/", 
                                   "regression_cptac3_CCRCC_tumor_PGDAC_MD_MAD_nonNA20.txt"),
                    data.table = F)

regression %>% nrow()
regression <- regression %>%
  mutate(GENE = KINASE) %>%
  mutate(SUB_GENE = SUBSTRATE)
reg_sig["phosphatase"] <- 0.05
regression <- adjust_regression_by_nonNA(regression = regression, reg_nonNA = 20, reg_sig = reg_sig)
regression %>% nrow()
regression <- annotate_ks_source(regression)


# input protein pair table ------------------------------------------------
# pair_tab_annotated <- fread(input = "./cptac2p/analysis_results/dependencies/tables/compile_protein_pair_table/protein_pair_table.txt", data.table = F)
pair_tab_annotated <- fread(input = paste0(resultD, 
                                           "dependencies/tables/compile_protein_pair_table/protein_pair_table_v2.txt"), data.table = F)

# set variables -----------------------------------------------------------
gene2check <- "VHL"
gene2check <- "PBRM1"
gene2check <- "SETD2"
gene2check <- "BAP1"
gene2check <- "PTEN" 
gene2check <- "MET" 
gene2check <- "TP53"
gene2check <- "HIF1A"
gene2check <- "EPAS1"
gene2check <- "SHC1"
gene2check <- "GAB1"
gene2check <- "CTTN"
gene2check <- "CTNND1"


# filter for VHL-interacting genes ----------------------------------------
pair_tab_annotated_gene <- pair_tab_annotated %>%
  filter(GENE == gene2check)

# select VHL related pairs ------------------------------------------------
regression2merge <- regression %>%
  filter(regulated == T) %>%
  filter(GENE %in% pair_tab_annotated_gene$SUB_GENE | SUB_GENE %in% pair_tab_annotated_gene$SUB_GENE) %>%
  filter(Cancer == "CCRCC") %>%
  mutate(phosphosite = paste0(SUB_GENE, ":", SUB_MOD_RSD))
 
regression2merge %>%
  select(GENE) %>%
  unique()

regression2merge_ex <- regression %>%
  filter(FDR_pho_kin < 0.1 & coef_sig) %>%
  filter(GENE %in% pair_tab_annotated_gene$SUB_GENE | SUB_GENE %in% pair_tab_annotated_gene$SUB_GENE) %>%
  filter(Cancer == "CCRCC") %>%
  mutate(phosphosite = paste0(SUB_GENE, ":", SUB_MOD_RSD))

regression2check <- regression %>%
  filter(GENE == gene2check | SUB_GENE == gene2check) %>%
  filter(Cancer == "CCRCC") %>%
  mutate(phosphosite = paste0(SUB_GENE, ":", SUB_MOD_RSD))

# input regulatory phosphosites -------------------------------------------
regulatory_sites <- fread(input = "./Ding_Lab/Projects_Current/CPTAC/pan3can_shared_data/Phospho_databases/PhosphositePlus/May_28_2019/Regulatory_sites", data.table = F)
regulatory_sites <- regulatory_sites %>%
  mutate(SUB_MOD_RSD = str_split_fixed(string = MOD_RSD, pattern = "-", n = 2)[,1])
regulatory_sites2merge <- regulatory_sites %>%
  select(GENE, SUB_MOD_RSD, ON_FUNCTION, ON_PROCESS, ON_PROT_INTERACT, ON_OTHER_INTERACT) %>%
  mutate(phosphosite = paste0(GENE, ":", SUB_MOD_RSD)) %>%
  unique

regression2merged <- merge(regression2merge_ex, regulatory_sites2merge, by = c("phosphosite"), all.x = T)

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


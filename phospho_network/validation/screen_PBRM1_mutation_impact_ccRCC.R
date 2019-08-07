# Yige Wu @ WashU 2019 July
## annotate the functions of substrate phosphosite for pairs associated with survival


# source ------------------------------------------------------------------
code_top_dir <- "~/Box/Ding_Lab/Projects_Current/PanCan_Phospho-signaling/phospho-signaling_analysis/"
path2phospho_network_shared <- paste0(code_top_dir, "phospho_network/phospho_network_shared.R")
source(path2phospho_network_shared)

# input mut_impact_tab result table --------------------------------------------------------
mut_impact_tab <- fread(input = paste0(ppnD, "genoalt/tables/test_ccRCC_SMG_mut_impact_global/PBRM1_SETD2_BAP1_mut_impact_proteome_RNA_cptac2p_cptac3_tab.txt"),
                    data.table = F)
# input protein pair table ------------------------------------------------
# pair_tab_annotated <- fread(input = "./cptac2p/analysis_results/dependencies/tables/compile_protein_pair_table/protein_pair_table.txt", data.table = F)
pair_tab_annotated <- fread(input = paste0(resultD, 
                                           "dependencies/tables/compile_protein_pair_table/protein_pair_table_v2.txt"), data.table = F)


# load KEGG and REACTOME --------------------------------------------------
load(paste0(dir2cptac2_retrospective , "Gene_family/2015-08-01_Gene_Set.RData"))

# set variables -----------------------------------------------------------
gene2check <- "PBRM1"

summary(mut_impact_tab$meddiff)
quantile(x = mut_impact_tab$meddiff, probs = 0.95)
quantile(x = mut_impact_tab$meddiff, probs = 0.05)

# filter for PBRM1 complex partners  ----------------------------------------
pair_tmp <- pair_tab_annotated %>%
  filter(GENE == gene2check) %>%
  filter(SUB_GENE.is_complex_partner)

mut_impact_tab2check <- mut_impact_tab %>%
  filter(GENE == gene2check & SUB_GENE %in% pair_tmp$SUB_GENE) %>%
  filter(fdr < 0.05) %>%
  arrange(fdr)


# filter for JAK-STAT pathway ---------------------------------------------
reactome_names <- data.frame(REACT_name = names(REACT))
react_name_tmp <- "198745	Signalling to STAT3"
react_name_tmp <- "1280215	Cytokine Signaling in Immune system"
react_name_tmp <- "168256	Immune System"
react_name_tmp <- "2262749	Cellular response to hypoxia"
react_name_tmp <- "912694	Regulation of IFNA signaling"
REACT[[react_name_tmp]]
"JAK1" %in% REACT[[react_name_tmp]]
"IL12" %in% REACT[[react_name_tmp]]

mut_impact_tab2check <- mut_impact_tab %>%
  filter(GENE == gene2check & SUB_GENE %in% REACT[[react_name_tmp]]) %>%
  filter(fdr < 0.05) %>%
  arrange(fdr)

# write out gene list -----------------------------------------------------
for (exp_type_tmp in c("RNA", "PHO", "PRO")) {
  for (variant_type_tmp in c("not_silent")) {
    tab2w <- mut_impact_tab %>%
      filter(affected_exp_type == exp_type_tmp) %>%
      filter(variant_class == variant_type_tmp) %>%
      filter(fdr < 0.05)
    genes2w <- unique(tab2w$SUB_GENE)
    write.table(x = genes2w, file = paste0(makeOutDir(), exp_type_tmp, "_", variant_type_tmp, "_PBRM1_DE_genes.txt"), quote = F, col.names = F, row.names = F)
    
    tab2w <- mut_impact_tab %>%
      filter(affected_exp_type == exp_type_tmp) %>%
      filter(variant_class == variant_type_tmp)
    genes2w <- unique(tab2w$SUB_GENE)
    write.table(x = genes2w, file = paste0(makeOutDir(), exp_type_tmp, "_", variant_type_tmp, "_PBRM1_background_genes.txt"), quote = F, col.names = F, row.names = F)
    
    tab2w <- mut_impact_tab %>%
      filter(affected_exp_type == exp_type_tmp) %>%
      filter(variant_class == variant_type_tmp) %>%
      filter(meddiff > 0) %>%
      filter(fdr < 0.05)
    genes2w <- unique(tab2w$SUB_GENE)
    write.table(x = genes2w, file = paste0(makeOutDir(), exp_type_tmp, "_", variant_type_tmp, "_PBRM1_DE_up_genes.txt"), quote = F, col.names = F, row.names = F)
    
    tab2w <- mut_impact_tab %>%
      filter(affected_exp_type == exp_type_tmp) %>%
      filter(variant_class == variant_type_tmp) %>%
      filter(meddiff > 0.5) %>%
      filter(fdr < 0.05)
    genes2w <- unique(tab2w$SUB_GENE)
    write.table(x = genes2w, file = paste0(makeOutDir(), exp_type_tmp, "_", variant_type_tmp, "_PBRM1_DE_up_genes.txt"), quote = F, col.names = F, row.names = F)
    
    
    tab2w <- mut_impact_tab %>%
      filter(affected_exp_type == exp_type_tmp) %>%
      filter(variant_class == variant_type_tmp) %>%
      filter(meddiff < 0) %>%
      filter(fdr < 0.05)
    genes2w <- unique(tab2w$SUB_GENE)
    write.table(x = genes2w, file = paste0(makeOutDir(), exp_type_tmp, "_", variant_type_tmp, "_PBRM1_DE_down_genes.txt"), quote = F, col.names = F, row.names = F)
    
    tab2w <- mut_impact_tab %>%
      filter(affected_exp_type == exp_type_tmp) %>%
      filter(variant_class == variant_type_tmp) %>%
      filter(meddiff < -0.42) %>%
      filter(fdr < 0.05)
    genes2w <- unique(tab2w$SUB_GENE)
    write.table(x = genes2w, file = paste0(makeOutDir(), exp_type_tmp, "_", variant_type_tmp, "_PBRM1_DE_down_capped_genes.txt"), quote = F, col.names = F, row.names = F)
    
    
  }
}

# # input regulatory phosphosites -------------------------------------------
# regulatory_sites <- fread(input = "./Ding_Lab/Projects_Current/CPTAC/pan3can_shared_data/Phospho_databases/PhosphositePlus/May_28_2019/Regulatory_sites", data.table = F)
# regulatory_sites <- regulatory_sites %>%
#   mutate(SUB_MOD_RSD = str_split_fixed(string = MOD_RSD, pattern = "-", n = 2)[,1])
# regulatory_sites2merge <- regulatory_sites %>%
#   select(GENE, SUB_MOD_RSD, ON_FUNCTION, ON_PROCESS, ON_PROT_INTERACT, ON_OTHER_INTERACT) %>%
#   mutate(phosphosite = paste0(GENE, ":", SUB_MOD_RSD)) %>%
#   unique
# 
# mut_impact_tab2merged <- merge(mut_impact_tab2merge_ex, regulatory_sites2merge, by = c("phosphosite"), all.x = T)
# 
# # input disease phosphosites ----------------------------------------------
# disease_sites <- fread(input = "./Ding_Lab/Projects_Current/CPTAC/pan3can_shared_data/Phospho_databases/PhosphositePlus/Mar_04_2019/Disease-associated_sites", data.table = F)
# disease_sites2merge <- disease_sites %>%
#   mutate(SUB_MOD_RSD = str_split_fixed(string = MOD_RSD, pattern = "-", n = 2)[,1]) %>%
#   mutate(phosphosite = paste0(GENE, ":", SUB_MOD_RSD)) %>%
#   select(phosphosite, DISEASE, ALTERATION, ORGANISM, PMIDs)
# mut_impact_tab2merged <- merge(mut_impact_tab2merged, disease_sites2merge, by = c("phosphosite"), all.x = T)
# 
# mut_impact_tab2merged %>%
#   select(pair) %>%
#   unique()
# 

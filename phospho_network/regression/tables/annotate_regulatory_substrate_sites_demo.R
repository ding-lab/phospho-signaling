# Yige Wu @ WashU June 2019
## calculate the number of proteins invovled in regulated kinase-substrate pairs

# source ------------------------------------------------------------------
baseD = "/Users/yigewu/Box\ Sync/"
setwd(baseD)
source("./cptac2p_analysis/phospho_network/phospho_network_shared.R")
source("./cptac2p_analysis/phospho_network/phospho_network_plotting.R")

# set variables -----------------------------------------------------------
omnipath_tab <- annotate_ks_source(omnipath_tab)
reg_sig <- c(0.05, 0.05); names(reg_sig) <- c("kinase", "phosphatase")
drug_genes <- fread(input = "./Ding_Lab/Projects_Current/CPTAC/pan3can_shared_data/reference_files/gene_drug_list/Premed_raw_databases/drugBank/drug_list.txt", data.table = F, col.names = "gene")

# input regression --------------------------------------------------------
regression <- fread(input = paste0(ppnD, "regression/tables/annotate_regression_with_mut_impact/", 
                                   "regression_cptac2p_cptac3_tumor_reg_nonNA20_mut_impact_cancer_specificity_annotated.txt"),
                    data.table = F)
regression %>% nrow()
regression <- regression %>%
  filter(pair_pro %in% omnipath_tab$pair_pro[omnipath_tab$is.direct]) %>%
  mutate(phosphosite2merge = paste0(SUB_GENE, "_", SUB_MOD_RSD))
regression %>% nrow()
regression <- adjust_regression_by_nonNA(regression = regression, reg_nonNA = 20, reg_sig = reg_sig)

regression_gene_tmp <- regression %>%
  # filter(SUB_GENE == "MYC") %>%
  filter(GENE == "MET") %>%
  filter(regulated == T) 


# input regulatory_sites --------------------------------------------------
regulatory_sites <- fread(input = "./Ding_Lab/Projects_Current/CPTAC/pan3can_shared_data/Phospho_databases/PhosphositePlus/May_28_2019/Regulatory_sites", data.table = F)
regulatory_sites <- regulatory_sites %>%
  mutate(SUB_MOD_RSD = str_split_fixed(string = MOD_RSD, pattern = "-", n = 2)[,1])

regulatory_sites %>%
  filter(GENE == "")

regulatory_sites2merge <- regulatory_sites %>%
  select(GENE, SUB_MOD_RSD, ON_FUNCTION, ON_PROCESS, ON_PROT_INTERACT, ON_OTHER_INTERACT, NOTES, PMIDs, ORGANISM) %>%
  mutate(phosphosite = paste0(GENE, "_", SUB_MOD_RSD)) %>%
  mutate(phosphosite2merge = toupper(phosphosite)) %>%
  unique
regulatory_sites2merge$phosphosite <- NULL
regression_regulatory_sites <- merge(regression %>%
                                       filter(regulated == T), 
                                     regulatory_sites2merge, by = c("phosphosite2merge"))
regression_regulatory_cancer_tmp <- regression_regulatory_sites %>%
  # filter(Cancer == "BRCA") %>%
  # filter(Cancer == "OV") %>%
  # filter(Cancer == "CCRCC") %>%
  filter(Cancer == "CO") %>%
  filter(!is.na(Source)) %>%
  filter(KINASE %in% drug_genes$gene) %>%
  filter(SELF == "trans")

# input disease sites -----------------------------------------------------
disease_sites <- fread(input = "./Ding_Lab/Projects_Current/CPTAC/pan3can_shared_data/Phospho_databases/PhosphositePlus/May_28_2019/Disease-associated_sites", data.table = F)
disease_sites <- disease_sites %>%
  mutate(SUB_MOD_RSD = str_split_fixed(string = MOD_RSD, pattern = "-", n = 2)[,1])
disease_sites %>% head()
disease_sites2merge <- disease_sites %>%
  # select(GENE, SUB_MOD_RSD, ON_FUNCTION, ON_PROCESS, ON_PROT_INTERACT, ON_OTHER_INTERACT, NOTES, PMIDs, ORGANISM) %>%
  mutate(phosphosite = paste0(GENE, "_", SUB_MOD_RSD)) %>%
  mutate(phosphosite2merge = toupper(phosphosite)) %>%
  unique
disease_sites2merge$phosphosite <- NULL 




regression_disease_sites <- merge(regression %>%
                                    filter(regulated == T), 
                                  disease_sites2merge, by = c("phosphosite2merge"))
regression_disease_cancer_tmp <- regression_disease_sites %>%
  # filter(Cancer == "BRCA") %>%
  # filter(Cancer == "OV") %>%
  filter(Cancer == "CCRCC") %>%
  filter(KINASE %in% drug_genes$gene) %>%
  filter(SELF == "trans")

# input pair table --------------------------------------------------------
regulated_kinase_substrate_major_pathways <- read_excel("cptac2p/manuscripts/supplementary_tables/regulated_kinase_substrate_major_pathways.xlsx", 
                                                        sheet = "v3")


# split the pairs to multiple rows---------------------------------------------------------
## keep cancer type, REACTOME ID, pathway name
## then split the kinases and substrates
pairs_tmp <- sapply(X = regulated_kinase_substrate_major_pathways$Representative_kinase_substrate_pairs, FUN = function(pairs) unlist(strsplit(x = pairs, split = "; "))) %>%
  unlist()
num_pairs_tmp <- sapply(X = regulated_kinase_substrate_major_pathways$Representative_kinase_substrate_pairs, FUN = function(pairs) length(unlist(strsplit(x = pairs, split = "; ")))) %>%
  unlist()

regulated_kinase_substrate_major_pathways_melt <- regulated_kinase_substrate_major_pathways[rep(1:nrow(regulated_kinase_substrate_major_pathways), num_pairs_tmp),]
regulated_kinase_substrate_major_pathways_melt$pair <- pairs_tmp
pairs_tmp <- str_split_fixed(string = pairs_tmp, pattern = ":", n = 3)
pairs_tmp <- data.frame(pairs_tmp); colnames(pairs_tmp) <- c("GENE", "SUB_GENE", "SUB_MOD_RSD")
pairs_tmp <- pairs_tmp %>%
  mutate(phosphosite2merge = paste0(SUB_GENE, "_", SUB_MOD_RSD))
regulated_kinase_substrate_major_pathways_melt$Representative_kinase_substrate_pairs <- NULL
regulated_kinase_substrate_major_pathways_melt$SUB_MOD_RSD <- NULL
regulated_kinase_substrate_major_pathways_melt <- cbind(regulated_kinase_substrate_major_pathways_melt, pairs_tmp)
regulated_kinase_substrate_major_pathways_melt <- merge(regulated_kinase_substrate_major_pathways_melt, regulatory_sites2merge, by = c("phosphosite2merge", "SUB_MOD_RSD") , all.x = T)
regulated_kinase_substrate_major_pathways_melt <- regulated_kinase_substrate_major_pathways_melt %>%
  arrange(Cancer_Type, phosphosite2merge) %>%
  select(Cancer_Type, Reactome_Pathway_Name, pair, ORGANISM, ON_FUNCTION, ON_PROCESS, ON_PROT_INTERACT, ON_OTHER_INTERACT, NOTES, PMIDs)
regulated_kinase_substrate_major_pathways_melt[regulated_kinase_substrate_major_pathways_melt == ""] <- NA
write.table(x = regulated_kinase_substrate_major_pathways_melt, file = paste0(makeOutDir(resultD = resultD), "regulated_kinase_substrate_major_pathways.txt"), quote = F, row.names = F, sep = "\t")

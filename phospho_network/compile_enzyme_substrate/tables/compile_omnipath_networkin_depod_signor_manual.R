# Yige Wu @ WashU 2018 Apr
# compile protein-level and site-level ks pairs from Omnipath + networkin + DEPOD


# source ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
library(biomaRt)

# set variables -----------------------------------------------------------

# inputs ------------------------------------------------------------------
ptms_site_pairs_sup <- read_csv(paste0(ppnD, "compile_enzyme_substrate/tables/compile_omnipath_networkin_depod_signor/omnipath_networkin_DEPOD_SignorNotSiteMapped.csv"))
ptms_site_pairs_sup <- data.frame(ptms_site_pairs_sup)

## clean enzyme type
ptms_site_pairs_sup$enzyme_type[grepl(pattern = "DEPOD", ptms_site_pairs_sup$Source)] <- "phosphatase"

ptms_site_pairs_sup$pair <- NULL
# manual addition of some phosphatase pairs -------------------------------------------
ptms_pairs2add <- readxl::read_excel(path = "./Ding_Lab/Projects_Current/PhosphoDrug/PhosphoDrug_shared_data/kinase_substrate_relations/manual_added_enzyme_substrate_pairs.xlsx")

ptms_pairs2add <- data.frame(GENE = ptms_pairs2add$GENE, SUB_GENE = ptms_pairs2add$SUB_GENE,SUB_MOD_RSD = ptms_pairs2add$SUB_MOD_RSD,
                            KINASE = ptms_pairs2add$GENE, KIN_ACC_ID = "NULL", KIN_ORGANISM = "human",
                            SUBSTRATE = ptms_pairs2add$SUB_GENE, SUB_GENE_ID = "NULL", SUB_ACC_ID = "NULL", SUB_ORGANISM = "human",
                            SITE_GRP_ID = "NULL", SITE_...7_AA = "NULL", Source = ptms_pairs2add$Source, networkin_score = Inf, enzyme_type = ptms_pairs2add$enzyme_type)
ptms_pairs2add$pair_pro <- paste0(ptms_pairs2add$GENE, ":", ptms_pairs2add$SUB_GENE)
ptms_site_pairs_sup <- rbind(ptms_pairs2add, ptms_site_pairs_sup)
write.table(x = ptms_site_pairs_sup, file = paste0(makeOutDir(resultD = resultD), "Omnipath_NetworKin_DEPOD_SignorNotSiteMapped_manualAdded.csv"), quote = F, row.names = F, sep = ",")


# Write table for pairs with SMGs as substrates ---------------------------
smg_as_sub_tab <- ptms_site_pairs_sup[ptms_site_pairs_sup$SUB_GENE %in% unlist(SMGs),]
smg_as_sub_tab <- smg_as_sub_tab[order(smg_as_sub_tab$SUB_GENE, smg_as_sub_tab$GENE),]
write.table(x = smg_as_sub_tab, file = paste0(makeOutDir(resultD = resultD), "Omnipath_NetworKin_DEPOD_SignorNotSiteMapped_manualAdded_SMG_as_substrate.csv"), quote = F, row.names = F, sep = ",")

unlist(SMGs)[!(unlist(SMGs) %in% smg_as_sub_tab$SUB_GENE)]

# Yige Wu @ WashU 2018 Aug
## generate a list of focal amplified/deleted genes for each cancer from scores.gistic


# source ------------------------------------------------------------------
source("Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R")

# set variables -----------------------------------------------------------
## set threshold for GISTIC2.0 q-value
gistic_qv_thres <- 0.25
## set threshold for the proportion of overlap between segment and gene to the gene length to be attribute the segment statistics to the gene
overlap_cutoff <- 0.25

# inputs ------------------------------------------------------------------
ptms_site_pairs_sup <- read_csv(paste0(ppnD, "compile_enzyme_substrate/tables/compile_omnipath_networkin_depod_signor_manual/Omnipath_NetworKin_DEPOD_SignorNotSiteMapped_manualAdded.csv"))
ptms_site_pairs_sup <- ptms_site_pairs_sup[ptms_site_pairs_sup$Source != "NetKIN" | (ptms_site_pairs_sup$Source == "NetKIN" & ptms_site_pairs_sup$networkin_score >= 5),]
ptms_site_pairs_sup <- data.frame(ptms_site_pairs_sup)
kinases <- unique(ptms_site_pairs_sup$GENE[ptms_site_pairs_sup$enzyme_type == "kinase"])
phosphatases <- unique(ptms_site_pairs_sup$GENE[ptms_site_pairs_sup$enzyme_type == "phosphatase"])

for (cancer in cancers_sort) {
  can_id <- tolower(substr(x = cancer, start = 1, stop = 2))

  for (cna_type in c("amp", "del")) {
    scores_gistic_gene_tab <- fread(input = paste0(cptac2pD, "copy_number/gistic/outputs/cptac2p.", can_id, "/somatic/run_gistic2/0.25/0.99/armpeel1/somatic_clean/", cna_type, "_genes.conf_99.txt"))
    genes <- unlist(scores_gistic_gene_tab[5:nrow(scores_gistic_gene_tab), 2:ncol(scores_gistic_gene_tab)], use.names = F)
    genes <- genes[!is.na(genes) & genes!= "" & genes != "NA"]
    
    print(genes[genes %in% c(kinases, phosphatases)])
    
  }
}
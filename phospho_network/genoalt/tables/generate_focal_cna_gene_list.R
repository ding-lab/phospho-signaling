# Yige Wu @ WashU 2018 Aug
## generate a list of focal amplified/deleted genes for each cancer from scores.gistic


# source ------------------------------------------------------------------
source("./cptac2p_analysis/phospho_network/phospho_network_shared.R")

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
  scores_gistic_gene_tab <- fread(input = paste0(cptac2pD, "copy_number/gistic/outputs/cptac2p.", can_id, "/somatic/run_gistic2/0.25/0.99/armpeel1/somatic_clean/scores.gistic.gene.txt"))
  ## calculate the proportion of overlap between segment and gene to the gene length
  scores_gistic_gene_tab$overlap2gene <- scores_gistic_gene_tab$length.overlap/scores_gistic_gene_tab$length.gene
  ## calcualte FDR
  scores_gistic_gene_tab$FDR <- 10^(-scores_gistic_gene_tab$`-log10(q-value)`)
  scores_gistic_gene_tab$`-log10(q-value)` <- NULL
  fn <- paste0(makeOutDir(resultD = resultD), cancer, "_scores.gistic.gene.txt")
  write.table(x = scores_gistic_gene_tab, file = fn, quote = F, row.names = F, sep = '\t')
  
  for (cna_type in c("Amp", "Del")) {
    ## filter table
    tab2w <- scores_gistic_gene_tab
    # print(nrow(tab2w))
    tab2w <- tab2w[tab2w$Type == cna_type,]
    # print(nrow(tab2w))
    tab2w <- tab2w[tab2w$FDR < gistic_qv_thres,]
    # print(nrow(tab2w))
    tab2w <- tab2w[tab2w$overlap2gene > overlap_cutoff,]
    # print(nrow(tab2w))
    # print(length(unique(tab2w$symbol.gene)))
    
    ## write filtered table,
    fn <- paste0(makeOutDir(resultD = resultD), cancer, "_", cna_type, "_gistic_qv_thres", gistic_qv_thres, "overlap_cutoff", overlap_cutoff, "_scores.gistic.gene.txt")
    write.table(x = tab2w, file = fn, quote = F, row.names = F, sep = '\t')
    tab2w_pk <- tab2w[tab2w$symbol.gene %in% kinases,]
    tab2w_pk$enzyme_type <- "kinase"
    
    tab2w_pp <- tab2w[tab2w$symbol.gene %in% phosphatases,]
    tab2w_pp$enzyme_type <- "phosphatase"
    
    tab2w <- rbind(tab2w_pp, tab2w_pk)
    fn <- paste0(makeOutDir(resultD = resultD), cancer, "_", cna_type, "_ptm_enzyme_gistic_qv_thres", gistic_qv_thres, "overlap_cutoff", overlap_cutoff, "_scores.gistic.gene.txt")
    write.table(x = tab2w, file = fn, quote = F, row.names = F, sep = '\t')
    print(table(unique(tab2w[, c("enzyme_type", "symbol.gene")])[, "enzyme_type"]))
  }
}
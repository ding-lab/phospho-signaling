# Yige Wu @ WashU 2018 Apr
# output enzyme-substrate pairs in hallmark pathways

source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
load("~/Box Sync/pan3can_shared_data/Gene_family/2015-08-01_Gene_Set.RData")
keywords <- c("adherens","Mismatch","EMT", "apoptosis","immunological","stromal","transmembrane","receptors","integrin",
              "TGFβ","LKB1","AMPK","TSC","mTOR","PI3K","Akt","Ras","MAPK","Notch","Wnt", "catenin", "cell cycle", "p53","RTK","erbb")


# get kegg pathways -------------------------------------------------------
all_keggs <- names(KEGG)
keggs <- NULL
for (key in keywords) {
  keggs <- c(keggs, all_keggs[grepl(key, all_keggs, ignore.case = T)])
}
keggs <- unique(keggs)


sup_tab <- read.table(file = paste0("./cptac2p/cptac_shared/analysis_results/phospho_network/kinase_activity/tables/integrate_enzyme_KSEA_diffexp_sub_FC_plus_cistrans_regression_plus_genoalt/", 
                                    "3can", "_KSEA_enzyme_diffexp_substrate_aa_reported.txt"), header = T)

for (sig in c(0.1)) {
  for (kegg in keggs[13]) {
    ## keep pairs whose enzyme and substrate are both in the samme pathway
    kegg_genes <- KEGG[[kegg]]
    kegg_tab <- sup_tab[sup_tab$GENE %in% kegg_genes & sup_tab$SUB_GENE %in% kegg_genes,]
    kegg_tab$reg_sig <- NA
    
    ## keep at most 1 in-significant edge between two proteins
    kegg_cis <- kegg_tab[!is.na(kegg_tab$FDR_pro_kin),]
    kegg_cis$reg_sig <- ifelse(kegg_cis$FDR_pro_kin.aa_reported < sig & kegg_cis$coef_pro_kin > 0, TRUE, FALSE)
    kegg_cis_sig <- kegg_cis[kegg_cis$reg_sig,]
    kegg_cis_nonsig <- kegg_cis[!kegg_cis$reg_sig,];kegg_cis_nonsig <- kegg_cis_nonsig[order(kegg_cis_nonsig$GENE, kegg_cis_nonsig$SUB_GENE, kegg_cis_nonsig$FDR_pro_kin.aa_reported, decreasing = F),]
    kegg_cis_nonsig_nodup <- kegg_cis_nonsig[!duplicated(kegg_cis_nonsig[, c("GENE", "SUB_GENE")]),]
    kegg_trans <- kegg_tab[!is.na(kegg_tab$FDR_pho_kin),]
    kegg_trans$reg_sig <- ifelse((kegg_trans$enzyme_type == "kinase" & kegg_trans$FDR_pho_kin.aa_reported < sig & kegg_trans$coef_pho_kin > 0) | (kegg_trans$enzyme_type == "phosphotase" & kegg_trans$FDR_pho_kin.aa_reported < sig), TRUE, FALSE)
    
    kegg_trans_sig <- kegg_trans[kegg_trans$reg_sig,]
    kegg_trans_nonsig <- kegg_trans[!kegg_trans$reg_sig,];kegg_trans_nonsig <- kegg_trans_nonsig[order(kegg_trans_nonsig$GENE, kegg_trans_nonsig$SUB_GENE, kegg_trans_nonsig$FDR_pho_kin.aa_reported, decreasing = F),]
    kegg_trans_nonsig_nodup <- kegg_trans_nonsig[!duplicated(kegg_trans_nonsig[, c("GENE", "SUB_GENE")]),]
    kegg_nonreg <- kegg_tab[is.na(kegg_tab$FDR_pro_kin) & is.na(kegg_tab$FDR_pho_kin),]
    
    kegg_tab4cys <- rbind(kegg_nonreg, kegg_cis_sig, kegg_trans_sig, kegg_cis_nonsig_nodup, kegg_trans_nonsig_nodup)

    ## add in identifier for phosphosites
    kegg_tab4cys$phosphosite <- paste0(kegg_tab4cys$SUB_GENE, ":", kegg_tab4cys$SUB_MOD_RSD)
    
    write.table(x = kegg_tab4cys[!is.na(kegg_tab4cys$substrate_log2FC),],
                file = paste0("./cptac2p/cptac_shared/analysis_results/phospho_network/signal_hub/tables/",
                              kegg, "_KSEA_enzyme_diffexp_substrate_aa_reported.txt"),
                row.names = F, col.names = T, quote = F, sep = "\t")
  }
  
}

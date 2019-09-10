# Yige Wu @ WashU 2018 Feb
# annotate regression kinase/substrate with a few pathways implicated in cancer
# based on : https://www.nature.com/articles/ncomms4887#s1
# The function space covered by the antibodies used in the RPPA analysis includes proliferation, DNA damage, polarity, vesicle function, EMT, invasiveness, hormone signalling, apoptosis, metabolism, immunological and stromal function as well as transmembrane receptors, integrin, TGFβ, LKB1/AMPK, TSC/mTOR, PI3K/Akt, Ras/MAPK, Hippo, Notch and Wnt/beta-catenin signalling. Thus, the function space encompasses major functional and signalling pathways of relevance to human cancer. 

# input -------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
resultDnow <- makeOutDir()
load("~/Box Sync/pan3can_shared_data/analysis_results/2015-08-01_Gene_Set.RData")
keywords <- c("adherens","Mismatch","EMT", "apoptosis","immunological","stromal","transmembrane","receptors","integrin",
              "TGFβ","LKB1","AMPK","TSC","mTOR","PI3K","Akt","Ras","MAPK","Notch","Wnt", "catenin", "cell cycle", "p53","RTK","erbb")


# get kegg pathways -------------------------------------------------------
all_keggs <- names(KEGG)
keggs <- NULL
for (key in keywords) {
  keggs <- c(keggs, all_keggs[grepl(key, all_keggs, ignore.case = T)])
}
keggs <- unique(keggs)


# get genes within regression table and annotate with above pathway --------
## get genes covered by regression kinase and substrates
for (protein in c("kinase")) {
  tn = paste0(resultD, "signal_hub/tables/overlap_diffexp_regression/", protein,"_substrate_regression_overlap_cptac2p_3can_diffPhosphositeAnnotated.txt")
  table_3can <- fread(tn, data.table = F)
  genes <- unique(c(table_3can$KINASE, table_3can$SUBSTRATE))
  genesinpath <- data.frame(matrix(data = NA, nrow = length(genes), ncol = length(keggs)))
  rownames(genesinpath) <- genes; colnames(genesinpath) <- keggs
  for (kegg in keggs) {
    genesinpath[, kegg] <- (genes %in% KEGG[[kegg]])
  }
  genesinpath$Gene <- genes
  table_3can <- markSigSiteCan(table_3can, sig_thres = sig, protein_type = protein)
  table_3can <- merge(table_3can, genesinpath, by.x = c("SUBSTRATE"), by.y = c("Gene"), all.x = T)
  table_3can$num_keggs <- rowSums(table_3can[, grepl("hsa", colnames(table_3can))])
  for (cancer in c("BRCA")) {
    table_can_sig <- table_3can[table_3can$Cancer==cancer & table_3can$fdr_sig & table_3can$coef_sig,]
  }
}

num_keggs_by_kinase <- table_can_sig %>%
  group_by(KINASE) %>%
  summarise(num_keggs_by_kinase = sum(num_keggs>0))
num_keggs_by_kinase <- data.frame(num_keggs_by_kinase)
num_keggs_by_kinase <- num_keggs_by_kinase[order(num_keggs_by_kinase$num_keggs_by_kinase, decreasing = T),]

multi_keggs_kinases <- as.vector(num_keggs_by_kinase$KINASE[num_keggs_by_kinase$num_keggs_by_kinase > 0])
for (kinase in multi_keggs_kinases) {
  tmp <- table_can_sig[table_can_sig$KINASE == kinase,]
  tmp_diff <- tmp[!is.na(tmp$diffexp_type),]
  tmp_diff_kegg <- tmp_diff[tmp_diff$num_keggs > 0,]
  if (nrow(tmp_diff_kegg) > 1) {
    print(tmp_diff_kegg)
  }
}


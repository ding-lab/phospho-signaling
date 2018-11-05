# Yige Wu @ WashU 2018 Oct
## parse outputs from multiple runs of mutational impact

# source ------------------------------------------------------------------
source("Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R")
source("./cptac2p_analysis/phospho_network/phospho_network_plotting.R")

library(dplyr)


# inputs ------------------------------------------------------------------
mut_cis_pro <- fread(input = paste0(ppnD, "genoalt/tables/test_mut_cna_impact_cis_pro_cptac3/UCEC_mut_cnv_tab_formatted.txt"), data.table = F)
mut_pho <- fread(input = paste0(ppnD, "genoalt/tables/test_mut_cna_impact_UCEC/mut_cnv_sig_cans.txt"), data.table = F)

## input enzyme-substrate table
ptms_site_pairs_sup <- read_csv(paste0(ppnD, "compile_enzyme_substrate/tables/compile_omnipath_networkin_depod_signor_manual/Omnipath_NetworKin_DEPOD_SignorNotSiteMapped_manualAdded.csv"))
ptms_site_pairs_sup <- ptms_site_pairs_sup[ptms_site_pairs_sup$Source == "NetKIN" | (ptms_site_pairs_sup$Source != "NetKIN" & ptms_site_pairs_sup$networkin_score >= 5),]
ptms_site_pairs_sup <- data.frame(ptms_site_pairs_sup)
enzyme_sub <- unique(ptms_site_pairs_sup[, c("GENE", "SUB_GENE")])

# adjust pvalue for protein events -----------------------------------------------------------
# fdr <- vector(mode = "numeric", length = nrow(mut_cis_pro))
# for (cancer in unique(mut_cis_pro$cancer)) {
#   rows_tmp <- which(mut_cis_pro$cancer == cancer & mut_cis_pro$GENE %in% SMGs[[cancer]])
#   fdr[rows_tmp] <- p.adjust(p = mut_cis_pro$p[rows_tmp], method = "fdr")
# }
# mut_cis_pro$fdr <- fdr

fdr <- vector(mode = "numeric", length = nrow(mut_cis_pro))
for (cancer in unique(mut_cis_pro$cancer)) {
  rows_tmp <- which(mut_cis_pro$cancer == cancer)
  fdr[rows_tmp] <- p.adjust(p = mut_cis_pro$p[rows_tmp], method = "fdr")
}
mut_cis_pro$fdr <- fdr

# adjust pvalue for phosphorylation events -----------------------------------------------------------
fdr <- vector(mode = "numeric", length = nrow(mut_pho))
for (cancer in unique(mut_pho$cancer)) {
  for (self in c("cis", "trans")) {
    rows_tmp <- which(mut_pho$cancer == cancer & mut_pho$SELF == self)
    fdr[rows_tmp] <- p.adjust(p = mut_pho$p[rows_tmp], method = "fdr")
  }
}
mut_pho$fdr <- fdr

# write out annotated phosphorylation events -------------------------------------------------------------------
mut_pho <- merge(mut_pho, mut_cis_pro[, c("GENE", "cancer", "p", "meddiff_bottom", "meddiff_upper", "meddiff", "num", "fdr")],
                 by = c("GENE", "cancer"), suffixes = c("", ".pro"), all.x = T)
write.table(x = mut_pho, file = paste0(makeOutDir(resultD = resultD), "cptac2p_mut_cnv.txt"), quote = F, sep = '\t', row.names = F)

# write out cis phosphorylation events ------------------------------------
#Among XX cis-regulation events significantly associated with its own mutations (P-value < 0.05), XX (XX%) of them are tested by at least 10 non-synonymous mutations, which indicate the limited statistical power
tab2w <- mut_pho
tab2w <- tab2w[tab2w$SELF == "cis",]
nrow(tab2w)
tab2w <- tab2w[tab2w$p < 0.05,]
nrow(tab2w)
316/3644
length(unique(tab2w$GENE))
nrow(tab2w[tab2w$fdr < 0.1 & tab2w$fdr > 0,])
tab2w[tab2w$fdr < 0.1 & tab2w$fdr > 0,]

## get the list of proteins whose protein abundance is associated with its own somatic mutations
tab2add <- mut_cis_pro[, c("GENE", "cancer", "p", "meddiff_bottom", "meddiff_upper", "meddiff", "num")]
tab2w <- merge(tab2w, tab2add, by = c("GENE", "cancer"), suffixes = c("", ".pro"), all.x = T)
nrow(tab2w)
tab2w <- tab2w[order(tab2w$num),]

num_mut2test <- 4:10
num_pair <- NULL
for (i in num_mut2test) {
  num_pair <- c(num_pair, nrow(tab2w[tab2w$num >= i,]))
}
c <- data.frame(num_mut = num_mut2test, num_pair = num_pair)

p <- ggplot()
p <- p + geom_point(data = tab2p, mapping = aes(x = num_mut, y = num_pair))
p <- p + theme_bw()
p
ggsave(filename = paste0(makeOutDir(resultD = resultD), "number_of_mutation_vs_sig_pairs.pdf"), width = 4, height = 4)
149/272
37/272

## Among XX events with available protein abundance data, XX (%) co-occurred with significant protein abundance changes while the rest XX are unique to the phosphorylation level
nrow(tab2w[!is.na(tab2w$p.pro),])
nrow(tab2w[!is.na(tab2w$p.pro) & tab2w$p.pro < 0.05,])
nrow(tab2w[!is.na(tab2w$p.pro) & tab2w$p.pro < 0.05,])/nrow(tab2w[!is.na(tab2w$p.pro),])

table(tab2w$cancer)
for (cancer in cancers_sort) {
  print(tab2w$GENE[tab2w$GENE %in% SMGs[[cancer]] & tab2w$cancer == cancer])
}

# examine trans phosphorylation events ------------------------------------
tab2w <- mut_pho[mut_pho$SELF == "trans",]
tab2w$pair_pro <- paste0(tab2w$GENE, ":", tab2w$SUB_GENE)
nrow(tab2w)
tab2w <- tab2w[tab2w$p < 0.05,]
nrow(tab2w)

nrow(tab2w[tab2w$pair_pro %in% ptms_site_pairs_sup$pair_pro,])

num_mut2test <- 4:10
num_pair <- NULL
for (i in num_mut2test) {
  num_pair <- c(num_pair, nrow(tab2w[tab2w$num >= i,]))
}
tab2p <- data.frame(num_mut = num_mut2test, num_pair = num_pair)

p <- ggplot()
p <- p + geom_point(data = tab2p, mapping = aes(x = num_mut, y = num_pair))
p <- p + theme_bw()
p
ggsave(filename = paste0(makeOutDir(resultD = resultD), "number_of_mutation_vs_sig_pairs_trans.pdf"), width = 4, height = 4)

nrow(tab2w[!is.na(tab2w$p.pro),])
nrow(tab2w[!is.na(tab2w$p.pro) & tab2w$p.pro < 0.05,])
nrow(tab2w[!is.na(tab2w$p.pro) & tab2w$p.pro < 0.05,])/nrow(tab2w[!is.na(tab2w$p.pro),])

reg_kin <- fread(input = paste0(ppnD, "regression/tables/", "change_regression_nonNA/", 
                                        "kinase", "_substrate_regression_", "cptac2p_3can", "_tumor",
                                        "_reg_nonNA", "20", ".txt"), data.table = F)
reg_kin <- markSigKS(regression = reg_kin, sig_thres = 0.1, enzyme_type = "kinase")
reg_kin$regulated <- (reg_kin$fdr_sig & reg_kin$coef_sig)
reg_kin <- reg_kin[reg_kin$regulated,]

reg_pp <- fread(input = paste0(ppnD, "regression/tables/", "change_regression_nonNA/", 
                               "phosphatase", "_substrate_regression_", "cptac2p_3can", "_tumor",
                               "_reg_nonNA", "5", ".txt"), data.table = F)
reg_pp <- markSigKS(regression = reg_pp, sig_thres = 0.1, enzyme_type = "phosphatase")
reg_pp$regulated <- (reg_pp$fdr_sig & reg_pp$coef_sig)
# reg_pp <- reg_pp[reg_pp$regulated,]



tab2w <- merge(tab2w, reg_kin[,  c("pair", "Cancer", "FDR_pho_kin")], by.x = c("pair", "cancer"), by.y = c("pair", "Cancer"), all.x = T)
nrow(unique(tab2w[!is.na(tab2w$FDR_pho_kin) & !is.na(tab2w$p.pro) & tab2w$p.pro > 0.05,]))
tab_trans_activity <- unique(tab2w[!is.na(tab2w$FDR_pho_kin) & !is.na(tab2w$p.pro) & tab2w$p.pro > 0.05,])

table(tab2w$cancer)




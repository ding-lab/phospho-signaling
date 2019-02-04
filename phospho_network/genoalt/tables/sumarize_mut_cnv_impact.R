# Yige Wu @ WashU 2018 Oct
## parse outputs from multiple runs of mutational impact

# source ------------------------------------------------------------------
source("Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R")
source("./cptac2p_analysis/phospho_network/phospho_network_plotting.R")
cptac_phase2process <- "cptac2p_3can"
# input and parse mutaional impact phosphorylation within enzyme-substrate pair ---------------------------
## mutational & CNA combined analysis
mut_cna_pho <- fread(input = paste0("./cptac2p/analysis_results/phospho_network/genoalt/tables/test_mut_cna_impact/mut_cnv_cans.txt"))
mut_cna_pho <- mut_cna_pho[mut_cna_pho$p > 0 & mut_cna_pho$num >=4 & mut_cna_pho$num_control >= 4,]
mut_cna_pho$pair_pro <- paste0(mut_cna_pho$GENE, ":", mut_cna_pho$SUB_GENE)
mut_cna_pho$pair <- paste0(mut_cna_pho$pair_pro, ":", mut_cna_pho$SUB_MOD_RSD)
mut_cna_pho$id <- paste0(mut_cna_pho$pair, ":", mut_cna_pho$cancer)

## mutational only analysis
mut_pho <- fread(input = paste0(ppnD, "genoalt/tables/test_mut_cna_impact_cptac2p/mut_cnv_sig_cans.txt"), data.table = F)
mut_pho$pair_pro <- paste0(mut_pho$GENE, ":", mut_pho$SUB_GENE)
mut_pho$pair <- paste0(mut_pho$pair_pro, ":", mut_cna_pho$SUB_MOD_RSD)
mut_pho$id <- paste0(mut_cna_pho$pair, ":", mut_cna_pho$cancer)

# mut_cis_pro <- fread(input = paste0(ppnD, "genoalt/tables/test_mut_cna_impact_cis_pro_cptac2p/mut_cnv_sig_cans.txt"), data.table = F)

# examine cis CNA effect -----------------------------------------
tab2x <- mut_cna_pho[mut_cna_pho$SELF == "cis" & mut_cna_pho$genoalt_type %in% c("shallow_amp", "shallow_del") & mut_cna_pho$p > 0 & mut_cna_pho$num >=4 & mut_cna_pho$num_control >= 4,]
table(tab2x$genoalt_type)
nrow(tab2x)
tab2x1 <- tab2x[tab2x$p < 0.05,]
nrow(tab2x1)
table(tab2x1[, c("GENE", "cancer")])

# examine trans CNA effect -----------------------------------------
tab2x <- mut_cna_pho[mut_cna_pho$SELF == "trans" & mut_cna_pho$genoalt_type %in% c("shallow_amp", "shallow_del") & mut_cna_pho$p > 0 & mut_cna_pho$num >=4 & mut_cna_pho$num_control >= 4,]
table(tab2x$genoalt_type)
nrow(tab2x)
tab2x1 <- tab2x[tab2x$p < 0.05,]
nrow(tab2x1)

# input regression result -------------------------------------------------
for (reg_nonNA in c(20)) {
  for (enzyme_type in c("kinase")) {
    regression <- fread(input = paste0(ppnD, "regression/tables/", "change_regression_nonNA/", 
                                       enzyme_type, "_substrate_regression_", cptac_phase2process, "_tumor",
                                       "_reg_nonNA", reg_nonNA, ".txt"), data.table = F)
    for (fdr_thres in c(0.1)) {
      regression <- markSigKS(regression = regression, sig_thres = fdr_thres, enzyme_type = enzyme_type)
      regression$regulated <- (regression$fdr_sig & regression$coef_sig)
      regression$id <- paste0(regression$pair, ":", regression$Cancer)
    }
  }
}


# examine overlap with cis regression results -----------------------------
length(regression$id[regression$SELF == "cis" & regression$regulated])
length(regression$id[regression$SELF == "cis" &  regression$regulated & regression$id %in% mut_cna_pho$id])
length(regression$id[regression$SELF == "cis" &  regression$regulated & regression$id %in% mut_cna_pho$id[mut_cna_pho$p < 0.05]])
cis_overlap <- regression$id[regression$SELF == "cis" &  regression$regulated & regression$id %in% mut_cna_pho$id[mut_cna_pho$p < 0.05]]
cis_overlap

# examine overlap with trans regression results -----------------------------
length(regression$id[regression$SELF == "trans" & regression$regulated])
length(regression$id[regression$SELF == "trans" &  regression$regulated & regression$id %in% mut_cna_pho$id])
trans_overlap1 <- regression$id[regression$SELF == "trans" &  regression$regulated & regression$id %in% mut_cna_pho$id[mut_cna_pho$p < 0.05]]
trans_overlap2 <- regression$id[regression$SELF == "trans" &  regression$regulated & regression$id %in% mut_pho_lm$id[mut_pho_lm$P_mut_kin < 0.05]]
intersect(trans_overlap1, trans_overlap2)

trans_overlap1


# examine how many trans pass FDR -----------------------------------------
fdr <- vector(mode = "numeric", length = nrow(mut_cna_pho))
for (cancer in unique(mut_cna_pho$cancer)) {
  for (self in c("cis", "trans")) {
    for (genoalt_type in unique(mut_cna_pho$genoalt_type)) {
      rows_tmp <- which(mut_cna_pho$cancer == cancer & mut_cna_pho$GENE %in% SMGs[[cancer]] & mut_cna_pho$SELF == self)
      fdr[rows_tmp] <- p.adjust(p = mut_cna_pho$p[rows_tmp], method = "fdr")
    }
  }
}
mut_cna_pho$fdr <- fdr
tab2x <- mut_cna_pho[mut_cna_pho$fdr < 0.1 & mut_cna_pho$fdr > 0,]
table(tab2x[, c("SELF", "genoalt_type")])
# genoalt_type
# SELF  mut shallow_amp shallow_del
# cis  10           2           7
tab2x[tab2x$SELF == "cis" & tab2x$genoalt_type == "mut",]
tab2x[tab2x$SELF == "cis" & tab2x$genoalt_type == "shallow_amp",]
tab2x[tab2x$SELF == "cis" & tab2x$genoalt_type == "shallow_del",]


# input and parse mutaional impact on protein within complex pair (BRCA+OV+CO+UCEC) ---------------------------
mut_pro_complex <- fread(input = paste0(ppnD, "genoalt/tables/test_mut_complex_pro/mut_cnv_sig_cans.txt"), data.table = F)
colnames(mut_pro_complex)[colnames(mut_pro_complex) == "geneA"] <- "GENE"
colnames(mut_pro_complex)[colnames(mut_pro_complex) == "geneB"] <- "SUB_GENE"
mut_cnv_cans <- unique(rbind(mut_pro_complex[,col2merge], mut_cnv_cans))

# examine CTNNB1 mutation affected pairs, espacially in CO ----------------------------------
tab2x <- mut_pro_complex[mut_pro_complex$geneA == "CTNNB1" & mut_pro_complex$p < 0.05,]

# examine other mutation affected CTNNB1 pairs, espacially in CO ----------------------------------
tab2x <- mut_pro_complex[mut_pro_complex$geneB == "CTNNB1" & mut_pro_complex$p < 0.05,]

# examine APC mutation affected pairs, espacially in CO ----------------------------------
tab2x <- mut_pro_complex[mut_pro_complex$geneA == "APC" & mut_pro_complex$p < 0.05,]

# examine other mutation affected CTNNB1 pairs, espacially in CO ----------------------------------
tab2x <- mut_pro_complex[mut_pro_complex$geneB == "APC" & mut_pro_complex$p < 0.05,]

# input and parse mutaional impact on phosphorylation within complex pair (BRCA+OV+CO+UCEC) ---------------------------
mut_cnv_cans_complex <- fread(input = paste0(ppnD, "genoalt/tables/test_mut_complex_cptac2p/mut_cnv_sig_cans.txt"), data.table = F)
colnames(mut_cnv_cans_complex)[colnames(mut_cnv_cans_complex) == "geneA"] <- "GENE"
colnames(mut_cnv_cans_complex)[colnames(mut_cnv_cans_complex) == "geneB"] <- "SUB_GENE"
mut_cnv_cans <- unique(rbind(mut_cnv_cans_complex[,col2merge], mut_cnv_cans))



# input and parse mutaional impact on protein within enzyme-substrate pair (BRCA+OV+CO+UCEC) ---------------------------
mut_pro_ks <- fread(input = paste0(ppnD, "genoalt/tables/test_mut_enzyme_substrate_pro/mut_cnv_sig_cans.txt"), data.table = F)
mut_cnv_cans <- unique(rbind(mut_pro_ks[,col2merge], mut_cnv_cans))
write.table(x = mut_cnv_cans, file = paste0(makeOutDir(resultD = resultD), "mut_cnv_cans.txt"), quote = F, sep = "\t", row.names = F)

# adjust pvalue for protein events -----------------------------------------------------------
fdr <- vector(mode = "numeric", length = nrow(mut_cis_pro))
for (cancer in unique(mut_cis_pro$cancer)) {
  rows_tmp <- which(mut_cis_pro$cancer == cancer & mut_cis_pro$GENE %in% SMGs[[cancer]])
  fdr[rows_tmp] <- p.adjust(p = mut_cis_pro$p[rows_tmp], method = "fdr")
}

mut_cis_pro$fdr <- fdr

# adjust pvalue for mutation impact phosphorylation events -----------------------------------------------------------
fdr <- vector(mode = "numeric", length = nrow(mut_pho))
for (cancer in unique(mut_pho$cancer)) {
  for (self in c("cis", "trans")) {
    rows_tmp <- which(mut_pho$cancer == cancer & mut_pho$GENE %in% SMGs[[cancer]] & mut_pho$SELF == self)
    fdr[rows_tmp] <- p.adjust(p = mut_pho$p[rows_tmp], method = "fdr")
  }
}
mut_pho$fdr <- fdr

# adjust pvalue for mutation impact phosphorylation events -----------------------------------------------------------
fdr <- vector(mode = "numeric", length = nrow(cna_pho))
for (cancer in unique(mut_pho$cancer)) {
  for (self in c("cis", "trans")) {
    rows_tmp <- which(mut_pho$cancer == cancer & mut_pho$GENE %in% SMGs[[cancer]] & mut_pho$SELF == self)
    fdr[rows_tmp] <- p.adjust(p = mut_pho$p[rows_tmp], method = "fdr")
  }
}
mut_pho$fdr <- fdr

# write out annotated phosphorylation events -------------------------------------------------------------------
mut_pho <- merge(mut_pho, mut_cis_pro[, c("GENE", "cancer", "p", "meddiff_bottom", "meddiff_upper", "meddiff", "num", "fdr")],
                 by = c("GENE", "cancer"), suffixes = c("", ".pro"), all.x = T)
write.table(x = mut_pho, file = paste0(makeOutDir(resultD = resultD), "cptac2p_mut_cnv.txt"), quote = F, sep = '\t', row.names = F)

# write out cis mutation impact phosphorylation events ------------------------------------
#Among XX cis-regulation events significantly associated with its own mutations (P-value < 0.05), XX (XX%) of them are tested by at least 10 non-synonymous mutations, which indicate the limited statistical power
tab2w <- mut_pho
tab2w <- tab2w[tab2w$SELF == "cis",]
nrow(tab2w)
tab2w <- tab2w[tab2w$p < 0.05,]
nrow(tab2w)
272/3768
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

# examine trans mutation phosphorylation events ------------------------------------
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

## check OV
tab2w <- mut_pho[mut_pho$cancer == "OV" & mut_pho$p < 0.05 & mut_pho$p > 0,]


tab2w <- merge(tab2w, reg_kin[,  c("pair", "Cancer", "FDR_pho_kin")], by.x = c("pair", "cancer"), by.y = c("pair", "Cancer"), all.x = T)
nrow(unique(tab2w[!is.na(tab2w$FDR_pho_kin) & !is.na(tab2w$p.pro) & tab2w$p.pro > 0.05,]))
tab_trans_activity <- unique(tab2w[!is.na(tab2w$FDR_pho_kin) & !is.na(tab2w$p.pro) & tab2w$p.pro > 0.05,])

table(tab2w$cancer)

# examine trans CNA phosphorylation events ------------------------------------
tab2w <- cna_pho[cna_pho$SELF == "trans",]
tab2w <- cna_pho[cna_pho$SELF == "trans" & cna_pho$cancer == "OV" & cna_pho$p < 0.05 & cna_pho$p > 0,]
tab2w <- cna_pho[cna_pho$SELF == "cis" & cna_pho$cancer == "OV" & cna_pho$p < 0.05 & cna_pho$p > 0,]



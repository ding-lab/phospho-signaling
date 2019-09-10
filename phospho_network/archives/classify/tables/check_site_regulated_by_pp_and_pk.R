# Yige Wu @ WashU 2018 Apr
# check phosphosites being regulated by more than 2 kinases
# check phosphosites being regulated by both kinase and phosphatase


# source ------------------------------------------------------------------
source("./cptac2p_analysis/phospho_network/phospho_network_shared.R")
library(dplyr)

# inputs ------------------------------------------------------------------
enzyme_type <- "kinase"
pk_tab <- fread(input = paste0("./cptac2p/analysis_results/phospho_network/regression/tables/regression_omnipath&newworkin_protein_level/", enzyme_type, "_substrate_regression_cptac2p_3can_tumor.txt"), data.table = F)
pp_tab <- fread(input = paste0("./cptac2p/analysis_results/phospho_network/regression/tables/regression_omnipath&newworkin&depod_protein_level/phosphatase_substrate_regression_cptac2p_3can_tumor.txt"), data.table = F)
ptms_sup_tab <- fread(input = "./cptac2p/analysis_results/phospho_network/compile_enzyme_substrate/tables/compile_omnipath_networkin_depod/omnipath_networkin_enzyme_substrate_site_level_union_source_cleaned_addDEPOD.csv", data.table = F)

# set variables -----------------------------------------------------------
fdr_pk <- 0.05
fdr_pp <- 0.1

# reformat --------------------------------------------------------------------
colnames(pk_tab)[1:2] <- c("GENE", "SUB_GENE")
colnames(pp_tab)[1:2] <- c("GENE", "SUB_GENE")
pk_tab$site <-  paste0(pk_tab$SUB_GENE, ":", pk_tab$SUB_MOD_RSD)
pp_tab$site <-  paste0(pp_tab$SUB_GENE, ":", pp_tab$SUB_MOD_RSD)

pk_tab <- markSigSiteCan(regression = pk_tab, sig_thres = fdr_pk, enzyme_type = "kinase")
pp_tab <- markSigSiteCan(regression = pp_tab, sig_thres = fdr_pp, enzyme_type = "phosphatase")
pk_tab$known_enzyme_phosphosite <- (pk_tab$pair %in% ptms_sup_tab$pair)
pp_tab$known_enzyme_phosphosite <- (pp_tab$pair %in% ptms_sup_tab$pair)

trans_pk_tab <- pk_tab[pk_tab$SELF == "trans",]
trans_pp_tab <- pp_tab[pp_tab$SELF == "trans",]

trans_pk_sig_tab <- trans_pk_tab[trans_pk_tab$fdr_sig & trans_pk_tab$coef_sig,]
trans_pp_sig_tab <- trans_pp_tab[trans_pp_tab$fdr_sig & trans_pp_tab$coef_sig,]

trans_pk_sig_multi <- data.frame(table(trans_pk_sig_tab[, c("site", "Cancer")]))
trans_pk_sig_multi <- trans_pk_sig_multi[trans_pk_sig_multi$Freq > 1,]

# check % of all phosphosites being trans-regulated -----------------------------
## number of sites examined in trans-regulation regression model 
trans_tab_can <- group_by(trans_tab, Cancer)
site_tab <- data.frame(summarise(trans_tab_can, num_sites_trans = length(unique(site))))
trans_pk_tab_can <- group_by(trans_pk_tab, Cancer)
tmp <- data.frame(summarise(trans_pk_tab_can, num_sites_trans_pk = length(unique(site))))
site_tab <- merge(site_tab, tmp)
trans_pp_tab_can <- group_by(trans_pp_tab, Cancer)
tmp <- data.frame(summarise(trans_pp_tab_can, num_sites_trans_pp = length(unique(site))))
site_tab <- merge(site_tab, tmp)

## number of sites trans-regulated by certain kinases
trans_pk_sig_tab_can <- group_by(trans_pk_sig_tab, Cancer)
tmp <- data.frame(summarise(trans_pk_sig_tab_can, num_sites_trans_pk_sig = length(unique(site))))
site_tab <- merge(site_tab, tmp)
site_tab$ratio_trans_pk_sig <- site_tab$num_sites_trans_pk_sig/site_tab$num_sites_trans_pk

## number of sites trans-regulated by multiple kinases
trans_pk_sig_multi_can <- group_by(trans_pk_sig_multi, Cancer)
tmp <- data.frame(summarise(trans_pk_sig_multi_can, num_sites_trans_pk_sig_multi = length(unique(site))))
site_tab <- merge(site_tab, tmp)
site_tab$ratio_trans_pk_sig_multi <- site_tab$num_sites_trans_pk_sig_multi/site_tab$num_sites_trans_pk_sig

## number of sites trans-regulated by certain phosphatase
trans_pp_sig_tab_can <- group_by(trans_pp_sig_tab, Cancer)
tmp <- data.frame(summarise(trans_pp_sig_tab_can, num_sites_trans_pp_sig = length(unique(site))))
site_tab <- merge(site_tab, tmp)
site_tab$ratio_trans_pp_sig <- site_tab$num_sites_trans_pp_sig/site_tab$num_sites_trans_pp

## site regulated by both a kinase and a phosphotases
# trans_pk_pp_sig_tab <- unique(merge(trans_pk_sig_tab[, c("SUB_GENE", "SUB_MOD_RSD", "site", "Cancer", "GENE")], 
#                                     trans_pp_sig_tab[, c("SUB_GENE", "SUB_MOD_RSD", "site", "Cancer", "GENE")], 
#                                     by = c("SUB_GENE", "SUB_MOD_RSD", "site", "Cancer"), suffixes = c(".k", ".p")))
trans_pk_pp_sig_tab <- unique(merge(trans_pk_sig_tab[, c("SUB_GENE", "SUB_MOD_RSD", "site", "Cancer", "GENE", "known_enzyme_phosphosite")],
             trans_pp_sig_tab[, c("SUB_GENE", "SUB_MOD_RSD", "site", "Cancer", "GENE", "known_enzyme_phosphosite")],
             by = c("SUB_GENE", "SUB_MOD_RSD", "site", "Cancer"), suffixes = c(".k", ".p")))
trans_pk_pp_sig_tab <- trans_pk_pp_sig_tab[order(trans_pk_pp_sig_tab$site, trans_pk_pp_sig_tab$Cancer),]

trans_pk_pp_sig_tab_can <- group_by(trans_pk_pp_sig_tab, Cancer)
tmp <- data.frame(summarise(trans_pk_pp_sig_tab_can, num_sites_trans_pk_pp_sig = length(unique(site))))
write.table(x = trans_pk_pp_sig_tab, file = paste0(makeOutDir(resultD = resultD), "trans_pk_pp_sig_tab.txt"), quote = F, col.names = T, row.names = F, sep = "\t")
## phospho-proteins regulated by both a kinase and a phosphatases
length(unique(trans_pk_pp_sig_tab$site))
unique(trans_pk_pp_sig_tab$GENE.p)


# check % of all phosphosites being cis-regulated -----------------------------
## number of sites examined in cis-regulation regression model 
cis_tab_can <- group_by(cis_tab, Cancer)
tmp <- data.frame(summarise(cis_tab_can, num_sites_cis = length(unique(site))))
site_tab <- merge(site_tab, tmp)

## number of sites cis-regulated by certain kinases
cis_pk_sig_tab <- cis_sig_tab
cis_pk_sig_tab_can <- group_by(cis_pk_sig_tab, Cancer)
tmp <- data.frame(summarise(cis_pk_sig_tab_can, num_sites_cis_pk_sig = length(unique(site))))
site_tab <- merge(site_tab, tmp)
site_tab$ratio_cis_pk_sig <- site_tab$num_sites_cis_pk_sig/site_tab$num_sites_cis


# check AKT1 --------------------------------------------------------------
trans_pk_sig_known_tab <- trans_pk_sig_tab[trans_pk_sig_tab$known_enzyme_phosphosite,]

## input mutations
somatic <- fread("./cptac2p/mutations/somatic/BRCA.maf", data.table = F)

## input phosphorylation
cancer <- "BRCA"
pho <- fread(input = paste0(cptac_sharedD, cancer,"/",prefix[cancer], "_PHO_formatted_normalized_", "noControl", ".txt"), data.table = F)
pho_head <- formatPhosphosite(phosphosite_vector = pho$Phosphosite, gene_vector = pho$Gene)
## check trans for AKT1
gene <- "AKT1"
cancer <- "BRCA"
trans_pk_sig_gene <- trans_pk_sig_tab[trans_pk_sig_tab$GENE == gene,]
trans_pk_sig_gene_can <- trans_pk_sig_gene[trans_pk_sig_gene$Cancer == cancer,]
trans_pk_sig_gene_known <- trans_pk_sig_gene[trans_pk_sig_gene$known_enzyme_phosphosite,]
trans_pk_sig_gene_known_BRCA <- trans_pk_sig_gene_known[trans_pk_sig_gene_known$Cancer == "BRCA",]
sub1s <- unique(trans_pk_sig_gene_can$SUB_GENE)
sub1s_e <- sub1s[sub1s %in% unique(sup_tab$GENE)]

## get AKT1 mutated samples
somatic_gene <- somatic[somatic$Hugo_Symbol == gene,]
rm(somatic)
mut_partIDs <- unique(str_split_fixed(string = somatic_gene$Tumor_Sample_Barcode, pattern = "_", 2)[,1])
sampIDs <- colnames(pho); sampIDs <- sampIDs[!(sampIDs %in% c("Gene", "Phosphosite"))]
clinical <- fread(input = paste0(cptac_sharedD, "Specimen_Data_20161005_Yige_20180307.txt"), data.table = F)
partIDs <- sampID2partID(sampleID_vector = sampIDs, sample_map = clinical)
mut_sampIDs <- sampIDs[partIDs %in% mut_partIDs]
unmut_sampIDs <- sampIDs[!(sampIDs %in% mut_sampIDs)]

sub1s_trans <- list()
for (sub1 in sub1s_e) {
  sub1_trans <- trans_pk_sig_tab[trans_pk_sig_tab$GENE == sub1 & trans_pk_sig_tab$Cancer == cancer,]
  sub1s_trans[[sub1]] <- sub1_trans
  
  for (i1 in (1:nrow(sub1_trans))) {
    sub2 <- as.vector(sub1_trans$SUB_GENE)[i1]
    rsd2 <- as.vector(sub1_trans$SUB_MOD_RSD)[i1]
    mut_pho <- t(pho[pho_head$SUBSTRATE == sub2 & pho_head$SUB_MOD_RSD == rsd2, mut_sampIDs]); mut_pho <- mut_pho[!is.na(mut_pho)]
    unmut_pho <- t(pho[pho_head$SUBSTRATE == sub2 & pho_head$SUB_MOD_RSD == rsd2, unmut_sampIDs]); unmut_pho[!(is.na(unmut_pho))]
    if (length(mut_pho) >=4 && length(unmut_pho) >=4) {
      stat <- wilcox.test(x = mut_pho, y = unmut_pho)
      if (stat$p.value < 0.2) {
        cat(paste0("Examine ", gene, " mutated vs unmuated samples at second substrate ", sub2, " ", rsd2, " for substrate ", sub1, "\n"))
        print(stat)
        
        sub2_trans <- trans_pk_sig_tab[trans_pk_sig_tab$GENE == sub2 & trans_pk_sig_tab$Cancer == cancer,]
        for (i2 in (1:nrow(sub2_trans))) {
          sub3 <- as.vector(sub2_trans$SUB_GENE)[i2]
          rsd3 <- as.vector(sub2_trans$SUB_MOD_RSD)[i2]
          mut_pho <- t(pho[pho_head$SUBSTRATE == sub3 & pho_head$SUB_MOD_RSD == rsd3, mut_sampIDs]); mut_pho <- mut_pho[!is.na(mut_pho)]
          unmut_pho <- t(pho[pho_head$SUBSTRATE == sub3 & pho_head$SUB_MOD_RSD == rsd3, unmut_sampIDs]); unmut_pho[!(is.na(unmut_pho))]
          
          if (length(mut_pho) >=4 && length(unmut_pho) >=4) {
            stat <- wilcox.test(x = mut_pho, y = unmut_pho)
            if (stat$p.value < 0.2) {
              cat(paste0("Examine ", gene, " mutated vs unmuated samples at third substrate ", sub3, " ", rsd3, " for second substrate ", sub2, "\n"))
              print(stat)
            }
          }
        }
        
      }
    }
  }
}


# check known ks pairs: experimental vs computational ---------------------
trans_pk_known <- unique(trans_pk_tab$pair[!is.na(trans_pk_tab$networkin_score)])
trans_pk_known_exp <- unique(trans_pk_tab$pair[is.infinite(trans_pk_tab$networkin_score) & !is.na(trans_pk_tab$networkin_score)])
trans_pk_known_com <- unique(trans_pk_tab$pair[!is.infinite(trans_pk_tab$networkin_score) & !is.na(trans_pk_tab$networkin_score)])
trans_pk_known_exp_com <- intersect(trans_pk_known_exp, trans_pk_known_com)

trans_pk_exp_tab <- trans_pk_tab[trans_pk_tab$pair %in% trans_pk_known_exp,]
trans_pk_exp_can <- group_by(trans_pk_exp_tab, Cancer)
site_tab <- data.frame(summarise(trans_pk_exp_can, num_sites_trans_pk_exp = length(unique(site))))

trans_pk_sig_exp_tab <- trans_pk_sig_tab[trans_pk_sig_tab$pair %in% trans_pk_known_exp,]
trans_pk_sig_exp_can <- group_by(trans_pk_sig_exp_tab, Cancer)
tmp <- data.frame(summarise(trans_pk_sig_exp_can, num_sites_trans_pk_sig_exp = length(unique(site))))
site_tab <- merge(site_tab, tmp)
site_tab$ratio_sites_trans_pk_sig_exp <- site_tab$num_sites_trans_pk_sig_exp/site_tab$num_sites_trans_pk_exp

trans_pk_com_tab <- trans_pk_tab[trans_pk_tab$pair %in% trans_pk_known_com,]
trans_pk_com_can <- group_by(trans_pk_com_tab, Cancer)
tmp <- data.frame(summarise(trans_pk_com_can, num_sites_trans_pk_com = length(unique(site))))
site_tab <- merge(site_tab, tmp)

trans_pk_sig_com_tab <- trans_pk_sig_tab[trans_pk_sig_tab$pair %in% trans_pk_known_com,]
trans_pk_sig_com_can <- group_by(trans_pk_sig_com_tab, Cancer)
tmp <- data.frame(summarise(trans_pk_sig_com_can, num_sites_trans_pk_sig_com = length(unique(site))))
site_tab <- merge(site_tab, tmp)
site_tab$ratio_sites_trans_pk_sig_com <- site_tab$num_sites_trans_pk_sig_com/site_tab$num_sites_trans_pk_com


# Yige Wu @ WashU 2018 Jul
# how many of mutaion impacted enzyme-substrate pairs are also regualted

# source ------------------------------------------------------------------
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_plotting.R')
source('/Users/yigewu/Box Sync/cptac2p_analysis/phospho_network/phospho_network_shared.R')
library(eulerr)


# set variables -----------------------------------------------------------
mut_sig_thres <- 0.1
cancers_sort <- c("BRCA", "OV", "CO")
reg_sig <- 0.05
diff_sig <- 0.2
diff_log2fc <- 1


# input super tables ---------------------------------------------------------
sup_cans_tab <- NULL
for (cancer in cancers_sort) {
  sup_tab <- fread(paste0(ppnD, "kinase_activity/tables/integrate_enzyme_KSEA_diffexp_sub_FC_plus_cistrans_regression_plus_genoalt/", 
                          cancer, "_KSEA_enzyme_diffexp_substrate.txt"), data.table = F)
  sup_tab$Cancer <- cancer
  sup_cans_tab <- rbind(sup_cans_tab, sup_tab)
}
sup_cans_tab_pk <- sup_cans_tab[sup_cans_tab$enzyme_type == "kinase",]
sup_cans_tab_pk <- markSigSiteCan(sup_cans_tab_pk, sig_thres = reg_sig, enzyme_type = "kinase")
sup_cans_tab_pp <- sup_cans_tab[sup_cans_tab$enzyme_type == "phosphotase",]
sup_cans_tab_pp <- markSigSiteCan(sup_cans_tab_pp, sig_thres = reg_sig, enzyme_type = "phosphotase")
sup_cans_tab <- rbind(sup_cans_tab_pk, sup_cans_tab_pp)
## annotate kinase substrate regulation
sup_cans_tab$regulated <- (sup_cans_tab$coef_sig & sup_cans_tab$fdr_sig)

# create data frame to plot -----------------------------------------------
sup_cans_mut_pho_tab <- sup_cans_tab[!is.na(sup_cans_tab$Mutated_Sample_Size) & sup_cans_tab$p_value < mut_sig_thres,]
sup_cans_mut_pho_tab$id <- paste0(sup_cans_mut_pho_tab$pair, ":", sup_cans_mut_pho_tab$Cancer)
sup_cans_mut_pho_reg_tab <- sup_cans_mut_pho_tab[sup_cans_mut_pho_tab$regulated,]

# venn for phosphosite impact --------------------------------------------------------------------
fn = paste(makeOutDir(resultD = resultD), 'regulated_pairs_with_mut_impacted_phosphosite_venn','.pdf',sep ="")
grid.newpage()
pdf(fn, height = 6, width = 6, useDingbats = FALSE)
all_pairs <- unique(sup_cans_mut_pho_tab$id)
dat <- data.frame(mut_pairs = ( all_pairs %in% sup_cans_mut_pho_tab$id),
                  mut_and_correlated_pairs = (all_pairs %in% sup_cans_mut_pho_reg_tab$id))

dat_euler <- euler(dat)
p <- plot(dat_euler, labels = list(fontsize = 20), quantities = list(fontsize = 20))
grid.draw(p)
dev.off()


# write table for the overlap ---------------------------------------------
sup_cans_mut_pho_reg_tab2p <- unique(sup_cans_mut_pho_reg_tab[, c("Cancer", "GENE", "SUB_GENE", "SUB_MOD_RSD", "Source", "FDR_pho_kin", "coef_pho_kin", "Mutated_Sample_Size", "Non_Mutated_Sample_Size","Fold_Change", "p_value", "id")])
sup_cans_mut_pho_reg_tab2p <- sup_cans_mut_pho_reg_tab2p[order(sup_cans_mut_pho_reg_tab2p$Source),]
sup_cans_mut_pho_reg_tab2p <- sup_cans_mut_pho_reg_tab2p[!duplicated(sup_cans_mut_pho_reg_tab2p$id),]
sup_cans_mut_pho_reg_tab2p$id <- NULL
write.table(x = sup_cans_mut_pho_reg_tab2p, file = paste0(makeOutDir(resultD = resultD), "overlap_mut_pho_regulated.txt"), quote = F, row.names = F, sep = '\t')
